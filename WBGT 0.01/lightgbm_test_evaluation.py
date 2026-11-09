
#!/usr/bin/env python
# coding: utf-8

"""
Final LightGBM residual-correction model and independent test evaluation.

1. Recreates the original 80/20 station-based split.
2. Calculates observed WBGT residuals.
3. Trains LightGBM to predict residuals using training stations.
4. Predicts residuals at held-out test stations.
5. Reconstructs predicted CIMIS WBGT.
6. Evaluates predicted CIMIS WBGT against observed CIMIS WBGT.
7. Saves predictions, metrics, model, and feature importance.
"""

import time
from pathlib import Path

import numpy as np
import pandas as pd
import joblib

from lightgbm import LGBMRegressor

from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import (
    mean_squared_error,
    mean_absolute_error,
    r2_score
)


# ============================================================
# 1. DEFINE PATHS
# ============================================================

print("\n" + "=" * 65)
print("STEP 1: DEFINING INPUT AND OUTPUT PATHS")
print("=" * 65)

file_path = Path(
    "CIMIS_WBGT_all_2000_2025_integrated_new.csv"
)

output_dir = Path("lightgbm_test_results")
output_dir.mkdir(parents=True, exist_ok=True)

if not file_path.exists():
    raise FileNotFoundError(
        f"Input CSV not found: {file_path.resolve()}"
    )

print(f"Input CSV: {file_path.resolve()}")
print(f"Output directory: {output_dir.resolve()}")


# ============================================================
# 2. DEFINE PREDICTORS
# ============================================================

print("\n" + "=" * 65)
print("STEP 2: DEFINING MODEL PREDICTORS")
print("=" * 65)

features = [
    # ERA5 WBGT is used as the residual reference,
    # NOT as a LightGBM predictor.
    "era5_wbgt",

    # Hourly meteorology
    "solar", "Tair", "relhum", "speed",
    "elev_m", "cza", "fdir", "pres",
    "zspeed", "dT",

    "urban",

    # ASTER terrain statistics
    "aster_elev_mean", "aster_elev_sd",
    "aster_elev_min", "aster_elev_max",
    "aster_elev_range", "aster_elev_median",
    "aster_elev_q25", "aster_elev_q75",
    "aster_elev_iqr",

    # NLCD land-cover variables
    "nlcd_urban", "nlcd_imp_frac", "nlcd_tree_frac",
    "nlcd_LC11_frac", "nlcd_LC21_frac",
    "nlcd_LC22_frac", "nlcd_LC23_frac",
    "nlcd_LC24_frac", "nlcd_LC31_frac",
    "nlcd_LC41_frac", "nlcd_LC42_frac",
    "nlcd_LC43_frac", "nlcd_LC52_frac",
    "nlcd_LC71_frac", "nlcd_LC81_frac",
    "nlcd_LC82_frac", "nlcd_LC90_frac",
    "nlcd_LC95_frac", "nlcd_developed_frac",
    "nlcd_forest_frac", "nlcd_shrub_frac",
    "nlcd_water_frac", "nlcd_crop_frac",
    "nlcd_wetland_frac",

    # Vegetation
    "modis_ndvi",

    # GridMET
    "gridmet_rmax", "gridmet_rmin", "gridmet_sph",
    "gridmet_srad", "gridmet_tmmn", "gridmet_tmmx",
    "gridmet_vs",

    # Daymet
    "daymet_vp", "daymet_tmin", "daymet_tmax",
    "daymet_srad",
]

# Exclude ERA5 WBGT from the LightGBM predictors.
# It is only used to calculate residuals and reconstruct WBGT.
all_features = [
    col for col in features
    if col != "era5_wbgt"
]

print(f"Number of LightGBM predictors: {len(all_features)}")
print("ERA5 WBGT is excluded from the predictor set.")


# ============================================================
# 3. LOAD DATA
# ============================================================

print("\n" + "=" * 65)
print("STEP 3: LOADING DATA")
print("=" * 65)

metadata_cols = [
    "station_id",
    "year",
    "month",
    "date_local",
    "datetime_utc",
    "lat",
    "lon",
    "wbgt",
    "era5_wbgt"
]

required_cols = list(dict.fromkeys(
    all_features + metadata_cols
))

print("Loading CSV in chunks of 500,000 rows...")

load_start = time.perf_counter()

chunks = []

for chunk_number, chunk in enumerate(
    pd.read_csv(
        file_path,
        usecols=required_cols,
        chunksize=500_000,
        dtype={"station_id": "string"}
    ),
    start=1
):

    print(f"Loading chunk {chunk_number}...")

    chunk["year"] = pd.to_numeric(
        chunk["year"],
        errors="coerce"
    )

    chunk = chunk.loc[
        chunk["year"].between(2000, 2025)
    ]

    if not chunk.empty:
        chunks.append(chunk)

if not chunks:
    raise ValueError("No observations found for 2000–2025.")

df = pd.concat(chunks, ignore_index=True)

load_time = time.perf_counter() - load_start

print("\nDATA LOADING COMPLETE")
print(f"Observations loaded: {len(df):,}")
print(f"Unique stations: {df['station_id'].nunique()}")
print(f"Loading time: {load_time / 60:.2f} minutes")


# ============================================================
# 4. CLEAN DATA AND CALCULATE OBSERVED RESIDUAL
# ============================================================

print("\n" + "=" * 65)
print("STEP 4: CLEANING DATA AND CALCULATING RESIDUAL")
print("=" * 65)

numeric_cols = ["wbgt", "era5_wbgt"] + all_features

print("Converting target, reference, and predictors to numeric...")

for col in numeric_cols:
    df[col] = pd.to_numeric(
        df[col],
        errors="coerce"
    )

print("Replacing infinite values with missing values...")

df[numeric_cols] = df[numeric_cols].replace(
    [np.inf, -np.inf],
    np.nan
)

print("Removing observations without observed WBGT, ERA5 WBGT, or station ID...")

df = df.dropna(
    subset=["wbgt", "era5_wbgt", "station_id"]
).copy()

# Calculate observed residual
df["wbgt_residual"] = (
    df["wbgt"] - df["era5_wbgt"]
)

df = df.reset_index(drop=True)

print("\nDATA CLEANING COMPLETE")
print(f"Observations after cleaning: {len(df):,}")
print(f"Unique stations: {df['station_id'].nunique()}")
print(f"Years represented: {df['year'].nunique()}")
print("Residual = CIMIS WBGT - ERA5 WBGT")


# ============================================================
# 5. RECREATE ORIGINAL 80/20 STATION-BASED SPLIT
# ============================================================

print("\n" + "=" * 65)
print("STEP 5: CREATING STATION-BASED TRAIN-TEST SPLIT")
print("=" * 65)

splitter = GroupShuffleSplit(
    n_splits=1,
    test_size=0.20,
    random_state=42
)

train_idx, test_idx = next(
    splitter.split(
        df,
        groups=df["station_id"]
    )
)

train = df.iloc[train_idx].copy()
test = df.iloc[test_idx].copy()

train_stations = set(train["station_id"].unique())
test_stations = set(test["station_id"].unique())

assert train_stations.isdisjoint(test_stations), (
    "ERROR: Training and testing stations overlap."
)

print("\nTRAIN-TEST SPLIT COMPLETE")
print(f"Training observations: {len(train):,}")
print(f"Testing observations: {len(test):,}")
print(f"Training stations: {len(train_stations)}")
print(f"Testing stations: {len(test_stations)}")
print("Station overlap: 0")


# ============================================================
# 6. TRAIN FINAL LIGHTGBM RESIDUAL MODEL
# ============================================================

print("\n" + "=" * 65)
print("STEP 6: TRAINING FINAL LIGHTGBM RESIDUAL MODEL")
print("=" * 65)

lgbm_params = {
    "objective": "regression",
    "n_estimators": 500,
    "learning_rate": 0.05,
    "num_leaves": 31,
    "max_depth": -1,
    "min_child_samples": 50,
    "subsample": 0.8,
    "subsample_freq": 1,
    "colsample_bytree": 0.8,
    "reg_lambda": 1.0,
    "random_state": 42,
    "n_jobs": -1,
    "verbosity": -1
}

X_train = train[all_features].replace(
    [np.inf, -np.inf], np.nan
)

y_train = train["wbgt_residual"].to_numpy()

print(f"Training observations: {len(X_train):,}")
print(f"Training stations: {train['station_id'].nunique()}")
print(f"Number of predictors: {len(all_features)}")

final_model = LGBMRegressor(**lgbm_params)

training_start = time.perf_counter()

final_model.fit(
    X_train,
    y_train
)

training_time = time.perf_counter() - training_start

print("\nFINAL RESIDUAL MODEL TRAINING COMPLETE")
print(f"Training time: {training_time / 60:.2f} minutes")


# ============================================================
# 7. PREDICT RESIDUALS AND RECONSTRUCT CIMIS WBGT
# ============================================================

print("\n" + "=" * 65)
print("STEP 7: PREDICTING CIMIS WBGT AT TEST STATIONS")
print("=" * 65)

X_test = test[all_features].replace(
    [np.inf, -np.inf], np.nan
)

# Observed CIMIS WBGT
observed_wbgt = test["wbgt"].to_numpy()

# ERA5 WBGT used only as the residual reference
era5_wbgt = test["era5_wbgt"].to_numpy()

# Observed residual
observed_residual = test["wbgt_residual"].to_numpy()

# Predict residual using LightGBM
predicted_residual = final_model.predict(X_test)

# Reconstruct predicted CIMIS WBGT
predicted_wbgt = (
    era5_wbgt + predicted_residual
)

print("Test predictions complete.")
print(f"Test observations: {len(test):,}")
print(f"Test stations: {test['station_id'].nunique()}")


# ============================================================
# 8. EVALUATE PREDICTED CIMIS WBGT
# ============================================================

print("\n" + "=" * 65)
print("STEP 8: EVALUATING PREDICTED CIMIS WBGT")
print("=" * 65)


def calculate_metrics(observed, predicted):

    return {
        "RMSE": np.sqrt(
            mean_squared_error(observed, predicted)
        ),
        "MAE": mean_absolute_error(
            observed, predicted
        ),
        "Bias": np.mean(predicted - observed),
        "R2": r2_score(observed, predicted)
    }


# Evaluate reconstructed WBGT directly against CIMIS
wbgt_metrics = calculate_metrics(
    observed_wbgt,
    predicted_wbgt
)

test_metrics = pd.DataFrame([
    {
        "Model": "LightGBM residual-corrected WBGT",
        **wbgt_metrics
    }
])

print("\nTEST-SET CIMIS WBGT PERFORMANCE:")
print(test_metrics.round(4).to_string(index=False))

# Residual prediction diagnostics
residual_metrics = calculate_metrics(
    observed_residual,
    predicted_residual
)

print("\nRESIDUAL PREDICTION DIAGNOSTICS:")
print(
    pd.DataFrame([residual_metrics])
    .round(4)
    .to_string(index=False)
)


# ============================================================
# 9. SAVE TEST-STATION PREDICTIONS
# ============================================================

print("\n" + "=" * 65)
print("STEP 9: SAVING TEST PREDICTIONS")
print("=" * 65)

metadata_to_save = [
    col for col in [
        "station_id",
        "year",
        "month",
        "date_local",
        "datetime_utc",
        "lat",
        "lon"
    ]
    if col in test.columns
]

test_predictions = test[metadata_to_save].copy()

test_predictions["wbgt_observed"] = observed_wbgt
test_predictions["era5_wbgt_reference"] = era5_wbgt
test_predictions["observed_residual"] = observed_residual
test_predictions["predicted_residual"] = predicted_residual
test_predictions["wbgt_predicted"] = predicted_wbgt

# Prediction error
test_predictions["prediction_error"] = (
    predicted_wbgt - observed_wbgt
)

prediction_path = (
    output_dir / "lightgbm_2000_2025_test_predictions.csv"
)

test_predictions.to_csv(
    prediction_path,
    index=False
)

print(f"Saved: {prediction_path}")


# ============================================================
# 10. SAVE PERFORMANCE METRICS
# ============================================================

print("\n" + "=" * 65)
print("STEP 10: SAVING PERFORMANCE METRICS")
print("=" * 65)

metrics_path = (
    output_dir / "lightgbm_2000_2025_test_metrics.csv"
)

test_metrics.to_csv(
    metrics_path,
    index=False
)

residual_metrics_path = (
    output_dir / "lightgbm_2000_2025_test_residual_metrics.csv"
)

pd.DataFrame([residual_metrics]).to_csv(
    residual_metrics_path,
    index=False
)

print(f"Saved: {metrics_path}")
print(f"Saved: {residual_metrics_path}")


# ============================================================
# 11. SAVE FINAL MODEL
# ============================================================

print("\n" + "=" * 65)
print("STEP 11: SAVING FINAL MODEL")
print("=" * 65)

model_path = (
    output_dir / "lightgbm_2000_2025_final_model.pkl"
)

joblib.dump(
    final_model,
    model_path
)

print(f"Saved: {model_path}")


# ============================================================
# 12. SAVE FEATURE IMPORTANCE
# ============================================================

print("\n" + "=" * 65)
print("STEP 12: SAVING FEATURE IMPORTANCE")
print("=" * 65)

importance_df = pd.DataFrame({
    "Feature": all_features,
    "Importance": final_model.feature_importances_
}).sort_values(
    "Importance",
    ascending=False
)

importance_path = (
    output_dir / "lightgbm_2000_2025_feature_importance.csv"
)

importance_df.to_csv(
    importance_path,
    index=False
)

print(f"Saved: {importance_path}")

print("\nTOP 15 FEATURES:")
print(importance_df.head(15).to_string(index=False))


# ============================================================
# 13. FINISH
# ============================================================

print("\n" + "=" * 65)
print("RESIDUAL-CORRECTION MODELING COMPLETE")
print("=" * 65)

print(f"Final training time: {training_time / 60:.2f} minutes")
print(f"Test observations: {len(test):,}")
print(f"Test stations: {test['station_id'].nunique()}")

print("\nFINAL TEST WBGT METRICS:")
for metric, value in wbgt_metrics.items():
    print(f"{metric}: {value:.4f}")

print(f"\nAll outputs saved in: {output_dir.resolve()}")
