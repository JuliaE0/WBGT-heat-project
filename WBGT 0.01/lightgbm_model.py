
#!/usr/bin/env python
# coding: utf-8

# ## Downscaling modeling with LightGBM

import time
from pathlib import Path

import pandas as pd
import numpy as np

from lightgbm import LGBMRegressor
from sklearn.model_selection import GroupShuffleSplit, GroupKFold
from sklearn.metrics import (
    mean_squared_error,
    mean_absolute_error,
    r2_score
)


# ============================================================
# 1. DEFINE FEATURES
# ============================================================

print("\n" + "=" * 65)
print("STEP 1: DEFINING MODEL FEATURES")
print("=" * 65)

features = [
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

all_features = [
    col for col in features
    if col != "era5_wbgt"
]

print(f"Total candidate predictors: {len(all_features)}")
print("Feature definition complete.")


# ============================================================
# 2. LOAD DATA
# ============================================================

print("\n" + "=" * 65)
print("STEP 2: LOADING CSV DATA")
print("=" * 65)

file_path = Path("CIMIS_WBGT_all_2000_2025_integrated_new.csv")

required_cols = list(dict.fromkeys(
    features + [
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
))

chunks = []

print(f"Input file: {file_path}")
print("Starting CSV loading in chunks of 500,000 rows...")

load_start = time.perf_counter()

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
        chunk["year"], errors="coerce"
    )

    chunk = chunk.loc[
        chunk["year"].between(2000, 2025)
    ]

    if not chunk.empty:
        chunks.append(chunk)

if not chunks:
    raise ValueError("No data found for years 2000–2025.")

df_full = pd.concat(chunks, ignore_index=True)

load_time = time.perf_counter() - load_start

print("\nCSV LOADING COMPLETE")
print(f"Loading time: {load_time / 60:.2f} minutes")
print(f"Total observations: {len(df_full):,}")
print("Years:", sorted(df_full["year"].dropna().unique()))
print(f"Number of stations: {df_full['station_id'].nunique()}")


# ============================================================
# 3. CLEAN DATA AND DEFINE RESIDUAL
# ============================================================

print("\n" + "=" * 65)
print("STEP 3: CLEANING DATA AND DEFINING WBGT RESIDUAL")
print("=" * 65)

print("Converting target variables and predictors to numeric...")

for col in ["wbgt", "era5_wbgt"]:
    df_full[col] = pd.to_numeric(
        df_full[col], errors="coerce"
    )

for col in all_features:
    df_full[col] = pd.to_numeric(
        df_full[col], errors="coerce"
    )

print("Replacing infinite values with missing values...")

numeric_cols = ["wbgt", "era5_wbgt"] + all_features

df_full[numeric_cols] = (
    df_full[numeric_cols]
    .replace([np.inf, -np.inf], np.nan)
)

print("Removing rows without target, baseline, or station ID...")

df_full = df_full.dropna(
    subset=["wbgt", "era5_wbgt", "station_id"]
).copy()

print("Defining WBGT residual...")

df_full["wbgt_residual"] = (
    df_full["wbgt"] - df_full["era5_wbgt"]
)

df_full = df_full.reset_index(drop=True)

print("\nRESIDUAL DEFINITION COMPLETE")
print("WBGT residual = CIMIS WBGT - ERA5 WBGT")
print(f"Observations after cleaning: {len(df_full):,}")
print(f"Unique stations: {df_full['station_id'].nunique()}")
print(f"Years represented: {df_full['year'].nunique()}")


# ============================================================
# 4. STATION-BASED TRAIN-TEST SPLIT
# ============================================================

print("\n" + "=" * 65)
print("STEP 4: CREATING STATION-BASED TRAIN-TEST SPLIT")
print("=" * 65)

print("Starting 80/20 station-based split...")

split_start = time.perf_counter()

splitter = GroupShuffleSplit(
    n_splits=1,
    test_size=0.20,
    random_state=42
)

train_idx, test_idx = next(
    splitter.split(
        df_full,
        groups=df_full["station_id"]
    )
)

train = df_full.iloc[train_idx].copy()
test = df_full.iloc[test_idx].copy()

split_time = time.perf_counter() - split_start

train_stations = set(train["station_id"].unique())
test_stations = set(test["station_id"].unique())

station_overlap = train_stations.intersection(test_stations)

assert not station_overlap, (
    f"Station overlap detected: {station_overlap}"
)

print("\nTRAIN-TEST SPLIT COMPLETE")
print(f"Split time: {split_time:.2f} seconds")
print(f"Training observations: {len(train):,}")
print(f"Testing observations: {len(test):,}")
print(f"Training stations: {len(train_stations)}")
print(f"Testing stations: {len(test_stations)}")
print(f"Station overlap: {station_overlap}")


# ============================================================
# 5. PREPARE LIGHTGBM CROSS-VALIDATION
# ============================================================

print("\n" + "=" * 65)
print("STEP 5: PREPARING 10-FOLD CROSS-VALIDATION")
print("=" * 65)

cv_data = train.copy()

y = cv_data["wbgt_residual"].to_numpy()
groups = cv_data["station_id"].to_numpy()

print(f"CV observations: {len(cv_data):,}")
print(f"CV training stations: {cv_data['station_id'].nunique()}")

assert cv_data["station_id"].nunique() >= 10, (
    "At least 10 training stations are required for 10-fold CV."
)

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

cv = GroupKFold(n_splits=10)

fold_splits = list(
    cv.split(
        np.zeros(len(cv_data)),
        y,
        groups
    )
)

print("10-fold station-based CV folds created.")
print("Each station is assigned to only one validation fold.")


# ============================================================
# 6. RUN 10-FOLD CROSS-VALIDATION
# ============================================================

print("\n" + "=" * 65)
print("STEP 6: STARTING 10-FOLD LIGHTGBM CROSS-VALIDATION")
print("=" * 65)

print("Three predictor sets will be evaluated in each fold.")
print("This may take some time for the full 2000–2025 dataset.\n")

cv_start = time.perf_counter()

results = []

for fold, (train_idx, val_idx) in enumerate(
    fold_splits, start=1
):

    fold_start = time.perf_counter()

    print("\n" + "-" * 65)
    print(f"STARTING FOLD {fold}/10")
    print("-" * 65)

    fold_train = cv_data.iloc[train_idx]
    fold_val = cv_data.iloc[val_idx]

    print(
        f"Training: {len(train_idx):,} observations, "
        f"{fold_train['station_id'].nunique()} stations"
    )

    print(
        f"Validation: {len(val_idx):,} observations, "
        f"{fold_val['station_id'].nunique()} stations"
    )

    # Calculate correlations using ONLY this fold's training data
    print("Calculating fold-specific predictor correlations...")

    corr = (
        fold_train[features + ["wbgt"]]
        .corr(numeric_only=True)["wbgt"]
        .drop("wbgt")
        .sort_values(key=abs, ascending=False)
    )

    features_selected = corr[
        (corr.abs() > 0.5) &
        (corr.index != "era5_wbgt")
    ].index.tolist()

    features_top = corr[
        (corr.abs() > 0.4) &
        (corr.index != "era5_wbgt")
    ].index.tolist()

    models = {
        "LightGBM - All predictors": all_features,
        "LightGBM - Selected predictors": features_selected,
        "LightGBM - Top predictors": features_top
    }

    print(f"Selected predictors: {len(features_selected)}")
    print(f"Top predictors: {len(features_top)}")

    y_train = fold_train["wbgt_residual"].to_numpy()
    y_val = fold_val["wbgt_residual"].to_numpy()

    wbgt_val = fold_val["wbgt"].to_numpy()
    era5_val = fold_val["era5_wbgt"].to_numpy()

    for model_name, feature_cols in models.items():

        print(f"\nRunning: {model_name}")

        if len(feature_cols) == 0:
            print(f"Skipping {model_name}: no predictors selected.")
            continue

        X_train = (
            fold_train[feature_cols]
            .replace([np.inf, -np.inf], np.nan)
        )

        X_val = (
            fold_val[feature_cols]
            .replace([np.inf, -np.inf], np.nan)
        )

        model = LGBMRegressor(**lgbm_params)

        model_start = time.perf_counter()

        model.fit(X_train, y_train)

        training_time = time.perf_counter() - model_start

        predicted_residual = model.predict(X_val)
        predicted_wbgt = era5_val + predicted_residual

        # Corrected WBGT metrics
        wbgt_rmse = np.sqrt(
            mean_squared_error(wbgt_val, predicted_wbgt)
        )
        wbgt_mae = mean_absolute_error(
            wbgt_val, predicted_wbgt
        )
        wbgt_bias = np.mean(predicted_wbgt - wbgt_val)
        wbgt_r2 = r2_score(wbgt_val, predicted_wbgt)

        # ERA5 baseline metrics
        era5_rmse = np.sqrt(
            mean_squared_error(wbgt_val, era5_val)
        )
        era5_mae = mean_absolute_error(wbgt_val, era5_val)
        era5_bias = np.mean(era5_val - wbgt_val)
        era5_r2 = r2_score(wbgt_val, era5_val)

        # Residual metrics
        residual_rmse = np.sqrt(
            mean_squared_error(y_val, predicted_residual)
        )
        residual_mae = mean_absolute_error(
            y_val, predicted_residual
        )
        residual_bias = np.mean(predicted_residual - y_val)
        residual_r2 = r2_score(y_val, predicted_residual)

        results.append({
            "Model": model_name,
            "Fold": fold,
            "Train_rows": len(train_idx),
            "Validation_rows": len(val_idx),
            "Train_stations": fold_train["station_id"].nunique(),
            "Validation_stations": fold_val["station_id"].nunique(),
            "Training_time_sec": training_time,
            "WBGT_RMSE": wbgt_rmse,
            "WBGT_MAE": wbgt_mae,
            "WBGT_Bias": wbgt_bias,
            "WBGT_R2": wbgt_r2,
            "ERA5_RMSE": era5_rmse,
            "ERA5_MAE": era5_mae,
            "ERA5_Bias": era5_bias,
            "ERA5_R2": era5_r2,
            "Residual_RMSE": residual_rmse,
            "Residual_MAE": residual_mae,
            "Residual_Bias": residual_bias,
            "Residual_R2": residual_r2
        })

        print(
            f"Completed {model_name} | "
            f"WBGT RMSE: {wbgt_rmse:.4f} °C | "
            f"WBGT R²: {wbgt_r2:.4f} | "
            f"Training time: {training_time:.2f} sec"
        )

    fold_time = time.perf_counter() - fold_start

    print(f"\nFOLD {fold}/10 COMPLETE")
    print(f"Fold elapsed time: {fold_time / 60:.2f} minutes")

cv_time = time.perf_counter() - cv_start

print("\n" + "=" * 65)
print("10-FOLD CROSS-VALIDATION COMPLETE")
print("=" * 65)
print(f"Total CV time: {cv_time / 60:.2f} minutes")

results_df = pd.DataFrame(results)

print(f"Total model-fold results: {len(results_df)}")
display(results_df)


# ============================================================
# 7. SUMMARIZE AND SAVE RESULTS
# ============================================================

print("\n" + "=" * 65)
print("STEP 7: SUMMARIZING AND SAVING RESULTS")
print("=" * 65)

metrics = [
    "WBGT_RMSE",
    "WBGT_MAE",
    "WBGT_Bias",
    "WBGT_R2",
    "ERA5_RMSE",
    "ERA5_MAE",
    "ERA5_Bias",
    "ERA5_R2",
    "Residual_RMSE",
    "Residual_MAE",
    "Residual_Bias",
    "Residual_R2",
    "Training_time_sec"
]

summary = (
    results_df
    .groupby("Model")[metrics]
    .agg(["mean", "std"])
)

print("\nCROSS-VALIDATION SUMMARY:")
display(summary.round(4))

print("\nSaving fold-level results...")

results_df.to_csv(
    "lightgbm_2000_2025_cv_fold_results.csv",
    index=False
)

print("Saving summary results...")

summary.to_csv(
    "lightgbm_2000_2025_cv_summary.csv"
)

print("\nRESULTS SAVED SUCCESSFULLY")
print(" - lightgbm_2000_2025_cv_fold_results.csv")
print(" - lightgbm_2000_2025_cv_summary.csv")

print("\n" + "=" * 65)
print("MODELING WORKFLOW COMPLETE")
print("=" * 65)

# Next steps:
# 1. Choose the best-performing predictor set.
# 2. Train the final LightGBM model using all 80% training stations.
# 3. Predict residuals for the untouched 20% test stations.
# 4. Add predicted residuals to ERA5 WBGT.
# 5. Compare corrected WBGT against observed CIMIS WBGT and the original ERA5 baseline.
