# -*- coding: utf-8 -*-
"""
ME4 Module – XGBoost 3-Year Mortality Prediction + SHAP Explanation
=====================================================================
Input:  ML_dataset_ME4_GeneSymbols_Survival.csv
        (exported by TCGA_COAD_READ_MYC_Analysis.R, section 17A)

Outputs:
  Fig_ML_ROC_3yr_ME4_600dpi.tiff   – ROC curve
  Fig_ML_SHAP_3yr_ME4_600dpi.tiff  – SHAP beeswarm summary

Features: ME4 WGCNA module genes (HGNC symbols) from full-cohort WGCNA.
Label:    death_3yr  (1 = died within 3 years, 0 = otherwise)
"""

# ── 0. Dependencies ────────────────────────────────────────────────────────────
# pip install xgboost shap scikit-learn pandas numpy matplotlib

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")          # non-interactive backend — safe on all platforms
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score, roc_curve, average_precision_score
from xgboost import XGBClassifier
import shap

# ── 1. Load data ───────────────────────────────────────────────────────────────
df = pd.read_csv("ML_dataset_ME4_GeneSymbols_Survival.csv", index_col=0)

print(f"Loaded dataset: {df.shape[0]} samples × {df.shape[1]} columns")
print("Last 5 columns:", df.columns[-5:].tolist())

# ── 2. Define 3-year mortality label ──────────────────────────────────────────
THRESHOLD_DAYS = 1095   # 3 years

df = df.dropna(subset=["OS_time", "OS_event"])

df["death_3yr"] = np.where(
    (df["OS_event"] == 1) & (df["OS_time"] <= THRESHOLD_DAYS),
    1, 0
)

print("\n3-year mortality label distribution:")
print(df["death_3yr"].value_counts())
print(f"Event rate: {df['death_3yr'].mean():.1%}")

# ── 3. Features and target ─────────────────────────────────────────────────────
NON_FEATURE_COLS = ["MYC_group", "MYC_label", "OS_time", "OS_event", "death_3yr"]

X = df.drop(columns=[c for c in NON_FEATURE_COLS if c in df.columns])
y = df["death_3yr"]

print(f"\nFeature matrix: {X.shape[0]} samples × {X.shape[1]} ME4 genes")

# ── 4. Stratified train / test split ──────────────────────────────────────────
X_train, X_test, y_train, y_test = train_test_split(
    X, y,
    test_size    = 0.25,
    random_state = 42,
    stratify     = y
)

print(f"\nTrain: {X_train.shape[0]}  Test: {X_test.shape[0]}")
print(f"Train event rate: {y_train.mean():.1%}  "
      f"Test event rate: {y_test.mean():.1%}")

# ── 5. XGBoost classifier ──────────────────────────────────────────────────────
# scale_pos_weight corrects for class imbalance (ratio of negatives to positives)
neg, pos = (y_train == 0).sum(), (y_train == 1).sum()

model = XGBClassifier(
    n_estimators      = 400,
    max_depth         = 4,
    learning_rate     = 0.05,
    subsample         = 0.8,
    colsample_bytree  = 0.8,
    scale_pos_weight  = neg / pos,      # handles imbalanced labels
    eval_metric       = "logloss",
    random_state      = 42,
    use_label_encoder = False
)

model.fit(X_train, y_train)

# ── 6. Evaluation ─────────────────────────────────────────────────────────────
y_pred_proba = model.predict_proba(X_test)[:, 1]

auc  = roc_auc_score(y_test, y_pred_proba)
auprc = average_precision_score(y_test, y_pred_proba)
fpr, tpr, _ = roc_curve(y_test, y_pred_proba)

# 5-fold cross-validated AUC on the full dataset for robustness check
cv_auc = cross_val_score(
    model, X, y,
    cv      = StratifiedKFold(n_splits=5, shuffle=True, random_state=42),
    scoring = "roc_auc"
)
print(f"\nTest  AUC:  {auc:.3f}")
print(f"Test  AUPRC: {auprc:.3f}")
print(f"5-fold CV AUC: {cv_auc.mean():.3f} ± {cv_auc.std():.3f}")

# ── 7. ROC curve figure ────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(fpr, tpr, color="#D73027", linewidth=2,
        label=f"XGBoost  AUC = {auc:.3f}\nAUPRC = {auprc:.3f}")
ax.plot([0, 1], [0, 1], linestyle="--", color="grey", linewidth=1)
ax.set_xlabel("False Positive Rate", fontsize=13)
ax.set_ylabel("True Positive Rate", fontsize=13)
ax.set_title("ROC – ME4 Module Predicting 3-Year Mortality (CRC)",
             fontsize=13, fontweight="bold")
ax.legend(fontsize=11, frameon=True)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig("Fig_ML_ROC_3yr_ME4_600dpi.tiff", dpi=600, bbox_inches="tight")
plt.close()
print("Saved: Fig_ML_ROC_3yr_ME4_600dpi.tiff")

# ── 8. SHAP explanation ────────────────────────────────────────────────────────
explainer   = shap.Explainer(model, X_train)
shap_values = explainer(X_test)

fig, ax = plt.subplots(figsize=(9, 7))
shap.plots.beeswarm(shap_values, max_display=20, show=False)
ax = plt.gca()
ax.set_title("SHAP Summary – ME4 Drivers of 3-Year Mortality (CRC)",
             fontsize=13, fontweight="bold", pad=12)
plt.tight_layout()
plt.savefig("Fig_ML_SHAP_3yr_ME4_600dpi.tiff", dpi=600, bbox_inches="tight")
plt.close()
print("Saved: Fig_ML_SHAP_3yr_ME4_600dpi.tiff")

print("\nDone.")
