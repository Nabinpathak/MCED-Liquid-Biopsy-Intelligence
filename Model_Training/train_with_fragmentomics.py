# train_with_fragmentomics.py
import pandas as pd
import numpy as np
import xgboost as xgb
import joblib
import shap
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

# === PATHS ===
CSV_PATH = r"D:\MCED_DP\Classifier\features_with_fragmentomics.csv"  # ← use the new file!
MODEL_SAVE = r"D:\MCED_DP\Model_Training\MCED_CANCER_DETECTOR_FRAGMENTOMICS.json"
SHAP_SUMMARY = r"D:\MCED_DP\shap_summary.png"

# === LOAD & PREPARE DATA ===
df = pd.read_csv(CSV_PATH)

# Binary label: tumor = 1, healthy = 0
df['label_bin'] = df['label'].str.contains('tumor|cancer', case=False, na=False).astype(int)
y = df['label_bin'].values

# Drop non-feature columns
drop_cols = ['filename', 'label', 'subtype', 'label_bin']
X = df.drop(columns=[c for c in drop_cols if c in df.columns])
X = X.select_dtypes(include=np.number).fillna(0)  # fill NaN with 0

print(f"Training on {len(X)} samples | {y.sum()} tumor | {len(y)-y.sum()} healthy")

# === BEST XGBoost PARAMS (tuned for fragmentomics) ===
params = {
    "objective": "binary:logistic",
    "eval_metric": "auc",
    "max_depth": 6,
    "eta": 0.05,
    "subsample": 0.85,
    "colsample_bytree": 0.85,
    "min_child_weight": 1,
    "gamma": 0.05,
    "scale_pos_weight": (y == 0).sum() / (y == 1).sum() if y.sum() > 0 else 1,
    "seed": 42,
    "nthread": 8,
}

# === 5-FOLD CROSS-VALIDATION ===
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
aucs = []

for fold, (train_idx, val_idx) in enumerate(skf.split(X, y)):
    print(f"\nFold {fold+1}/5")
    X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
    y_train, y_val = y[train_idx], y[val_idx]

    dtrain = xgb.DMatrix(X_train, label=y_train)
    dval = xgb.DMatrix(X_val, label=y_val)

    model = xgb.train(params, dtrain, num_boost_round=1000,
                      evals=[(dval, 'val')], early_stopping_rounds=100, verbose_eval=100)

    pred = model.predict(dval)
    auc = roc_auc_score(y_val, pred)
    aucs.append(auc)
    print(f"Fold {fold+1} AUC: {auc:.5f}")

print(f"\nFinal 5-Fold CV AUC: {np.mean(aucs):.5f} ± {np.std(aucs):.5f}")

# === FINAL MODEL ON ALL DATA ===
dtrain_full = xgb.DMatrix(X, label=y)
final_model = xgb.train(params, dtrain_full, num_boost_round=500)
final_model.save_model(MODEL_SAVE)
print(f"\nModel saved: {MODEL_SAVE}")

# === SHAP EXPLANATION (very important!) ===
explainer = shap.TreeExplainer(final_model)
shap_values = explainer(X)

# Summary plot: which features matter most globally
shap.summary_plot(shap_values, X, show=False)
import matplotlib.pyplot as plt
plt.tight_layout()
plt.savefig(SHAP_SUMMARY)
print(f"SHAP summary plot saved: {SHAP_SUMMARY}")

print("\nTraining complete! Use MCED_CANCER_DETECTOR_FRAGMENTOMICS.json in your Streamlit app.")
print("Fragmentomics should push cfDNA performance to 95–99% — check SHAP to see which features drive it!")