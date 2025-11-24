# app.py - FINAL WORKING VERSION (November 21, 2025)
import streamlit as st
import pandas as pd
import numpy as np
import xgboost as xgb
import os

# === PAGE CONFIG ===
st.set_page_config(page_title="MCED Cancer Detector", page_icon="🧬", layout="centered")
st.title("🧬 MCED Cancer Detector")
st.markdown("### Upload a bedGraph file → Instant Cancer vs Healthy prediction")

# === LOAD MODEL (JSON - NEVER BREAKS) ===
@st.cache_resource
def load_model():
    model_path = "MCED_CANCER_DETECTOR.json"
    if not os.path.exists(model_path):
        st.error(f"Model not found: {model_path}\nPlace it in the same folder as app.py")
        return None
    return xgb.Booster(model_file=model_path)

model = load_model()
if model is None:
    st.stop()

# === FEATURE COLUMNS (exact match from training) ===
FEATURE_COLUMNS = [
    'n_CpG', 'global_mean', 'global_median', 'global_std', 'pct_hyper_0.8', 'pct_hypo_0.2',
    'chr1_count', 'chr1_mean', 'chr2_count', 'chr2_mean', 'chr3_count', 'chr3_mean',
    'chr4_count', 'chr4_mean', 'chr5_count', 'chr5_mean', 'chr6_count', 'chr6_mean',
    'chr7_count', 'chr7_mean', 'chr8_count', 'chr8_mean', 'chr9_count', 'chr9_mean',
    'chr10_count', 'chr10_mean', 'chr11_count', 'chr11_mean', 'chr12_count', 'chr12_mean',
    'chr13_count', 'chr13_mean', 'chr14_count', 'chr14_mean', 'chr15_count', 'chr15_mean',
    'chr16_count', 'chr16_mean', 'chr17_count', 'chr17_mean', 'chr18_count', 'chr18_mean',
    'chr19_count', 'chr19_mean', 'chr20_count', 'chr20_mean', 'chr21_count', 'chr21_mean',
    'chr22_count', 'chr22_mean', 'chrX_count', 'chrX_mean', 'chrY_count', 'chrY_mean'
]

# === EXTRACT FEATURES FROM BEDGRAPH ===
def extract_features_from_bedgraph(file_content):
    lines = file_content.decode("utf-8").splitlines()
    global_vals = []
    global_weights = []
    chrom_data = {f"chr{i}": {"vals": [], "w": []} for i in list(range(1,23)) + ["X", "Y"]}

    for line in lines:
        if not line.strip() or line.startswith(("track", "#")):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 4:
            continue
        chrom, start, end, methyl = parts[0], parts[1], parts[2], parts[3]
        try:
            methyl = float(methyl)
            if methyl > 1: methyl /= 100
        except:
            continue
        weight = 1

        global_vals.append(methyl)
        global_weights.append(weight)
        if chrom in chrom_data:
            chrom_data[chrom]["vals"].append(methyl)
            chrom_data[chrom]["w"].append(weight)

    if not global_vals:
        return None

    def w_mean(v, w):
        return np.average(v, weights=w) if sum(w) > 0 else np.nan

    feats = {}
    total_w = sum(global_weights)
    feats["n_CpG"] = total_w
    feats["global_mean"] = w_mean(global_vals, global_weights)
    feats["global_median"] = np.median(global_vals)
    feats["global_std"] = np.sqrt(np.average((np.array(global_vals) - feats["global_mean"])**2, weights=global_weights)) if total_w > 1 else 0
    arr = np.array(global_vals)
    w = np.array(global_weights)
    feats["pct_hyper_0.8"] = np.sum((arr >= 0.8) * w) / total_w if total_w > 0 else 0
    feats["pct_hypo_0.2"] = np.sum((arr <= 0.2) * w) / total_w if total_w > 0 else 0

    for c in chrom_data:
        prefix = c.replace("chr", "") if c not in ["chrX", "chrY"] else c[3:]
        count = sum(chrom_data[c]["w"])
        feats[f"chr{prefix}_count"] = count
        feats[f"chr{prefix}_mean"] = w_mean(chrom_data[c]["vals"], chrom_data[c]["w"])

    return pd.DataFrame([feats])[FEATURE_COLUMNS]

# === UPLOAD & PREDICT ===
uploaded_file = st.file_uploader("Choose a bedGraph file", type=["bedGraph", "bedgraph", "txt"])

if uploaded_file is not None:
    with st.spinner("Extracting features & predicting..."):
        features_df = extract_features_from_bedgraph(uploaded_file.getvalue())
        
    if features_df is None or features_df.isna().all().all():
        st.error("No valid CpG data found in file. Check format.")
    else:
        # PREDICTION WITH JSON MODEL (NO MORE ERRORS!)
        dtest = xgb.DMatrix(features_df)
        cancer_prob = float(model.predict(dtest)[0])
        healthy_prob = 1 - cancer_prob

        col1, col2 = st.columns(2)
        with col1:
            st.metric("Healthy", f"{healthy_prob:.1%}")
        with col2:
            st.metric("Cancer", f"{cancer_prob:.1%}", delta=f"{cancer_prob - 0.5:+.1%}")

        if cancer_prob > 0.5:
            st.error(f"**CANCER DETECTED** — {cancer_prob:.2%} confidence")
            st.balloons()
        else:
            st.success(f"**HEALTHY** — {healthy_prob:.2%} confidence")
            st.balloons()

        st.success("Prediction complete!")