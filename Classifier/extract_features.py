# extract_features_FINAL.py 
import os
import numpy as np
import pandas as pd
from tqdm import tqdm

LABEL_CSV = r"D:\MCED_DP\mced_dataset\labelled_data.xlsx"
HEALTHY_FOLDER = r"D:\MCED_DP\mced_dataset\healthy_CpG"
TUMOR_FOLDER = r"D:\MCED_DP\mced_dataset\tumor_CpG"
OUT_FEATURES = r"D:\MCED_DP\features_dataset_FINAL.csv"

CHRS = [f"chr{i}" for i in range(1,23)] + ["chrX","chrY"]

def read_labels(path):
    df = pd.read_excel(path) if path.lower().endswith(".xlsx") else pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    return df

def parse_bedgraph(filepath):
    with open(filepath, "r", errors="ignore") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("track", "#")):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except:
                continue

            # Try to get coverage if present (6-column Bismark format)
            if len(parts) >= 6:
                try:
                    cov = int(parts[4])
                    meth_reads = int(parts[5])
                    if cov == 0: continue
                    methyl = meth_reads / cov
                    weight = cov
                except:
                    methyl = float(parts[3])
                    if methyl > 1: methyl /= 100
                    weight = 1
            else:
                # 4-column format → no coverage info
                methyl = float(parts[3])
                if methyl > 1: methyl /= 100
                weight = 1

            yield chrom, methyl, weight

def compute_features(filepath):
    chrom_vals = {c: [] for c in CHRS}
    chrom_weights = {c: [] for c in CHRS}
    global_vals = []
    global_weights = []

    for chrom, methyl, weight in parse_bedgraph(filepath):
        global_vals.append(methyl)
        global_weights.append(weight)
        if chrom in chrom_vals:
            chrom_vals[chrom].append(methyl)
            chrom_weights[chrom].append(weight)

    def wmean(arr, w):
        if not arr: return np.nan
        a, w = np.array(arr), np.array(w)
        return np.average(a, weights=w)

    def wcount(w):
        return int(sum(w)) if w else 0

    feats = {}
    # Global
    feats["n_CpG"] = wcount(global_weights)
    feats["global_mean"] = wmean(global_vals, global_weights)
    feats["global_median"] = np.median(global_vals) if global_vals else np.nan
    feats["global_std"] = np.sqrt(np.average((np.array(global_vals) - feats["global_mean"])**2, weights=global_weights)) if global_vals else np.nan

    if feats["n_CpG"] > 0:
        arr = np.array(global_vals)
        w = np.array(global_weights)
        feats["pct_hyper_0.8"] = np.sum((arr >= 0.8) * w) / w.sum()
        feats["pct_hypo_0.2"] = np.sum((arr <= 0.2) * w) / w.sum()

    # Per chromosome
    for c in CHRS:
        prefix = c.replace("chr", "")
        feats[f"chr{prefix}_count"] = wcount(chrom_weights[c])
        feats[f"chr{prefix}_mean"] = wmean(chrom_vals[c], chrom_weights[c])

    return feats

def main():
    df_labels = read_labels(LABEL_CSV)
    df_labels['filename'] = df_labels['filename'].astype(str).str.strip()

    rows = []
    for _, row in tqdm(df_labels.iterrows(), total=len(df_labels)):
        fname = row['filename']
        paths = [
            os.path.join(HEALTHY_FOLDER, fname),
            os.path.join(TUMOR_FOLDER, fname),
            os.path.join(HEALTHY_FOLDER, fname.replace(".bedGraph", ".bedgraph")),
            os.path.join(TUMOR_FOLDER, fname.replace(".bedGraph", ".bedgraph")),
        ]
        fp = next((p for p in paths if os.path.exists(p)), None)
        if not fp:
            print(f"Missing: {fname}")
            continue

        feats = compute_features(fp)
        feats['filename'] = fname
        feats['label'] = row.get('label', '')
        feats['subtype'] = row.get('subtype', '')
        rows.append(feats)

    pd.DataFrame(rows).to_csv(OUT_FEATURES, index=False)
    print(f"\nDONE! Real features saved to:\n{OUT_FEATURES}")

if __name__ == "__main__":
    main()
