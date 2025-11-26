# Multi-Cancer Early Detection (MCED) using cfDNA Methylation  
A Machine Learning Pipeline for Early-Stage Tumor Detection from Blood

This repository contains my research project where I build a complete end-to-end workflow for detecting cancer from **cell-free DNA (cfDNA)** methylation patterns. It includes data preprocessing, feature extraction, model training, evaluation, and a prediction pipeline that works on new patient samples.  
This project is still in progress, and I am improving the dataset size and accuracy over time.
<p align="center"> <img src="https://img.shields.io/badge/Python-3.10%2B-blue"> <img src="https://img.shields.io/badge/ML-XGBoost-green"> <img src="https://img.shields.io/badge/Data-Methylation%20CpG-orange"> <img src="https://img.shields.io/badge/Status-Research%20Prototype-brightgreen"> </p>

## Background

DNA methylation is one of the strongest biomarkers for cancer. Tumor cells release fragmented DNA into the bloodstream (cfDNA), and their methylation pattern becomes different compared to healthy cells.  
By analyzing genome-wide methylation from cfDNA, it is possible to detect cancer early — even before imaging.

This project uses **bedGraph methylation profiles** (4-column format: `chr, start, end, methyl_fraction`) extracted from tumor and healthy samples.

##  Project Directory Structure

```
MCED-Liquid-Biopsy-Intelligence/
│
├── raw_bedgraph/               # Raw bedGraph files (healthy + tumor)
│   ├── *.bedGraph
│   └── *.bedgraph
│
├── processed_bedgraph/         # Standardized 4-column cleaned bedGraph files
│   ├── tumor_*.bedGraph
│   └── healthy_*.bedGraph
│
├── features/                   # Extracted feature vectors
│   └── features_dataset.csv    # Final dataset used for training
│
├── scripts/                    # All core ML and processing scripts
│   ├── extract_features_corrected.py
│   ├── train_xgb.py
│   ├── inference.py
│   └── visualize_progress.py
│
├── models/                     # Trained ML models
│   └── MCED_CANCER_DETECTOR.joblib
│
├── requirements.txt            # Python dependencies
└── README.md                   # This documentation
```



## Dataset Information

The dataset used for this project is large (20–100 GB) and cannot be uploaded directly to GitHub. For now, the dataset is stored locally.
Where:  
- **chr1** — chromosome  
- **10468–10469** — genomic interval (1 bp CpG site)  
- **100 / 83** — methylation percentage or fraction  
- **6 / 5** — number of reads covering the CpG  
- **0 / 1** — methylated vs unmethylated read count


CpG position = chr1_10468
methylation = m_i
coverage     = c_i

You can download the raw cfDNA dataset from the official website of National Center for Biotechnology Information(NCBI).(https://www.ncbi.nlm.nih.gov/)
```
D:\MCED-Liquid-Biopsy-Intelligence\mced_dataset
├── healthy_CpG
├── tumor_CpG
└── labelled_data.xlsx
```
## Features Extracted

Each sample is converted to a fixed-length feature vector including:

Let each CpG site be indexed by i = 1, 2, …, N, with:
*	m_i — methylation fraction at site i (range 0–1)
*	c_i — coverage (number of reads) at site i
*	chr_i — chromosome of CpG i
*	Only CpG sites with c_i ≥ 10 are included

Each CpG site contains:
*	m_i = methylation value (0–1)
*	c_i = coverage (reads)

### Coverage-Weighted Global Mean Methylation
$$
\mu = \frac{\sum_{i=1}^{N} c_i m_i}{\sum_{i=1}^{N} c_i}
$$

### Coverage-Weighted Standard Deviation(Captures genome-wide methylation variability.)
$$
\sigma = \sqrt{
\frac{\sum_{i=1}^{N} c_i (m_i - \mu)^2}
     {\sum_{i=1}^{N} c_i}
}
$$

### Hyper-methylated Fraction (m ≥ 0.8)

$$
Pct_{hyper} =
\frac{\sum_{i=1}^{N} c_i \cdot \mathbf{1}(m_i \ge 0.8)}
     {\sum_{i=1}^{N} c_i}
$$

### Hypo-methylated Fraction (m ≤ 0.2)(Cancer cfDNA typically shows global hypomethylation.)

$$
Pct_{hypo} =
\frac{\sum_{i=1}^{N} c_i \cdot \mathbf{1}(m_i \le 0.2)}
     {\sum_{i=1}^{N} c_i}
$$

### Chromosome-wise Mean Methylation(Where 𝐶𝑘 is the set of CpGs on chromosome 𝑘.)

$$
\mu_k =
\frac{\sum_{i \in C_k} c_i m_i}
     {\sum_{i \in C_k} c_i}
$$

### Final Feature Vector
$$
F =
[
\mu,\ \sigma,\ Pct_{hyper},\ Pct_{hypo},\
\{\mu_k\}_{k=1}^{24},\ \{Count_k\}_{k=1}^{24},\ CpG_{eff}
]
$$

These features are stored in `features_dataset.csv` and used for model training.

## Model Training

I train an **XGBoost classifier** using 5-fold cross-validation.

**What the model learns:**  
Patterns in methylation that separate tumor DNA from healthy DNA.

**Training Output Example (my latest run):**

- Train samples: 16 (20gb raw dataset) 
- 5-fold CV AUC: **1.0000 ± 0.0000**  
- Sensitivity @ 99% specificity: **100%**  
- Top features:
  - `n_CpG`
  - `global_mean`
Model saved as: models/MCED_CANCER_DETECTOR.joblib

[Healthy vs Tumor]
<p align="center">
  <img src="images/healthy_vs_tumour.png" width="600">
</p>

[Hypomethylation Comparison]
<p align="center">
  <img src="images/hypomethylation_comparison.png" width="600">
</p>



To run prediction on a new sample: python predict_sample.py --input path/to/sample.bedgraph

The script will:

1. extract features from the uploaded methylation file  
2. load the trained model  
3. output cancer probability (0–1)  
4. give recommended clinical threshold  
5. show top contributing features  

## Current Status

-  Built full data pipeline  
-  Feature extraction validated  
-  Binary cancer detection model trained  
-  Achieved perfect AUC on initial dataset  
-  Ready for multi-class expansion  

## Next Steps (Planned)

- Increase dataset from ~20 GB → **60–80 GB** (multi-cancer from GEO + EGA)  
- Add **tissue-of-origin classifier** (multiclass XGBoost)  
- Implement a **Streamlit web app** for real sample uploads  
- Apply for **Research Assistant** positions at Harvard / MIT  
- Perform deeper biological analysis on DMRs  

## Why I Built This

I have a strong interest in Machine Learning for healthcare. My goal is to work in a research group and contribute to developing early detection methods for cancer using cfDNA.  
This project is my attempt to show I can build a complete pipeline independently — from raw sequencing data to an actual working cancer detector.

It also helped me understand:

- how methylation behaves in solid tumors  
- how cfDNA leaks tumor signals into blood  
- how machine learning can pick up these subtle methylation patterns  

There may be some mistakes in the code as I am still learning, but I keep improving things as dataset grows.

## Contact

If you want to discuss or collaborate feel free to reach out.  
**Nabin Pathak**  
Computer Engineering Graduate  
Machine Learning Research Enthusiast  
Kathmandu, Nepal  
**Gmail:** nabinpathak520@gmail.com

## Acknowledgements

This project uses public methylation datasets from GEO and other open sources.  



