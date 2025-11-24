# Multi-Cancer Early Detection (MCED) using cfDNA Methylation  
A Machine Learning Pipeline for Early-Stage Tumor Detection from Blood

This repository contains my research project where I build a complete end-to-end workflow for detecting cancer from **cell-free DNA (cfDNA)** methylation patterns. It includes data preprocessing, feature extraction, model training, evaluation, and a prediction pipeline that works on new patient samples.  
This project is still in progress, and I am improving the dataset size and accuracy over time.

## Background

DNA methylation is one of the strongest biomarkers for cancer. Tumor cells release fragmented DNA into the bloodstream (cfDNA), and their methylation pattern becomes different compared to healthy cells.  
By analyzing genome-wide methylation from cfDNA, it is possible to detect cancer early — even before imaging.

This project uses **bedGraph methylation profiles** (4-column format: `chr, start, end, methyl_fraction`) extracted from tumor and healthy samples.

## Project Structure

MCED-Liquid-Biopsy-Intelligence/
│
├── raw_bedgraph/ # raw bedGraph files (tumor + healthy)
├── processed_bedgraph/ # standardized 4-column files
├── features/ # extracted features for each sample
│ └── features_dataset.csv
├── scripts/ # feature extractor, training, inference
├── models/ # trained ML models (*.joblib)
└── README.md # this file

## Features Extracted

Each sample is converted to a fixed-length feature vector including:

- global mean methylation  
- global standard deviation  
- percent hyper-methylated CpGs  
- percent hypo-methylated CpGs  
- chromosome-wise methylation  
- CpG count (coverage weighted)

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
Nabin Pathak
Computer Engineering Graduate
Machine Learning Research Enthusiast
Kathmandu, Nepal
Gmail: nabinpathak520@gmail.com

## Acknowledgements

This project uses public methylation datasets from GEO and other open sources.  
Thanks to the researchers who made their cfDNA data accessible for scientific progress.


