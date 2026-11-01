# Multi-Cancer Early Detection (MCED) using cfDNA Methylation  
A Machine Learning Pipeline for Early-Stage Tumor Detection from Blood

This repository contains my research project where I build a complete end-to-end workflow for detecting cancer from **cell-free DNA (cfDNA)** methylation patterns. It includes data preprocessing, feature extraction, model training, evaluation, and a prediction pipeline that works on new patient samples.  
This project is still in progress, and I am improving the dataset size and accuracy over time.
<p align="center"> <img src="https://img.shields.io/badge/Python-3.10%2B-blue"> <img src="https://img.shields.io/badge/ML-XGBoost-green"> <img src="https://img.shields.io/badge/Data-Methylation%20CpG-orange"> <img src="https://img.shields.io/badge/Status-Research%20Prototype-brightgreen"> </p>

## Background

DNA methylation is one of the strongest biomarkers for cancer. Tumor cells release fragmented DNA into the bloodstream (cfDNA), and their methylation pattern becomes different compared to healthy cells.Detecting such modifications early can provide inexpensive cancer screening, even before symptoms or imaging results.  However, cfDNA signals are loud, sparse, and extremely changeable, making reliable early detection problematic.By analyzing genome-wide methylation from cfDNA, it is possible to detect cancer early — even before imaging.

<p align="center">
  <img src="cfDNA_Methylation.png" width="420">
</p>

<p align="center"><em>Figure 1: Overview of circulating cell-free DNA methylation and its potential clinical application[1]</em></p>

Tumors typically exhibit global hypomethylation and focal promoter hypermethylation, both of which can be captured through cfDNA sequencing. Several major studies, including the GRAIL/CCGA investigations, demonstrate that targeted methylation patterns in cfDNA can identify more than 50 cancer types with high specificity and can also predict the tissue of origin. However, multiple researchers highlight technical challenges such as coverage heterogeneity, fragment-size bias, and dataset batch effects, making cfDNA methylation a powerful but complex domain.

DNA methylation has been repeatedly validated as one of the most reliable and stable biomarkers for detecting cancer. Tumors undergo widespread epigenetic disruption, including decreased methylation across large genomic regions and increased methylation at promoter sites, and these aberrant signatures remain detectable even after DNA fragments enter the bloodstream (Vincenza et al., 2021). Because of this stability, cfDNA methylation enables non-invasive early detection, often before imaging or clinical symptoms appear.[1,2]

Recent systematic reviews on population-scale screening (Ros et al., 2024) support the idea that cfDNA methylation can reveal dozens of cancer types with excellent specificity. Furthermore, advanced machine-learning approaches—such as those described by Janet et al. (2024)—have shown that integrating sparse genome-wide CpG signals enables accurate detection of early-stage cancers across multiple tissues. Policy-oriented studies (Mussab et al., 2024) also suggest that MCED tests can complement existing screening workflows and significantly improve early diagnostic pathways.[3,4]

Technological progress in cfDNA sequencing has accelerated the field even further. THUNDER (Gao et al., 2023) demonstrated that low-input cfDNA methylation sequencing can robustly distinguish tumor-derived fragments from healthy DNA. TOTEM (Xiong et al., 2024) introduced a powerful approach capable of predicting the tissue of origin using methylation markers, even when coverage is extremely sparse. Real-world evidence from over 100,000 clinical tests (Matrana et al., 2025) shows that cfDNA methylation-based MCED assays can be deployed safely and reproducibly across diverse populations, offering consistent diagnostic utility.[5,6,7]

These research findings show that genome-wide cfDNA methylation—when combined with well-engineered computational pipelines—forms a promising yet technically demanding foundation for multi-cancer early detection. This scientific landscape directly motivates the design and development of the MCED model implemented in this project.
The findings indicate that, when paired with reliable computational pipelines, genome-wide cfDNA methylation offers a promising but technically challenging basis for early cancer diagnosis, which directly inspired the model created in this effort.

<p align="center">
  <img src="MCED_PIPELINES.png" width="250">
</p>

<p align="center"><em>Figure 2: Working Mechanism of MCED</em></p>


Cancer changes DNA methylation across the genome, creating global hypomethylation and localized hypermethylation. Tumor-derived DNA fragments (cfDNA) circulate in blood and keep these aberrant methylation patterns.  Using bisulfite sequencing, methylation at millions of CpG sites may be measured and summarized into global and chromosome-level characteristics.  These parameters capture crucial cancer-specific methylation changes and allow an ML classifier to discriminate tumor cfDNA from healthy cfDNA with high accuracy.  

This project implements the pipeline: cfDNA methylation extraction → feature computation → XGBoost classification for early cancer diagnosis.

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
│   ├── extract_features.py
│   ├── train_xgb.py
│   ├── visualize_progress.py
│   
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

# MCED-Liquid-Biopsy-Intelligence

Exploratory computational research on cell-free DNA methylation for multi-cancer early detection (MCED).

> [!NOTE]
> This is an independent research project. It is exploratory and has not been clinically validated.

## Contents

- [Dataset preparation and preprocessing](#dataset-preparation-and-preprocessing)
  - [Tumor dataset](#1-tumor-dataset)
  - [Healthy dataset](#2-healthy-dataset)
  - [Standardized data structure](#3-standardized-data-structure)
  - [Output organization](#4-output-organization)

## Dataset preparation and preprocessing

The project uses separate tumor and healthy DNA methylation datasets. The original data are kept locally and are not in this repository because of their size.

### 1. Tumor dataset

| | |
|---|---|
| Source | NCBI Gene Expression Omnibus (GEO) |
| Accession | [GSE212391](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212391) |
| Study | Towards the DNA methylation haplotype map of 11 common solid cancers |
| Data type | Whole-genome bisulfite sequencing (WGBS) |
| Samples | 110 fresh-frozen tumor tissue samples, 11 solid cancer types |

The downloaded files are converted from bedGraph into a standard four-column format.

<details>
<summary><b>Tumor preprocessing steps</b></summary>

1. Read input files with a `.bedGraph` or `.bedGraph.gz` extension.
2. Skip metadata and comment lines that begin with `track` or `#`.
3. Drop records with fewer than six columns.
4. Take the chromosome, genomic start and genomic end coordinates.
5. Read coverage from the fifth column and the methylated-read count from the sixth.
6. Calculate the methylation fraction:

   ```text
   Methylation Fraction = Methylated Reads / Coverage
   ```

7. Set the fraction to zero when coverage is zero.
8. Write the records to normalized four-column bedGraph files, rounding the fraction to six decimal places.

</details>

### 2. Healthy dataset

The healthy methylation BED files are stored separately in the `healthy_CpG` directory. A Python script using pandas does the following:

1. Reads each `.bed` file as tab-separated with no header.
2. Keeps the first four columns (columns 0 to 3).
3. Leaves the values as they are, with no numerical transformation.
4. Writes the result as `.bedgraph` files in the `healthy_bedgraph` directory.

This step standardizes the file structure. It does not calculate or normalize methylation fractions.

> [!WARNING]
> The accession of the original healthy dataset and the meaning of its fourth column still need to be checked against the source.

### 3. Standardized data structure

Both workflows produce a four-column representation:

| Column | Description |
|---|---|
| Chromosome | Genomic chromosome identifier |
| Start | Genomic start coordinate |
| End | Genomic end coordinate |
| Value | Tumor: calculated methylation fraction. Healthy: original fourth-column value |

> [!IMPORTANT]
> The fourth column has a different origin in each workflow. Before combining the two datasets, confirm what the healthy value represents and whether it is comparable to the tumor methylation fraction.

### 4. Output organization

```text
mced_dataset/
├── healthy_CpG/         # Original healthy BED files
├── healthy_bedgraph/    # Standardized healthy bedGraph files
├── tumor_CpG/           # Original tumor bedGraph files
└── tumor_processed/     # Processed tumor bedGraph files
```

These files are the input for the later feature extraction and machine-learning experiments.
## Features Extracted

The structure utilized in contemporary MCED (Multi-Cancer Early Detection) systems is adopted in this project:

-	 **Epigenetic divergence**: Systematic hyper/hypomethylation causes tumor methylomes to diverge from normal cells.

-	 **cfDNA shedding model**: Sensitive pattern extraction is necessary because tumor-derived fragments mingle with vast amounts of healthy DNA.

-	 **Coverage-aware feature engineering**: Low-coverage  Weighting by coverage and filtering enhances signal quality since CpG sites add noise.

-	 **Machine learning classification**: Extracted methylation statistics (global, chromosome-wise, hyper/hypo ratios) operate as features for an ML classifier such as XGBoost.
  
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


Mu represents the coverage-weighted global mean methylation across all CpG sites
### Coverage-Weighted Standard Deviation
$$
\sigma = \sqrt{
\frac{\sum_{i=1}^{N} c_i (m_i - \mu)^2}
     {\sum_{i=1}^{N} c_i}
}
$$


Sigma measures the coverage-weighted variability of methylation levels across the genome.
### Hyper-methylated Fraction (m ≥ 0.8)

$$
Pct_{hyper} =
\frac{\sum_{i=1}^{N} c_i \cdot \mathbf{1}(m_i \ge 0.8)}
     {\sum_{i=1}^{N} c_i}
$$


The hyper-methylated fraction is the proportion of CpG coverage coming from sites with methylation levels of at least eighty percent.
### Hypo-methylated Fraction (m ≤ 0.2)

$$
Pct_{hypo} =
\frac{\sum_{i=1}^{N} c_i \cdot \mathbf{1}(m_i \le 0.2)}
     {\sum_{i=1}^{N} c_i}
$$


The hypomethylated percentage is calculated as the total coverage of CpG sites with methylation less than or equal to twenty percent, divided by the total CpG coverage.
### Chromosome-wise Mean Methylation(Where 𝐶𝑘 is the set of CpGs on chromosome 𝑘.)

$$
\mu_k =
\frac{\sum_{i \in C_k} c_i m_i}
     {\sum_{i \in C_k} c_i}
$$


Mu k represents the coverage-weighted mean methylation level for chromosome k.
### Feature Vector for Methylation
$$
F =
[
\mu,\ \sigma,\ Pct_{hyper},\ Pct_{hypo},\
\{\mu_k\}_{k=1}^{24},\ \{Count_k\}_{k=1}^{24},\ CpG_{eff}
]
$$


The feature vector includes global mean and variance of methylation, hyper- and hypo-methylation percentages, chromosome-wise mean methylation and CpG counts across all twenty-four chromosomes, and the effective number of CpG sites.

### Updated: Fragmentomics Integration 

We enhanced the feature set with **fragmentomics** to capture cancer-specific cfDNA fragmentation patterns:
- Cancer cfDNA tends to have **shorter** (<120 bp) and more **irregular** fragments due to abnormal apoptosis and open chromatin.
- These are calculated directly from **inter-CpG distances** in bedGraph files (no additional data needed).

This brings the pipeline in line with 2025 state-of-the-art approaches (e.g., GRAIL Galleri, DELFI, SPOT-MAS).

### Feature Vector (F)

The model uses the following comprehensive feature vector **F** for each sample:

$$
F = \left[ 
\mu, \sigma, P_{ct_{hyper}}, P_{ct_{hypo}}, 
\{\mu_k, C_k\}_{k=1}^{24}, 
C_pG_{eff},
\mu_f, m_f, P_{short}, \sigma_f, CV_f, H_f 
\right]
$$

Where:
 **Fragmentomics** (captures cancer-specific cfDNA fragmentation):
  - $\mu_f$: Mean inter-CpG fragment length
  - $m_f$: Median inter-CpG fragment length
  - $P_{short}$: Percentage of fragments < 120 bp (strong cancer signal)
  - $\sigma_f$: Standard deviation of fragment lengths
  - $CV_f$: Coefficient of variation (std / mean)
  - $H_f$: Shannon entropy of fragment length distribution (measures irregularity)
- **Chromosome-wise** (k = 1 to 22, X, Y — 24 chromosomes):
  - $\mu_k$: Mean methylation on chromosome k
  - $C_k$: Weighted CpG count on chromosome k
- **Effective coverage**:
  - $C_pG_{eff}$: Total effective (coverage-weighted) CpG count
  

These features are computed in a coverage-weighted manner to reduce noise in low-depth cfDNA samples and are stored in `features_with_fragmentomics.csv` and used for model training.

## Model Training

I train an **XGBoost classifier** using 5-fold cross-validation.

**What the model learns:**

Patterns in methylation, coverage, and genomic features that may distinguish tumor-associated cfDNA signals from healthy cfDNA.

**Training Output Example (my latest run):**

- **Train samples:** 16 (20 GB raw dataset)
- **5-fold CV AUC:** **1.0000 ± 0.0000**
- **Top features:**
  - `n_CpG`
  - `global_mean`
- **Model:** XGBoost
- **Model saved as:** `models/MCED_CANCER_DETECTOR.joblib`

> **Note:** This is an initial result from a small cohort and should not be interpreted as evidence of clinical performance or generalization. Further experiments with larger and independent cohorts are planned to evaluate the robustness of the observed signal.

### Healthy vs Tumor

<p align="left">
  <img src="healthy vs tumour.png" width="500">
</p>

### Hypomethylation Comparison

<p align="left">
  <img src="hypomethylation comparison.png" width="500">
</p>


The script will:

1. Extract features from the input methylation file
2. Load the trained model
3. Output the model's predicted cancer probability (0–1)
4. Report the prediction threshold used by the model
5. Show the top contributing features

> **Note:** The prediction pipeline is intended for research and experimentation, not clinical diagnosis.

---

## Current Status

- Built full data pipeline
- Feature extraction validated
- Binary cancer classification model trained
- Obtained an initial **5-fold CV AUC of 1.0000** on a 16-sample cohort
- Identified methylation and coverage-related features for further investigation
- Pipeline is ready to scale to larger and more diverse datasets

---

## Next Steps (Planned)

- Increase dataset from ~20 GB → **60–80 GB** using multi-cancer data from GEO + EGA
- Perform **independent and cross-cohort validation**
- Investigate model robustness across different sequencing coverage and dataset characteristics
- Perform **feature-ablation experiments** to determine the contribution of different feature groups
- Investigate whether `n_CpG` and coverage-related features reflect biological signal or technical/sample-quality effects
- Add a **tissue-of-origin classifier** after establishing robust binary cancer-vs-control performance
- Perform deeper biological analysis of DMRs
- Explore interpretable machine-learning methods to connect model predictions with biologically meaningful methylation patterns

---

## Why I Built This

I have a strong interest in **Machine Learning for healthcare and computational biology**. My goal is to develop machine-learning methods that can contribute to early cancer detection using molecular data such as cfDNA.

This project is my attempt to build a complete research pipeline independently — from processing raw sequencing-derived data and designing molecular features to training a machine-learning model and investigating which features contribute to its predictions.

It also helped me explore:

- How methylation patterns in cfDNA may differ between tumor-associated and healthy samples
- How cfDNA can carry molecular signals associated with tumors
- How sequencing coverage and genomic features influence sample-level representations
- How machine learning can identify combinations of molecular features associated with cancer
- How computational models can be used to investigate biological questions

There is still much for me to learn about the biological foundations of cancer genomics, and I consider this project an ongoing research effort.

Rather than treating the initial results as a final conclusion, I am using them to identify better experiments, improve the methodology, and test whether the observed signals remain consistent across larger and more heterogeneous datasets.
## Contact

If you want to discuss or collaborate feel free to reach out.  
**Nabin Pathak**  
Computer Engineering Graduate  
Machine Learning Research Enthusiast  
Kathmandu, Nepal  
**Gmail:** nabinpathak520@gmail.com

### REFERENCES

1. Vincenza C, Judy H, Yasutaka Y, Sheng-Y, Himisha B. *Epigenetics in prostate cancer: clinical implications*. doi: 10.21037/tau-20-1339, 2021.

2. Ros W, Sarah N, Melissa H, Yiwen L, Claire K, Gary R, Rachel C, Sofia D.  
   *Multi-cancer early detection tests for general population screening: a systematic literature review.*  
   https://doi.org/10.1101/2024.02.14.24302576 (2024)

3. Janet V, David G, Alex G, Eric A, Jordan J.  
   *A multi-cancer early detection blood test using machine learning detects early-stage cancers.*  
   https://doi.org/10.1038/s41698-024-00568-z (2024)

4. Mussab F, Hadi A, Stephen Q, Özge K, Jon E, Maarten J.  
   *Integrating Multi-Cancer Early Detection (MCED) Tests with Standard Cancer Screening.*  
   October 2024.

5. Gao, Q., Lin, Y. P., Li, B. S., Li, G. Q., Dong, L. Q., Shen, B. Y., Lou, W. H., Wu, W. C., Ge, D., Zhu, Q. L...  
   *Unintrusive multi-cancer detection by circulating cell-free DNA methylation sequencing (THUNDER).*  
   February 2023.

6. Xiong, D., Han, T., Li, Y., Hong, Y., Li, S., Li, X., Tao, W., Huang, Y. S., Chen, W., Li, C.  
   *TOTEM: A multi-cancer detection and localization approach using circulating tumor DNA methylation markers.* (2024)

7. Matrana, M., Shukla, V., Kingsbury, D., Poliak, M., Lipton, J., McMillin, M., Malinow, L. B., Venn, O.,  
   Beausang, J. F., Stanley, G., Hubbell, E., Kurtzman, K. N., Venstrom, J. M., Shaknovich, R., Westgate, C.  
   *Real-world data and clinical experience from over 100,000 multi-cancer early detection tests.* (2025)



  



