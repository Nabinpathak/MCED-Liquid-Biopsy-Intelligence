<h1 align="center">MCED-Liquid-Biopsy-Intelligence</h1>

<p align="center">
  Exploratory computational research on cell-free DNA methylation for multi-cancer early detection
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Python-3.10%2B-blue">
  <img src="https://img.shields.io/badge/ML-XGBoost-green">
  <img src="https://img.shields.io/badge/Data-Methylation%20CpG-orange">
  <img src="https://img.shields.io/badge/Status-Research%20Prototype-brightgreen">
  <img src="https://img.shields.io/badge/Clinical%20validation-None-red">
</p>

> [!IMPORTANT]
> This is an independent research prototype. It has **not** been clinically validated, and its results come from 16 samples. Nothing here should be used for diagnosis.

## Contents

- [Abstract](#abstract)
- [At a glance](#at-a-glance)
- [Background](#background)
- [Research questions](#research-questions)
- [Pipeline](#pipeline)
- [Repository structure](#repository-structure)
- [Data](#data)
- [Features](#features)
- [Current results](#current-results)
- [Planned validation protocol](#planned-validation-protocol)
- [Known limitations](#known-limitations)
- [Roadmap](#roadmap)
- [Reproducibility](#reproducibility)
- [Related work](#related-work)
- [Citation](#citation)
- [Author and contact](#author-and-contact)

## Abstract

Methylation patterns in cell-free DNA (cfDNA) can carry a signal from tumors, and published tests reach high specificity. Their sensitivity for early-stage cancer is much lower, so the open problem is how to represent a weak tumor signal so that it can be learned and tested honestly. This repository is my attempt to build a full pipeline for that problem on my own: from raw bedGraph methylation profiles to sample-level features and an XGBoost classifier. The current model was trained on 16 samples and reached a cross-validated AUC of 1.0000, which I treat as a warning sign and not as a result, because the top features depend on coverage and the tumor and healthy samples come from different sources. The next stage of the project is to test whether the signal survives strict validation.

## At a glance

| | |
|---|---|
| Goal | Learn and test sample-level representations of cfDNA methylation |
| Input | Four-column bedGraph files: `chr, start, end, methyl_fraction` |
| Features | Coverage-weighted global, chromosome-level and hyper/hypomethylation statistics, plus inter-CpG spacing features |
| Model | XGBoost, 5-fold cross-validation |
| Training set | 16 samples (about 20 GB of raw data) |
| Validation | Internal cross-validation only. No independent cohort yet |
| Clinical claims | None |

## Background

Tumors typically show global hypomethylation together with focal hypermethylation at promoters, and part of that signal circulates in blood as cfDNA [11]. Published methylation tests report high specificity but low early-stage sensitivity:

| Study | Participants | Specificity | Stage I sensitivity |
|---|---|---|---|
| Liu et al., 2020 [1] | 6,689 | 99.3% | 18% across all cancer types |
| Klein et al., 2021 [2] | 4,077 | 99.5% | 16.8% |

A systematic review of MCED tests for population screening found specificity above 96% for every test, lower and more variable sensitivity for early-stage cancers, and no completed randomized trials [3]. In the prospective SPOT-MAS study of 9,024 adults without symptoms, specificity was 99.71% but the positive predictive value was 39.53%, and the sensitivity estimate came from only 17 true positives [4]. A test can therefore have excellent specificity and still send many healthy people to follow-up.

<p align="center">
  <img src="cfDNA_Methylation.png" width="420">
</p>

<p align="center"><em>Figure 1: Overview of circulating cell-free DNA methylation and its potential clinical application [11]</em></p>

## Research questions

The project is organized around three questions that the literature leaves open. I have not answered any of them yet.

1. Does combining feature types help? In the 2023 SPOT-MAS study, the best single feature (end motifs) reached an AUC of 0.90 in discovery cross-validation, and a stacked ensemble reached 0.95 in validation, but one model given all features at once scored 0.88 [5]. Whether extra signals help depends on how they are combined, and I did not find a study that changes only that choice on identical samples and folds.
2. Do reported estimates survive strict checks? These are patient-level splits, label shuffling and trivial baselines such as age or sequencing depth alone.
3. Is stage the right axis? Jamshidi et al. found that tumor allele fraction predicted classifier performance better than clinical stage or tumor type [6].

## Pipeline

<p align="center">
  <img src="MCED_PIPELINES.png" width="250">
</p>

<p align="center"><em>Figure 2: Working mechanism of MCED</em></p>

```text
bedGraph methylation profiles
        │
        ▼
standardization (4 columns)
        │
        ▼
feature extraction (coverage-weighted)
        │
        ▼
XGBoost classifier (5-fold CV)
        │
        ▼
predicted cancer probability
```

## Repository structure

```text
MCED-Liquid-Biopsy-Intelligence/
│
├── raw_bedgraph/               # Raw bedGraph files (healthy + tumor)
├── processed_bedgraph/         # Standardized 4-column bedGraph files
├── features/
│   └── features_dataset.csv    # Final dataset used for training
├── scripts/
│   ├── extract_features.py
│   ├── train_xgb.py
│   └── visualize_progress.py
├── models/
│   └── MCED_CANCER_DETECTOR.joblib
├── requirements.txt
└── README.md
```

## Data

The data are large (20 to 100 GB) and are not in this repository. They are stored locally. The project uses separate tumor and healthy datasets.

<details>
<summary><b>Raw record format</b></summary>

Each line of a raw methylation file describes one CpG site:

- `chr1`: chromosome
- `10468–10469`: genomic interval (1 bp CpG site)
- `100 / 83`: methylation percentage or fraction
- `6 / 5`: number of reads covering the CpG
- `0 / 1`: methylated vs unmethylated read count

In the formulas below, a CpG site is written as:

```text
CpG position = chr1_10468
methylation  = m_i
coverage     = c_i
```

</details>

### Tumor dataset

| | |
|---|---|
| Source | NCBI Gene Expression Omnibus (GEO) |
| Accession | [GSE212391](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212391) |
| Study | Towards the DNA methylation haplotype map of 11 common solid cancers |
| Data type | Whole-genome bisulfite sequencing (WGBS) |
| Samples | 110 fresh-frozen tumor tissue samples, 11 solid cancer types |

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

### Healthy dataset

The healthy methylation BED files are stored separately in the `healthy_CpG` directory. A Python script using pandas reads each `.bed` file as tab-separated with no header, keeps the first four columns (0 to 3) unchanged, and writes `.bedgraph` files to `healthy_bedgraph`. It does not calculate or normalize methylation fractions.

> [!WARNING]
> The accession of the original healthy dataset and the meaning of its fourth column still need to be checked against the source.

### Standardized structure

| Column | Description |
|---|---|
| Chromosome | Genomic chromosome identifier |
| Start | Genomic start coordinate |
| End | Genomic end coordinate |
| Value | Tumor: calculated methylation fraction. Healthy: original fourth-column value |

> [!IMPORTANT]
> The fourth column has a different origin in each workflow. Before combining the two datasets, confirm what the healthy value represents and whether it is comparable to the tumor methylation fraction.

<details>
<summary><b>Local data layout</b></summary>

```text
mced_dataset/
├── healthy_CpG/         # Original healthy BED files
├── healthy_bedgraph/    # Standardized healthy bedGraph files
├── tumor_CpG/           # Original tumor bedGraph files
└── tumor_processed/     # Processed tumor bedGraph files
```

</details>

## Features

Each sample is converted to a fixed-length vector. Let CpG sites be indexed by i = 1, 2, …, N, with methylation fraction `m_i` (0 to 1), coverage `c_i` and chromosome `chr_i`. Only sites with c_i ≥ 10 are used. All statistics are coverage-weighted to reduce noise in low-depth samples.

| Feature | Definition |
|---|---|
| Global mean methylation (μ) | Coverage-weighted mean across all CpG sites |
| Standard deviation (σ) | Coverage-weighted variability of methylation |
| Hypermethylated fraction | Share of total coverage from sites with m ≥ 0.8 |
| Hypomethylated fraction | Share of total coverage from sites with m ≤ 0.2 |
| Chromosome-wise mean (μ_k) | Coverage-weighted mean on chromosome k, 24 chromosomes |
| Chromosome-wise count | Coverage-weighted CpG count on chromosome k |
| Effective CpG count | Total coverage-weighted number of CpG sites |
| Inter-CpG spacing | Mean, median, share below 120 bp, standard deviation, coefficient of variation and entropy of distances between consecutive CpG observations |

<details>
<summary><b>Formulas</b></summary>

$$
\mu = \frac{\sum_{i=1}^{N} c_i m_i}{\sum_{i=1}^{N} c_i}
\qquad
\sigma = \sqrt{\frac{\sum_{i=1}^{N} c_i (m_i - \mu)^2}{\sum_{i=1}^{N} c_i}}
$$

$$
Pct_{hyper} = \frac{\sum_{i=1}^{N} c_i \cdot \mathbf{1}(m_i \ge 0.8)}{\sum_{i=1}^{N} c_i}
\qquad
Pct_{hypo} = \frac{\sum_{i=1}^{N} c_i \cdot \mathbf{1}(m_i \le 0.2)}{\sum_{i=1}^{N} c_i}
$$

$$
\mu_k = \frac{\sum_{i \in C_k} c_i m_i}{\sum_{i \in C_k} c_i}
$$

where C_k is the set of CpG sites on chromosome k.

$$
F = \left[
\mu,\ \sigma,\ Pct_{hyper},\ Pct_{hypo},\
\{\mu_k\}_{k=1}^{24},\ \{Count_k\}_{k=1}^{24},\ CpG_{eff},\
\mu_d,\ m_d,\ P_{d<120},\ \sigma_d,\ CV_d,\ H_d
\right]
$$

The last six terms are the inter-CpG spacing statistics.

</details>

> [!IMPORTANT]
> Inter-CpG distance is **not** cfDNA fragment length. It depends on which CpG sites were covered and how they are distributed across the genome, not on where a DNA fragment starts and ends. These are methylation-spacing features, not fragmentomics. True fragmentomics needs fragment-level information such as paired-end coordinates or insert sizes, fragment ends, end motifs or nucleosomal patterns, which methylation calls alone cannot supply [7, 8].

## Current results

| | |
|---|---|
| Training samples | 16 (about 20 GB raw data) |
| 5-fold CV AUC | 1.0000 ± 0.0000 |
| Top features | `n_CpG`, `global_mean` |
| Saved model | `models/MCED_CANCER_DETECTOR.joblib` |

> [!CAUTION]
> A perfect AUC on 16 samples is not evidence of clinical performance or generalization. Both top features depend on how many CpG sites were measured and on overall methylation level, and the tumor and healthy samples come from different datasets. The score may therefore reflect differences between data sources or sequencing coverage, not cancer biology. Ruling this out is the first task in the [roadmap](#roadmap).

<p align="center">
  <img src="healthy%20vs%20tumour.png" width="420">
  <img src="hypomethylation%20comparison.png" width="420">
</p>

<p align="center"><em>Left: healthy vs tumor. Right: hypomethylation comparison.</em></p>

## Planned validation protocol

The next experiment is small and its rules are fixed before any result is seen.

1. Data audit. Find a public cohort with paired-end reads, a cancer or healthy label, and one identifier per patient. If none exists, the experiment stops.
2. Splits. Repeated nested cross-validation with patient-level folds. Scaling, feature selection and imputation are fitted inside the training folds only.
3. Models on identical folds. Metadata-only baseline (age, sex, sequencing depth), methylation-only, spacing- or fragment-only, concatenated features, and stacked models.
4. Controls. At least 1,000 label permutations per model, and a model that tries to predict the data source or batch from the features.
5. Metrics. AUC and sensitivity at fixed specificity (95% and 98%) with bootstrap confidence intervals, plus calibration.

| Result | Interpretation |
|---|---|
| Observed AUC inside the label-shuffle distribution | No detectable signal. Performance is not reported. |
| Metadata-only baseline within 0.03 AUC of the best model | The model is not beating trivial baselines. The signal is treated as confounded. |
| Data source predictable from features well above chance | Features carry provenance, and any gain may be a batch effect. |
| Combined-model gain smaller than the spread across folds | Combining features is not shown to help. Reported as a negative result. |

A clean negative result is an acceptable outcome. This experiment can show whether a feature set adds information under controlled conditions. It cannot show clinical value, generalization across populations, or tissue of origin.

## Known limitations

| Limitation | Why it matters |
|---|---|
| 16 training samples | Estimates are unstable and cross-validation can be overly optimistic |
| Tumor and healthy data come from different sources | A classifier can learn the source instead of cancer |
| Healthy dataset not yet verified | The meaning of its fourth column is unconfirmed |
| Top features depend on coverage | Technical differences may drive the score |
| Inter-CpG distance is not fragment length | The spacing features are not fragmentomics |
| No independent cohort | Generalization is untested |
| No stage or tissue-of-origin analysis | Early-stage performance is unknown |

## Roadmap

- [ ] Check whether `n_CpG` and coverage-related features reflect biology or technical effects
- [ ] Verify the healthy dataset accession and the meaning of its fourth column
- [ ] Run the validation protocol above on a public paired-end cohort
- [ ] Run independent and cross-cohort validation
- [ ] Run feature-ablation experiments to find which feature groups contribute
- [ ] Test robustness across sequencing coverage and dataset characteristics
- [ ] Increase the dataset from about 20 GB to 60–80 GB using multi-cancer data from GEO and EGA
- [ ] Add a tissue-of-origin classifier once binary performance is established
- [ ] Analyze differentially methylated regions (DMRs) in more depth
- [ ] Explore interpretable methods that link predictions to methylation patterns

## Reproducibility

Requirements: Python 3.10 or later and the packages in `requirements.txt`.

```bash
git clone https://github.com/Nabinpathak/MCED-Liquid-Biopsy-Intelligence.git
cd MCED-Liquid-Biopsy-Intelligence
pip install -r requirements.txt
python scripts/extract_features.py --help
python scripts/train_xgb.py --help
```

The raw data are not included. Download the tumor data from [GSE212391](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212391) and preprocess them as described in [Data](#data).

The prediction script extracts features from an input methylation file, loads the trained model, and outputs a cancer probability between 0 and 1 with the threshold used and the top contributing features.

> [!NOTE]
> The prediction pipeline is for research and experimentation, not clinical diagnosis.

## Related work

1. Liu MC, Oxnard GR, Klein EA, et al. Sensitive and specific multi-cancer detection and localization using methylation signatures in cell-free DNA. *Ann Oncol*. 2020;31(6):745–759. doi:10.1016/j.annonc.2020.02.011
2. Klein EA, Richards D, Cohn A, et al. Clinical validation of a targeted methylation-based multi-cancer early detection test using an independent validation set. *Ann Oncol*. 2021;32(9):1167–1177.
3. Wade R, Nevitt S, Harden M, et al. Multi-cancer early detection tests for general population screening: a systematic literature review. medRxiv. 2024. doi:10.1101/2024.02.14.24302576
4. Nguyen LHD, Nguyen THH, et al. Prospective validation study: a non-invasive circulating tumor DNA-based assay for simultaneous early detection of multiple cancers in asymptomatic adults. *BMC Med*. 2025;23:90. doi:10.1186/s12916-025-03929-y
5. Nguyen VTC, et al. Multimodal analysis of methylomics and fragmentomics in plasma cell-free DNA for multi-cancer early detection and localization. *eLife*. 2023;12:89083. doi:10.7554/eLife.89083
6. Jamshidi A, et al. Evaluation of cell-free DNA approaches for multi-cancer early detection. *Cancer Cell*. 2022;40(12):1537–1549. doi:10.1016/j.ccell.2022.10.022
7. Mouliere F, Chandrananda D, Piskorz AM, et al. Enhanced detection of circulating tumor DNA by fragment size analysis. *Sci Transl Med*. 2018;10(466). doi:10.1126/scitranslmed.aat4921
8. Cristiano S, Leal A, Phallen J, et al. Genome-wide cell-free DNA fragmentation in patients with cancer. *Nature*. 2019;570:385–389. doi:10.1038/s41586-019-1272-6
9. Gao Q, Lin YP, et al. Unintrusive multi-cancer detection by circulating cell-free DNA methylation sequencing (THUNDER): development and independent validation studies. *Ann Oncol*. 2023;34(5):486–495.
10. Esfahani MS, Hamilton EG, Mehrmohamadi M, et al. Inferring gene expression from cell-free DNA fragmentation profiles. *Nat Biotechnol*. 2022;40:585–597. doi:10.1038/s41587-022-01222-4
11. Conteduca V, Hess J, Yamada Y, Ku SY, Beltran H. Epigenetics in prostate cancer: clinical implications. *Transl Androl Urol*. 2021. doi:10.21037/tau-20-1339

## Citation

If you refer to this project, please cite it as:

```bibtex
@misc{pathak_mced_liquid_biopsy,
  author = {Pathak, Nabin},
  title  = {MCED-Liquid-Biopsy-Intelligence},
  year   = {2026},
  note   = {Independent research project},
  url    = {https://github.com/Nabinpathak/MCED-Liquid-Biopsy-Intelligence}
}
```

## Author and contact

**Nabin Pathak**, Kathmandu, Nepal
Machine learning research, computational biology and cancer genomics
Email: nabinpathak520@gmail.com

I built this project to learn how to design and test molecular features for early cancer detection. There is still much for me to learn about the biology of cancer genomics, and I treat the initial results as a starting point for better experiments. To discuss the project or collaborate, write to me at the email above.
