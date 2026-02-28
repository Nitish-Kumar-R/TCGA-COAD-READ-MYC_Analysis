# TCGA-COAD-READ-MYC_Analysis
# TCGA COAD/READ – MYC-Centered Integrative Analysis

An end-to-end reproducible pipeline for MYC-driven colorectal cancer (CRC) transcriptomics using TCGA COAD + READ RNA-seq data (STAR counts), including differential expression, WGCNA, survival analysis, GSEA, ceRNA network reconstruction, and XGBoost-based mortality prediction with SHAP interpretability.

---

## Repository Structure

```
TCGA-COAD-READ-MYC-Analysis/
│
├── TCGA_COAD_READ_MYC_Analysis.R     # Main R pipeline (Sections 00–18)
├── ml10_ME4_XGBoost_SHAP.py          # ME4 XGBoost 3-year mortality + SHAP
│
├── figures/                          # 15 publication-ready TIFF figures (600 DPI)
├── tables/                           # 12 CSV result tables
├── results/
│   ├── ml/                           # ML model outputs (ROC, SHAP TIFFs)
│   └── ceRNA/                        # Cytoscape-ready edge/node tables
│
└── data/
    └── external/                     # Place miRTarBase & ENCORI files here (see §18)
```

---

## Pipeline Overview

| Section | Description | Key Output |
|---------|-------------|------------|
| 00 | Setup, packages, helpers | — |
| 01 | TCGA COAD+READ download via TCGAbiolinks | `counts`, `coldata` |
| 02 | Sample annotation & QC | Tumor / Normal split |
| 03 | Gene filtering (≥10 counts in ≥20% samples) | Filtered count matrix |
| 04 | DESeq2 – Tumor vs Normal | `tables/01_DGE_Tumor_vs_Normal.csv` |
| 05 | PCA | `figures/Fig_01_PCA_Tumor_vs_Normal.tiff` |
| 06 | Heatmap – Top 50 DEGs (Tumor vs Normal) | `figures/Fig_02_Heatmap_Tumor_vs_Normal.tiff` |
| 07 | MYC stratification (median split) | `coldata$MYC_group` |
| 08 | DESeq2 – MYC-High vs MYC-Low + apeglm LFC shrinkage | `tables/02_DGE_MYC_High_vs_Low.csv` |
| 09 | Heatmap – Top 30 DEGs (MYC-High vs Low) | `figures/Fig_03_Heatmap_MYC_High_vs_Low.tiff` |
| 10 | Volcano plots – Tumor vs Normal & MYC High/Low | `figures/Fig_04_*, Fig_05_*` |
| 11 | GO enrichment – MYC-High upregulated genes | `tables/03_GO_MYC_High_Up.csv` |
| 12 | WGCNA – Tumor samples | Module dendrogram, TOM |
| 13 | MYC-correlated module (ME4) & hub genes | `tables/04_WGCNA_ME4_hubs.csv` |
| 14 | GO enrichment – ME4 hub genes | `tables/05_GO_ME4_hubs.csv` |
| 15 | GSEA – MYC + Hallmark gene sets | `tables/06_GSEA_results.csv` |
| 16 | Survival – ME4 score KM + Cox regression | `tables/07_Cox_ME4_*.csv` |
| 17 | ML dataset export (leakage-free split) | `ML_dataset_ME4_GeneSymbols_Survival.csv` |
| 18 | ceRNA network – miRTarBase + ENCORI | `tables/11_ceRNA_network_edges.csv` |

---

## Machine Learning Module (`ml_ME4_XGBoost_SHAP.py`)

Uses the ME4 WGCNA module genes to predict **3-year mortality** in CRC patients.

**Inputs:**
- `ML_dataset_ME4_GeneSymbols_Survival.csv` — exported by Section 17A of the R pipeline

**Model:**
- XGBoost classifier with `scale_pos_weight` for class imbalance
- 75/25 stratified train/test split
- 5-fold cross-validated AUC reported

**Outputs:**
- `results/ml/Fig_ML_ROC_3yr_ME4_600dpi.tiff` — ROC curve (AUC + AUPRC)
- `results/ml/Fig_ML_SHAP_3yr_ME4_600dpi.tiff` — SHAP beeswarm (top 20 features)

---

## Requirements

### R (≥ 4.2)

Packages are auto-installed on first run. Key dependencies:

**CRAN:** `ggplot2`, `ggrepel`, `dplyr`, `stringr`, `data.table`, `tibble`, `survival`, `survminer`, `scales`

**Bioconductor:** `TCGAbiolinks`, `SummarizedExperiment`, `DESeq2`, `clusterProfiler`, `org.Hs.eg.db`, `AnnotationDbi`, `msigdbr`, `enrichplot`, `WGCNA`, `pheatmap`

### Python (≥ 3.9)

```bash
pip install xgboost shap scikit-learn pandas numpy matplotlib
```

---

## Usage

### Step 1 – Run the R pipeline

```r
# In R or RStudio, from the repo root:
source("TCGA_COAD_READ_MYC_Analysis.R")
```

> **Note:** Section 01 downloads ~5 GB of TCGA data via `TCGAbiolinks`. This may take 20–60 minutes depending on network speed. Downloaded files are cached locally by GDCdownload.

### Step 2 – Run the Python ML module

```bash
# Ensure ML_dataset_ME4_GeneSymbols_Survival.csv is in the working directory
python ml10_ME4_XGBoost_SHAP.py
```

### Step 3 (Optional) – ceRNA network (Section 18)

Download the required external files and place them in `data/external/`:

| File | Source |
|------|--------|
| `miRTarBase_MTI.csv` | https://mirtarbase.cuhk.edu.cn |
| `ENCORI_miRNA_lncRNA.tsv` | https://rnasysu.com/encori |

Then re-run Section 18 of the R script. Output tables are ready for direct import into **Cytoscape**.

---

## Outputs Summary

### Figures (`figures/`)

| File | Description |
|------|-------------|
| `Fig_01_PCA_Tumor_vs_Normal.tiff` | PCA – all samples |
| `Fig_02_Heatmap_Tumor_vs_Normal.tiff` | Heatmap – top 50 DEGs |
| `Fig_03_Heatmap_MYC_High_vs_Low.tiff` | Heatmap – top 30 MYC DEGs |
| `Fig_04_Volcano_Tumor_vs_Normal.tiff` | Volcano – Tumor vs Normal |
| `Fig_05_Volcano_MYC_High_vs_Low.tiff` | Volcano – MYC High vs Low |
| `Fig_06_GO_MYC_High_Up.tiff` | GO dotplot – MYC-High up |
| `Fig_07_WGCNA_Dendrogram.tiff` | WGCNA module dendrogram |
| `Fig_08_Module_Trait_Heatmap.tiff` | Module–trait correlation |
| `Fig_09_ME4_Hub_Network.tiff` | ME4 hub gene network |
| `Fig_10_GO_ME4_Hubs.tiff` | GO dotplot – ME4 hubs |
| `Fig_11_GSEA_MYC_Hallmark.tiff` | GSEA enrichment plot |
| `Fig_12_KM_ME4_Survival.tiff` | Kaplan–Meier – ME4 score |
| `Fig_13_ML_ClassBalance.tiff` | ML split class balance |

### Tables (`tables/`)

| File | Description |
|------|-------------|
| `01_DGE_Tumor_vs_Normal.csv` | DESeq2 results – Tumor vs Normal |
| `02_DGE_MYC_High_vs_Low.csv` | DESeq2 results – MYC High vs Low |
| `03_GO_MYC_High_Up.csv` | GO enrichment – MYC-High up |
| `04_WGCNA_ME4_hubs.csv` | ME4 hub genes |
| `05_GO_ME4_hubs.csv` | GO enrichment – ME4 hubs |
| `06_GSEA_results.csv` | GSEA results |
| `07_Cox_ME4_continuous.csv` | Cox regression – ME4 continuous |
| `08_Cox_ME4_stage_adjusted.csv` | Cox regression – stage-adjusted |
| `09_ML_train_ME4.csv` | ML training set (leakage-free) |
| `10_ML_test_ME4.csv` | ML test set |
| `11_ceRNA_network_edges.csv` | ceRNA edges (Cytoscape) |
| `12_ceRNA_network_nodes.csv` | ceRNA nodes (Cytoscape) |

---

## Key Design Decisions

**WGCNA `cor()` swap:** `WGCNA::blockwiseModules()` internally calls `WGCNA::cor()`, which accepts extra arguments that `stats::cor()` does not. The script temporarily sets `cor <- WGCNA::cor` immediately before each `blockwiseModules()` call and restores `cor <- stats::cor` immediately after to prevent conflicts with all other correlation calls in the pipeline.

**Leakage-free ML split (Section 17B):** WGCNA is run independently on the training set only. The resulting module gene list is then applied to the test set, ensuring no test-set information leaks into feature selection.

**Class imbalance:** XGBoost uses `scale_pos_weight = negatives / positives` to correct for imbalanced 3-year mortality labels.

---

## Citation

If you use this pipeline, please cite the relevant tools:

- Love MI et al. *Genome Biology* 2014 (DESeq2)
- Langfelder P & Horvath S. *BMC Bioinformatics* 2008 (WGCNA)
- Yu G et al. *OMICS* 2012 (clusterProfiler)
- Chen T & Guestrin C. *KDD* 2016 (XGBoost)
- Lundberg SM & Lee SI. *NeurIPS* 2017 (SHAP)
- Colaprico A et al. *Nucleic Acids Research* 2016 (TCGAbiolinks)

---

## License

MIT License. See `LICENSE` for details.
