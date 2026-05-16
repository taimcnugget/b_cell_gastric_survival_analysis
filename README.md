# TCGA-STAD B Cell Survival Analysis
## Overview 
This project validates transcriptional findings from a companion scRNA-seq analysis of B cell states in gastric cancer 
([bcell-gastric-tme](https://github.com/taimcnugget/bcell-gastric-tme))  using bulk RNA-seq and clinical outcome data from The Cancer Genome Atlas Stomach Adenocarcinoma (TCGA-STAD) cohort.

B cell maturation suppression was identified across disease stages in single-cell analysis. This project asks whether those transcriptional programs associate with patient survival outcomes in an independent clinical cohort of 390 TCGA-STAD patients.

## Background 
Tumor-infiltrating B cells play emerging roles in anti-tumor immunity, yet their functional states and clinical relevance in gastric cancer remain poorly characterized. Using differentially expressed genes derived  from scRNA-seq analysis of memory B cells and plasma cells across healthy, primary, and metastatic gastric cancer conditions, we scored 390 TCGA-STAD patients using single-sample GSEA (ssGSEA) and assessed associations with overall survival. Immune cell deconvolution via CIBERSORTx was applied to estimate B cell infiltration fractions and further investigate the  relationship between B cell abundance and clinical outcomes.

## Data Sources
- **Bulk RNA-seq:** TCGA-STAD STAR-Counts (GDC Data Portal, n=448 files)
- **Clinical data:** GDC API (overall survival, vital status, tumor stage, gender)
- **B cell signatures:** DEGs from GSE163558 scRNA-seq analysis
  ([bcell-gastric-tme](https://github.com/taimcnugget/bcell-gastric-tme))

## Methods & Notebooks
### Notebook 01: Data Prep
- Downloaded 448 TCGA-STAD STAR-Counts TSV files via GDC client
- Built expression matrix from raw counts using data.table
- Pulled clinical data via GDC REST API
- Mapped file UUIDs to TCGA patient barcodes
- Deduplicated samples (390 unique patients retained)
- Built TPM matrix for CIBERSORTx input
- Built ENSEMBL → HUGO gene symbol mapping

### Notebook 02: Normalization
- Variance stabilizing transformation (VST) via DESeq2
- Filtered lowly expressed genes (≥10 counts in ≥10 samples)
- Converted ENSEMBL IDs to HUGO gene symbols
- Removed duplicate gene symbols (251 duplicates dropped)
- Final VST matrix: 39,896 genes × 390 samples

### Notebook 03: Survival Analysis
- Derived memory B cell and plasma cell signatures from scRNA-seq DEGs
- Scored 390 patients using ssGSEA (GSVA package, ssgseaParam)
- Stratified patients by median signature score
- Kaplan-Meier survival curves with log-rank test (survminer)
- Cox proportional hazards modeling — univariate and multivariate
- Sex-stratified sensitivity analysis

### Notebook 04: CIBERSORTx Deconvolution **(in progress)**
- TPM matrix uploaded to CIBERSORTx portal (LM22 signature matrix)
- B cell naive, memory, and plasma cell fractions extracted per patient
- Merged with clinical survival data
- Kaplan-Meier and Cox survival analysis on B cell infiltration fractions

## Results
### Cohort 
390 TCGA-STAD patients with matched bulk RNA-seq and clinical data were included in this survival analysis. In this cohort, there were 161 death events (41% event rate). In total, patients were stratified by median overall survival time of 424 days(~ 1 day and 2 months). This is an unfortuate outcome for this type of cancer and this analysis aims to assess whether there are signatures in B cells that can be used to asssess survival outcomes in stomach adeocarcinsoma patients. 

### Signature Analysis 
B cell transcriptional signatures were derived from differentially expresssed genes (DEG) identified in the GSE163558 scRNA-seq analysis (see bcell_gastric_tme):
    * **Memory B cell signature:** 301 upregulated genes (~84% of total DEG genes) in primary tumor vs. healthy memory B cells (avg_log2FC > 0.5, p_val_adj_global < 0.05) (total DEG = 367)
    * **Plasma cell signature:** 531 downregulated genes (~91% of total DEG genes) in primary tumor vs. healthy plasma cells (avg_log2FC > 0.5, p_val_adj_global < 0.05) (total DEG = 581)

### ssGEA Scoring 
Single-sample GSEA (ssGSEA) was applied to VST-normalized bulk RNA-seq expression data to score each patient on memory B cell and plasma cell transcriptional activity. Patients were stratified into high and low groups by median signature score.

### Survival Analysis
Kaplan-Meier survival curves and Cox proportional hazards models were used to assess association between signature scores and overall survival.

| Analysis | Memory B Cell | Plasma Cell |
| --------- | ----------- | ------------ |
| Kaplan-Meier (median split) | p = 0.24 | p = 0.33 |
| Cox (continuous score) | p = 0.60 | p = 0.30 | 
| Cox (combined model) | p = 0.60 | ----- | 

Neither memory B cell nor plasma cell transcriptional signatures significantly associated with overall survival in TCGA-STAD bulk RNA-seq across Kaplan-Meier stratification or Cox proportional hazards modeling.

**Note:** The non-sig result is consistnet with known limitations of applying sc-derived transcriptional signatures to bulk RNA-seq data, where immune cell signals are diluted by the dominant epithelial and stromal compartments. hese findings motivate immune cell fraction estimation via CIBERSORTx deconvolution (Notebook 04, in progress), which isolates B cell-specific signals prior to survival analysis.

### Requirements

```r
# CRAN
install.packages(c("tidyverse", "data.table", "survival",  "survminer", "httr", "jsonlite"))

# Bioconductor
BiocManager::install(c("DESeq2", "GSVA", "BiocParallel"))
```

### Related Projects
- [bcell-gastric-tme](https://github.com/taimcnugget/bcell-gastric-tme)