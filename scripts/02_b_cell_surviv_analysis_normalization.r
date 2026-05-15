# =============================================================
# Script: 02_bcell_surviv_analysis_normalization.r
# Purpose: Normalize the expression data from patients 
#          that had accurate patient records. This will be used
#          to build the survival analysis with previously
#          determined B cell markers.
#
# Dataset: Expression data (STAR Counts) were retrived for all 390 
#          patients with matching clinical data - which can be
#          found in notebook 1.
#
# Input:   01_b_cell_survival_expression_filtered.rds
#          01_b_cell_survival_sample_table.rds
#          01_b_cell_survival_stad_gene_map.rds
# Output:  02_b_cell_survival_vst_matrix.rds
# Author:  Tailynn McCarty
# Date:    May 2026
# =============================================================

# --- Dependencies
# R packages required to run this script:
# DESeq2, BiocParallel (optional)
# BiocParallel was used to speed up VST, but it is optional

# --- Directories
output_directory <- "/kaggle/working"

# --- Import files
gene_map <- readRDS("/kaggle/input/datasets/taimcnugget/b-cell-gastric-cancer-survival-analysis/01_b_cell_survival_stad_gene_map.rds")
expr_filtered <- readRDS("/kaggle/input/datasets/taimcnugget/b-cell-gastric-cancer-survival-analysis/01_b_cell_survival_expression_filtered.rds")
sample_table <- readRDS("/kaggle/input/datasets/taimcnugget/b-cell-gastric-cancer-survival-analysis/01_b_cell_survival_sample_table.rds")

# --- Parameters
smallest_gene_count <- 10
smallest_group_size <- 10

# --- Install packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "BiocParallel"))

# --- Download libraries
library(DESeq2)
library(BiocParallel)

# --- Start parallel cores
register(MulticoreParam(workers = 3))

# --- Normalize expression data
dds <- DESeqDataSetFromMatrix(countData = round(expr_filtered),
                              colData = sample_table,
                              design = ~ 1
)

# -- Remove lowly expressed genes
keep <- rowSums(counts(dds) >= smallest_gene_count) >= smallest_group_size 
dds <- dds[keep,]
cat("Genes after filtering:", nrow(dds), "\n")

# --- Run VST
vst_matrix <- assay(vst(dds, blind = TRUE))
cat("VST matrix dimensions:", dim(vst_matrix), "\n")

# --- Map ENSEMBL to gene symbols (needed for next notebook)
rownames(vst_matrix) <- sub("\\..*", "", rownames(vst_matrix))
matched <- gene_map[match(rownames(vst_matrix), gene_map$gene_id_clean), ]
valid <- !is.na(matched$gene_name) & matched$gene_name != ""
rownames(vst_matrix)[valid] <- matched$gene_name[valid]

# --- Save results
saveRDS(vst_matrix, file.path(output_directory, "02_b_cell_survival_vst_matrix.rds"))
cat("VST complete\n")

message("Saved: 02_b_cell_survival_vst_matrix.rds")