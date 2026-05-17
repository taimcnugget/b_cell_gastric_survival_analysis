# =============================================================
# Script: 03_bcell_survival_scoring_analysis.r
# Purpose: Using memory B cells= and plasma cell DEGs from 
#          previous scRNA-seq analysis from gastric cancer  
#          (taimcnugget/bcell-gastric-tme) to score patients
#          survival outcomes. 
#
# Dataset: DEGs from scRNA-seq (GSE163558) and VST matrix (notebook 02)
#
# Input:   05_bcell_DEG_significant.csv (taimcnugget/bcell-gastric-tme)
#          02_b_cell_survival_vst_matrix.rds
# Output:  03_b_cell_survival_scoring_analysis.rds
#          03_b_cell_survival_tpm.txt
#          figures/03_b_cell_memory_survival_plot.png
#          figures/03_b_cell_plasma_survival_plot.png
# Author:  Tailynn McCarty
# Date:    May 2026
# =============================================================

# --- Dependencies
# R packages required to run this script:
# GSVA, tidyverse, survival, survminer

# --- Directories
output_directory <- "/kaggle/working"
figure_directory <- file.path(output_directory, "figures")

dir.create("figures")

# --- Import files
sample_table <- readRDS("/kaggle/input/datasets/taimcnugget/b-cell-gastric-cancer-survival-analysis/01_b_cell_survival_sample_table.rds")
vst_matrix <- readRDS("/kaggle/input/datasets/taimcnugget/b-cell-gastric-cancer-survival-analysis/02_b_cell_survival_vst_matrix.rds")
deg <- read.csv("/kaggle/input/datasets/taimcnugget/b-cell-gastric-cancer-survival-analysis/05_bcell_DEG_significant.csv")
tmp_matrix <- readRDS("/kaggle/input/datasets/taimcnugget/b-cell-gastric-cancer-survival-analysis/01_b_cell_survival_tpm_matrix.rds")

# --- Install packages
BiocManager::install("GSVA", update = FALSE, ask = FALSE)

# --- Load libraries
library(tidyverse)
library(GSVA)
library(survival)
library(survminer)

# --- Parameters
avg_log2FC_threshold <- 0.5 
p_val_adj_global_threshold <- 0.05

# --- Filter DEG file for upregulated genes in memory B cells in primary TME vs healthy
# Note: DEGs were pre-filtered to have avg_log2FC > 0.5 and avg_log2FC < 0.5 and p_val_adj_global < 0.05
#       but in case this is not true for other files imported, the filters are included in code below
memory_sig <- deg %>%
  filter(
    subtype == "Memory B cells",
    condition1 == "Primary",
    condition2 == "Healthy",
    avg_log2FC > avg_log2FC_threshold, #upregulated
    p_val_adj_global < p_val_adj_global_threshold
  ) %>%
  mutate(gene = trimws(gene)) %>%
  arrange(desc(avg_log2FC))

cat("Memory B cell primary signature genes:", nrow(memory_sig), "\n")
cat("Genes present in VST matrix:", sum(memory_sig$gene %in% rownames(vst_matrix)), "\n")

# --- Filter DEG file for downregulated genes in plasma cells in primary TME vs healthy
# Note: DEGs were pre-filtered to have avg_log2FC > 0.5 and avg_log2FC < 0.5 and p_val_adj_global < 0.05
#       but in case this is not true for other files imported, the filters are included in code below

plasma_sig <- deg %>%
  filter(
    subtype == "Plasma cells",
    condition1 == "Primary",
    condition2 == "Healthy",
    avg_log2FC < -avg_log2FC_threshold, #downregulated
    p_val_adj_global < p_val_adj_global_threshold
  ) %>%
  mutate(gene = trimws(gene)) %>%
  arrange(avg_log2FC)

cat("Plasma cell primary signature genes:", nrow(plasma_sig), "\n")
cat("Genes present in VST matrix:", sum(plasma_sig$gene %in% rownames(vst_matrix)), "\n")

# --- Removing duplicates to run
sum(duplicated(rownames(vst_matrix)))

vst_matrix <- vst_matrix[!duplicated(rownames(vst_matrix)), ]
cat("VST matrix after deduplication:", dim(vst_matrix), "\n")

# --- Scoring 
memory_genes <- memory_sig$gene
plasma_genes <- plasma_sig$gene

memory_genes <- memory_genes[memory_genes %in% rownames(vst_matrix)]
plasma_genes <- plasma_genes[plasma_genes %in% rownames(vst_matrix)]

cat("Final memory B signature:", length(memory_genes), "genes\n")
cat("Final plasma signature:", length(plasma_genes), "genes\n")

gsva_params <- ssgseaParam(
    exprData = vst_matrix,
    geneSets = list(
        memory_signature = memory_genes,
        plasma_signature = plasma_genes
    ) 
)

gsva_result <- gsva(gsva_params, verbose = TRUE)

cat("GSVA dimensions:", dim(gsva_result), "\n")

# --- Update sample table with scores
sample_table_scored <- sample_table %>%
  mutate(
    memory_score = as.numeric(gsva_result["memory_signature", file_uuid]),
    plasma_score = as.numeric(gsva_result["plasma_signature", file_uuid])
)

# --- Stratification (Kaplan-Meier)
sample_table_scored <- sample_table_scored %>%
    mutate(
        memory_group = ifelse(memory_score >= median(memory_score), "High", "Low"),
        plasma_group = ifelse(plasma_score >= median(plasma_score), "High", "Low")
    )

# - Memory B cell survival fit
fit_memory <- survfit(
    Surv(os_time, os_event) ~ memory_group,
    data = sample_table_scored
)


# - Plasma B cell survival fit
fit_plasma <- survfit(
    Surv(os_time, os_event) ~ plasma_group,
    data = sample_table_scored
)

png(filename = file.path(figure_directory, "03_b_cell_memory_survival_plot.png"),
    width = 800,
    height = 800,
    res = 150
)

ggsurvplot(
  fit_memory,
  data = sample_table_scored,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("#e08214", "#8073ac"),
  title = "Overall Survival by Memory B Cell Signature",
  xlab = "Time (days)",
  ylab = "Survival Probability",
  legend.labs = c("High", "Low"),
  ggtheme = theme_classic()
)
dev.off()

png(filename = file.path(figure_directory, "03_b_cell_plasma_survival_plot.png"),
    width = 800,
    height = 800,
    res = 150
)

ggsurvplot(
  fit_plasma,
  data = sample_table_scored,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("#e08214", "#8073ac"),
  title = "Overall Survival by Plasma Cell Signature",
  xlab = "Time (days)",
  ylab = "Survival Probability",
  legend.labs = c("High", "Low"),
  ggtheme = theme_classic()
)

dev.off()

# --- Cox models
# -- Memory
cox_memory <- coxph(Surv(os_time, os_event) ~ memory_score, data = sample_table_scored)
summary(cox_memory)

# -- Plasma
cox_plasma <- coxph(Surv(os_time, os_event) ~ plasma_score, data = sample_table_scored)
summary(cox_plasma)

# -- Combined
cox_combined <- coxph(Surv(os_time, os_event) ~ memory_score + plasma_score, data = sample_table_scored)
summary(cox_combined)

# --- Filter TMP file for deconvo
tpm_filtered <- tmp_matrix[, sample_table$file_uuid]
cat("TPM filtered dimensions:", dim(tpm_filtered), "\n")

colnames(tpm_filtered) <- sample_table$submitter_id[match(colnames(tpm_filtered), sample_table$file_uuid)]

colnames(tpm_filtered) <- gsub("\\.", "-", colnames(tpm_filtered))

write.table(tpm_filtered,
            file.path(output_directory, "03_b_cell_survival_tpm.txt"),
            sep = "\t",
            quote = FALSE, 
            col.names = NA
)
cat("deconvo input file exported\n")

# --- Save files
saveRDS(sample_table_scored, file.path(output_directory, "03_b_cell_survival_scoring_analysis.rds"))


message("Saved: 03_b_cell_survival_scoring_analysis.rds")