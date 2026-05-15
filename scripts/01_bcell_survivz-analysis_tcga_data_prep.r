# =============================================================
# Script: 01_bcell_surviv_analysis_tcga_data_prep.R
# Purpose: Download the STAD dataset from The Cancer Genome Atlas
#          Stomach Adenocarcinoma to use for survival analysis with
#          previously generated B cell analysis from gastric cancer
#          patients.
#
# Dataset: The Cancer Genome Atlas - Stomach Adenocarcinoma cohort data
#
# Input:   tcga-stad.zip
# Output:  01_b_cell_survival_stad_gene_map.rds
#          01_b_cell_survival_stad_expression_matrix.rds
#          01_b_cell_survival_stad_clinical_cohort.rds
#          01_b_cell_survival_sample_table.rds
#          01_b_cell_survival_expression_filtered.rds
#          01_b_cell_survival_tpm_matrix.rds
#
# Author:  Tailynn McCarty
# Date:    May 2026
# =============================================================

# --- Dependencies
# R packages required to run this script:
# tidyverse, data.table, httr, jsonlite

# --- Directories
output_directory <- "/kaggle/working"

# ---- Load libraries
library(tidyverse)
library(data.table)
library(httr)
library(jsonlite)

# --- Import files 
# --- Note, this step is added here because there was 
#     an issue with importing the files directly from the TCGA using TCAbiolinks
#     The files were downloaded locally via manifest and zipped using bash commands

files <- list.files(
  "/kaggle/input/datasets/taimcnugget/tcga-stad",
  pattern = "rna_seq.augmented_star_gene_counts.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Total TSV files found:", length(files), "\n")

# --- Create expression matrix
file_read <- function (f) {
    dt <- fread(f, skip =1, header = TRUE)
    dt <- dt[grepl("^ENSG", gene_id)]
    dt <- dt[, .(gene_id, unstranded)]
    sample_id <- basename(dirname(f))
    setnames(dt, "unstranded", sample_id)
    return(dt)
}

all_files <- lapply(files, file_read)

expression_matrix <- Reduce(
    function(a,b) merge(a,b, by = "gene_id"),
    all_files)

expression_matrix <- expression_matrix %>%
  column_to_rownames("gene_id")

expression_matrix <- as.matrix(expression_matrix)

cat("Dimensions:", dim(expression_matrix), "\n")

# --- Build ENSEMBL to gene symbol mapping
gene_map <- fread(files[1], skip = 1, header = TRUE) %>%
    filter(grepl("^ENSG", gene_id)) %>%
    select(gene_id, gene_name) %>%
    mutate(
        gene_id_clean = sub("\\..*", "", gene_id),
        gene_name = trimws(gene_name)
)
cat("Gene map saved:", nrow(gene_map), "genes\n")

# --- Pull TCGA-STAD clinical data via GDC API
demo_response <- GET(
  "https://api.gdc.cancer.gov/cases",
  query = list(
    filters = '{"op":"in","content":{"field":"project.project_id","value":["TCGA-STAD"]}}',
    fields = "submitter_id,demographic.vital_status,demographic.days_to_death,demographic.gender",
    size = 500,
    format = "JSON"
  )
)

demo_raw <- fromJSON(content(demo_response, "text"))
demo <- demo_raw$data$hits

# Pull diagnoses separately
diag_response <- GET(
  "https://api.gdc.cancer.gov/cases",
  query = list(
    filters = '{"op":"in","content":{"field":"project.project_id","value":["TCGA-STAD"]}}',
    fields = "submitter_id,diagnoses.days_to_last_follow_up,diagnoses.ajcc_pathologic_stage",
    size = 500,
    format = "JSON"
  )
)

diag_raw <- fromJSON(content(diag_response, "text"))
diag <- diag_raw$data$hits

# Check both came back correctly
str(demo$demographic, max.level = 1)
str(diag$diagnoses, max.level = 1)

extraction <- function(list_col, field) {
    sapply(list_col, function(x) {
        if(is.data.frame(x) && field %in% colnames(x)) x[[field]][1] else NA
    })
}

clinical_clean <- data.frame(
    case_id = demo$id,
    submitter_id = demo$submitter_id,
    vital_status = demo$demographic$vital_status,
    days_to_death = demo$demographic$days_to_death,
    days_to_last_followup = extraction(diag$diagnoses, "days_to_last_follow_up"),
    stage = extraction(diag$diagnoses, "ajcc_pathologic_stage"),
    gender = demo$demographic$gender,
  stringsAsFactors = FALSE
)

# --- Build survival columns
clinical_final <- clinical_clean %>%
    mutate(
        os_time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_followup),
        os_event = ifelse(trimws(vital_status) == "Dead", 1, 0),
        remove_record = (os_event == 1 & os_time == 0)
    ) %>%
    filter(!remove_record, !is.na(os_time), os_time > 0)

cat("Clinical cohort:", nrow(clinical_final), "patients\n")

# --- Map expression data (UUID) to patient records (TCGA)
response_map <- GET(
    "https://api.gdc.cancer.gov/files",
    query = list(
        filters = '{"op":"and", "content":[
        {"op":"in", "content":{"field":"cases.project.project_id","value":["TCGA-STAD"]}},
        {"op":"in", "content":{"field":"data_type","value":["Gene Expression Quantification"]}},
        {"op":"in", "content":{"field":"analysis.workflow_type","value":["STAR - Counts"]}}
        ]}',
        fields = "file_id,cases.submitter_id",
        size = 500, 
        format ="JSON"
    )
)

mapping_raw <- fromJSON(content(response_map, "text"))
mapping <- mapping_raw$data$hits

mapping_clean <- data.frame(
    file_uuid = mapping$file_id,
    submitter_id = extraction(mapping$cases, "submitter_id"),
    stringsAsFactors = FALSE
)

cat("Mapping rows:", nrow(mapping_clean), "\n")
cat("Overlap with expression matrix:", 
    sum(mapping_clean$file_uuid %in% colnames(expression_matrix)), "\n")

sample_table <- mapping_clean %>%
    inner_join(clinical_final, by = "submitter_id")

cat("Samples with expression + clinical:", nrow(sample_table), "\n")
expr_filtered <- expression_matrix[, sample_table$file_uuid]
cat("Expression matrix dimensions:", dim(expr_filtered), "\n")

# --- Remove duplicates found after mapping
sample_table$total_counts <- colSums(expr_filtered)[sample_table$file_uuid]

# Keep one sample per patient - highest total counts wins
sample_table_cleaned <- sample_table %>%
  group_by(submitter_id) %>%
  slice_max(total_counts, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("After deduplication:", nrow(sample_table_cleaned), "\n")

# Re-filter expression matrix
expr_filtered <- expr_filtered[, sample_table_cleaned$file_uuid]
cat("Expression matrix dimensions:", dim(expr_filtered), "\n")

# --- TPM matrix for later analysis
tpm_read <- function (f){
    dt <- fread(f, skip =1, header = TRUE)
    dt <- dt[grepl("^ENSG", gene_id)]
    dt <- dt[, .(gene_id, gene_name, tpm_unstranded)]
    sample_id <- basename(dirname(f))
    setnames(dt, "tpm_unstranded", sample_id)
    return(dt)
}

all_tpm <- lapply(files, tpm_read)
tpm_matrix <- Reduce(function(a,b) merge(a,b, by = c("gene_id", "gene_name")), all_tpm)

# - EcoTyper needs gene symbols
tpm_matrix <- tpm_matrix %>%
    select(-gene_id) %>%
    filter(trimws(gene_name) != "") %>%
    distinct(gene_name, .keep_all = TRUE) %>%
    column_to_rownames("gene_name")

cat("TMP matrix dimensions:", dim(tpm_matrix), "\n")

# --- Save files
saveRDS(gene_map, file.path(output_directory, "01_b_cell_survival_stad_gene_map.rds"))
saveRDS(expression_matrix, file.path(output_directory, "01_b_cell_survival_stad_expression_matrix.rds"))
saveRDS(clinical_final, file.path(output_directory, "01_b_cell_survival_stad_clinical_cohort.rds"))
saveRDS(sample_table_cleaned, file.path(output_directory, "01_b_cell_survival_sample_table.rds"))
saveRDS(expr_filtered, file.path(output_directory, "01_b_cell_survival_expression_filtered.rds"))
saveRDS(tpm_matrix, file.path(output_directory, "01_b_cell_survival_tpm_matrix.rds"))



message("Saved:
        01_b_cell_survival_stad_gene_map.rds
        01_b_cell_survival_stad_expression_matrix.rds, 
        01_b_cell_survival_stad_clinical_cohort.rds, 
        01_b_cell_survival_sample_table.rds, 
        01_b_cell_survival_expression_filtered.rds,
        01_b_cell_survival_tpm_matrix.rds"
)