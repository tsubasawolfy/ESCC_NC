# KM plots with risk tables for RStudio manual export.
# This copy keeps the original grouping logic and does not save PDF/PNG files.

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(survival)
library(survminer)

ensure_tcga_workdir <- function() {
  required_file <- file.path(
    "esca_tcga_pan_can_atlas_2018",
    "data_mrna_seq_v2_rsem_zscores_ref_normal_samples.txt"
  )

  if (file.exists(required_file)) {
    return(invisible(getwd()))
  }

  script_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (!is.na(script_file)) {
    candidate_dir <- normalizePath(file.path(dirname(script_file), "..", "..", "TCGA"), mustWork = FALSE)
    if (file.exists(file.path(candidate_dir, required_file))) {
      setwd(candidate_dir)
      return(invisible(getwd()))
    }
  }

  stop(
    "Set the R working directory to the TCGA data folder, then rerun this script.",
    call. = FALSE
  )
}

read_cbioportal_table <- function(path, skip_lines = 4) {
  raw_lines <- readLines(path)
  data_lines <- raw_lines[-seq_len(skip_lines)]
  temp_file <- tempfile(fileext = ".txt")
  on.exit(unlink(temp_file), add = TRUE)
  writeLines(data_lines, temp_file)

  read.table(
    temp_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = ""
  )
}

make_km_plot <- function(fit, data, title, ylab, legend.labs) {
  ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    risk.table = TRUE,
    risk.table.title = "Number at risk",
    risk.table.height = 0.30,
    tables.theme = theme_cleantable(),
    title = title,
    xlab = "Time in Months",
    ylab = ylab,
    legend.labs = legend.labs
  )
}

ensure_tcga_workdir()

maf_data <- read.table(
  "esca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_normal_samples.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

gene_list <- read_excel("Extract_excel/EHMT2.xlsx")

metadata <- read_cbioportal_table("esca_tcga_pan_can_atlas_2018/data_clinical_sample.txt") %>%
  filter(CANCER_TYPE_DETAILED == "Esophageal Squamous Cell Carcinoma") %>%
  select(SAMPLE_ID, PATIENT_ID, CANCER_TYPE_DETAILED) %>%
  rename(Sample_Identifier = SAMPLE_ID, Patient_Identifier = PATIENT_ID)

clinical_data <- read_excel("esca_tcga_pan_can_atlas_2018/data_clinical_patient_mod.xlsx")

filtered_data <- maf_data %>%
  filter(Hugo_Symbol %in% gene_list$Gene)

colnames(filtered_data)[3:ncol(filtered_data)] <- str_replace_all(
  colnames(filtered_data)[3:ncol(filtered_data)],
  "\\.",
  "-"
)

existing_columns <- intersect(metadata$Sample_Identifier, colnames(filtered_data)[3:ncol(filtered_data)])

filtered_data2 <- filtered_data %>%
  select(Hugo_Symbol, all_of(existing_columns))

filtered_clinical_data <- clinical_data %>%
  filter(PATIENT_ID %in% metadata$Patient_Identifier) %>%
  left_join(metadata, by = c("PATIENT_ID" = "Patient_Identifier"))

transposed_data <- filtered_data2 %>%
  pivot_longer(-Hugo_Symbol, names_to = "Sample_Identifier", values_to = "Value") %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = Value)

merged_data <- left_join(filtered_clinical_data, transposed_data, by = "Sample_Identifier")
write.csv(merged_data, "merged_data_EHMT2_MYC_clinical.csv", row.names = FALSE)

#### Overall survival ####
median_cutoff <- median(merged_data$EHMT2, na.rm = TRUE)
zero_cutoff <- 0

merged_data <- merged_data %>%
  mutate(
    EHMT2_median_group = ifelse(EHMT2 >= median_cutoff, "High", "Low"),
    EHMT2_zero_group = ifelse(EHMT2 >= zero_cutoff, "High", "Low"),
    OS_STATUS_binary = ifelse(OS_STATUS == "1:DECEASED", 1, 0)
  )

surv_object_median <- Surv(time = merged_data$OS_MONTHS, event = merged_data$OS_STATUS_binary)
fit_median <- survfit(surv_object_median ~ EHMT2_median_group, data = merged_data)
plot_median <- make_km_plot(
  fit_median,
  merged_data,
  "Survival Analysis (Median Cut-off)",
  "Survival Probability",
  c("Low EHMT2", "High EHMT2")
)
print(plot_median)

myc_median_cutoff <- median(merged_data$MYC, na.rm = TRUE)
merged_data <- merged_data %>%
  mutate(MYC_median_group = ifelse(MYC >= myc_median_cutoff, "High", "Low"))

surv_object_myc <- Surv(time = merged_data$OS_MONTHS, event = merged_data$OS_STATUS_binary)
fit_myc <- survfit(surv_object_myc ~ MYC_median_group, data = merged_data)
plot_myc <- make_km_plot(
  fit_myc,
  merged_data,
  "Survival Analysis (MYC Median Cut-off)",
  "Survival Probability",
  c("Low MYC", "High MYC")
)
print(plot_myc)

merged_data <- merged_data %>%
  mutate(EHMT2_MYC_group = case_when(
    EHMT2_median_group == "High" & MYC_median_group == "High" ~ "High EHMT2 & High MYC",
    EHMT2_median_group == "Low" & MYC_median_group == "Low" ~ "Low EHMT2 & Low MYC",
    TRUE ~ "Other"
  ))

merged_data_combined <- merged_data %>%
  filter(EHMT2_MYC_group %in% c("High EHMT2 & High MYC", "Low EHMT2 & Low MYC"))
write.csv(merged_data_combined, "merged_data_EHMT2_MYC_CLinical_group.csv", row.names = FALSE)

surv_object_combined <- Surv(time = merged_data_combined$OS_MONTHS, event = merged_data_combined$OS_STATUS_binary)
fit_combined <- survfit(surv_object_combined ~ EHMT2_MYC_group, data = merged_data_combined)
plot_combined <- make_km_plot(
  fit_combined,
  merged_data_combined,
  "Survival Analysis (Combined EHMT2 and MYC Median Cut-offs)",
  "Survival Probability",
  c("Low EHMT2 & Low MYC", "High EHMT2 & High MYC")
)
print(plot_combined)

#### PFS ####
merged_data <- merged_data %>%
  mutate(PFS_STATUS_binary = ifelse(PFS_STATUS == "1:PROGRESSION", 1, 0))

pfs_object_ehmt2 <- Surv(time = merged_data$PFS_MONTHS, event = merged_data$PFS_STATUS_binary)
fit_pfs_ehmt2 <- survfit(pfs_object_ehmt2 ~ EHMT2_median_group, data = merged_data)
plot_pfs_ehmt2 <- make_km_plot(
  fit_pfs_ehmt2,
  merged_data,
  "Progression-Free Survival Analysis (EHMT2 Median Cut-off)",
  "Progression-Free Survival Probability",
  c("Low EHMT2", "High EHMT2")
)
print(plot_pfs_ehmt2)

pfs_object_myc <- Surv(time = merged_data$PFS_MONTHS, event = merged_data$PFS_STATUS_binary)
fit_pfs_myc <- survfit(pfs_object_myc ~ MYC_median_group, data = merged_data)
plot_pfs_myc <- make_km_plot(
  fit_pfs_myc,
  merged_data,
  "Progression-Free Survival Analysis (MYC Median Cut-off)",
  "Progression-Free Survival Probability",
  c("Low MYC", "High MYC")
)
print(plot_pfs_myc)

merged_data_combined <- merged_data %>%
  mutate(EHMT2_MYC_group = case_when(
    EHMT2_median_group == "High" & MYC_median_group == "High" ~ "High EHMT2 & High MYC",
    EHMT2_median_group == "Low" & MYC_median_group == "Low" ~ "Low EHMT2 & Low MYC",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(EHMT2_MYC_group))

pfs_object_combined <- Surv(time = merged_data_combined$PFS_MONTHS, event = merged_data_combined$PFS_STATUS_binary)
fit_pfs_combined <- survfit(pfs_object_combined ~ EHMT2_MYC_group, data = merged_data_combined)
plot_pfs_combined <- make_km_plot(
  fit_pfs_combined,
  merged_data_combined,
  "Progression-Free Survival Analysis (Combined EHMT2 and MYC Median Cut-offs)",
  "Progression-Free Survival Probability",
  c("Low EHMT2 & Low MYC", "High EHMT2 & High MYC")
)
print(plot_pfs_combined)

km_plots <- list(
  OS_EHMT2 = plot_median,
  OS_MYC = plot_myc,
  OS_EHMT2_MYC = plot_combined,
  PFS_EHMT2 = plot_pfs_ehmt2,
  PFS_MYC = plot_pfs_myc,
  PFS_EHMT2_MYC = plot_pfs_combined
)