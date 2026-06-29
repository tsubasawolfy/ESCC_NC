# DepMap niclosamide sensitivity and ESCC RNA-expression correlation analysis.
# Outputs are written to Data/DepMap and Figures/DepMap relative to the current
# working directory.

suppressPackageStartupMessages({
  library(depmap)
  library(ggplot2)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
})

required_packages <- c("depmap", "ggplot2", "stringr", "dplyr", "tidyr", "purrr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages before running: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

data_dir <- file.path("Data", "DepMap")
figure_dir <- file.path("Figures", "DepMap")
rds_dir <- file.path("RDS", "DepMap")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

read_or_fetch_depmap <- function(file_name, fetch_fun) {
  path <- file.path(rds_dir, file_name)
  if (file.exists(path)) {
    return(readRDS(path))
  }
  object <- fetch_fun()
  saveRDS(object, path)
  object
}

TPM <- read_or_fetch_depmap("TPM_data.rds", depmap::depmap_TPM)
metadata <- read_or_fetch_depmap("metadata.rds", depmap::depmap_metadata)
drug_sensitivity <- read_or_fetch_depmap("drug_sensitivity.rds", depmap::depmap_drug_sensitivity)

drug_subset <- drug_sensitivity %>%
  filter(str_detect(compound, fixed("BRD-K35960502-001-20")))

metadata_subset <- metadata %>%
  filter(primary_disease == "Esophageal Cancer", lineage_subtype == "esophagus_squamous")

merged_data <- inner_join(drug_subset, metadata_subset, by = "depmap_id")
merged_tpm_data <- inner_join(merged_data, TPM, by = "depmap_id") %>%
  select(!duplicated(names(.)))

run_cor_tests <- function(dependency, expression) {
  ok <- is.finite(dependency) & is.finite(expression)
  dependency <- dependency[ok]
  expression <- expression[ok]

  if (length(dependency) < 3 || sd(dependency) == 0 || sd(expression) == 0) {
    return(tibble(
      n = length(dependency),
      Spearman_Correlation = NA_real_,
      Spearman_P_Value = NA_real_,
      Pearson_Correlation = NA_real_,
      Pearson_P_Value = NA_real_
    ))
  }

  spearman <- suppressWarnings(cor.test(dependency, expression, method = "spearman", exact = FALSE))
  pearson <- suppressWarnings(cor.test(dependency, expression, method = "pearson"))

  tibble(
    n = length(dependency),
    Spearman_Correlation = unname(spearman$estimate),
    Spearman_P_Value = spearman$p.value,
    Pearson_Correlation = unname(pearson$estimate),
    Pearson_P_Value = pearson$p.value
  )
}

correlation_results <- merged_tpm_data %>%
  group_by(gene_name) %>%
  summarise(run_cor_tests(dependency, rna_expression), .groups = "drop") %>%
  mutate(
    Spearman_FDR = p.adjust(Spearman_P_Value, method = "BH"),
    Pearson_FDR = p.adjust(Pearson_P_Value, method = "BH")
  )

significant_results <- correlation_results %>%
  filter(
    !is.na(Spearman_Correlation),
    Spearman_P_Value < 0.05,
    !is.na(Pearson_Correlation),
    Pearson_P_Value < 0.05
  )

write_csv(correlation_results, file.path(data_dir, "Niclosamide_ESCC_allgene_correlations.csv"))
write_csv(significant_results, file.path(data_dir, "Niclosamide_ESCC_allgene_significant_results.csv"))

data_sorted <- merged_data %>%
  arrange(desc(dependency))

niclosamide_bar_plot <- ggplot(data_sorted, aes(x = reorder(cell_line_name, dependency), y = dependency)) +
  geom_col(fill = "lightblue", color = "black", alpha = 0.8) +
  labs(x = "ESCC cell lines", y = "Niclosamide dependency") +
  coord_flip() +
  theme_bw() +
  theme(
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_blank()
  )

ggsave(
  file.path(figure_dir, "Niclosamide_ESCC_dependency_barplot.png"),
  niclosamide_bar_plot,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600
)
print(niclosamide_bar_plot)

significant_plot_data <- significant_results %>%
  mutate(
    Spearman_neg_log10_pval = -log10(Spearman_P_Value),
    Pearson_neg_log10_pval = -log10(Pearson_P_Value)
  )

find_thresholds <- function(correlations) {
  c(
    max(correlations[correlations < 0], na.rm = TRUE),
    min(correlations[correlations > 0], na.rm = TRUE)
  )
}

if (nrow(significant_plot_data) > 0) {
  spearman_thresholds <- find_thresholds(significant_plot_data$Spearman_Correlation)
  pearson_thresholds <- find_thresholds(significant_plot_data$Pearson_Correlation)

  spearman_volcano <- ggplot(significant_plot_data, aes(x = Spearman_Correlation, y = Spearman_neg_log10_pval)) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = spearman_thresholds, linetype = "dashed", color = "blue") +
    labs(x = "Spearman correlation coefficient", y = "-log10(Spearman p-value)") +
    theme_minimal()

  pearson_volcano <- ggplot(significant_plot_data, aes(x = Pearson_Correlation, y = Pearson_neg_log10_pval)) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = pearson_thresholds, linetype = "dashed", color = "blue") +
    labs(x = "Pearson correlation coefficient", y = "-log10(Pearson p-value)") +
    theme_minimal()

  ggsave(file.path(figure_dir, "Spearman_significant_gene_volcano.png"), spearman_volcano, width = 6, height = 5, dpi = 600)
  ggsave(file.path(figure_dir, "Pearson_significant_gene_volcano.png"), pearson_volcano, width = 6, height = 5, dpi = 600)
  print(spearman_volcano)
  print(pearson_volcano)
}

genes_to_plot <- c("EHMT2", "MEN1")

for (gene_name_i in genes_to_plot) {
  tpm_subset <- TPM %>% filter(gene_name == gene_name_i)

  gene_data <- inner_join(merged_data, tpm_subset, by = "depmap_id") %>%
    select(!duplicated(names(.))) %>%
    mutate(cell_line_label = str_split(cell_line.x, "_", simplify = TRUE)[, 1])

  spearman <- cor(gene_data$dependency, gene_data$rna_expression, method = "spearman", use = "complete.obs")
  pearson <- cor(gene_data$dependency, gene_data$rna_expression, method = "pearson", use = "complete.obs")

  scatter_plot <- ggplot(gene_data, aes(x = dependency, y = rna_expression)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    geom_text(aes(label = cell_line_label), hjust = 1.2, vjust = 0.5) +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, label = paste("Spearman:", round(spearman, 3))) +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 3, label = paste("Pearson:", round(pearson, 3))) +
    labs(
      x = "Niclosamide dependency",
      y = "RNA expression log2(TPM + 1)",
      title = paste("Niclosamide dependency vs RNA expression in ESCC:", gene_name_i)
    ) +
    theme_bw()

  ggsave(file.path(figure_dir, paste0(gene_name_i, "_dependency_expression_scatter.png")), scatter_plot, width = 6, height = 5, dpi = 600)
  print(scatter_plot)
}

writeLines(capture.output(sessionInfo()), file.path(data_dir, "sessionInfo_DepMap.txt"))