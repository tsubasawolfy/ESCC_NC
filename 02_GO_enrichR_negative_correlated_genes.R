# enrichR GO enrichment for genes negatively correlated with niclosamide
# dependency in the DepMap ESCC analysis.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(enrichR)
})

required_packages <- c("dplyr", "readr", "enrichR")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages before running: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

input_csv <- file.path("Data", "DepMap", "Niclosamide_ESCC_allgene_significant_results.csv")
output_dir <- file.path("Data", "DepMap", "GO")
negative_genes_file <- file.path(output_dir, "negative_pearson_genes.csv")
methylation_results_file <- file.path(output_dir, "Methylation_GO_results.csv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_csv)) {
  stop("Missing input file: ", input_csv, ". Run 01_DepMap_niclosamide_ESCC_correlation.R first.", call. = FALSE)
}

significant_results <- read_csv(input_csv, show_col_types = FALSE)
required_columns <- c("gene_name", "Pearson_Correlation")
missing_columns <- setdiff(required_columns, names(significant_results))
if (length(missing_columns)) {
  stop("Input file is missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

negative_pearson_results <- significant_results %>%
  filter(Pearson_Correlation < 0)

write_csv(negative_pearson_results, negative_genes_file)

negative_genes <- unique(na.omit(negative_pearson_results$gene_name))
if (!length(negative_genes)) {
  stop("No negative Pearson-correlated genes were found for enrichR analysis.", call. = FALSE)
}

ontologies <- c("GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023")

results <- list()
for (ontology in ontologies) {
  message("Running enrichR: ", ontology)
  results[[ontology]] <- enrichr(gene = negative_genes, database = ontology)
  write_csv(as.data.frame(results[[ontology]]), file.path(output_dir, paste0("enrichR_result_", ontology, "_Negative.csv")))
}

methylation_go_terms <- lapply(ontologies, function(ontology) {
  result_file <- file.path(output_dir, paste0("enrichR_result_", ontology, "_Negative.csv"))
  result_df <- read_csv(result_file, show_col_types = FALSE)

  term_column <- paste0(ontology, ".Term")
  adjusted_pvalue_column <- paste0(ontology, ".Adjusted.P.value")
  odds_ratio_column <- paste0(ontology, ".Odds.Ratio")
  genes_column <- paste0(ontology, ".Genes")

  required_result_columns <- c(term_column, adjusted_pvalue_column, odds_ratio_column, genes_column)
  if (!all(required_result_columns %in% names(result_df))) {
    return(tibble())
  }

  result_df %>%
    filter(grepl("Methylation", .data[[term_column]], ignore.case = TRUE)) %>%
    transmute(
      ontology = ontology,
      Term = .data[[term_column]],
      Adjusted.P.value = .data[[adjusted_pvalue_column]],
      Odds.Ratio = .data[[odds_ratio_column]],
      Genes = .data[[genes_column]]
    )
}) %>%
  bind_rows()

print(methylation_go_terms)
write_csv(methylation_go_terms, methylation_results_file)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_enrichR.txt"))