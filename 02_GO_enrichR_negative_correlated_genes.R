# GO enrichment for the historical dual nominal-p selection and negative Pearson direction.
# Offline: DEPMAP_GO_CACHE_DIR must contain three result CSVs and negative_pearson_genes.csv.
# Live: enrichR is optional; each fresh query is archived with its exact input gene list.
if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', '.UTF-8'))
read_go_csv <- function(path) {
  if (!file.exists(path)) stop("Missing CSV: ", path, call. = FALSE)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}
validate_query <- function(d) {
  required <- c("gene_name", "Pearson_Correlation", "Spearman_Correlation", "Pearson_P_Value", "Spearman_P_Value")
  missing <- setdiff(required, names(d))
  if (length(missing)) stop("Input missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!nrow(d) || anyNA(d$gene_name) || any(!nzchar(trimws(d$gene_name))) || anyDuplicated(d$gene_name))
    stop("Input requires nonempty, unique gene symbols.", call. = FALSE)
  for (column in required[-1]) {
    if (!is.numeric(d[[column]]) || any(!is.finite(d[[column]]))) stop("Invalid numeric input: ", column, call. = FALSE)
    limit <- if (grepl("Correlation", column)) c(-1, 1) else c(0, 1)
    if (any(d[[column]] < limit[1] | d[[column]] > limit[2])) stop("Out-of-range input: ", column, call. = FALSE)
  }
  if (any(d$Pearson_P_Value >= .05 | d$Spearman_P_Value >= .05))
    stop("Significant-results input must satisfy BOTH nominal Pearson and Spearman p < 0.05.", call. = FALSE)
  out <- d[d$Pearson_Correlation < 0, , drop = FALSE]
  if (!nrow(out)) stop("No negative Pearson-correlated genes found.", call. = FALSE)
  out
}
normalize_go <- function(d, ontology, genes) {
  prefix <- paste0(ontology, ".")
  names(d) <- sub(paste0("^", gsub("\\.", "\\\\.", prefix)), "", names(d))
  required <- c("Term", "Overlap", "P.value", "Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes")
  missing <- setdiff(required, names(d))
  if (length(missing)) stop("Malformed GO result for ", ontology, ": missing ", paste(missing, collapse = ", "), call. = FALSE)
  for (column in c("P.value", "Adjusted.P.value", "Odds.Ratio", "Combined.Score")) {
    if (!is.numeric(d[[column]]) || any(!is.finite(d[[column]]))) stop("Invalid GO numeric column: ", column, call. = FALSE)
  }
  for (column in c("P.value", "Adjusted.P.value")) {
    if (any(d[[column]] < 0 | d[[column]] > 1)) stop("Invalid GO p value: ", column, call. = FALSE)
  }
  if (anyNA(d$Term) || any(!nzchar(d$Term)) || anyNA(d$Genes) || any(!nzchar(d$Genes))) stop("Missing GO term or overlap genes.", call. = FALSE)
  members <- unique(unlist(strsplit(d$Genes, ";", fixed = TRUE)))
  if (length(setdiff(toupper(members), toupper(genes)))) stop("GO overlap genes are not a subset of the query; cache provenance mismatch.", call. = FALSE)
  # Remove historical write.csv row-number columns, retain all enrichment statistics.
  d <- d[, !names(d) %in% c("", "X"), drop = FALSE]
  d
}
run_go <- function() {
  root <- Sys.getenv("DEPMAP_OUTPUT_DIR", file.path("Data", "DepMap"))
  output <- file.path(root, "GO")
  input <- file.path(root, "Niclosamide_ESCC_allgene_significant_results.csv")
  if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE))) stop("GO output directory is not empty; choose a fresh DEPMAP_OUTPUT_DIR.", call. = FALSE)
  negative <- validate_query(read_go_csv(input))
  genes <- sort(unique(negative$gene_name))
  ontologies <- c("GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023")
  cache <- Sys.getenv("DEPMAP_GO_CACHE_DIR", "")
  results <- list()
  if (nzchar(cache)) {
    cached_input <- read_go_csv(file.path(cache, "negative_pearson_genes.csv"))
    cached_query <- validate_query(cached_input)
    if (nrow(cached_query) != nrow(cached_input)) stop("Cached query list contains nonnegative Pearson genes.", call. = FALSE)
    if (!identical(sort(unique(cached_query$gene_name)), genes)) stop("Cached negative gene list does not match current query.", call. = FALSE)
    for (ontology in ontologies) results[[ontology]] <- normalize_go(read_go_csv(file.path(cache, paste0("enrichR_result_", ontology, "_Negative.csv"))), ontology, genes)
    metadata_file <- file.path(cache, "GO_query_metadata.txt")
    verified_metadata <- character()
    if (file.exists(metadata_file)) {
      metadata <- readLines(metadata_file, warn = FALSE)
      get_meta <- function(key) {
        matched <- metadata[startsWith(metadata, paste0(key, "="))]
        if (length(matched) != 1L) stop("Malformed GO query metadata: ", key, call. = FALSE)
        substring(matched, nchar(key) + 2L)
      }
      if (get_meta("mode") != "live_official_enrichr_api" || get_meta("query_readback") != "verified" || get_meta("query_gene_count") != as.character(length(genes))) stop("Invalid GO query metadata.", call. = FALSE)
      query_file <- file.path(cache, "negative_pearson_genes.csv")
      if (get_meta("query_md5") != unname(tools::md5sum(query_file))) stop("GO query metadata fingerprint mismatch.", call. = FALSE)
      for (ontology in ontologies) {
        name <- paste0("enrichR_result_", ontology, "_Negative.csv")
        if (get_meta(paste0("result_md5:", name)) != unname(tools::md5sum(file.path(cache, name)))) stop("GO result metadata fingerprint mismatch.", call. = FALSE)
      }
      verified_metadata <- c("original_query_metadata=recorded fresh official API query; readback verified after case normalization; input/result fingerprints validated", paste0("query_metadata:", metadata))
    }
    source_mode <- "offline_archive_replay"
    archive <- normalizePath(cache, winslash = "/", mustWork = TRUE)
  } else {
    if (!requireNamespace("enrichR", quietly = TRUE)) stop("Live GO queries require enrichR. For offline replay set DEPMAP_GO_CACHE_DIR.", call. = FALSE)
    for (ontology in ontologies) {
      message("Querying Enrichr: ", ontology)
      response <- enrichR::enrichr(genes, ontology)
      if (!is.list(response) || !is.data.frame(response[[ontology]])) stop("Enrichr returned no ontology data frame: ", ontology, call. = FALSE)
      results[[ontology]] <- normalize_go(response[[ontology]], ontology, genes)
    }
    source_mode <- "live_enrichr_query"
    archive <- tempfile(paste0("live_", format(Sys.time(), "%Y%m%dT%H%M%S"), "_"), tmpdir = file.path(output, "archives"))
  }
  # All schemas and provenance are validated before any output is changed.
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  if (!nzchar(cache)) dir.create(archive, recursive = TRUE, showWarnings = FALSE)
  write.csv(negative, file.path(output, "negative_pearson_genes.csv"), row.names = FALSE)
  combined <- do.call(rbind, lapply(ontologies, function(o) data.frame(ontology = rep(o, nrow(results[[o]])), results[[o]], check.names = FALSE)))
  for (ontology in ontologies) {
    filename <- paste0("enrichR_result_", ontology, "_Negative.csv")
    write.csv(results[[ontology]], file.path(output, filename), row.names = FALSE)
    if (!nzchar(cache)) write.csv(results[[ontology]], file.path(archive, filename), row.names = FALSE)
  }
  write.csv(combined[combined$Adjusted.P.value < .05, , drop = FALSE], file.path(output, "GO_FDR_significant.csv"), row.names = FALSE)
  methylation <- combined[grepl("methylation", combined$Term, ignore.case = TRUE), , drop = FALSE]
  methylation$FDR_significant <- methylation$Adjusted.P.value < .05
  write.csv(methylation, file.path(output, "Methylation_GO_results.csv"), row.names = FALSE)
  write.csv(methylation[methylation$FDR_significant, , drop = FALSE], file.path(output, "Methylation_GO_FDR_significant.csv"), row.names = FALSE)
  ehmt2 <- vapply(strsplit(combined$Genes, ";", fixed = TRUE), function(x) "EHMT2" %in% x, logical(1))
  membership <- combined[ehmt2, , drop = FALSE]
  membership$FDR_significant <- membership$Adjusted.P.value < .05
  write.csv(membership, file.path(output, "EHMT2_GO_terms.csv"), row.names = FALSE)
  provenance <- c(paste0("mode=", source_mode), paste0("archive=", archive), paste0("input=", normalizePath(input, winslash = "/")), paste0("input_md5=", tools::md5sum(input)), paste0("query_gene_count=", length(genes)), "selection=both nominal p < 0.05; Pearson correlation < 0", "universe=Enrichr library default; no custom tested-gene universe supplied or verified", "FDR=Enrichr Adjusted.P.value per library; threshold < 0.05", "methylation=text matches; EHMT2=exact overlap gene membership; neither implies significant enrichment")
  if (nzchar(cache)) {
    cache_files <- file.path(cache, c("negative_pearson_genes.csv", paste0("enrichR_result_", ontologies, "_Negative.csv")))
    provenance <- c(provenance, if (length(verified_metadata)) verified_metadata else "original_query_metadata=unverified; companion list match and overlap subset checks only", paste0("cache_md5:", basename(cache_files), "=", unname(tools::md5sum(cache_files))))
  } else provenance <- c(provenance, "original_query_metadata=exact list archived by this run")
  output_files <- list.files(output, pattern = "[.]csv$", full.names = TRUE)
  provenance <- c(provenance, paste0("output_md5:", basename(output_files), "=", unname(tools::md5sum(output_files))))
  writeLines(provenance, file.path(output, "GO_provenance.txt"))
  if (!nzchar(cache)) {
    write.csv(negative, file.path(archive, "negative_pearson_genes.csv"), row.names = FALSE)
    writeLines(provenance, file.path(archive, "GO_provenance.txt"))
  }
  writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo_enrichR.txt"))
  message(length(genes), " query genes; ", nrow(methylation), " methylation text matches; ", sum(methylation$FDR_significant), " with adjusted p < 0.05.")
  invisible(results)
}
if (Sys.getenv("DEPMAP_TEST_MODE") != "1") run_go()
