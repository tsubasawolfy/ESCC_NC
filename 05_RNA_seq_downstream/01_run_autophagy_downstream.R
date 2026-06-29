ensure_packages <- function(cran = character(), bioc = character()) {
  missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  missing_packages <- c(missing_cran, missing_bioc)
  if (length(missing_packages)) {
    stop(
      "Install required packages before running this publication script: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
}

ensure_packages(
  cran = c("ggplot2", "ggrepel", "dplyr", "readr", "stringr", "pheatmap", "tidyr", "tibble", "msigdbr"),
  bioc = c("clusterProfiler", "ReactomePA", "reactome.db", "org.Hs.eg.db", "enrichplot", "fgsea", "ComplexHeatmap")
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(readr)
  library(stringr)
  library(pheatmap)
})

project_dir <- Sys.getenv("AUTOPHAGY_PROJECT_DIR", unset = "F:/RNA_Seq/Autophagy")
analysis_dir <- file.path(project_dir, "Analysis")
script_dir <- file.path(analysis_dir, "Script")
input_dir <- file.path(project_dir, "downstream_inputs", "nfcore_star_salmon_20260616")
tpm_file <- file.path(input_dir, "expression_matrices", "salmon.merged.gene_tpm.tsv")
metadata_file <- file.path(input_dir, "metadata", "sample_metadata.tsv")
gtf_file <- file.path(project_dir, "workflow", "references", "gencode_v49", "gencode.v49.primary_assembly.annotation.gtf.gz")

dirs <- c(
  "01_differential_expression",
  "02_log2fc_scatter",
  "03_pathway_enrichment",
  "04_gsea",
  "05_heatmaps",
  "06_maplots",
  "07_cancer_signature_analysis",
  "08_immune_gene_analysis",
  "09_emt_invasion_score",
  "validation",
  "versions"
)
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(script_dir, showWarnings = FALSE, recursive = TRUE)

clean_existing_outputs <- tolower(Sys.getenv("AUTOPHAGY_CLEAN_OUTPUTS", "false")) %in% c("1", "true", "yes")
if (clean_existing_outputs) {
  for (path in file.path(analysis_dir, dirs)) {
    if (dir.exists(path)) {
      old_files <- list.files(path, all.files = TRUE, no.. = TRUE, full.names = TRUE)
      if (length(old_files)) unlink(old_files, recursive = TRUE, force = TRUE)
    }
  }
}
invisible(vapply(file.path(analysis_dir, dirs), dir.create, logical(1), showWarnings = FALSE, recursive = TRUE))

samples <- c("29-81T-Con", "30-81T-NC", "31-OE21-Con", "32-OE21-NC", "33-KYSE-Con", "34-KYSE-NC")
comparisons <- tibble::tribble(
  ~comparison, ~label, ~cell_line, ~control, ~nc,
  "ComparisonA", "NC_Treatment", "81T", "29-81T-Con", "30-81T-NC",
  "ComparisonB", "NC_Treatment", "OE21", "31-OE21-Con", "32-OE21-NC",
  "ComparisonC", "NC_Treatment", "KYSE", "33-KYSE-Con", "34-KYSE-NC"
)
pseudocount <- 0.1
log2fc_cutoff <- 1.5

write_skip <- function(path, reason) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  writeLines(reason, file.path(path, "SKIPPED.txt"))
}

save_plot <- function(plot, filename, width = 7, height = 5) {
  ggsave(filename, plot, width = width, height = height, units = "in", dpi = 300)
}

extract_attr <- function(attrs, key) {
  m <- regexec(paste0(key, ' "([^"]+)"'), attrs, perl = TRUE)
  hits <- regmatches(attrs, m)
  vapply(hits, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1))
}

read_gencode_gene_annotation <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  out <- list()
  idx <- 1L
  repeat {
    chunk <- readLines(con, n = 50000, warn = FALSE)
    if (!length(chunk)) break
    gene_lines <- chunk[grepl("\tgene\t", chunk, fixed = TRUE)]
    if (!length(gene_lines)) next
    fields <- strsplit(gene_lines, "\t", fixed = TRUE)
    attrs <- vapply(fields, function(x) x[[9]], character(1))
    out[[idx]] <- tibble::tibble(
      gene_id = extract_attr(attrs, "gene_id"),
      gene_id_base = sub("\\..*$", "", gene_id),
      gencode_gene_name = extract_attr(attrs, "gene_name"),
      gene_type = extract_attr(attrs, "gene_type")
    )
    idx <- idx + 1L
  }
  bind_rows(out) %>% distinct(gene_id_base, .keep_all = TRUE)
}

required_files <- c(tpm_file, metadata_file, gtf_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

metadata <- read.delim(metadata_file, check.names = FALSE)
tpm_raw <- read.delim(tpm_file, check.names = FALSE)
missing_samples <- setdiff(samples, names(tpm_raw))
if (length(missing_samples)) stop("Missing TPM sample columns: ", paste(missing_samples, collapse = ", "))
if (!all(samples %in% metadata$sample)) stop("Sample metadata does not contain all expected samples.")

annotation <- read_gencode_gene_annotation(gtf_file)
genes <- tpm_raw %>%
  mutate(gene_id_base = sub("\\..*$", "", gene_id)) %>%
  inner_join(annotation, by = "gene_id_base") %>%
  rename(gene_id = gene_id.x, gencode_gene_id = gene_id.y) %>%
  filter(gene_type == "protein_coding") %>%
  mutate(
    gene_name = ifelse(is.na(gene_name) | gene_name == "", gencode_gene_name, gene_name),
    gene_label = ifelse(is.na(gene_name) | gene_name == "", gene_id_base, gene_name),
    gene_label_unique = ifelse(duplicated(gene_label) | duplicated(gene_label, fromLast = TRUE),
                               paste0(gene_label, " (", gene_id_base, ")"),
                               gene_label)
  ) %>%
  distinct(gene_id_base, .keep_all = TRUE)

ehmt2_ids <- genes %>%
  filter(toupper(gene_name) %in% c("EHMT2", "G9A")) %>%
  pull(gene_id_base)

log_expr <- as.matrix(genes[, samples])
rownames(log_expr) <- genes$gene_id_base
log_expr <- log2(log_expr + pseudocount)

calc_comparison <- function(comp_row) {
  control <- comp_row$control
  nc <- comp_row$nc
  genes %>%
    transmute(
      gene_id,
      gene_id_base,
      gene_name,
      gene_label,
      gene_label_unique,
      gene_type,
      comparison = comp_row$comparison,
      cell_line = comp_row$cell_line,
      control_sample = control,
      nc_sample = nc,
      control_tpm = .data[[control]],
      nc_tpm = .data[[nc]],
      mean_tpm = (control_tpm + nc_tpm) / 2,
      pass_expression_filter = control_tpm >= 1 | nc_tpm >= 1,
      log2FC = log2((nc_tpm + pseudocount) / (control_tpm + pseudocount)),
      direction = case_when(
        pass_expression_filter & log2FC >= log2fc_cutoff ~ "up_in_NC",
        pass_expression_filter & log2FC <= -log2fc_cutoff ~ "down_in_NC",
        TRUE ~ "not_DE"
      ),
      is_EHMT2_G9a = gene_id_base %in% ehmt2_ids
    )
}

de_all <- bind_rows(lapply(seq_len(nrow(comparisons)), function(i) calc_comparison(comparisons[i, ])))
de_dir <- file.path(analysis_dir, "01_differential_expression")
readr::write_tsv(de_all, file.path(de_dir, "all_comparisons_full_log2fc.tsv"))

gene_sets <- list()
for (comp in comparisons$comparison) {
  tab <- de_all %>% filter(comparison == comp)
  readr::write_tsv(tab, file.path(de_dir, paste0(comp, "_full_log2fc.tsv")))
  up <- tab %>% filter(direction == "up_in_NC")
  down <- tab %>% filter(direction == "down_in_NC")
  readr::write_tsv(up, file.path(de_dir, paste0(comp, "_up_in_NC_log2FC_ge_", log2fc_cutoff, ".tsv")))
  readr::write_tsv(down, file.path(de_dir, paste0(comp, "_down_in_NC_log2FC_le_minus_", log2fc_cutoff, ".tsv")))
  gene_sets[[paste0(comp, "_up")]] <- up$gene_id_base
  gene_sets[[paste0(comp, "_down")]] <- down$gene_id_base
}

up_overlap_ids <- Reduce(intersect, gene_sets[paste0(comparisons$comparison, "_up")])
down_overlap_ids <- Reduce(intersect, gene_sets[paste0(comparisons$comparison, "_down")])
overlap_table <- de_all %>%
  filter(gene_id_base %in% c(up_overlap_ids, down_overlap_ids)) %>%
  mutate(overlap_group = ifelse(gene_id_base %in% up_overlap_ids, "reciprocal_up_in_NC", "reciprocal_down_in_NC")) %>%
  dplyr::select(overlap_group, comparison, gene_id, gene_id_base, gene_name, log2FC, control_tpm, nc_tpm, direction)
readr::write_tsv(overlap_table, file.path(de_dir, "reciprocal_overlap_long.tsv"))
readr::write_tsv(de_all %>% filter(gene_id_base %in% up_overlap_ids), file.path(de_dir, "reciprocal_overlap_up_in_NC.tsv"))
readr::write_tsv(de_all %>% filter(gene_id_base %in% down_overlap_ids), file.path(de_dir, "reciprocal_overlap_down_in_NC.tsv"))

readme_de <- c(
  "# Exploratory TPM Fold-Change Differential Expression",
  "",
  "This analysis uses one sample per condition and therefore does not run DESeq2, p-value testing, or FDR filtering.",
  "",
  paste0("Formula: log2FC = log2((NC TPM + ", pseudocount, ") / (Control TPM + ", pseudocount, "))"),
  paste0("Expression filter: TPM >= 1 in at least one sample in the comparison."),
  paste0("DE threshold: absolute log2FC >= ", log2fc_cutoff, "."),
  "Positive log2FC means higher expression after Niclosamide treatment.",
  "Only protein-coding genes from local GENCODE v49 annotation are analyzed."
)
writeLines(readme_de, file.path(de_dir, "README.md"))

scatter_dir <- file.path(analysis_dir, "02_log2fc_scatter")
fc_wide <- de_all %>%
  dplyr::select(gene_id_base, gene_name, gene_label, comparison, log2FC, pass_expression_filter, is_EHMT2_G9a) %>%
  tidyr::pivot_wider(names_from = comparison, values_from = c(log2FC, pass_expression_filter))
scatter_pairs <- list(
  c("ComparisonA", "ComparisonB"),
  c("ComparisonB", "ComparisonC"),
  c("ComparisonA", "ComparisonC")
)
for (pair in scatter_pairs) {
  xcol <- paste0("log2FC_", pair[[1]])
  ycol <- paste0("log2FC_", pair[[2]])
  pdat <- fc_wide %>%
    mutate(any_large_fc = abs(.data[[xcol]]) >= log2fc_cutoff | abs(.data[[ycol]]) >= log2fc_cutoff)
  eh <- pdat %>% filter(is_EHMT2_G9a)
  p <- ggplot(pdat, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_point(aes(color = any_large_fc), alpha = 0.55, size = 1.1) +
    scale_color_manual(values = c("FALSE" = "grey65", "TRUE" = "#c7372f"), name = paste0("|log2FC| >= ", log2fc_cutoff)) +
    geom_text_repel(data = eh, aes(label = "EHMT2"), size = 3.4, color = "black", min.segment.length = 0) +
    labs(x = paste(pair[[1]], "log2FC (NC / Control)"), y = paste(pair[[2]], "log2FC (NC / Control)")) +
    theme_bw(base_size = 11)
  save_plot(p, file.path(scatter_dir, paste0(pair[[1]], "_vs_", pair[[2]], "_log2fc_scatter.png")), 6.5, 5.5)
}

maplot_dir <- file.path(analysis_dir, "06_maplots")
for (comp in comparisons$comparison) {
  pdat <- de_all %>% filter(comparison == comp)
  eh <- pdat %>% filter(is_EHMT2_G9a)
  p <- ggplot(pdat, aes(x = log10(mean_tpm + pseudocount), y = log2FC)) +
    geom_hline(yintercept = c(-log2fc_cutoff, 0, log2fc_cutoff), color = c("grey80", "grey45", "grey80"), linewidth = c(0.3, 0.4, 0.3)) +
    geom_point(aes(color = abs(log2FC) >= log2fc_cutoff & pass_expression_filter), alpha = 0.6, size = 1.1) +
    scale_color_manual(values = c("FALSE" = "grey65", "TRUE" = "#c7372f"), name = paste0("|log2FC| >= ", log2fc_cutoff)) +
    geom_text_repel(data = eh, aes(label = "EHMT2"), size = 3.4, color = "black", min.segment.length = 0) +
    labs(x = "log10(mean TPM + 0.1)", y = "log2FC (NC / Control)", title = comp) +
    theme_bw(base_size = 11)
  save_plot(p, file.path(maplot_dir, paste0(comp, "_MA_plot.png")), 6.5, 5.5)
}

write_labeled_publication_maplots <- function(de_all, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  label_genes <- c("EHMT2", "MYC", "MAP1LC3B", "SQSTM1")
  comparison_titles <- c(
    "ComparisonA" = "CE81T2-NC vs CE81T2-Control",
    "ComparisonB" = "OE21-NC vs OE21-Control",
    "ComparisonC" = "KYSE270-NC vs KYSE270-Control"
  )
  label_rows <- de_all %>%
    filter(gene_name %in% label_genes) %>%
    dplyr::select(comparison, gene_id, gene_name, control_sample, nc_sample, control_tpm, nc_tpm, mean_tpm, log2FC, direction)
  readr::write_tsv(label_rows, file.path(out_dir, "labeled_gene_log2fc_values.tsv"))

  for (comp in comparisons$comparison) {
    pdat <- de_all %>% filter(comparison == comp)
    labels <- pdat %>%
      filter(gene_name %in% label_genes) %>%
      mutate(label = gene_name)
    p <- ggplot(pdat, aes(x = log10(mean_tpm + pseudocount), y = log2FC)) +
      geom_hline(yintercept = c(-log2fc_cutoff, 0, log2fc_cutoff), color = c("grey80", "grey45", "grey80"), linewidth = c(0.3, 0.4, 0.3)) +
      geom_point(aes(color = abs(log2FC) >= log2fc_cutoff & pass_expression_filter), alpha = 0.6, size = 1.1) +
      scale_color_manual(
        values = c("FALSE" = "grey65", "TRUE" = "#c7372f"),
        labels = c("FALSE" = "No", "TRUE" = "Yes"),
        name = paste0("|log2FC| >= ", log2fc_cutoff)
      ) +
      geom_point(data = labels, aes(x = log10(mean_tpm + pseudocount), y = log2FC), color = "black", fill = "#ffd166", shape = 21, size = 2.4, stroke = 0.45) +
      geom_text_repel(
        data = labels,
        aes(label = label),
        size = 3.6,
        color = "black",
        min.segment.length = 0,
        box.padding = 0.4,
        point.padding = 0.35,
        max.overlaps = Inf
      ) +
      labs(
        x = "log10(mean TPM + 0.1)",
        y = "log2FC (NC / Control)",
        title = unname(comparison_titles[comp])
      ) +
      theme_bw(base_size = 11) +
      theme(
        legend.position = c(0.985, 0.985),
        legend.justification = c(1, 1),
        legend.direction = "vertical",
        legend.background = element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.85), color = "grey70", linewidth = 0.25),
        legend.key = element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0)),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 15, face = "bold")
      )
    ggsave(
      file.path(out_dir, paste0(comp, "_MA_plot_labeled_EHMT2_MYC_MAP1LC3B_SQSTM1.png")),
      p,
      width = 7.2,
      height = 5.8,
      units = "in",
      dpi = 300
    )
  }
}

write_labeled_publication_maplots(
  de_all,
  file.path(maplot_dir, "labeled_EHMT2_MYC_MAP1LC3B_SQSTM1")
)

plot_heatmap <- function(ids, filename, title, max_rows = 100) {
  ids <- unique(ids)
  ids <- ids[ids %in% rownames(log_expr)]
  if (!length(ids)) {
    writeLines("No genes available for this heatmap.", sub("\\.png$", "_SKIPPED.txt", filename))
    return(invisible(FALSE))
  }
  if (length(ids) > max_rows) {
    fc_rank <- de_all %>%
      filter(gene_id_base %in% ids) %>%
      group_by(gene_id_base) %>%
      summarize(rank_abs_fc = max(abs(log2FC), na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(rank_abs_fc))
    keep <- head(fc_rank$gene_id_base, max_rows)
    if (any(ids %in% ehmt2_ids)) keep <- unique(c(ehmt2_ids[ehmt2_ids %in% ids], keep))
    ids <- head(keep, max_rows)
  }
  mat <- log_expr[ids, samples, drop = FALSE]
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[is.na(mat_scaled)] <- 0
  labels <- genes$gene_label[match(rownames(mat_scaled), genes$gene_id_base)]
  if (nrow(mat_scaled) > 45) labels[!(rownames(mat_scaled) %in% ehmt2_ids)] <- ""
  rownames(mat_scaled) <- labels
  ann_col <- metadata %>%
    filter(sample %in% samples) %>%
    dplyr::select(sample, cell_line, condition) %>%
    tibble::column_to_rownames("sample")
  pheatmap::pheatmap(
    mat_scaled,
    filename = filename,
    width = 7,
    height = max(4, min(12, 1.8 + 0.12 * nrow(mat_scaled))),
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    labels_row = labels,
    annotation_col = ann_col,
    main = title
  )
  invisible(TRUE)
}

heatmap_dir <- file.path(analysis_dir, "05_heatmaps")
top_ids <- de_all %>%
  filter(pass_expression_filter) %>%
  group_by(gene_id_base) %>%
  summarize(max_abs_log2FC = max(abs(log2FC), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_abs_log2FC)) %>%
  slice_head(n = 100) %>%
  pull(gene_id_base)
if (length(intersect(ehmt2_ids, genes$gene_id_base))) top_ids <- unique(c(intersect(ehmt2_ids, genes$gene_id_base), top_ids))
plot_heatmap(top_ids, file.path(heatmap_dir, "top_DE_genes_heatmap.png"), "Top exploratory DE genes")
for (nm in names(gene_sets)) {
  plot_heatmap(gene_sets[[nm]], file.path(heatmap_dir, paste0(nm, "_heatmap.png")), nm)
}
plot_heatmap(up_overlap_ids, file.path(heatmap_dir, "reciprocal_overlap_up_in_NC_heatmap.png"), "Reciprocal up in NC")
plot_heatmap(down_overlap_ids, file.path(heatmap_dir, "reciprocal_overlap_down_in_NC_heatmap.png"), "Reciprocal down in NC")

package_available <- function(pkg) requireNamespace(pkg, quietly = TRUE)
installed_status <- tibble::tibble(
  package = c("BiocManager", "clusterProfiler", "ReactomePA", "reactome.db", "org.Hs.eg.db", "msigdbr", "enrichplot", "fgsea", "pheatmap", "ComplexHeatmap", "ggplot2", "ggrepel", "dplyr", "readr", "tidyr", "stringr"),
  available = vapply(package, package_available, logical(1)),
  version = vapply(package, function(p) if (package_available(p)) as.character(utils::packageVersion(p)) else NA_character_, character(1))
)
readr::write_tsv(installed_status, file.path(analysis_dir, "versions", "package_versions.tsv"))

run_enrichment <- package_available("clusterProfiler") && package_available("org.Hs.eg.db")
if (run_enrichment) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })
  universe_map <- suppressMessages(clusterProfiler::bitr(unique(genes$gene_id_base), fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db))
  universe_entrez <- unique(universe_map$ENTREZID)
  comparison_universe_maps <- setNames(vector("list", nrow(comparisons)), comparisons$comparison)
  for (comp in comparisons$comparison) {
    expressed_ids <- de_all %>%
      filter(comparison == comp, pass_expression_filter) %>%
      pull(gene_id_base) %>%
      unique()
    comparison_universe_maps[[comp]] <- suppressMessages(
      clusterProfiler::bitr(expressed_ids, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
    )
  }
} else {
  write_skip(file.path(analysis_dir, "03_pathway_enrichment"), "Skipped pathway enrichment because clusterProfiler or org.Hs.eg.db is unavailable.")
  write_skip(file.path(analysis_dir, "04_gsea"), "Skipped GSEA because clusterProfiler or org.Hs.eg.db is unavailable.")
}

bar_enrich <- function(df, xcol, title) {
  top <- df %>%
    arrange(.data[[xcol]]) %>%
    slice_head(n = 20) %>%
    mutate(
      neg_log10_value = -log10(.data[[xcol]]),
      Description = factor(Description, levels = rev(Description))
    )
  ggplot(top, aes(x = neg_log10_value, y = Description)) +
    geom_col(fill = "#33658a") +
    labs(x = paste0("-log10(", xcol, ")"), y = NULL, title = title) +
    theme_bw(base_size = 10)
}

save_enrich_result <- function(res, out_dir, prefix) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  raw_dir <- file.path(out_dir, "raw_p_lt_0.05")
  padj_dir <- file.path(out_dir, "padj_lt_0.05")
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(padj_dir, showWarnings = FALSE, recursive = TRUE)
  df <- as.data.frame(res)
  readr::write_tsv(df, file.path(out_dir, paste0(prefix, "_all.tsv")))
  if (!nrow(df)) {
    writeLines("No enrichment terms returned.", file.path(out_dir, paste0(prefix, "_NO_TERMS.txt")))
    return(invisible(FALSE))
  }
  raw <- df %>% filter(pvalue < 0.05)
  padj <- df %>% filter(p.adjust < 0.05)
  readr::write_tsv(raw, file.path(raw_dir, paste0(prefix, "_raw_p_lt_0.05.tsv")))
  readr::write_tsv(padj, file.path(padj_dir, paste0(prefix, "_padj_lt_0.05.tsv")))
  if (nrow(raw)) save_plot(bar_enrich(raw, "pvalue", paste(prefix, "raw p < 0.05")), file.path(raw_dir, paste0(prefix, "_barplot_raw_p.png")), 7, 5)
  if (nrow(padj)) save_plot(bar_enrich(padj, "p.adjust", paste(prefix, "padj < 0.05")), file.path(padj_dir, paste0(prefix, "_barplot_padj.png")), 7, 5)
  invisible(TRUE)
}

bar_gsea <- function(df, title) {
  top <- df %>%
    arrange(pvalue) %>%
    slice_head(n = 20) %>%
    mutate(
      Description_wrapped = stringr::str_wrap(Description, width = 48),
      Description_wrapped = factor(Description_wrapped, levels = Description_wrapped)
    )
  ggplot(top, aes(x = NES, y = Description_wrapped, fill = NES > 0)) +
    geom_vline(xintercept = 0, color = "grey45", linewidth = 0.35) +
    geom_col() +
    scale_fill_manual(values = c("FALSE" = "#5b8cbe", "TRUE" = "#c7372f"), labels = c("Negative", "Positive"), name = "NES direction") +
    labs(x = "NES", y = NULL, title = title) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 11),
      plot.margin = margin(8, 18, 8, 8)
    )
}

save_gsea_result <- function(res, out_dir, prefix) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  raw_dir <- file.path(out_dir, "raw_p_lt_0.05")
  padj_dir <- file.path(out_dir, "padj_lt_0.05")
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(padj_dir, showWarnings = FALSE, recursive = TRUE)
  df <- as.data.frame(res)
  readr::write_tsv(df, file.path(out_dir, paste0(prefix, "_all.tsv")))
  if (!nrow(df)) {
    writeLines("No GSEA terms returned.", file.path(out_dir, paste0(prefix, "_NO_TERMS.txt")))
    return(invisible(FALSE))
  }
  raw <- df %>% filter(pvalue < 0.05)
  padj <- df %>% filter(p.adjust < 0.05)
  readr::write_tsv(raw, file.path(raw_dir, paste0(prefix, "_raw_p_lt_0.05.tsv")))
  readr::write_tsv(padj, file.path(padj_dir, paste0(prefix, "_padj_lt_0.05.tsv")))
  if (nrow(raw)) save_plot(
    bar_gsea(raw, paste(prefix, "raw p < 0.05")),
    file.path(raw_dir, paste0(prefix, "_NES_barplot_raw_p.png")),
    13,
    max(6, 2.5 + 0.28 * min(20, nrow(raw)))
  )
  if (nrow(padj)) save_plot(
    bar_gsea(padj, paste(prefix, "padj < 0.05")),
    file.path(padj_dir, paste0(prefix, "_NES_barplot_padj.png")),
    13,
    max(6, 2.5 + 0.28 * min(20, nrow(padj)))
  )
  invisible(TRUE)
}

if (run_enrichment) {
  ora_root <- file.path(analysis_dir, "03_pathway_enrichment")
  for (set_name in names(gene_sets)) {
    comp_name <- sub("_(up|down)$", "", set_name)
    comparison_universe_entrez <- unique(comparison_universe_maps[[comp_name]]$ENTREZID)
    ids <- gene_sets[[set_name]]
    mapped <- suppressMessages(clusterProfiler::bitr(ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
    entrez <- unique(mapped$ENTREZID)
    if (length(entrez) < 5) {
      write_skip(file.path(ora_root, set_name), paste("Skipped", set_name, "because fewer than five genes mapped to Entrez IDs."))
      next
    }
    kegg_res <- tryCatch(
      clusterProfiler::enrichKEGG(gene = entrez, organism = "hsa", universe = comparison_universe_entrez, pvalueCutoff = 1, qvalueCutoff = 1),
      error = function(e) e
    )
    if (inherits(kegg_res, "error")) {
      write_skip(file.path(ora_root, set_name, "KEGG"), paste("KEGG ORA failed:", conditionMessage(kegg_res)))
    } else {
      save_enrich_result(kegg_res, file.path(ora_root, set_name, "KEGG"), paste0(set_name, "_KEGG_ORA"))
    }
    if (package_available("ReactomePA")) {
      reactome_res <- tryCatch(
        ReactomePA::enrichPathway(gene = entrez, organism = "human", universe = comparison_universe_entrez, pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE),
        error = function(e) e
      )
      if (inherits(reactome_res, "error")) {
        write_skip(file.path(ora_root, set_name, "Reactome"), paste("Reactome ORA failed:", conditionMessage(reactome_res)))
      } else {
        save_enrich_result(reactome_res, file.path(ora_root, set_name, "Reactome"), paste0(set_name, "_Reactome_ORA"))
      }
    } else {
      write_skip(file.path(ora_root, set_name, "Reactome"), "Skipped Reactome ORA because ReactomePA is not installed.")
    }
  }

  gsea_root <- file.path(analysis_dir, "04_gsea")
  for (comp in comparisons$comparison) {
    ranked_df <- de_all %>%
      filter(comparison == comp, pass_expression_filter) %>%
      dplyr::select(gene_id_base, log2FC) %>%
      inner_join(universe_map %>% dplyr::select(gene_id_base = ENSEMBL, ENTREZID), by = "gene_id_base") %>%
      group_by(ENTREZID) %>%
      summarize(log2FC = log2FC[which.max(abs(log2FC))], .groups = "drop") %>%
      arrange(desc(log2FC))
    rank_vec <- ranked_df$log2FC
    names(rank_vec) <- ranked_df$ENTREZID
    kegg_gsea <- tryCatch(
      suppressMessages(clusterProfiler::gseKEGG(geneList = rank_vec, organism = "hsa", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, verbose = FALSE)),
      error = function(e) e
    )
    if (inherits(kegg_gsea, "error")) {
      write_skip(file.path(gsea_root, comp, "KEGG"), paste("KEGG GSEA failed:", conditionMessage(kegg_gsea)))
    } else {
      save_gsea_result(kegg_gsea, file.path(gsea_root, comp, "KEGG"), paste0(comp, "_KEGG_GSEA"))
    }
    if (package_available("ReactomePA")) {
      reactome_gsea <- tryCatch(
        suppressMessages(ReactomePA::gsePathway(geneList = rank_vec, organism = "human", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, verbose = FALSE)),
        error = function(e) e
      )
      if (inherits(reactome_gsea, "error")) {
        write_skip(file.path(gsea_root, comp, "Reactome"), paste("Reactome GSEA failed:", conditionMessage(reactome_gsea)))
      } else {
        save_gsea_result(reactome_gsea, file.path(gsea_root, comp, "Reactome"), paste0(comp, "_Reactome_GSEA"))
      }
    } else {
      write_skip(file.path(gsea_root, comp, "Reactome"), "Skipped Reactome GSEA because ReactomePA is not installed.")
    }
  }
}

get_hallmark_sets <- function() {
  if (!package_available("msigdbr")) return(NULL)
  res <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "H"), error = function(e2) NULL)
  )
  if (is.null(res) || !all(c("gs_name", "gene_symbol") %in% names(res))) return(NULL)
  res %>% dplyr::select(gs_name, gene_symbol) %>% distinct()
}

format_hallmark_label <- function(x) {
  out <- x %>%
    stringr::str_remove("^HALLMARK_") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_title()
  replacements <- c(
    "Dna" = "DNA", "E2f" = "E2F", "G2m" = "G2M", "Myc" = "MYC",
    "P53" = "P53", "Pi3k" = "PI3K", "Akt" = "AKT", "Mtor" = "MTOR",
    "Mtorc1" = "MTORC1", "Kras" = "KRAS", "Wnt" = "WNT", "Tgf" = "TGF",
    "Il6" = "IL6", "Jak" = "JAK", "Stat3" = "STAT3", "Il2" = "IL2",
    "Stat5" = "STAT5", "Tnfa" = "TNFA", "Nfkb" = "NFKB", "Dn" = "DN",
    "Up" = "UP", "Atr" = "ATR", "Atf4" = "ATF4"
  )
  for (from in names(replacements)) {
    out <- stringr::str_replace_all(out, paste0("\\b", from, "\\b"), replacements[[from]])
  }
  out
}

clean_sample_label <- function(x) stringr::str_remove(x, "^\\d+-")

score_gene_sets <- function(set_tbl, selected_sets, out_dir, prefix) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  expr_by_symbol <- log_expr
  rownames(expr_by_symbol) <- make.unique(genes$gene_name[match(rownames(log_expr), genes$gene_id_base)])
  score_rows <- list()
  used_sets <- intersect(selected_sets, unique(set_tbl$gs_name))
  for (set_name in used_sets) {
    members <- unique(set_tbl$gene_symbol[set_tbl$gs_name == set_name])
    present <- intersect(members, rownames(expr_by_symbol))
    if (length(present) < 3) next
    mat <- expr_by_symbol[present, samples, drop = FALSE]
    z <- t(scale(t(mat)))
    z[is.na(z)] <- 0
    score <- colMeans(z)
    score_rows[[set_name]] <- tibble::tibble(panel = set_name, sample = names(score), score = as.numeric(score), genes_used = length(present))
  }
  scores <- bind_rows(score_rows) %>%
    left_join(metadata %>% dplyr::select(sample, cell_line, condition), by = "sample") %>%
    mutate(
      panel_label = format_hallmark_label(panel),
      sample_label = factor(clean_sample_label(sample), levels = clean_sample_label(samples))
    )
  readr::write_tsv(scores, file.path(out_dir, paste0(prefix, "_sample_scores.tsv")))
  if (!nrow(scores)) {
    writeLines("No selected gene sets had at least three genes present.", file.path(out_dir, paste0(prefix, "_NO_SCORES.txt")))
    return(invisible(FALSE))
  }
  delta <- scores %>%
    dplyr::select(panel, panel_label, cell_line, condition, score) %>%
    tidyr::pivot_wider(names_from = condition, values_from = score) %>%
    mutate(delta_NC_minus_Control = NC - Con)
  readr::write_tsv(delta, file.path(out_dir, paste0(prefix, "_delta_scores.tsv")))
  p1 <- ggplot(scores, aes(x = sample_label, y = score, fill = condition)) +
    geom_col() +
    facet_wrap(~panel_label, scales = "free_y") +
    labs(x = NULL, y = "Mean row-z score") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p1, file.path(out_dir, paste0(prefix, "_sample_score_plot.png")), 10, 7)
  p2 <- ggplot(delta, aes(x = cell_line, y = delta_NC_minus_Control, fill = delta_NC_minus_Control > 0)) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.3) +
    geom_col() +
    facet_wrap(~panel_label, scales = "free_y") +
    scale_fill_manual(values = c("FALSE" = "#5b8cbe", "TRUE" = "#c7372f"), guide = "none") +
    labs(x = NULL, y = "Delta score (NC - Control)") +
    theme_bw(base_size = 9)
  save_plot(p2, file.path(out_dir, paste0(prefix, "_delta_plot.png")), 10, 7)
  score_mat <- scores %>%
    dplyr::select(panel, sample, score) %>%
    tidyr::pivot_wider(names_from = sample, values_from = score) %>%
    tibble::column_to_rownames("panel") %>%
    as.matrix()
  pheatmap::pheatmap(
    score_mat,
    filename = file.path(out_dir, paste0(prefix, "_score_heatmap.png")),
    width = 7,
    height = max(4, 0.3 * nrow(score_mat) + 2),
    cluster_rows = nrow(score_mat) > 1,
    cluster_cols = ncol(score_mat) > 1
  )
  invisible(TRUE)
}

expression_by_symbol <- function() {
  as_tibble(log_expr, rownames = "gene_id_base") %>%
    left_join(genes %>% dplyr::select(gene_id_base, gene_name), by = "gene_id_base") %>%
    filter(!is.na(gene_name), gene_name != "") %>%
    dplyr::select(gene_name, all_of(samples)) %>%
    group_by(gene_name) %>%
    summarize(across(all_of(samples), mean), .groups = "drop") %>%
    tibble::column_to_rownames("gene_name") %>%
    as.matrix()
}

plot_gene_set_heatmap <- function(gene_symbols, out_file, title, show_gene_names = FALSE, label_fontsize = 7, height_per_gene = 0.08) {
  expr_symbol <- expression_by_symbol()
  present <- intersect(unique(gene_symbols), rownames(expr_symbol))
  if (length(present) < 2) {
    writeLines("Fewer than two genes from this set are present in the RNA-seq matrix.", sub("\\.png$", "_SKIPPED.txt", out_file))
    return(invisible(FALSE))
  }
  mat <- expr_symbol[present, samples, drop = FALSE]
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[is.na(mat_scaled)] <- 0
  row_order <- order(apply(mat_scaled, 1, function(x) max(abs(x), na.rm = TRUE)), decreasing = TRUE)
  mat_scaled <- mat_scaled[row_order, , drop = FALSE]
  ann_col <- metadata %>%
    filter(sample %in% samples) %>%
    mutate(sample_label = clean_sample_label(sample)) %>%
    dplyr::select(sample, cell_line, condition) %>%
    tibble::column_to_rownames("sample")
  pheatmap::pheatmap(
    mat_scaled,
    filename = out_file,
    width = 7,
    height = max(6, min(28, 2.5 + height_per_gene * nrow(mat_scaled))),
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_rownames = show_gene_names || nrow(mat_scaled) <= 60,
    fontsize_row = label_fontsize,
    labels_col = clean_sample_label(samples),
    annotation_col = ann_col,
    main = title
  )
  invisible(TRUE)
}

plot_gene_set_delta_heatmap <- function(gene_symbols, out_file, title, show_gene_names = TRUE, label_fontsize = 5.5, height_per_gene = 0.18) {
  expr_symbol <- expression_by_symbol()
  present <- intersect(unique(gene_symbols), rownames(expr_symbol))
  if (length(present) < 2) {
    writeLines("Fewer than two genes from this set are present in the RNA-seq matrix.", sub("\\.png$", "_SKIPPED.txt", out_file))
    return(invisible(FALSE))
  }
  mat <- expr_symbol[present, samples, drop = FALSE]
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[is.na(mat_scaled)] <- 0
  delta_mat <- sapply(seq_len(nrow(comparisons)), function(i) {
    comp <- comparisons[i, ]
    mat_scaled[, comp$nc] - mat_scaled[, comp$control]
  })
  colnames(delta_mat) <- comparisons$cell_line
  row_order <- order(apply(delta_mat, 1, function(x) max(abs(x), na.rm = TRUE)), decreasing = TRUE)
  delta_mat <- delta_mat[row_order, , drop = FALSE]
  pheatmap::pheatmap(
    delta_mat,
    filename = out_file,
    width = 5.8,
    height = max(6, min(28, 2.5 + height_per_gene * nrow(delta_mat))),
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_rownames = show_gene_names,
    fontsize_row = label_fontsize,
    main = title
  )
  invisible(TRUE)
}

write_gene_set_gene_tables <- function(set_tbl, selected_sets, out_dir, prefix) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  expr_symbol <- expression_by_symbol()
  for (set_name in selected_sets) {
    all_genes <- set_tbl %>%
      filter(gs_name == set_name) %>%
      transmute(gene_symbol = gene_symbol) %>%
      distinct() %>%
      arrange(gene_symbol)
    present_genes <- all_genes %>%
      mutate(present_in_matrix = gene_symbol %in% rownames(expr_symbol)) %>%
      filter(present_in_matrix)
    file_stub <- stringr::str_remove(set_name, "^HALLMARK_") %>% stringr::str_to_lower()
    readr::write_tsv(all_genes, file.path(out_dir, paste0(prefix, "_", file_stub, "_all_msigdb_genes.tsv")))
    readr::write_tsv(present_genes, file.path(out_dir, paste0(prefix, "_", file_stub, "_present_genes.tsv")))
    plot_gene_set_heatmap(
      present_genes$gene_symbol,
      file.path(out_dir, paste0(prefix, "_", file_stub, "_row_z_heatmap.png")),
      format_hallmark_label(set_name),
      show_gene_names = FALSE
    )
    top_present <- rownames(expression_by_symbol())[rownames(expression_by_symbol()) %in% present_genes$gene_symbol]
    if (length(top_present) > 60) {
      mat <- expression_by_symbol()[top_present, samples, drop = FALSE]
      z <- t(scale(t(mat)))
      z[is.na(z)] <- 0
      top_present <- names(sort(apply(z, 1, function(x) max(abs(x), na.rm = TRUE)), decreasing = TRUE))[seq_len(60)]
    }
    plot_gene_set_heatmap(
      top_present,
      file.path(out_dir, paste0(prefix, "_", file_stub, "_top60_labeled_row_z_heatmap.png")),
      paste(format_hallmark_label(set_name), "top altered genes"),
      show_gene_names = TRUE,
      label_fontsize = 5.5,
      height_per_gene = 0.18
    )
    plot_gene_set_delta_heatmap(
      top_present,
      file.path(out_dir, paste0(prefix, "_", file_stub, "_top60_delta_NC_minus_Control_heatmap.png")),
      paste(format_hallmark_label(set_name), "top altered genes delta"),
      show_gene_names = TRUE,
      label_fontsize = 5.5,
      height_per_gene = 0.18
    )
  }
}

write_myc_targets_readme <- function(out_dir) {
  readme <- c(
    "# MYC Targets Analysis",
    "",
    "## Gene Set Source",
    "",
    "MYC Targets V1 and MYC Targets V2 come from MSigDB Hallmark gene sets:",
    "",
    "- `HALLMARK_MYC_TARGETS_V1`",
    "- `HALLMARK_MYC_TARGETS_V2`",
    "",
    "The `*_all_msigdb_genes.tsv` files list all MSigDB genes in each set.",
    "The `*_present_genes.tsv` files list the subset found in the protein-coding RNA-seq TPM matrix.",
    "",
    "## Heatmap Values",
    "",
    "Heatmaps use `log2(TPM + 0.1)` values.",
    "Each gene is row-z scaled across the six samples, so colors show relative high or low expression of that gene across samples, not absolute TPM.",
    "",
    "Sample order is fixed as:",
    "",
    "`81T-Con`, `81T-NC`, `OE21-Con`, `OE21-NC`, `KYSE-Con`, `KYSE-NC`.",
    "",
    "## Top Altered Gene Criterion",
    "",
    "For the labeled `top60` heatmaps, genes are ranked by their strongest relative alteration across the six samples:",
    "",
    "`top_alteration_score = max(abs(row_z_score))`",
    "",
    "The top 60 genes with the largest `top_alteration_score` are shown with gene symbols.",
    "This criterion selects genes with the most extreme sample-to-sample variation within the MYC target set.",
    "It is not based on p-values, FDR, or replicate-aware differential expression because this experiment has one sample per condition.",
    "",
    "The full MYC V1/V2 heatmaps include all present genes from each MYC target set.",
    "",
    "## Delta Heatmaps",
    "",
    "The `*_top60_delta_NC_minus_Control_heatmap.png` files use the same top-60 genes.",
    "Each value is the difference in row-z expression score between paired NC and control samples within a cell line:",
    "",
    "`delta_alteration_score = row_z(NC) - row_z(Control)`",
    "",
    "The x-axis is reduced to the three cell lines:",
    "",
    "`81T`, `OE21`, `KYSE`.",
    "",
    "Positive values indicate higher relative expression in NC than control for that gene within the cell line; negative values indicate lower relative expression in NC."
  )
  writeLines(readme, file.path(out_dir, "README_MYC_targets_analysis.md"))
}

write_stat3_related_analysis <- function(hallmark, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  all_msig <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens"), error = function(e) NULL)
  if (!is.null(all_msig)) {
    name_cols <- intersect(c("gs_collection", "gs_subcollection", "gs_collection_name", "gs_name", "gs_description"), names(all_msig))
    stat3_sets <- all_msig %>%
      distinct(across(all_of(name_cols))) %>%
      filter(if_any(any_of(c("gs_name", "gs_description")), ~ grepl("STAT3|JAK", .x, ignore.case = TRUE))) %>%
      arrange(gs_name)
    readr::write_tsv(stat3_sets, file.path(out_dir, "msigdb_STAT3_JAK_related_sets.tsv"))
  }
  stat3_set <- "HALLMARK_IL6_JAK_STAT3_SIGNALING"
  stat3_genes <- hallmark %>% filter(gs_name == stat3_set) %>% pull(gene_symbol) %>% unique()
  write_gene_set_gene_tables(hallmark, stat3_set, out_dir, "stat3_downstream")
  score_gene_sets(hallmark, stat3_set, out_dir, "stat3_downstream")
  plot_gene_set_heatmap(
    stat3_genes,
    file.path(out_dir, "stat3_downstream_IL6_JAK_STAT3_row_z_heatmap.png"),
    "IL6 JAK STAT3 Signaling",
    show_gene_names = TRUE
  )
}

write_selected_hallmark_delta_figure <- function(out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  delta_file <- file.path(analysis_dir, "07_cancer_signature_analysis", "hallmark_cancer_delta_scores.tsv")
  if (!file.exists(delta_file)) {
    writeLines("Missing hallmark_cancer_delta_scores.tsv; run hallmark cancer scoring first.", file.path(out_dir, "SKIPPED.txt"))
    return(invisible(FALSE))
  }
  selected_sets <- c(
    "HALLMARK_APOPTOSIS",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_HEDGEHOG_SIGNALING",
    "HALLMARK_HYPOXIA",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MYC_TARGETS_V2",
    "HALLMARK_NOTCH_SIGNALING",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
  )
  cell_line_labels <- c("81T" = "CE81T2", "OE21" = "OE21", "KYSE" = "KYSE270")
  selected_labels <- setNames(format_hallmark_label(selected_sets), selected_sets)
  delta <- read.delim(delta_file, check.names = FALSE) %>%
    filter(panel %in% selected_sets) %>%
    mutate(
      panel_label = unname(selected_labels[panel]),
      cell_line_label = factor(unname(cell_line_labels[cell_line]), levels = c("CE81T2", "OE21", "KYSE270")),
      regulation = ifelse(delta_NC_minus_Control >= 0, "Upregulated", "Downregulated")
    ) %>%
    filter(!is.na(cell_line_label))
  facet_levels <- c(
    paste("Upregulated", sort(unique(delta$panel_label[delta$regulation == "Upregulated"])), sep = " - "),
    paste("Downregulated", sort(unique(delta$panel_label[delta$regulation == "Downregulated"])), sep = " - ")
  )
  plot_data <- delta %>%
    mutate(
      facet_label = factor(paste(regulation, panel_label, sep = " - "), levels = facet_levels),
      regulation = factor(regulation, levels = c("Upregulated", "Downregulated"))
    )
  readr::write_tsv(
    plot_data %>% dplyr::select(panel, panel_label, cell_line, cell_line_label, Con, NC, delta_NC_minus_Control, regulation),
    file.path(out_dir, "selected_hallmark_cancer_delta_values.tsv")
  )
  p <- ggplot(plot_data, aes(x = cell_line_label, y = delta_NC_minus_Control, fill = cell_line_label)) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35) +
    geom_col(width = 0.72) +
    facet_wrap(~facet_label, ncol = 4, scales = "free_y") +
    scale_fill_manual(values = c("CE81T2" = "#6f9ed4", "OE21" = "#d95f02", "KYSE270" = "#1b9e77"), name = "Cell line") +
    labs(x = NULL, y = "Delta score (NC - Control)") +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 8, face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.major.x = element_blank()
    )
  ggsave(
    file.path(out_dir, "selected_hallmark_cancer_delta_up_down_600dpi.png"),
    p,
    width = 13,
    height = 13,
    units = "in",
    dpi = 600
  )
  p_large_yaxis <- p +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 9, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 9, angle = 35, hjust = 1),
      panel.grid.major.x = element_blank()
    )
  ggsave(
    file.path(out_dir, "selected_hallmark_cancer_delta_up_down_600dpi_large_yaxis.png"),
    p_large_yaxis,
    width = 13,
    height = 13,
    units = "in",
    dpi = 600
  )
  invisible(TRUE)
}

hallmark <- get_hallmark_sets()
cancer_sets <- c(
  "HALLMARK_P53_PATHWAY", "HALLMARK_HYPOXIA", "HALLMARK_ANGIOGENESIS", "HALLMARK_GLYCOLYSIS",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_APOPTOSIS",
  "HALLMARK_DNA_REPAIR", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_KRAS_SIGNALING_DN", "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_NOTCH_SIGNALING", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_HEDGEHOG_SIGNALING"
)
immune_sets <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_COMPLEMENT", "HALLMARK_ALLOGRAFT_REJECTION"
)

if (is.null(hallmark)) {
  write_skip(file.path(analysis_dir, "07_cancer_signature_analysis"), "Skipped Hallmark cancer signature scoring because msigdbr is not installed or could not return Hallmark gene sets.")
  write_skip(file.path(analysis_dir, "08_immune_gene_analysis"), "Skipped Hallmark immune/cytokine signature scoring because msigdbr is not installed or could not return Hallmark gene sets.")
  write_skip(file.path(analysis_dir, "09_emt_invasion_score"), "Skipped EMT scoring because msigdbr is not installed or could not return Hallmark gene sets.")
} else {
  score_gene_sets(hallmark, cancer_sets, file.path(analysis_dir, "07_cancer_signature_analysis"), "hallmark_cancer")
  write_selected_hallmark_delta_figure(file.path(analysis_dir, "07_cancer_signature_analysis", "selected_hallmark_delta_figure"))
  write_gene_set_gene_tables(
    hallmark,
    c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"),
    file.path(analysis_dir, "07_cancer_signature_analysis", "MYC_targets_analysis"),
    "MYC"
  )
  write_myc_targets_readme(file.path(analysis_dir, "07_cancer_signature_analysis", "MYC_targets_analysis"))
  write_stat3_related_analysis(
    hallmark,
    file.path(analysis_dir, "07_cancer_signature_analysis", "STAT3_related_analysis")
  )
  score_gene_sets(hallmark, immune_sets, file.path(analysis_dir, "08_immune_gene_analysis"), "hallmark_immune_cytokine")
  emt_dir <- file.path(analysis_dir, "09_emt_invasion_score")
  emt_genes <- hallmark %>% filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% pull(gene_symbol) %>% unique()
  expr_by_symbol <- log_expr
  rownames(expr_by_symbol) <- make.unique(genes$gene_name[match(rownames(log_expr), genes$gene_id_base)])
  present_emt <- intersect(emt_genes, rownames(expr_by_symbol))
  readr::write_tsv(tibble::tibble(gene_symbol = present_emt), file.path(emt_dir, "emt_genes_present.tsv"))
  if (length(present_emt) < 3) {
    write_skip(emt_dir, "Skipped EMT scoring because fewer than three EMT genes are present.")
  } else {
    emt_mat <- expr_by_symbol[present_emt, samples, drop = FALSE]
    emt_z <- t(scale(t(emt_mat)))
    emt_z[is.na(emt_z)] <- 0
    emt_score <- tibble::tibble(sample = samples, EMT_score = colMeans(emt_z)) %>%
      left_join(metadata %>% dplyr::select(sample, cell_line, condition), by = "sample")
    readr::write_tsv(emt_score, file.path(emt_dir, "emt_sample_scores.tsv"))
    emt_delta <- emt_score %>%
      dplyr::select(cell_line, condition, EMT_score) %>%
      tidyr::pivot_wider(names_from = condition, values_from = EMT_score) %>%
      mutate(delta_NC_minus_Control = NC - Con)
    readr::write_tsv(emt_delta, file.path(emt_dir, "emt_delta_scores.tsv"))
    p1 <- ggplot(emt_score, aes(x = sample, y = EMT_score, fill = condition)) +
      geom_col() +
      labs(x = NULL, y = "EMT mean row-z score") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(p1, file.path(emt_dir, "emt_individual_sample_score_plot.png"), 7, 5)
    p2 <- ggplot(emt_delta, aes(x = cell_line, y = delta_NC_minus_Control, fill = delta_NC_minus_Control > 0)) +
      geom_hline(yintercept = 0, color = "grey45", linewidth = 0.3) +
      geom_col() +
      scale_fill_manual(values = c("FALSE" = "#5b8cbe", "TRUE" = "#c7372f"), guide = "none") +
      labs(x = NULL, y = "EMT delta score (NC - Control)") +
      theme_bw(base_size = 10)
    save_plot(p2, file.path(emt_dir, "emt_comparison_delta_plot.png"), 6, 4)
    pheatmap::pheatmap(emt_z, filename = file.path(emt_dir, "full_EMT_gene_heatmap.png"), width = 7, height = 10, show_rownames = nrow(emt_z) <= 80)
    delta_mat <- sapply(seq_len(nrow(comparisons)), function(i) {
      comp <- comparisons[i, ]
      emt_mat[, comp$nc] - emt_mat[, comp$control]
    })
    colnames(delta_mat) <- comparisons$cell_line
    pheatmap::pheatmap(delta_mat, filename = file.path(emt_dir, "EMT_gene_delta_heatmap.png"), width = 5, height = 10, show_rownames = nrow(delta_mat) <= 80)
  }
}

direction_check <- isTRUE(all.equal(
  de_all$log2FC,
  log2((de_all$nc_tpm + pseudocount) / (de_all$control_tpm + pseudocount)),
  tolerance = 1e-12,
  check.attributes = FALSE
))
transcript_input_check <- basename(tpm_file) == "salmon.merged.gene_tpm.tsv" &&
  !grepl("transcript", normalizePath(tpm_file, winslash = "/", mustWork = FALSE), ignore.case = TRUE)

validation <- tibble::tibble(
  check = c(
    "all_expected_samples_in_tpm",
    "all_expected_samples_in_metadata",
    "protein_coding_filter_applied",
    "EHMT2_present_after_filter",
    "no_statistical_de_columns",
    "positive_log2fc_is_NC_over_control",
    "transcript_level_files_not_used"
  ),
  passed = c(
    all(samples %in% names(tpm_raw)),
    all(samples %in% metadata$sample),
    all(genes$gene_type == "protein_coding"),
    length(ehmt2_ids) > 0,
    !any(c("pvalue", "padj", "FDR", "p.adjust") %in% names(de_all)),
    direction_check,
    transcript_input_check
  ),
  detail = c(
    paste(samples, collapse = ", "),
    paste(metadata$sample, collapse = ", "),
    paste(nrow(genes), "protein-coding genes retained"),
    paste(genes$gene_name[match(ehmt2_ids, genes$gene_id_base)], collapse = ", "),
    "DE output uses fold-change only.",
    "log2FC = log2((NC TPM + 0.1) / (Control TPM + 0.1))",
    "Only salmon.merged.gene_tpm.tsv is read."
  )
)
readr::write_tsv(validation, file.path(analysis_dir, "validation", "validation_checks.tsv"))
if (any(!validation$passed)) stop("One or more validation checks failed. See validation_checks.tsv.")

summary_lines <- c(
  "# Autophagy RNA-seq Downstream Analysis Summary",
  "",
  "## Inputs",
  paste0("- TPM matrix: ", tpm_file),
  paste0("- Sample metadata: ", metadata_file),
  paste0("- GENCODE annotation: ", gtf_file),
  "",
  "## Samples and Comparisons",
  "- ComparisonA: 30-81T-NC / 29-81T-Con",
  "- ComparisonB: 32-OE21-NC / 31-OE21-Con",
  "- ComparisonC: 34-KYSE-NC / 33-KYSE-Con",
  "",
  "## Exploratory Fold-Change Settings",
  paste0("- TPM pseudocount: ", pseudocount),
  "- Direction: NC / Control, so positive log2FC means higher expression after Niclosamide treatment.",
  "- Protein-coding genes only, using local GENCODE v49 gene_type annotation.",
  paste0("- Expression filter: TPM >= 1 in at least one sample per comparison."),
  paste0("- Exploratory DE threshold: absolute log2FC >= ", log2fc_cutoff),
  "- No DESeq2, p-values, or FDR are used because each condition has one sample.",
  "",
  "## Output Folders",
  "- 01_differential_expression: fold-change tables, filtered lists, reciprocal overlaps.",
  "- 02_log2fc_scatter: comparison-pair log2FC scatter plots.",
  "- 03_pathway_enrichment: KEGG/Reactome ORA outputs or skip notes.",
  "- 04_gsea: KEGG/Reactome GSEA outputs or skip notes.",
  "- 05_heatmaps: top DE, comparison-specific, and reciprocal-overlap heatmaps.",
  "- 06_maplots: MA plots for each comparison.",
  "- 07_cancer_signature_analysis: Hallmark cancer panel scores or skip note.",
  "- 08_immune_gene_analysis: Hallmark immune/cytokine panel scores or skip note.",
  "- 09_emt_invasion_score: Hallmark EMT scores or skip note.",
  "",
  "## EHMT2/G9a Labeling",
  if (length(ehmt2_ids)) "- EHMT2 is present after filtering and is labeled on gene-level scatter/MA figures." else "- EHMT2/G9a was not present after filtering.",
  "",
  "## Package Versions",
  paste(capture.output(print(installed_status, n = Inf)), collapse = "\n")
)
writeLines(summary_lines, file.path(analysis_dir, "README_summary.md"))

cat("Autophagy downstream analysis completed.\n")
cat("Analysis directory:", analysis_dir, "\n")
