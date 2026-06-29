suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

project_dir <- Sys.getenv("AUTOPHAGY_PROJECT_DIR", unset = "F:/RNA_Seq/Autophagy")
analysis_dir <- file.path(project_dir, "Analysis")
de_file <- file.path(analysis_dir, "01_differential_expression", "all_comparisons_full_log2fc.tsv")
out_dir <- file.path(analysis_dir, "06_maplots", "highlighted_EHMT2_MYC_MAP1LC3B_SQSTM1_600px_panel")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

highlight_genes <- c("EHMT2", "MAP1LC3B", "MYC", "SQSTM1")
pseudocount <- 0.1
log2fc_cutoff <- 1.5

comparison_labels <- c(
  "ComparisonA" = "CE81T2",
  "ComparisonB" = "OE21",
  "ComparisonC" = "KYSE270"
)

gene_colors <- c(
  "EHMT2" = "#7B3294",
  "MAP1LC3B" = "#008837",
  "MYC" = "#E66101",
  "SQSTM1" = "#2166AC"
)

de_all <- readr::read_tsv(de_file, show_col_types = FALSE) %>%
  mutate(
    comparison_label = factor(unname(comparison_labels[comparison]), levels = unname(comparison_labels)),
    log10_mean_tpm = log10(mean_tpm + pseudocount),
    large_fc = abs(log2FC) >= log2fc_cutoff & pass_expression_filter,
    large_fc_label = factor(ifelse(large_fc, "Yes", "No"), levels = c("No", "Yes"))
  )

highlight_rows <- de_all %>%
  filter(gene_name %in% highlight_genes) %>%
  mutate(gene_name = factor(gene_name, levels = highlight_genes)) %>%
  select(
    comparison, comparison_label, gene_id, gene_name, control_sample, nc_sample,
    control_tpm, nc_tpm, mean_tpm, log2FC, direction, large_fc_label
  )

readr::write_tsv(highlight_rows, file.path(out_dir, "highlighted_gene_log2fc_values.tsv"))

p <- ggplot(de_all, aes(x = log10_mean_tpm, y = log2FC)) +
  geom_hline(yintercept = -log2fc_cutoff, color = "grey80", linewidth = 0.25) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_hline(yintercept = log2fc_cutoff, color = "grey80", linewidth = 0.25) +
  geom_point(
    aes(color = large_fc_label),
    alpha = 0.42,
    size = 0.35
  ) +
  geom_point(
    data = highlight_rows,
    aes(x = log10(mean_tpm + pseudocount), y = log2FC, fill = gene_name),
    inherit.aes = FALSE,
    shape = 21,
    color = "black",
    size = 2.3,
    stroke = 0.35
  ) +
  facet_wrap(~comparison_label, nrow = 1) +
  scale_color_manual(
    values = c("No" = "grey72", "Yes" = "#C7372F"),
    name = paste0("|log2FC| >= ", log2fc_cutoff)
  ) +
  scale_fill_manual(
    values = gene_colors,
    name = "Highlighted gene"
  ) +
  labs(
    x = "log10(mean TPM + 0.1)",
    y = "log2FC (NC / Control)"
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 2.0, alpha = 0.9)),
    fill = guide_legend(order = 2, override.aes = list(size = 3.0))
  ) +
  theme_bw(base_size = 8) +
  theme(
    strip.background = element_rect(fill = "grey92", color = "grey55", linewidth = 0.25),
    strip.text = element_text(size = 9, face = "bold"),
    axis.title = element_text(size = 8.5, face = "bold"),
    axis.text = element_text(size = 7.5, color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing.x = grid::unit(4, "pt"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7.5),
    legend.key.size = grid::unit(8, "pt"),
    legend.spacing.y = grid::unit(1, "pt"),
    plot.margin = margin(4, 5, 3, 4)
  )

ggsave(
  file.path(out_dir, "MA_plot_highlighted_genes_600px_panel.png"),
  p,
  width = 6.25,
  height = 2.5,
  units = "in",
  dpi = 600
)

ggsave(
  file.path(out_dir, "MA_plot_highlighted_genes_600px_panel.pdf"),
  p,
  width = 6.25,
  height = 2.5,
  units = "in",
  device = cairo_pdf
)

cat("Wrote highlighted-gene 600 px MA panel to:", out_dir, "\n")
