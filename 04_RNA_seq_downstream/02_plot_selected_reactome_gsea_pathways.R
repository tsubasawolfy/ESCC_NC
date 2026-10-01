library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
library(stringr)

project_dir <- Sys.getenv("AUTOPHAGY_PROJECT_DIR", unset = "F:/RNA_Seq/Autophagy")
analysis_dir <- file.path(project_dir, "Analysis")
input_file <- file.path(analysis_dir, "04_gsea", "selected_cell_line_specific_reactome_gsea_pathways.csv")
out_dir <- file.path(analysis_dir, "04_gsea", "selected_cell_line_specific_reactome_gsea_pathways")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

selected <- readr::read_csv(input_file, show_col_types = FALSE)

cell_line_labels <- c(
  "81T" = "CE81T2",
  "CE81T2" = "CE81T2",
  "KYSE" = "KYSE270",
  "KYSE270" = "KYSE270",
  "OE21" = "OE21"
)

format_pathway <- function(x) {
  x %>%
    str_replace_all("rRNA", "rRNA") %>%
    str_replace_all("tRNA", "tRNA") %>%
    str_wrap(width = 34)
}

plot_data <- selected %>%
  mutate(
    cell_line_label = recode(cell_line, !!!cell_line_labels, .default = cell_line),
    cell_line_label = factor(cell_line_label, levels = c("CE81T2", "OE21", "KYSE270")),
    pathway_label = format_pathway(reactome_pathway),
    direction_label = ifelse(NES >= 0, "Up in NC", "Down in NC"),
    neg_log10_padj = -log10(p_adjust),
    pathway_order = reorder(pathway_label, NES)
  ) %>%
  arrange(cell_line_label, NES)

readr::write_tsv(plot_data, file.path(out_dir, "selected_cell_line_specific_reactome_gsea_plot_data.tsv"))

cell_line_colors <- c(
  "CE81T2" = "#4C78A8",
  "OE21" = "#D95F02",
  "KYSE270" = "#1B9E77"
)

p <- ggplot(plot_data, aes(x = NES, y = pathway_order)) +
  geom_vline(xintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_segment(
    aes(x = 0, xend = NES, yend = pathway_order, color = cell_line_label),
    linewidth = 4.8,
    lineend = "round"
  ) +
  geom_point(
    aes(size = neg_log10_padj, fill = direction_label),
    shape = 21,
    color = "black",
    stroke = 0.35
  ) +
  facet_grid(cell_line_label ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = cell_line_colors, guide = "none") +
  scale_fill_manual(
    values = c("Up in NC" = "#B2182B", "Down in NC" = "#2166AC"),
    name = "GSEA direction"
  ) +
  scale_size_continuous(
    range = c(2.8, 6.0),
    name = expression(-log[10]("FDR"))
  ) +
  scale_x_continuous(
    limits = c(-3.1, 2.8),
    breaks = seq(-3, 3, 1),
    expand = expansion(mult = c(0.03, 0.04))
  ) +
  labs(
    x = "Reactome GSEA normalized enrichment score (NES)",
    y = NULL,
    title = "Selected cell-line-specific Reactome GSEA pathways",
    subtitle = "Niclosamide-treated samples compared with matched controls"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 10.5),
    strip.background = element_rect(fill = "grey92", color = "grey55"),
    strip.text.y = element_text(size = 12, face = "bold", angle = 0),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

ggsave(
  file.path(out_dir, "selected_cell_line_specific_reactome_gsea_NES_plot.png"),
  p,
  width = 10,
  height = 7.2,
  units = "in",
  dpi = 600
)

ggsave(
  file.path(out_dir, "selected_cell_line_specific_reactome_gsea_NES_plot.pdf"),
  p,
  width = 10,
  height = 7.2,
  units = "in",
  device = cairo_pdf
)

cat("Wrote selected Reactome GSEA plot to:", out_dir, "\n")

