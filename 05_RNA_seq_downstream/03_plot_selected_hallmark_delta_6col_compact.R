suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(patchwork)
})

project_dir <- Sys.getenv("AUTOPHAGY_PROJECT_DIR", unset = "F:/RNA_Seq/Autophagy")
analysis_dir <- file.path(project_dir, "Analysis")
input_file <- file.path(
  analysis_dir,
  "07_cancer_signature_analysis",
  "selected_hallmark_delta_figure",
  "selected_hallmark_cancer_delta_values.tsv"
)
output_file <- file.path(
  analysis_dir,
  "07_cancer_signature_analysis",
  "selected_hallmark_delta_figure",
  "selected_hallmark_cancer_delta_up_down_600dpi_6col_compact.png"
)

plot_data <- readr::read_tsv(input_file, show_col_types = FALSE) %>%
  mutate(
    cell_line_label = factor(cell_line_label, levels = c("CE81T2", "OE21", "KYSE270")),
    regulation = factor(regulation, levels = c("Upregulated", "Downregulated"))
  )

plot_data <- plot_data %>%
  group_by(regulation) %>%
  mutate(panel_label = factor(panel_label, levels = sort(unique(panel_label)))) %>%
  ungroup()

plot_one_regulation <- function(df, title_text, show_legend = FALSE) {
  ggplot(df, aes(x = cell_line_label, y = delta_NC_minus_Control, fill = cell_line_label)) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35) +
    geom_col(width = 0.72) +
    facet_wrap(~panel_label, ncol = 6, scales = "free_y") +
    scale_fill_manual(
      values = c("CE81T2" = "#6f9ed4", "OE21" = "#d95f02", "KYSE270" = "#1b9e77"),
      name = "Cell line"
    ) +
    labs(x = NULL, y = "Delta score (NC - Control)", title = title_text) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      plot.title = element_text(size = 13, face = "bold", hjust = 0),
      strip.text = element_text(size = 8.2, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 8)),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 8, angle = 35, hjust = 1),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(5.5, 5.5, 3, 5.5)
    )
}

up_plot <- plot_one_regulation(
  plot_data %>% filter(regulation == "Upregulated"),
  "Upregulated",
  show_legend = FALSE
)

down_plot <- plot_one_regulation(
  plot_data %>% filter(regulation == "Downregulated"),
  "Downregulated",
  show_legend = TRUE
)

p <- up_plot / down_plot + plot_layout(heights = c(1, 1.05), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  output_file,
  p,
  width = 12,
  height = 6.6,
  units = "in",
  dpi = 600
)

cat("Wrote:", output_file, "\n")
