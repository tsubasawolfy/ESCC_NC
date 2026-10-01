get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

code_dir <- get_script_dir()

scripts <- c(
  "01_run_autophagy_downstream.R",
  "02_plot_selected_reactome_gsea_pathways.R",
  "03_plot_selected_hallmark_delta_6col_compact.R",
  "04_plot_MA_highlighted_genes_600dpi_panel.R"
)

for (script in scripts) {
  script_path <- file.path(code_dir, script)
  if (!file.exists(script_path)) {
    stop("Missing review script: ", script_path)
  }
  message("Running: ", script)
  source(script_path, local = FALSE)
}

message("NC downstream review code completed.")
