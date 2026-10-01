# NC Downstream RNA-seq Review Code

This folder contains the R code used to generate the downstream niclosamide RNA-seq results included in the manuscript figures.

## Code Files

1. `01_run_autophagy_downstream.R`
   - Main downstream analysis script.
   - Reads gene-level Salmon TPM output, applies GENCODE v49 protein-coding filtering, calculates exploratory log2FC values, runs GSEA, and generates Hallmark cancer signature scores.
   - Uses `log2FC = log2((NC TPM + 0.1) / (Control TPM + 0.1))`.
   - Positive log2FC indicates higher expression after niclosamide treatment.

2. `02_plot_selected_reactome_gsea_pathways.R`
   - Generates the selected cell-line-specific Reactome GSEA pathway plot from `selected_cell_line_specific_reactome_gsea_pathways.csv`.
   - Used for the Reactome pathway panel emphasizing cornified-envelope/keratinocyte programs and OE21 mitotic/RNA-processing programs.

3. `03_plot_selected_hallmark_delta_6col_compact.R`
   - Generates the compact selected Hallmark cancer delta figure.
   - Uses six columns, separates upregulated and downregulated signatures, and plots delta scores as `NC - Control`.

4. `04_plot_MA_highlighted_genes_600dpi_panel.R`
   - Generates the combined three-cell-line MA plot panel.
   - Highlights `EHMT2`, `MAP1LC3B`, `MYC`, and `SQSTM1` by color.
   - Exports a 600 DPI PNG and PDF.

5. `00_run_NC_downstream_review_code.R`
   - Convenience runner that executes scripts 1-4 in order.
   - The main workflow can take several minutes.

## Supporting Files

- `Materials_and_Methods_RNAseq_downstream.md`: publication-style methods summary for the Results-used analyses.
- `package_versions.tsv`: R package versions recorded from the analysis environment.
- `code_file_manifest_sha256.csv`: file sizes and SHA256 checksums for the review-code folder.

## Inputs Expected by the Scripts

The scripts use absolute project paths under:

```text
F:/RNA_Seq/Autophagy
```

Key input files:

```text
F:/RNA_Seq/Autophagy/downstream_inputs/nfcore_star_salmon_20260616/expression_matrices/salmon.merged.gene_tpm.tsv
F:/RNA_Seq/Autophagy/downstream_inputs/nfcore_star_salmon_20260616/metadata/sample_metadata.tsv
F:/RNA_Seq/Autophagy/workflow/references/gencode_v49/gencode.v49.primary_assembly.annotation.gtf.gz
F:/RNA_Seq/Autophagy/Analysis/04_gsea/selected_cell_line_specific_reactome_gsea_pathways.csv
```

## Analyses Not Included in Manuscript Results

The following generated analyses were not used in the manuscript Results and are not emphasized in this review code summary: pairwise log2FC scatter plots, over-representation pathway enrichment, exploratory heatmaps, immune Hallmark scoring, and EMT/invasion scoring.
