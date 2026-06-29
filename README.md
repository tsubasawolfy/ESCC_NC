# Revised Publication Scripts

This folder contains the publication-facing code set for the manuscript:

**Niclosamide Suppresses Esophageal Squamous Cell Carcinoma in Association with EHMT2/G9a Downregulation and Reduced MYC/E2F Proliferative Output**

For the GitHub-facing overview, see the repository-level `README.md`.

## Contents and Manuscript Mapping

| Path | Purpose | Manuscript output |
| --- | --- | --- |
| `01_DepMap_niclosamide_ESCC_correlation.R` | DepMap ESCC niclosamide dependency and RNA-expression correlation analysis | Figure 1A-C |
| `02_GO_enrichR_negative_correlated_genes.R` | enrichR GO analysis for negatively correlated genes and methylation-related prioritization | Supplementary Figure 1 |
| `03_TCGA_EHMT2_MYC_KM_risk_table_RStudio.R` | TCGA ESCC Kaplan-Meier plots for EHMT2, MYC, and combined EHMT2/MYC expression | Figure 5, Supplementary Figure 2 |
| `04_RNA_seq_upstream/` | nf-core/rnaseq STAR-Salmon workflow archive, sample sheet, run records, parameters, trace, and software versions | Supplementary RNA-seq methods |
| `05_RNA_seq_downstream/` | Exploratory TPM/log2FC, Reactome GSEA, Hallmark score, and highlighted-gene MA plot scripts | Figure 4A-C, Supplementary Table 1 |
| `code_manifest_sha256.csv` | File size and SHA256 manifest for this revised folder | Code provenance |

## Suggested Execution

Run from the repository root unless a subfolder README says otherwise:

```r
source("Script/revised/01_DepMap_niclosamide_ESCC_correlation.R")
source("Script/revised/02_GO_enrichR_negative_correlated_genes.R")
source("Script/revised/03_TCGA_EHMT2_MYC_KM_risk_table_RStudio.R")
```

For RNA-seq downstream:

```r
Sys.setenv(AUTOPHAGY_PROJECT_DIR = "F:/RNA_Seq/Autophagy")
source("Script/revised/05_RNA_seq_downstream/00_run_NC_downstream_review_code.R")
```

## Main Revisions for Publication Release

- Removed automatic package installation from revised RNA-seq downstream code.
- Added explicit package checks and clearer failure messages.
- Added stable output folders for DepMap/enrichR outputs.
- Added session information output for revised DepMap/enrichR analyses.
- Retained original TCGA KM grouping logic while adding risk-table plotting for RStudio export.
- Preserved upstream RNA-seq run records as provenance.

## Validation

R parse checks passed for all revised R scripts. Python syntax checks passed for RNA-seq upstream helper scripts. The revised-code manifest is recorded in `code_manifest_sha256.csv`.