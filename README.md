# Revised Publication Scripts

This folder contains the publication-facing revised script set. Originals remain one level above in `Script/`.

## Contents

- `01_DepMap_niclosamide_ESCC_correlation.R`: DepMap ESCC niclosamide dependency and RNA-expression correlation analysis.
- `02_GO_enrichR_negative_correlated_genes.R`: enrichR GO analysis for negatively correlated genes from the DepMap analysis.
- `03_TCGA_EHMT2_MYC_KM_risk_table_RStudio.R`: TCGA ESCC KM plots with risk tables for manual RStudio export.
- `04_RNA_seq_upstream/`: nf-core/rnaseq upstream workflow archive and run records.
- `05_RNA_seq_downstream/`: downstream RNA-seq analysis and manuscript plotting scripts.
- `code_manifest_sha256.csv`: file size and SHA256 manifest for this revised folder.

## Suggested Execution

Run from the project root unless a subfolder README says otherwise:

```r
source("Script/revised/01_DepMap_niclosamide_ESCC_correlation.R")
source("Script/revised/02_GO_enrichR_negative_correlated_genes.R")
source("Script/revised/03_TCGA_EHMT2_MYC_KM_risk_table_RStudio.R")
```

For RNA-seq downstream, set the project root if it differs from the original location:

```r
Sys.setenv(AUTOPHAGY_PROJECT_DIR = "F:/RNA_Seq/Autophagy")
source("Script/revised/05_RNA_seq_downstream/00_run_NC_downstream_review_code.R")
```

## Main Revisions

- Removed automatic package installation from revised RNA-seq downstream code.
- Added explicit package checks and clearer failure messages.
- Added stable output folders for DepMap/enrichR outputs.
- Added session information output for revised DepMap/enrichR analyses.
- Retained original TCGA KM grouping logic while adding risk-table plotting for RStudio export.

## Validation

R parse checks passed for all revised R scripts. Python bytecode checks passed for RNA-seq upstream helper scripts.