# Materials and Methods: Downstream RNA-seq Analysis

## RNA-seq Expression Matrix and Gene Annotation

Downstream RNA-seq analysis was performed using the gene-level TPM matrix generated from the STAR-Salmon workflow (`salmon.merged.gene_tpm.tsv`). Transcript-level abundance files were not used for the downstream analyses. The six analyzed samples represented matched control and niclosamide-treated conditions in three esophageal squamous cell carcinoma cell lines: CE81T2/81T (`29-81T-Con` and `30-81T-NC`), OE21 (`31-OE21-Con` and `32-OE21-NC`), and KYSE270 (`33-KYSE-Con` and `34-KYSE-NC`).

Gene annotations were derived from the local GENCODE v49 primary assembly GTF file (`gencode.v49.primary_assembly.annotation.gtf.gz`). Gene identifiers were matched after removing Ensembl version suffixes, and only genes annotated as `gene_type == "protein_coding"` were retained for analysis. Gene symbols from the TPM matrix were used when available; otherwise, GENCODE gene names were used.

## Exploratory Gene-level Fold-change Analysis

Because each condition contained one sample per cell line, the analysis was treated as exploratory and replicate-aware differential expression testing was not performed. DESeq2, gene-level P values, and gene-level FDR values were not calculated. For each cell line, niclosamide-treated and control TPM values were compared using a pseudocount of 0.1:

```text
log2FC = log2((TPM_NC + 0.1) / (TPM_Control + 0.1))
```

Positive log2FC values therefore indicate higher expression after niclosamide treatment. Genes were retained for comparison-level fold-change summaries if TPM was at least 1 in either the control or niclosamide-treated sample. Genes with absolute log2FC at least 1.5 were highlighted as strongly altered in the exploratory analysis. The three comparisons were CE81T2/81T niclosamide versus control, OE21 niclosamide versus control, and KYSE270 niclosamide versus control.

MA plots were generated from log10(mean TPM + 0.1) and log2FC values for each comparison. For manuscript visualization, EHMT2, MAP1LC3B, MYC, and SQSTM1 were highlighted as colored points on a combined three-panel MA plot. The plot was generated at 600 DPI for downstream figure assembly.

## Gene Set Enrichment Analysis

Gene set enrichment analysis (GSEA) was performed using full ranked gene lists from each comparison. Genes were ranked by exploratory log2FC, with positive ranks corresponding to higher expression after niclosamide treatment. Ensembl gene identifiers were converted to Entrez identifiers using `org.Hs.eg.db`; when multiple entries mapped to the same Entrez identifier, the entry with the largest absolute log2FC was retained. Ranked vectors were analyzed using `clusterProfiler::gseKEGG` for KEGG pathways and `ReactomePA::gsePathway` for Reactome pathways, with gene set size limits of 10 to 500 genes and `pvalueCutoff = 1` to retain complete ranked enrichment output for review.

Reactome GSEA results were used for the manuscript-focused selected pathway visualization. Selected cell-line-specific pathways were curated from the Reactome GSEA output based on normalized enrichment score (NES), adjusted P value, rank, and biological interpretability. The selected pathways emphasized cornified-envelope and keratinocyte differentiation programs in KYSE270 and CE81T2/81T, and cell-cycle, mitotic, rRNA-processing, and tRNA-processing pathways in OE21. NES values were displayed with pathway direction after niclosamide treatment, and adjusted P values from the GSEA output were used as the enrichment significance metric.

## Hallmark Cancer Signature Scoring

MSigDB Hallmark gene sets were retrieved for Homo sapiens using `msigdbr`. Cancer-relevant Hallmark pathways were scored using the protein-coding TPM matrix transformed as log2(TPM + 0.1). For each Hallmark gene set, genes present in the expression matrix were extracted and row-scaled across the six samples. A sample-level pathway score was calculated as the mean row-z score across the genes present in the pathway. For each cell line, the niclosamide-induced signature alteration was calculated as:

```text
Delta score = Hallmark score_NC - Hallmark score_Control
```

Positive delta scores indicate higher relative pathway activity after niclosamide treatment, whereas negative delta scores indicate lower relative pathway activity. Selected Hallmark cancer signatures for manuscript visualization included Apoptosis, E2F Targets, G2M Checkpoint, Hedgehog Signaling, Hypoxia, MYC Targets V1, MYC Targets V2, Notch Signaling, Oxidative Phosphorylation, P53 Pathway, and Unfolded Protein Response. Selected signatures were plotted as niclosamide-control delta scores across CE81T2, OE21, and KYSE270, with upregulated and downregulated signatures displayed separately.

## Software

Analyses were performed in R 4.4.3. The main packages used for the Results-relevant analyses were `clusterProfiler` 4.12.6, `ReactomePA` 1.48.0, `reactome.db` 1.88.0, `org.Hs.eg.db` 3.19.1, `msigdbr` 26.1.0, `fgsea` 1.30.0, `ggplot2` 4.0.3, `ggrepel` 0.9.8, `dplyr` 1.2.1, `readr` 2.2.0, `tidyr` 1.3.2, and `stringr` 1.6.0.

