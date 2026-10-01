# NC Upstream RNA-seq Code

This folder contains the code and run records used to process the Autophagy/NC paired-end FASTQ files into expression matrices for publication review.

## Folder Contents

```text
workflow_scripts/          Shell and Python scripts used to run nf-core/rnaseq STAR + Salmon.
config_and_samplesheet/    Samplesheet and Nextflow resource-limit config used for the run.
qc_and_handoff_scripts/    Output validation, run-envelope, and downstream handoff scripts.
run_records/               Non-executable provenance files from the completed run.
```

## Execution Order

The upstream processing was run from WSL with:

```bash
cd /mnt/f/RNA_Seq/Autophagy/workflow
bash run_full_pipeline_wsl.sh
```

For review, the copied equivalent is:

```bash
cd /mnt/f/RNA_Seq/Autophagy/NC_upstream_code
bash workflow_scripts/00_run_upstream_pipeline.sh
```

The wrapper calls the same pipeline launcher and validation steps. It assumes the original project folder structure and FASTQ paths are present under `/mnt/f/RNA_Seq/Autophagy`.

## Main Pipeline

The analysis used:

```text
nf-core/rnaseq 3.26.0
Nextflow 25.10.0
STAR + Salmon route
GENCODE human Release 49 GRCh38 FASTA/GTF
fastp trimming with --trim_poly_g --length_required 20
paired-end FASTQ input
strandedness=auto, resolved as reverse
```

Existing GENCODE v49 GRCh38 STAR and Salmon indexes were reused from:

```text
/mnt/e/rnaseq_project/workflow/results_star_salmon/genome/index/star
/mnt/e/rnaseq_project/workflow/results_star_salmon/genome/index/salmon
```

## Main Outputs

The final expression matrices are archived in:

```text
F:\RNA_Seq\Autophagy\downstream_inputs\nfcore_star_salmon_20260616\expression_matrices
```

For differential expression, use:

```text
salmon.merged.gene_counts.tsv
```

For expression visualization and sanity checks, use:

```text
salmon.merged.gene_tpm.tsv
```

## Provenance

The completed run validation is copied to:

```text
run_records/check_outputs.txt
```

The run manifest, nf-core parameters, final Nextflow trace, and software versions are copied to `run_records/`.

