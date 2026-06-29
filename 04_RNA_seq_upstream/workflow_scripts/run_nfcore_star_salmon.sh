#!/usr/bin/env bash
set -euo pipefail

# Primary modern route for transcript/gene TPM and counts.
# Run from this directory in Linux/WSL2/HPC with Nextflow and Docker/Singularity/Apptainer.

export PATH="/root/bin:${PATH}"

PROFILE="${PROFILE:-docker}"
OUTDIR="${OUTDIR:-/mnt/f/RNA_Seq/Autophagy/results_star_salmon}"
WORKDIR="${WORKDIR:-/root/rnaseq_work_nfcore_autophagy}"
REPORT_DIR="${REPORT_DIR:-${OUTDIR}/pipeline_info}"
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"
EXTRA_NXF_CONFIG="${EXTRA_NXF_CONFIG:-rnaseq_resource_limits.config}"

# Use a matched FASTA + GTF from the same provider/release.
# Recommended: download pinned GENCODE v49 GRCh38 FASTA + GTF into references/gencode_v49/.
FASTA="${FASTA:-references/gencode_v49/GRCh38.primary_assembly.genome.fa.gz}"
GTF="${GTF:-references/gencode_v49/gencode.annotation.gtf.gz}"
STAR_INDEX="${STAR_INDEX:-}"
SALMON_INDEX="${SALMON_INDEX:-}"
DEFAULT_STAR_INDEX="/mnt/e/rnaseq_project/workflow/results_star_salmon/genome/index/star"
DEFAULT_SALMON_INDEX="/mnt/e/rnaseq_project/workflow/results_star_salmon/genome/index/salmon"

if [[ -z "${STAR_INDEX}" && -d "${DEFAULT_STAR_INDEX}" ]]; then
  STAR_INDEX="${DEFAULT_STAR_INDEX}"
fi
if [[ -z "${SALMON_INDEX}" && -d "${DEFAULT_SALMON_INDEX}" ]]; then
  SALMON_INDEX="${DEFAULT_SALMON_INDEX}"
fi

mkdir -p "${REPORT_DIR}"

config_args=()
if [[ -n "${EXTRA_NXF_CONFIG}" && -f "${EXTRA_NXF_CONFIG}" ]]; then
  config_args=(-c "${EXTRA_NXF_CONFIG}")
fi

index_args=()
if [[ -n "${STAR_INDEX}" ]]; then
  echo "Reusing STAR index: ${STAR_INDEX}"
  index_args+=(--star_index "${STAR_INDEX}")
fi
if [[ -n "${SALMON_INDEX}" ]]; then
  echo "Reusing Salmon index: ${SALMON_INDEX}"
  index_args+=(--salmon_index "${SALMON_INDEX}")
fi

nextflow run nf-core/rnaseq -r 3.26.0 \
  "${config_args[@]}" \
  --input samplesheet.csv \
  --outdir "${OUTDIR}" \
  -work-dir "${WORKDIR}" \
  --aligner star_salmon \
  --fasta "${FASTA}" \
  --gtf "${GTF}" \
  "${index_args[@]}" \
  --gencode \
  --trimmer fastp \
  --extra_fastp_args "--trim_poly_g --length_required 20" \
  --seq_platform ILLUMINA \
  --save_reference \
  --save_trimmed \
  --save_align_intermeds \
  -with-report "${REPORT_DIR}/nextflow_report_${RUN_TAG}.html" \
  -with-timeline "${REPORT_DIR}/nextflow_timeline_${RUN_TAG}.html" \
  -with-trace "${REPORT_DIR}/nextflow_trace_${RUN_TAG}.txt" \
  -profile "${PROFILE}" \
  -resume
