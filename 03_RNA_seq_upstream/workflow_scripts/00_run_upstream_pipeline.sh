#!/usr/bin/env bash
set -euo pipefail

# Review wrapper for the Autophagy/NC upstream FASTQ-to-expression workflow.
# This script delegates to the project-local pipeline launcher that was used
# for the completed run.

PROJECT_WORKFLOW="/mnt/f/RNA_Seq/Autophagy/workflow"

cd "${PROJECT_WORKFLOW}"

export OUTDIR="${OUTDIR:-/mnt/f/RNA_Seq/Autophagy/results_star_salmon}"
export WORKDIR="${WORKDIR:-/root/rnaseq_work_nfcore_autophagy}"
export LOG_DIR="${LOG_DIR:-/mnt/f/RNA_Seq/Autophagy/logs}"
export PROFILE="${PROFILE:-docker}"

bash run_full_pipeline_wsl.sh

