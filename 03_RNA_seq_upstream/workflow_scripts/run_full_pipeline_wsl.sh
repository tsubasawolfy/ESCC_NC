#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

LOG_DIR="${LOG_DIR:-/mnt/f/RNA_Seq/Autophagy/logs}"
OUTDIR="${OUTDIR:-/mnt/f/RNA_Seq/Autophagy/results_star_salmon}"
WORKDIR="${WORKDIR:-/root/rnaseq_work_nfcore_autophagy}"
PROFILE="${PROFILE:-docker}"

export OUTDIR WORKDIR PROFILE

mkdir -p "$LOG_DIR"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
STDOUT_LOG="${LOG_DIR}/nfcore_star_salmon_${TIMESTAMP}.out.log"
STDERR_LOG="${LOG_DIR}/nfcore_star_salmon_${TIMESTAMP}.err.log"

exec > >(tee -a "$STDOUT_LOG") 2> >(tee -a "$STDERR_LOG" >&2)

export PATH="/root/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
export NXF_HOME="/root/.nextflow_rnaseq"

echo "started_at=$(date --iso-8601=seconds)"
echo "workflow_path=$PWD"
echo "stdout_log=$STDOUT_LOG"
echo "stderr_log=$STDERR_LOG"

java -version
nextflow -version
if ! docker info >/dev/null 2>&1; then
  service docker start >/dev/null 2>&1 || true
  sleep 5
fi
docker info --format '{{.ServerVersion}} {{.OSType}} {{.OperatingSystem}}'
python3 validate_reference.py

bash run_nfcore_star_salmon.sh

mkdir -p "${OUTDIR}/validation"
python3 check_outputs.py | tee "${OUTDIR}/validation/check_outputs.txt"
python3 make_run_envelope.py
echo "finished_at=$(date --iso-8601=seconds)"
