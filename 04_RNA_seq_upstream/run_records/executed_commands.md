# Executed Upstream Commands

The upstream RNA-seq processing was run in WSL from the Autophagy workflow folder.

## Runtime Check

```bash
cd /mnt/f/RNA_Seq/Autophagy/workflow
python3 -m py_compile check_outputs.py make_run_envelope.py validate_reference.py
python3 validate_reference.py
```

## Full Pipeline Run

```bash
cd /mnt/f/RNA_Seq/Autophagy/workflow
bash run_full_pipeline_wsl.sh
```

The launcher executed nf-core/rnaseq with these key parameters:

```bash
nextflow run nf-core/rnaseq -r 3.26.0 \
  -c rnaseq_resource_limits.config \
  --input samplesheet.csv \
  --outdir /mnt/f/RNA_Seq/Autophagy/results_star_salmon \
  -work-dir /root/rnaseq_work_nfcore_autophagy \
  --aligner star_salmon \
  --fasta references/gencode_v49/GRCh38.primary_assembly.genome.fa.gz \
  --gtf references/gencode_v49/gencode.annotation.gtf.gz \
  --star_index /mnt/e/rnaseq_project/workflow/results_star_salmon/genome/index/star \
  --salmon_index /mnt/e/rnaseq_project/workflow/results_star_salmon/genome/index/salmon \
  --gencode \
  --trimmer fastp \
  --extra_fastp_args "--trim_poly_g --length_required 20" \
  --seq_platform ILLUMINA \
  --save_reference \
  --save_trimmed \
  --save_align_intermeds \
  -profile docker \
  -resume
```

## Output Validation

```bash
OUTDIR=/mnt/f/RNA_Seq/Autophagy/results_star_salmon python3 check_outputs.py
OUTDIR=/mnt/f/RNA_Seq/Autophagy/results_star_salmon python3 make_run_envelope.py
```

## Downstream Handoff Creation

```powershell
powershell -ExecutionPolicy Bypass -File F:\RNA_Seq\Autophagy\NC_upstream_code\qc_and_handoff_scripts\create_downstream_handoff.ps1
```

