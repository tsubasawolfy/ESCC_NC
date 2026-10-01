$ErrorActionPreference = "Stop"

$ProjectRoot = "F:\RNA_Seq\Autophagy"
$RunLabel = "nfcore_star_salmon_20260616"
$Base = Join-Path $ProjectRoot "downstream_inputs\$RunLabel"

$Dirs = @(
    "$Base\expression_matrices",
    "$Base\metadata",
    "$Base\qc_tables",
    "$Base\reports",
    "$Base\run_manifests"
)
New-Item -ItemType Directory -Force -Path $Dirs | Out-Null

$StarSalmon = Join-Path $ProjectRoot "results_star_salmon\star_salmon"
$MatrixFiles = @(
    "salmon.merged.gene_counts.tsv",
    "salmon.merged.gene_counts_scaled.tsv",
    "salmon.merged.gene_counts_length_scaled.tsv",
    "salmon.merged.gene_tpm.tsv",
    "salmon.merged.gene_lengths.tsv",
    "salmon.merged.gene.SummarizedExperiment.rds",
    "salmon.merged.transcript_counts.tsv",
    "salmon.merged.transcript_tpm.tsv",
    "salmon.merged.transcript_lengths.tsv",
    "salmon.merged.transcript.SummarizedExperiment.rds",
    "salmon.merged.tx2gene.tsv",
    "salmon.merged.tx2gene_augmented.tsv"
)
foreach ($File in $MatrixFiles) {
    Copy-Item -LiteralPath (Join-Path $StarSalmon $File) -Destination "$Base\expression_matrices\$File" -Force
}

Copy-Item -LiteralPath "$ProjectRoot\workflow\samplesheet.csv" -Destination "$Base\metadata\samplesheet.input.csv" -Force
Copy-Item -LiteralPath "$ProjectRoot\results_star_salmon\samplesheets\samplesheet_with_bams.csv" -Destination "$Base\metadata\samplesheet_with_bams.csv" -Force

$SampleRows = Import-Csv -Path "$ProjectRoot\results_star_salmon\samplesheets\samplesheet_with_bams.csv"
$Metadata = foreach ($Row in $SampleRows) {
    $Parts = $Row.sample -split "-"
    [PSCustomObject]@{
        sample = $Row.sample
        batch_number = $Parts[0]
        cell_line = $Parts[1]
        condition = $Parts[2]
        strandedness = $Row.strandedness
        fastq_1 = $Row.fastq_1
        fastq_2 = $Row.fastq_2
        percent_mapped = $Row.percent_mapped
    }
}
$Metadata | Export-Csv -Path "$Base\metadata\sample_metadata.tsv" -Delimiter "`t" -NoTypeInformation

$QCSrc = Join-Path $ProjectRoot "results_star_salmon\multiqc\star_salmon\multiqc_report_data"
$QCFiles = @(
    "multiqc_general_stats.txt",
    "star_summary_table.txt",
    "multiqc_strand_check_summary_table.txt",
    "fastp_filtered_reads_plot.txt",
    "multiqc_fastp.txt",
    "multiqc_star.txt",
    "multiqc_picard_dups.txt",
    "multiqc_rseqc_bam_stat.txt",
    "multiqc_rseqc_infer_experiment.txt",
    "multiqc_rseqc_read_distribution.txt",
    "multiqc_samtools_flagstat.txt",
    "multiqc_samtools_stats.txt"
)
foreach ($File in $QCFiles) {
    Copy-Item -LiteralPath (Join-Path $QCSrc $File) -Destination "$Base\qc_tables\$File" -Force
}

Copy-Item -LiteralPath "$ProjectRoot\results_star_salmon\multiqc\star_salmon\multiqc_report.html" -Destination "$Base\reports\multiqc_report.html" -Force
Copy-Item -LiteralPath "$ProjectRoot\results_star_salmon\summary.md" -Destination "$Base\reports\run_summary.md" -Force
Copy-Item -LiteralPath "$ProjectRoot\results_star_salmon\validation\check_outputs.txt" -Destination "$Base\reports\check_outputs.txt" -Force

$PipelineInfo = Join-Path $ProjectRoot "results_star_salmon\pipeline_info"
$ManifestFiles = @(
    "$ProjectRoot\results_star_salmon\run_manifest.json",
    "$ProjectRoot\results_star_salmon\artifact_index.json",
    "$ProjectRoot\results_star_salmon\artifact_index.tsv",
    "$ProjectRoot\results_star_salmon\resources\resource_plan.json",
    "$ProjectRoot\results_star_salmon\resources\resource_manifest.tsv",
    "$ProjectRoot\results_star_salmon\resources\resource_readiness.md",
    "$PipelineInfo\nextflow_trace_20260616_160625.txt",
    "$PipelineInfo\nextflow_report_20260616_160625.html",
    "$PipelineInfo\nextflow_timeline_20260616_160625.html",
    "$PipelineInfo\pipeline_dag_2026-06-16_16-06-30.html",
    "$PipelineInfo\params_2026-06-16_16-06-48.json",
    "$PipelineInfo\nf_core_rnaseq_software_mqc_versions.yml"
)
foreach ($Source in $ManifestFiles) {
    if (Test-Path -LiteralPath $Source) {
        Copy-Item -LiteralPath $Source -Destination "$Base\run_manifests\$(Split-Path $Source -Leaf)" -Force
    }
}

Write-Host "Downstream handoff refreshed: $Base"

