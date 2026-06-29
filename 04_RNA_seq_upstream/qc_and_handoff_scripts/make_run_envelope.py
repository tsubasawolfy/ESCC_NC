#!/usr/bin/env python3
"""Create a compact run envelope for the nf-core/rnaseq output folder."""

from __future__ import annotations

import csv
import json
import os
from datetime import datetime, timezone
from pathlib import Path


OUTDIR = Path(os.environ.get("OUTDIR", "results_star_salmon"))
INPUT_SAMPLESHEET = Path(os.environ.get("SAMPLESHEET", "samplesheet.csv"))
REFERENCE_DIR = Path("references/gencode_v49")


def rel(path: Path) -> str:
    try:
        return path.relative_to(Path.cwd()).as_posix()
    except ValueError:
        return path.as_posix()


def add_existing(artifacts: list[dict[str, object]], category: str, pattern: str) -> None:
    for path in sorted(OUTDIR.glob(pattern)):
        if path.is_file():
            artifacts.append(
                {
                    "category": category,
                    "path": rel(path),
                    "bytes": path.stat().st_size,
                }
            )


def read_checksums() -> list[dict[str, str]]:
    checksum_file = REFERENCE_DIR / "checksums.sha256"
    checksums: list[dict[str, str]] = []
    if not checksum_file.exists():
        return checksums
    for line in checksum_file.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) >= 2:
            checksums.append({"sha256": parts[0], "file": parts[1]})
    return checksums


def read_validated_samplesheet() -> list[dict[str, str]]:
    path = OUTDIR / "samplesheets/samplesheet_with_bams.csv"
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def read_input_samplesheet() -> list[dict[str, str]]:
    if not INPUT_SAMPLESHEET.exists():
        return []
    with INPUT_SAMPLESHEET.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def sample_names(rows: list[dict[str, str]]) -> list[str]:
    return [row.get("sample", "").strip() for row in rows if row.get("sample", "").strip()]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if not OUTDIR.exists():
        raise SystemExit(f"Missing output directory: {OUTDIR}")

    artifacts: list[dict[str, object]] = []
    add_existing(artifacts, "multiqc", "multiqc/star_salmon/multiqc_report.html")
    add_existing(artifacts, "software_versions", "pipeline_info/*software*mqc_versions.yml")
    add_existing(artifacts, "params", "pipeline_info/params_*.json")
    add_existing(artifacts, "nextflow_report", "pipeline_info/nextflow_report*.html")
    add_existing(artifacts, "nextflow_timeline", "pipeline_info/nextflow_timeline*.html")
    add_existing(artifacts, "nextflow_trace", "pipeline_info/nextflow_trace*.txt")
    add_existing(artifacts, "samplesheet", "samplesheets/*.csv")
    add_existing(artifacts, "expression_matrix", "star_salmon/salmon.merged.*")
    add_existing(artifacts, "clean_export_manifest", "star_salmon/_clean_exports/99_manifests/*.tsv")

    samples = read_validated_samplesheet()
    input_samples = read_input_samplesheet()
    sample_list = sample_names(samples) or sample_names(input_samples)

    for sample in sample_list:
        add_existing(artifacts, "salmon_quant", f"star_salmon/{sample}/quant.sf")
        add_existing(artifacts, "salmon_gene_quant", f"star_salmon/{sample}/quant.genes.sf")
        add_existing(artifacts, "genome_bam", f"star_salmon/{sample}.sorted.bam")
        add_existing(artifacts, "genome_bam_index", f"star_salmon/{sample}.sorted.bam.bai")

    artifacts.sort(key=lambda x: (str(x["category"]), str(x["path"])))

    reference_files = [
        REFERENCE_DIR / "GRCh38.primary_assembly.genome.fa.gz",
        REFERENCE_DIR / "gencode.annotation.gtf.gz",
    ]
    resource_rows = []
    for path in reference_files:
        resource_rows.append(
            {
                "kind": "reference",
                "path": rel(path),
                "exists": path.exists(),
                "bytes": path.stat().st_size if path.exists() else 0,
            }
        )
    for path in [
        OUTDIR / "genome/index/star",
        OUTDIR / "genome/index/salmon",
    ]:
        resource_rows.append(
            {
                "kind": "index",
                "path": rel(path),
                "exists": path.exists(),
                "bytes": "",
            }
        )

    validation_dir = OUTDIR / "validation"
    resources_dir = OUTDIR / "resources"
    validation_dir.mkdir(parents=True, exist_ok=True)
    resources_dir.mkdir(parents=True, exist_ok=True)

    write_tsv(OUTDIR / "artifact_index.tsv", artifacts, ["category", "path", "bytes"])
    (OUTDIR / "artifact_index.json").write_text(
        json.dumps(artifacts, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_tsv(resources_dir / "resource_manifest.tsv", resource_rows, ["kind", "path", "exists", "bytes"])

    missing_resources = [row for row in resource_rows if not row["exists"]]
    resource_plan = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "status": "ready" if not missing_resources else "missing_resources",
        "missing": missing_resources,
        "resources": resource_rows,
        "reference_checksums": read_checksums(),
    }
    (resources_dir / "resource_plan.json").write_text(
        json.dumps(resource_plan, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    readiness_lines = [
        "# Resource Readiness",
        "",
        f"Status: `{resource_plan['status']}`",
        "",
        "| Kind | Path | Exists | Bytes |",
        "|---|---|---:|---:|",
    ]
    for row in resource_rows:
        readiness_lines.append(f"| {row['kind']} | `{row['path']}` | {row['exists']} | {row['bytes']} |")
    (resources_dir / "resource_readiness.md").write_text("\n".join(readiness_lines) + "\n", encoding="utf-8")

    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "workflow": "nf-core/rnaseq",
        "workflow_version": "3.26.0",
        "route": "STAR + Salmon",
        "outdir": rel(OUTDIR),
        "samples": sample_names(samples) or sample_list,
        "validated_samplesheet": rel(OUTDIR / "samplesheets/samplesheet_with_bams.csv"),
        "reference": {
            "release": "GENCODE human Release 49 GRCh38",
            "fasta": rel(REFERENCE_DIR / "GRCh38.primary_assembly.genome.fa.gz"),
            "gtf": rel(REFERENCE_DIR / "gencode.annotation.gtf.gz"),
            "checksums": read_checksums(),
        },
        "artifact_index": rel(OUTDIR / "artifact_index.json"),
        "resource_plan": rel(resources_dir / "resource_plan.json"),
    }
    (OUTDIR / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    summary_lines = [
        "# RNA-seq Run Summary",
        "",
        "Workflow: nf-core/rnaseq 3.26.0",
        "Route: STAR + Salmon",
        "Reference: GENCODE human Release 49 GRCh38",
        "",
        "## Key Outputs",
        "",
        "- MultiQC: `multiqc/star_salmon/multiqc_report.html`",
        "- Gene TPM: `star_salmon/salmon.merged.gene_tpm.tsv`",
        "- Gene counts: `star_salmon/salmon.merged.gene_counts.tsv`",
        "- Transcript TPM: `star_salmon/salmon.merged.transcript_tpm.tsv`",
        "- Transcript counts: `star_salmon/salmon.merged.transcript_counts.tsv`",
        "",
        "## Samples",
        "",
    ]
    if samples:
        summary_lines.extend(f"- {row.get('sample', '')}: strandedness={row.get('strandedness', '')}, percent_mapped={row.get('percent_mapped', '')}" for row in samples)
    else:
        summary_lines.extend(f"- {sample}" for sample in sample_list)
    summary_lines.extend(
        [
            "",
            "## Run Envelope",
            "",
            "- `run_manifest.json`",
            "- `artifact_index.json`",
            "- `artifact_index.tsv`",
            "- `resources/resource_plan.json`",
            "- `resources/resource_manifest.tsv`",
            "- `resources/resource_readiness.md`",
        ]
    )
    (OUTDIR / "summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    print(f"run_manifest={OUTDIR / 'run_manifest.json'}")
    print(f"artifact_index={OUTDIR / 'artifact_index.tsv'}")
    print(f"resource_plan={resources_dir / 'resource_plan.json'}")


if __name__ == "__main__":
    main()
