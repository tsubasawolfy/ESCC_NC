#!/usr/bin/env python3
"""Check expected nf-core/rnaseq STAR + Salmon output files."""

from __future__ import annotations

import csv
import os
from pathlib import Path


OUTDIR = Path(os.environ.get("OUTDIR", "results_star_salmon"))
INPUT_SAMPLESHEET = Path(os.environ.get("SAMPLESHEET", "samplesheet.csv"))


def any_file(pattern: str) -> bool:
    return any(path.is_file() for path in OUTDIR.rglob(pattern))


def read_samples_from_csv(path: Path) -> list[str]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [
            row["sample"].strip()
            for row in csv.DictReader(handle)
            if row.get("sample", "").strip()
        ]


def expected_samples() -> list[str]:
    validated = OUTDIR / "samplesheets/samplesheet_with_bams.csv"
    return read_samples_from_csv(validated) or read_samples_from_csv(INPUT_SAMPLESHEET)


def main() -> None:
    if not OUTDIR.exists():
        raise SystemExit(f"Missing output directory: {OUTDIR}")

    files = list(OUTDIR.rglob("*"))
    names = [path.as_posix() for path in files]

    checks = {
        "multiqc_report": any(path.name == "multiqc_report.html" for path in files),
        "nextflow_report": any_file("nextflow_report*.html"),
        "nextflow_timeline": any_file("nextflow_timeline*.html"),
        "nextflow_trace": any_file("nextflow_trace*.txt"),
        "salmon_outputs": any("salmon" in path.lower() for path in names),
        "star_outputs": any("star" in path.lower() for path in names),
    }

    missing_checks = [name for name, ok in checks.items() if not ok]
    if missing_checks:
        raise SystemExit("Missing expected outputs: " + ", ".join(missing_checks))

    samples = expected_samples()
    if not samples:
        raise SystemExit(f"No expected samples found in {INPUT_SAMPLESHEET}")

    missing_samples = [sample for sample in sorted(samples) if not any(sample in path for path in names)]
    if missing_samples:
        raise SystemExit("Missing sample names in outputs: " + ", ".join(missing_samples))

    print("outputs_ok")
    print("samples=" + ",".join(samples))
    for name in sorted(checks):
        print(f"{name}=ok")


if __name__ == "__main__":
    main()
