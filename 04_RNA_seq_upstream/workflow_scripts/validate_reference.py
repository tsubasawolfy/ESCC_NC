#!/usr/bin/env python3
"""Validate the pinned GENCODE RNA-seq reference bundle."""

from __future__ import annotations

import gzip
from pathlib import Path


REFERENCE_DIR = Path("references/gencode_v49")
FASTA = REFERENCE_DIR / "GRCh38.primary_assembly.genome.fa.gz"
GTF = REFERENCE_DIR / "gencode.annotation.gtf.gz"
CHECKSUMS = REFERENCE_DIR / "checksums.sha256"


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def first_fasta_contigs(path: Path, limit: int = 200) -> set[str]:
    contigs: set[str] = set()
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                contigs.add(line[1:].split()[0])
                if len(contigs) >= limit:
                    break
    return contigs


def first_gtf_contigs(path: Path, limit: int = 200) -> set[str]:
    contigs: set[str] = set()
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 9:
                contigs.add(fields[0])
                if len(contigs) >= limit:
                    break
    return contigs


def main() -> None:
    missing = [str(path) for path in (FASTA, GTF, CHECKSUMS) if not path.exists()]
    if missing:
        raise SystemExit("Missing reference files:\n" + "\n".join(missing))

    fasta_contigs = first_fasta_contigs(FASTA)
    gtf_contigs = first_gtf_contigs(GTF)
    required = {"chr1", "chr2", "chrM"}
    missing_required = sorted(required - fasta_contigs - gtf_contigs)
    shared_required = sorted(required & fasta_contigs & gtf_contigs)

    if missing_required or len(shared_required) != len(required):
        raise SystemExit(
            "Reference contig validation failed.\n"
            f"Required in both FASTA and GTF: {sorted(required)}\n"
            f"Found in both: {shared_required}\n"
            f"FASTA sample: {sorted(list(fasta_contigs))[:10]}\n"
            f"GTF sample: {sorted(list(gtf_contigs))[:10]}"
        )

    print("reference_ok")
    print(f"fasta={FASTA}")
    print(f"gtf={GTF}")
    print(f"shared_required={','.join(shared_required)}")


if __name__ == "__main__":
    main()
