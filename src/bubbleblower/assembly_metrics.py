"""Assembly metrics from minimap2 alignments to reference genomes.

A misassembly is a contig whose alignments of at least ``min_align`` bases
hit more than one reference. Genome fraction is the share of reference bases
covered by those alignments. Mismatches per 100 kbp use the ``NM`` tag.
"""

from __future__ import annotations

import subprocess
from pathlib import Path


def n50(lengths: list[int]) -> int:
    """Contig N50. Empty input returns 0."""
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2.0
    acc = 0
    for length in ordered:
        acc += length
        if acc >= half:
            return length
    return ordered[-1]


def fasta_lengths(path: str | Path) -> list[int]:
    """Sequence lengths in a FASTA file."""
    lengths: list[int] = []
    current = 0
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if current:
                lengths.append(current)
            current = 0
        else:
            current += len(line.strip())
    if current:
        lengths.append(current)
    return lengths


def _reference_lengths(reference_fasta: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name = ""
    current = 0
    for line in reference_fasta.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name:
                lengths[name] = current
            name = line[1:].split()[0]
            current = 0
        else:
            current += len(line.strip())
    if name:
        lengths[name] = current
    return lengths


def _covered_length(spans: list[tuple[int, int]]) -> int:
    if not spans:
        return 0
    ordered = sorted(spans)
    total = 0
    start, end = ordered[0]
    for left, right in ordered[1:]:
        if left <= end:
            end = max(end, right)
        else:
            total += end - start
            start, end = left, right
    return total + end - start


def quast_like(
    assembly: str | Path,
    references: str | Path,
    *,
    minimap2: str = "minimap2",
    min_align: int = 500,
) -> dict[str, float]:
    """Score one assembly against a multi-FASTA reference.

    ``minimap2`` is the executable. The function raises if the aligner fails.
    """
    assembly = Path(assembly)
    references = Path(references)
    lengths = fasta_lengths(assembly)
    ref_lengths = _reference_lengths(references)
    proc = subprocess.run(
        [minimap2, "-x", "asm5", "-c", "--secondary=no", str(references), str(assembly)],
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or "minimap2 failed")
    intervals: dict[str, list[tuple[int, int]]] = {name: [] for name in ref_lengths}
    by_query: dict[str, set[str]] = {}
    mismatch = 0
    aligned = 0
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        query, ref = fields[0], fields[5]
        qstart, qend = int(fields[2]), int(fields[3])
        rstart, rend = sorted((int(fields[7]), int(fields[8])))
        span = rend - rstart
        if span < min_align or qend - qstart < min_align:
            continue
        if ref not in intervals:
            continue
        by_query.setdefault(query, set()).add(ref)
        aligned += span
        intervals[ref].append((rstart, rend))
        for tag in fields[12:]:
            if tag.startswith("NM:i:"):
                mismatch += int(tag.split(":")[-1])
    ref_total = sum(ref_lengths.values())
    covered_bases = sum(_covered_length(spans) for spans in intervals.values())
    misassemblies = sum(1 for refs in by_query.values() if len(refs) > 1)
    return {
        "n50": float(n50(lengths)),
        "n_contigs": float(len(lengths)),
        "total_length": float(sum(lengths)),
        "genome_fraction": covered_bases / ref_total if ref_total else 0.0,
        "misassemblies": float(misassemblies),
        "mismatches_per_100kbp": (100000.0 * mismatch / aligned) if aligned else 0.0,
        "duplication": (aligned / ref_total) if ref_total else 0.0,
    }
