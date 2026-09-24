"""Split a contig where the supporting read colour changes.

Read identifiers carry an accession, as in the InSilicoSeq headers. A k-mer
is painted with the accession that contains it most often. A contig that
switches accession is emitted as separate pieces and those pieces are not
joined back together. Short flickers are absorbed into the neighbouring run.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path


def _genome_id(read_id: str) -> str:
    token = read_id.split("_", 2)
    if len(token) < 2:
        raise ValueError(f"read id has no accession: {read_id}")
    return f"{token[0]}_{token[1]}"


def _kmers_of(sequence: str, k: int) -> set[str]:
    if len(sequence) < k:
        return set()
    return {sequence[index : index + k] for index in range(len(sequence) - k + 1)}


def _read_kmer_counts(fastq_paths: list[Path], needed: set[str], k: int) -> dict[str, Counter[str]]:
    hits: dict[str, Counter[str]] = defaultdict(Counter)
    for path in fastq_paths:
        lines = path.read_text(encoding="utf-8").splitlines()
        for offset in range(0, len(lines), 4):
            genome = _genome_id(lines[offset][1:].split()[0])
            sequence = lines[offset + 1].strip().upper()
            seen: set[str] = set()
            if len(sequence) < k:
                continue
            for index in range(len(sequence) - k + 1):
                kmer = sequence[index : index + k]
                if kmer not in needed or kmer in seen:
                    continue
                seen.add(kmer)
                hits[kmer][genome] += 1
    return hits


def _pieces(sequence: str, hits: dict[str, Counter[str]], k: int, min_run: int) -> list[tuple[str, str | None]]:
    if len(sequence) < k:
        return [(sequence, None)]
    labels: list[str | None] = []
    for index in range(len(sequence) - k + 1):
        counts = hits.get(sequence[index : index + k])
        labels.append(None if not counts else counts.most_common(1)[0][0])
    runs: list[tuple[int, int, str | None]] = []
    current = labels[0]
    start = 0
    for index, label in enumerate(labels[1:], start=1):
        if label != current:
            runs.append((start, index, current))
            current = label
            start = index
    runs.append((start, len(labels), current))
    changed = True
    while changed and len(runs) > 1:
        changed = False
        merged: list[tuple[int, int, str | None]] = []
        for run in runs:
            if run[1] - run[0] < min_run and merged:
                left, _right, label = merged[-1]
                merged[-1] = (left, run[1], label)
                changed = True
            else:
                merged.append(run)
        runs = merged
    pieces: list[tuple[str, str | None]] = []
    for start, end, label in runs:
        stop = min(len(sequence), end + k - 1)
        pieces.append((sequence[start:stop], label))
    return pieces


def break_read_colour_chimeras(
    records: list[tuple[str, str]],
    fastq_paths: list[str | Path],
    *,
    k: int = 21,
    min_run: int = 40,
    min_piece: int = 80,
) -> list[tuple[str, str]]:
    """Return contigs, splitting any whose read colour changes along the sequence."""
    if k < 1:
        raise ValueError("k must be positive")
    needed: set[str] = set()
    for _name, sequence in records:
        needed |= _kmers_of(sequence, k)
    if not needed:
        return list(records)
    hits = _read_kmer_counts([Path(path) for path in fastq_paths], needed, k)
    if not hits:
        raise ValueError("no read k-mers hit the contigs; colours were not invented")
    broken: list[tuple[str, str]] = []
    for name, sequence in records:
        parts = _pieces(sequence, hits, k, min_run)
        genomes = {label for _piece, label in parts if label}
        if len(genomes) < 2 or len(parts) < 2:
            broken.append((name, sequence))
            continue
        for index, (piece, _label) in enumerate(parts):
            if len(piece) >= min_piece:
                broken.append((f"{name}_{index}", piece))
    return broken
