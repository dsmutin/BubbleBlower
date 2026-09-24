"""Colour unitigs from reads whose identifiers name a genome.

Read ids are expected to start with an accession such as ``GCF_001549955``.
A unitig receives a colour when at least ``min_depth`` of its k-mers are seen
in that genome's reads. Coverage already stored on the graph is not replaced.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from bubbleblower.graph import AssemblyGraph


def genome_id(read_id: str) -> str:
    """Return the accession token at the start of an InSilicoSeq read id."""
    token = read_id.split("_", 2)
    if len(token) < 2:
        raise ValueError(f"read id has no accession: {read_id}")
    return f"{token[0]}_{token[1]}"


def _kmers(sequence: str, k: int) -> list[str]:
    if len(sequence) < k:
        return []
    return [sequence[index : index + k] for index in range(len(sequence) - k + 1)]


def colour_from_fastq(
    graph: AssemblyGraph,
    fastq_paths: list[str | Path],
    *,
    k: int = 31,
    min_depth: int = 2,
) -> AssemblyGraph:
    """Attach genome colours from exact k-mer hits. The input graph is copied."""
    if k < 1:
        raise ValueError("k must be positive")
    updated = graph.copy()
    index: dict[str, list[str]] = {}
    for unitig in updated.cdbg.unitigs:
        for kmer in _kmers(unitig.sequence, k):
            index.setdefault(kmer, []).append(unitig.unitig_id)
    hits: dict[str, Counter[str]] = {unitig.unitig_id: Counter() for unitig in updated.cdbg.unitigs}
    for path in fastq_paths:
        header: str | None = None
        for line_number, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines()):
            if line_number % 4 == 0:
                header = raw[1:].split()[0]
                continue
            if line_number % 4 != 1 or header is None:
                continue
            genome = genome_id(header)
            seen: set[str] = set()
            for kmer in _kmers(raw.strip().upper(), k):
                for unitig_id in index.get(kmer, ()):
                    if unitig_id in seen:
                        continue
                    seen.add(unitig_id)
                    hits[unitig_id][genome] += 1
    genomes = sorted({genome for counts in hits.values() for genome in counts})
    if not genomes:
        raise ValueError("no read k-mers hit the graph; colours were not invented")
    palette = {genome: index for index, genome in enumerate(genomes)}
    updated.cdbg.colors = [
        {"color_id": str(palette[genome]), "namespace": "genome", "value": genome} for genome in genomes
    ]
    for unitig in updated.cdbg.unitigs:
        kept = [palette[genome] for genome, depth in hits[unitig.unitig_id].items() if depth >= min_depth]
        unitig.color_ids = sorted(kept)
        for row in updated.cdbg.mapping:
            if row.unitig_id == unitig.unitig_id:
                row.color_ids = list(unitig.color_ids)
    for link in updated.cdbg.links:
        source = set(updated.unitig(link.source).color_ids)
        target = set(updated.unitig(link.target).color_ids)
        shared = source & target
        link.color_ids = sorted(shared if shared else source | target)
    updated.coverage_source = graph.coverage_source
    updated.validate()
    return updated
