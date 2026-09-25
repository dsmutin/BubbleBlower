"""Colour unitigs from reads whose identifiers name a genome.

Read colouring is implemented in MetaMetro. This module expands the ToCUMG
to a CFA, calls ``colour_by_read_accessions``, and copies colour sets back
onto the existing unitigs. Coverage already stored on the graph is not
replaced.
"""

from __future__ import annotations

from pathlib import Path

from bubbleblower.graph import AssemblyGraph
from metametro.contracts.assembly import read_fastq
from metametro.contracts.colouring import accession_from_read_id, colour_by_read_accessions
from metametro.converters.cdbg_to_cfa import cdbg_to_cfa
from metametro.errors import ContractError
from metametro.formats.cfa.validator import parse_color_set


def genome_id(read_id: str) -> str:
    """Return the accession token at the start of an InSilicoSeq read id."""
    try:
        return accession_from_read_id(read_id)
    except ContractError as exc:
        raise ValueError(str(exc)) from exc


def colour_from_fastq(
    graph: AssemblyGraph,
    fastq_paths: list[str | Path],
    *,
    k: int = 31,
    min_depth: int = 2,
) -> AssemblyGraph:
    """Attach genome colours from MetaMetro read-accession colouring. The input graph is copied."""
    if k < 1:
        raise ValueError("k must be positive")
    reads: list[tuple[str, str, str]] = []
    for path in fastq_paths:
        for read_id, sequence in read_fastq(path):
            reads.append((read_id, genome_id(read_id), sequence))
    if not reads:
        raise ValueError("no reads in the FASTQ files")
    cfa = cdbg_to_cfa(graph.cdbg)
    metadata = dict(cfa.metadata)
    if not isinstance(metadata.get("k"), int) or isinstance(metadata.get("k"), bool) or int(metadata.get("k") or 0) <= 0:
        metadata["k"] = k
        cfa.metadata = metadata
    coloured = colour_by_read_accessions(
        cfa,
        reads,
        min_vertex_depth=min_depth,
        min_edge_kmer_density=max(1, min_depth),
        operation="replace",
    )
    by_node = {row["node_id"]: parse_color_set(row.get("color_set", "")) for row in coloured.nodes}
    by_edge = {row["edge_id"]: parse_color_set(row.get("color_set", "")) for row in coloured.edges}
    updated = graph.copy()
    updated.cdbg.colors = [dict(row) for row in coloured.colors or []]
    for mapping in updated.cdbg.mapping:
        mapping.color_ids = list(by_node.get(mapping.cfa_node_id, []))
    for unitig in updated.cdbg.unitigs:
        members: set[int] = set()
        for node_id in unitig.members:
            members.update(by_node.get(node_id, []))
        unitig.color_ids = sorted(members)
        unitig.internal_edge_colors = [list(group) for group in unitig.internal_edge_colors]
    for link in updated.cdbg.links:
        link.color_ids = list(by_edge.get(link.link_id, []))
    updated.validate()
    return updated
