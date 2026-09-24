"""Walk a resolved graph into contig sequences.

A contig is a path whose internal unitigs have in-degree 1 and out-degree 1.
Overlap stored on the outgoing link is consumed. A missing overlap is an error.
"""

from __future__ import annotations

from pathlib import Path

from bubbleblower.graph import AssemblyGraph


def contig_sequences(graph: AssemblyGraph) -> list[tuple[str, str]]:
    """Return ``(contig_id, sequence)`` in stable order."""
    outgoing: dict[str, list] = {unitig.unitig_id: [] for unitig in graph.cdbg.unitigs}
    indeg = {unitig.unitig_id: 0 for unitig in graph.cdbg.unitigs}
    for link in graph.cdbg.links:
        outgoing[link.source].append(link)
        indeg[link.target] += 1
    starts = [
        unitig.unitig_id
        for unitig in graph.cdbg.unitigs
        if indeg[unitig.unitig_id] != 1 or len(outgoing[unitig.unitig_id]) != 1
    ]
    if not starts:
        starts = [graph.cdbg.unitigs[0].unitig_id] if graph.cdbg.unitigs else []
    seen: set[str] = set()
    records: list[tuple[str, str]] = []
    for start in sorted(starts):
        if start in seen:
            continue
        node = start
        sequence = graph.unitig(node).sequence
        path = [node]
        seen.add(node)
        while len(outgoing[node]) == 1 and indeg[outgoing[node][0].target] == 1:
            link = outgoing[node][0]
            nxt = link.target
            if nxt in seen:
                break
            if link.overlap is None:
                raise ValueError(f"link {link.link_id} has no overlap; refusing to invent a join")
            extra = graph.unitig(nxt).sequence[link.overlap :]
            sequence += extra
            path.append(nxt)
            seen.add(nxt)
            node = nxt
        records.append((f"contig_{len(records) + 1:05d}", sequence))
    for unitig in graph.cdbg.unitigs:
        if unitig.unitig_id not in seen:
            records.append((f"contig_{len(records) + 1:05d}", unitig.sequence))
    return records


def write_fasta(records: list[tuple[str, str]], path: str | Path) -> None:
    """Write contig records as FASTA."""
    lines = []
    for name, sequence in records:
        lines.append(f">{name}")
        lines.append(sequence)
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
