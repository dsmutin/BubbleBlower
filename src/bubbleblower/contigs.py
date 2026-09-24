"""Walk a resolved graph into contig sequences.

A contig is a path whose internal oriented unitigs have in-degree 1 and
out-degree 1. Link orientation selects the strand. Overlap stored on the
outgoing link is consumed. A missing overlap is an error.
"""

from __future__ import annotations

from pathlib import Path

from bubbleblower.graph import AssemblyGraph

_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _oriented(sequence: str, orientation: str) -> str:
    if orientation == "+":
        return sequence
    if orientation == "-":
        return sequence.translate(_RC)[::-1]
    raise ValueError(f"orientation must be + or -, got {orientation}")


def contig_sequences(graph: AssemblyGraph) -> list[tuple[str, str]]:
    """Return ``(contig_id, sequence)`` in stable order."""
    sequences = {unitig.unitig_id: unitig.sequence for unitig in graph.cdbg.unitigs}
    outgoing: dict[tuple[str, str], list] = {(unitig_id, strand): [] for unitig_id in sequences for strand in "+-"}
    indeg = {(unitig_id, strand): 0 for unitig_id in sequences for strand in "+-"}
    for link in graph.cdbg.links:
        if len(link.orientation) != 2:
            raise ValueError(f"link {link.link_id} orientation must be two characters")
        source_strand, target_strand = link.orientation
        outgoing[(link.source, source_strand)].append(link)
        indeg[(link.target, target_strand)] += 1
    starts = [
        (unitig_id, strand)
        for unitig_id in sorted(sequences)
        for strand in "+-"
        if indeg[(unitig_id, strand)] != 1 or len(outgoing[(unitig_id, strand)]) != 1
    ]
    if not starts and sequences:
        starts = [(sorted(sequences)[0], "+")]
    seen: set[str] = set()
    records: list[tuple[str, str]] = []

    def append_record(sequence: str) -> None:
        records.append((f"contig_{len(records) + 1:05d}", sequence))

    for start_node, start_strand in starts:
        if start_node in seen:
            continue
        node, strand = start_node, start_strand
        sequence = _oriented(sequences[node], strand)
        seen.add(node)
        while len(outgoing[(node, strand)]) == 1:
            link = outgoing[(node, strand)][0]
            nxt = link.target
            nxt_strand = link.orientation[1]
            if indeg[(nxt, nxt_strand)] != 1 or nxt in seen:
                break
            if link.overlap is None:
                raise ValueError(f"link {link.link_id} has no overlap; refusing to invent a join")
            piece = _oriented(sequences[nxt], nxt_strand)
            if link.overlap > len(piece):
                raise ValueError(f"overlap on {link.link_id} is longer than the oriented target")
            sequence += piece[link.overlap :]
            seen.add(nxt)
            node, strand = nxt, nxt_strand
        append_record(sequence)
    for unitig_id in sorted(sequences):
        if unitig_id not in seen:
            append_record(sequences[unitig_id])
    return records


def write_fasta(records: list[tuple[str, str]], path: str | Path) -> None:
    """Write contig records as FASTA."""
    lines = []
    for name, sequence in records:
        lines.append(f">{name}")
        lines.append(sequence)
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
