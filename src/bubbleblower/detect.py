"""Simple bubble detection on a totally coloured graph.

A bubble is two or more internally disjoint paths that leave one oriented
unitig and enter one oriented unitig. Internal unitigs have in-degree 1 and
out-degree 1. Each link is an arc between oriented unitigs (``orientation``,
default ``++``). A graph that stores both strands of an adjacency (MEGAHIT)
yields every bubble once per strand; the copy with the smaller ``bubble_id``
is kept. Superbubbles are left for a later contract.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass

from bubbleblower.graph import AssemblyGraph, sequence_hash


_FLIP_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

Oriented = tuple[str, str]


@dataclass(frozen=True)
class Branch:
    """One path from the source exclusive to the sink exclusive.

    ``strands`` gives the strand (``+`` or ``-``) of each unitig in ``path``.
    """

    path: tuple[str, ...]
    links: tuple[str, ...]
    strands: tuple[str, ...] = ()


@dataclass(frozen=True)
class Bubble:
    """One detected bubble. ``bubble_id`` is a sequence signature over oriented sequences."""

    bubble_id: str
    source: str
    sink: str
    branches: tuple[Branch, ...]
    parent_bubble_id: str | None = None
    source_strand: str = "+"
    sink_strand: str = "+"


def _oriented_sequence(graph: AssemblyGraph, node: Oriented) -> str:
    sequence = graph.unitig(node[0]).sequence
    return sequence if node[1] == "+" else sequence.translate(_FLIP_COMPLEMENT)[::-1]


def _arcs(graph: AssemblyGraph):
    outgoing: dict[Oriented, list[tuple[str, Oriented]]] = {}
    indeg: dict[Oriented, int] = {}
    outdeg: dict[Oriented, int] = {}
    for link in graph.cdbg.links:
        orientation = link.orientation or "++"
        source, target = (link.source, orientation[0]), (link.target, orientation[1])
        outgoing.setdefault(source, []).append((link.link_id, target))
        outdeg[source] = outdeg.get(source, 0) + 1
        indeg[target] = indeg.get(target, 0) + 1
    for arcs in outgoing.values():
        arcs.sort()
    return outgoing, indeg, outdeg


def _walk(link_id: str, node: Oriented, source: Oriented, outgoing, indeg, outdeg, limit: int):
    """Follow a non-branching oriented path. The far node is the candidate sink."""
    path: list[Oriented] = []
    links = [link_id]
    seen = {source}
    while True:
        if node in seen:
            return None
        if indeg.get(node, 0) != 1 or outdeg.get(node, 0) != 1:
            return tuple(path), tuple(links), node
        if len(path) >= limit:
            return None
        seen.add(node)
        path.append(node)
        links.append(outgoing[node][0][0])
        node = outgoing[node][0][1]


def _signature(graph: AssemblyGraph, source: Oriented, sink: Oriented, branches) -> str:
    tokens = [sequence_hash(_oriented_sequence(graph, source)), sequence_hash(_oriented_sequence(graph, sink))]
    for path, links in branches:
        text = "".join(_oriented_sequence(graph, node) for node in path) if path else links[0]
        tokens.append(sequence_hash(text))
    return hashlib.sha256("|".join(tokens).encode("ascii")).hexdigest()[:16]


def detect_bubbles(graph: AssemblyGraph, *, max_internal: int = 32) -> list[Bubble]:
    """Return simple bubbles, one per strand pair. Order is stable across identical graphs."""
    outgoing, indeg, outdeg = _arcs(graph)
    found: dict[tuple, Bubble] = {}
    for source in sorted(outgoing):
        arcs = outgoing[source]
        if len(arcs) < 2:
            continue
        grouped: dict[Oriented, list[tuple[tuple[Oriented, ...], tuple[str, ...]]]] = {}
        for link_id, node in arcs:
            walked = _walk(link_id, node, source, outgoing, indeg, outdeg, max_internal)
            if walked is None:
                continue
            path, link_ids, sink = walked
            if sink == source:
                continue
            grouped.setdefault(sink, []).append((path, link_ids))
        for sink, branches in grouped.items():
            unique = list(dict.fromkeys(sorted(branches, key=lambda item: item[1])))
            if len(unique) < 2:
                continue
            internal = [node[0] for path, _ in unique for node in path]
            if len(internal) != len(set(internal)) or {source[0], sink[0]} & set(internal):
                continue
            bubble = Bubble(
                bubble_id=_signature(graph, source, sink, unique),
                source=source[0],
                sink=sink[0],
                branches=tuple(
                    Branch(tuple(node[0] for node in path), links, tuple(node[1] for node in path))
                    for path, links in unique
                ),
                source_strand=source[1],
                sink_strand=sink[1],
            )
            key = (
                frozenset((source[0], sink[0])),
                frozenset(frozenset(branch.path) for branch in bubble.branches),
            )
            if key not in found or bubble.bubble_id < found[key].bubble_id:
                found[key] = bubble
    return sorted(found.values(), key=lambda bubble: bubble.bubble_id)
