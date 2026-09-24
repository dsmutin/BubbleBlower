"""Simple bubble detection on a totally coloured graph.

A bubble is two or more internally disjoint paths that leave one unitig and
enter one unitig. Internal unitigs have in-degree 1 and out-degree 1.
Superbubbles are left for a later contract.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass

from bubbleblower.graph import AssemblyGraph, sequence_hash


@dataclass(frozen=True)
class Branch:
    """One path from the source exclusive to the sink exclusive."""

    path: tuple[str, ...]
    links: tuple[str, ...]


@dataclass(frozen=True)
class Bubble:
    """One detected bubble. ``bubble_id`` is a sequence signature."""

    bubble_id: str
    source: str
    sink: str
    branches: tuple[Branch, ...]
    parent_bubble_id: str | None = None


def _adjacency(graph: AssemblyGraph) -> dict[str, list]:
    outgoing: dict[str, list] = {unitig.unitig_id: [] for unitig in graph.cdbg.unitigs}
    for link in graph.cdbg.links:
        outgoing[link.source].append(link)
    for links in outgoing.values():
        links.sort(key=lambda item: item.link_id)
    return outgoing


def _degrees(graph: AssemblyGraph) -> tuple[dict[str, int], dict[str, int]]:
    indeg = {unitig.unitig_id: 0 for unitig in graph.cdbg.unitigs}
    outdeg = {unitig.unitig_id: 0 for unitig in graph.cdbg.unitigs}
    for link in graph.cdbg.links:
        outdeg[link.source] += 1
        indeg[link.target] += 1
    return indeg, outdeg


def _walk(start_link, source: str, outgoing: dict[str, list], indeg: dict[str, int], outdeg: dict[str, int], limit: int):
    """Follow a non-branching path. The far node is the candidate sink."""
    path: list[str] = []
    links = [start_link.link_id]
    seen = {source}
    node = start_link.target
    while True:
        if node in seen:
            return None
        if indeg[node] != 1 or outdeg[node] != 1:
            return tuple(path), tuple(links), node
        if len(path) >= limit:
            return None
        seen.add(node)
        path.append(node)
        nxt = outgoing[node]
        if len(nxt) != 1:
            return tuple(path), tuple(links), node
        links.append(nxt[0].link_id)
        node = nxt[0].target


def _signature(graph: AssemblyGraph, source: str, sink: str, branches: list[Branch]) -> str:
    source_hash = sequence_hash(graph.unitig(source).sequence)
    sink_hash = sequence_hash(graph.unitig(sink).sequence)
    branch_tokens = []
    for branch in branches:
        if branch.path:
            sequences = [graph.unitig(unitig_id).sequence for unitig_id in branch.path]
        else:
            sequences = [graph.link(branch.links[0]).link_id]
        branch_tokens.append(sequence_hash("".join(sequences)))
    payload = "|".join([source_hash, sink_hash, *branch_tokens])
    return hashlib.sha256(payload.encode("ascii")).hexdigest()[:16]


def detect_bubbles(graph: AssemblyGraph, *, max_internal: int = 32) -> list[Bubble]:
    """Return simple bubbles. Order is stable across identical graphs."""
    outgoing = _adjacency(graph)
    indeg, outdeg = _degrees(graph)
    grouped: dict[tuple[str, str], list[Branch]] = {}
    for source, links in outgoing.items():
        if len(links) < 2:
            continue
        for link in links:
            walked = _walk(link, source, outgoing, indeg, outdeg, max_internal)
            if walked is None:
                continue
            path, link_ids, sink = walked
            if sink == source:
                continue
            grouped.setdefault((source, sink), []).append(Branch(path, link_ids))
    bubbles: list[Bubble] = []
    for (source, sink), branches in grouped.items():
        branches = sorted(branches, key=lambda branch: branch.links)
        unique: list[Branch] = []
        seen_links: set[tuple[str, ...]] = set()
        for branch in branches:
            if branch.links in seen_links:
                continue
            seen_links.add(branch.links)
            unique.append(branch)
        if len(unique) < 2:
            continue
        internals: list[str] = []
        overlap = False
        for branch in unique:
            for unitig_id in branch.path:
                if unitig_id in internals:
                    overlap = True
                internals.append(unitig_id)
        if overlap:
            continue
        signature = _signature(graph, source, sink, unique)
        bubbles.append(Bubble(bubble_id=signature, source=source, sink=sink, branches=tuple(unique)))
    bubbles.sort(key=lambda bubble: bubble.bubble_id)
    return bubbles
