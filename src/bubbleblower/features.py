"""Immutable bubble features. Extraction does not edit the graph."""

from __future__ import annotations

import math
from dataclasses import dataclass

from bubbleblower.detect import Branch, Bubble
from bubbleblower.graph import AssemblyGraph


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def _median(values: list[float]) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def _gc(sequence: str) -> float:
    if not sequence:
        return 0.0
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


@dataclass(frozen=True)
class BranchFeatures:
    """Coverage, sequence, and colour summary of one branch."""

    path: tuple[str, ...]
    links: tuple[str, ...]
    length: int
    mean_coverage: float
    median_coverage: float
    variance: float
    minimum: float
    maximum: float
    gc: float
    colours: tuple[int, ...]


@dataclass(frozen=True)
class BubbleFeatures:
    """Feature record ``D_B``. Optional evidence stays None when it was not supplied."""

    bubble_id: str
    branches: tuple[BranchFeatures, ...]
    coverage_ratio: float
    source_degree: int
    sink_degree: int
    n_branches: int
    colour_intersection: tuple[int, ...]
    colour_union: tuple[int, ...]
    colour_entropy: float
    colours_equal: bool
    colours_disjoint: bool
    colour_consistent: bool
    linkage: bool | None = None
    low_kmer: bool | None = None


def _branch_coverage(graph: AssemblyGraph, branch: Branch) -> list[float]:
    if branch.path:
        return [graph.node_coverage[unitig_id] for unitig_id in branch.path]
    return [graph.link_coverage[link_id] for link_id in branch.links]


def _entropy(colour_sets: list[tuple[int, ...]]) -> float:
    if not colour_sets:
        return 0.0
    counts: dict[int, int] = {}
    for colours in colour_sets:
        for colour in colours:
            counts[colour] = counts.get(colour, 0) + 1
    total = sum(counts.values())
    if total == 0:
        return 0.0
    value = 0.0
    for count in counts.values():
        probability = count / total
        value -= probability * math.log2(probability)
    return value


def branch_features(graph: AssemblyGraph, branch: Branch) -> BranchFeatures:
    """Summarise one branch. Raises KeyError if coverage was not provided."""
    coverages = _branch_coverage(graph, branch)
    sequences = [graph.unitig(unitig_id).sequence for unitig_id in branch.path]
    colour_set: set[int] = set()
    for unitig_id in branch.path:
        colour_set.update(graph.unitig(unitig_id).color_ids)
    for link_id in branch.links:
        colour_set.update(graph.link(link_id).color_ids)
    mean = _mean(coverages)
    variance = _mean([(value - mean) ** 2 for value in coverages]) if coverages else 0.0
    return BranchFeatures(
        path=branch.path,
        links=branch.links,
        length=sum(len(sequence) for sequence in sequences),
        mean_coverage=mean,
        median_coverage=_median(coverages),
        variance=variance,
        minimum=min(coverages) if coverages else 0.0,
        maximum=max(coverages) if coverages else 0.0,
        gc=_gc("".join(sequences)),
        colours=tuple(sorted(colour_set)),
    )


def extract_features(
    graph: AssemblyGraph,
    bubble: Bubble,
    *,
    linkage: bool | None = None,
    low_kmer: bool | None = None,
) -> BubbleFeatures:
    """Build ``D_B``. ``linkage`` and ``low_kmer`` are used only when the caller has them."""
    branches = tuple(branch_features(graph, branch) for branch in bubble.branches)
    means = [branch.mean_coverage for branch in branches]
    positive = [value for value in means if value > 0]
    ratio = (min(positive) / max(positive)) if len(positive) >= 2 else 1.0
    sets = [set(branch.colours) for branch in branches]
    intersection = set.intersection(*sets) if sets else set()
    union = set.union(*sets) if sets else set()
    outdeg: dict[str, int] = {}
    indeg: dict[str, int] = {}
    for link in graph.cdbg.links:
        outdeg[link.source] = outdeg.get(link.source, 0) + 1
        indeg[link.target] = indeg.get(link.target, 0) + 1
    equal = len({branch.colours for branch in branches}) == 1
    disjoint = len(branches) >= 2 and set.intersection(*sets) == set() and all(sets)
    source_colours = set(graph.unitig(bubble.source).color_ids)
    consistent = all(set(branch.colours) <= source_colours for branch in branches if branch.colours)
    return BubbleFeatures(
        bubble_id=bubble.bubble_id,
        branches=branches,
        coverage_ratio=ratio,
        source_degree=outdeg.get(bubble.source, 0),
        sink_degree=indeg.get(bubble.sink, 0),
        n_branches=len(branches),
        colour_intersection=tuple(sorted(intersection)),
        colour_union=tuple(sorted(union)),
        colour_entropy=_entropy([branch.colours for branch in branches]),
        colours_equal=equal,
        colours_disjoint=disjoint,
        colour_consistent=consistent,
        linkage=linkage,
        low_kmer=low_kmer,
    )
