"""Global graph score.

Coverage of instances should look like a small mixture. Complexity is penalised
so extra instances are not free.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from bubbleblower.detect import detect_bubbles
from bubbleblower.graph import AssemblyGraph


@dataclass(frozen=True)
class Score:
    """Weighted score and its parts. Higher is better."""

    total: float
    coverage: float
    flow: float
    topology: float
    complexity: float
    n_bubbles: int
    n_instances: int


def _means(values: list[float], k: int) -> list[float]:
    ordered = sorted(values)
    means = [ordered[min(len(ordered) - 1, int((index + 0.5) * len(ordered) / k))] for index in range(k)]
    for _ in range(12):
        clusters: list[list[float]] = [[] for _ in range(k)]
        for value in ordered:
            nearest = min(range(k), key=lambda index: abs(value - means[index]))
            clusters[nearest].append(value)
        means = [sum(cluster) / len(cluster) if cluster else means[index] for index, cluster in enumerate(clusters)]
    return means


def _bic(values: list[float], k: int) -> float:
    n = len(values)
    means = _means(values, k)
    variance = 0.0
    for value in values:
        center = min(means, key=lambda item: abs(value - item))
        variance += (value - center) ** 2
    variance = variance / n + 1.0
    nll = 0.0
    for value in values:
        center = min(means, key=lambda item: abs(value - item))
        nll += 0.5 * math.log(2.0 * math.pi * variance) + (value - center) ** 2 / (2.0 * variance)
    params = k + 1
    return nll + 0.5 * params * math.log(n)


def coverage_score(values: list[float], k_max: int = 3) -> float:
    """Negative BIC of the best Gaussian mixture with ``K <= k_max``."""
    if not values:
        return 0.0
    k_limit = min(k_max, len(values))
    return -min(_bic(values, k) for k in range(1, k_limit + 1))


def flow_penalty(graph: AssemblyGraph) -> float:
    """Squared mismatch between a unitig's coverage and its outgoing links."""
    outgoing: dict[str, float] = {}
    incoming: dict[str, float] = {}
    for link in graph.cdbg.links:
        outgoing[link.source] = outgoing.get(link.source, 0.0) + graph.link_coverage[link.link_id]
        incoming[link.target] = incoming.get(link.target, 0.0) + graph.link_coverage[link.link_id]
    penalty = 0.0
    for unitig_id, coverage in graph.node_coverage.items():
        scale = abs(coverage) + 1.0
        if unitig_id in outgoing:
            penalty += (coverage - outgoing[unitig_id]) ** 2 / scale
        if unitig_id in incoming:
            penalty += (coverage - incoming[unitig_id]) ** 2 / scale
    return penalty


def score_graph(graph: AssemblyGraph, *, k_max: int = 3) -> Score:
    """``S(G)``. Read support is omitted until read evidence is attached."""
    bubbles = detect_bubbles(graph)
    n_instances = len(graph.cdbg.unitigs)
    n_links = len(graph.cdbg.links)
    coverage = coverage_score(list(graph.node_coverage.values()), k_max=k_max)
    flow = -flow_penalty(graph)
    topology = -0.05 * len(bubbles)
    complexity = -(0.15 * n_instances + 0.02 * n_links)
    total = coverage + 0.35 * flow + topology + complexity
    return Score(
        total=total,
        coverage=coverage,
        flow=flow,
        topology=topology,
        complexity=complexity,
        n_bubbles=len(bubbles),
        n_instances=n_instances,
    )
