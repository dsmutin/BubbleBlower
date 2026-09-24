"""Global graph score.

Coverage of remaining instances may look like a small mixture. That term is
intentionally weak: fitting variance by deleting real strain mass is the
failure the contract forbids. Unexplained flow and stranded colours dominate.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from bubbleblower.graph import AssemblyGraph

# Sequencing-error leftover is a small fraction of parent coverage.
_ERROR_FRACTION = 0.15


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
    colour: float = 0.0


def _means(values: list[float], k: int) -> list[float]:
    """Deterministic k-means centres. Initialisation is the sorted quantiles."""
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


def _incident(graph: AssemblyGraph) -> tuple[dict[str, list], dict[str, list]]:
    outgoing: dict[str, list] = {unitig.unitig_id: [] for unitig in graph.cdbg.unitigs}
    incoming: dict[str, list] = {unitig.unitig_id: [] for unitig in graph.cdbg.unitigs}
    for link in graph.cdbg.links:
        outgoing[link.source].append(link)
        incoming[link.target].append(link)
    return outgoing, incoming


def residual_likelihood(graph: AssemblyGraph) -> float:
    """Score leftover flow after comparing a unitig to its incident links.

    A leftover no larger than ``_ERROR_FRACTION`` of the unitig is treated as
    sequencing error (small bonus). A larger leftover is unexplained biological
    mass and is subtracted in full.
    """
    outgoing, incoming = _incident(graph)
    score = 0.0
    for unitig_id, coverage in graph.node_coverage.items():
        for links in (outgoing[unitig_id], incoming[unitig_id]):
            if not links:
                continue
            observed = sum(graph.link_coverage[link.link_id] for link in links)
            residual = abs(coverage - observed)
            budget = _ERROR_FRACTION * coverage + 1.0
            if residual <= 1e-9:
                score += 0.25
            elif residual <= budget:
                score += 1.0
            else:
                score -= residual
    return score


def colour_consistency(graph: AssemblyGraph) -> float:
    """Penalise colours that sit on a unitig but on no incident link.

    Popping the only branch of a taxon leaves that taxon's colour stranded on
    the source and sink. An error branch shares the parent colour, so popping
    it does not strand a colour.
    """
    outgoing, incoming = _incident(graph)
    penalty = 0.0
    for unitig in graph.cdbg.unitigs:
        node_colours = set(unitig.color_ids)
        if not node_colours:
            continue
        coverage = graph.node_coverage[unitig.unitig_id]
        for links in (outgoing[unitig.unitig_id], incoming[unitig.unitig_id]):
            if not links:
                continue
            link_colours: set[int] = set()
            for link in links:
                link_colours.update(link.color_ids)
            missing = node_colours - link_colours
            if missing:
                penalty += coverage * len(missing)
    return -penalty


def score_graph(graph: AssemblyGraph, *, k_max: int = 3) -> Score:
    """``S(G)``. Read support is omitted until read evidence is attached.

    The coverage mixture is down-weighted so deleting a low-abundance strain
    cannot improve the total. Colour and residual terms decide pops.
    """
    outgoing, _incoming = _incident(graph)
    n_instances = len(graph.cdbg.unitigs)
    n_links = len(graph.cdbg.links)
    n_forks = sum(1 for links in outgoing.values() if len(links) >= 2)
    coverage = coverage_score(list(graph.node_coverage.values()), k_max=k_max)
    flow = residual_likelihood(graph)
    colour = colour_consistency(graph)
    topology = -0.05 * n_forks
    complexity = -(0.15 * n_instances + 0.02 * n_links)
    total = 0.05 * coverage + flow + colour + topology + complexity
    return Score(
        total=total,
        coverage=coverage,
        flow=flow,
        topology=topology,
        complexity=complexity,
        n_bubbles=n_forks,
        n_instances=n_instances,
        colour=colour,
    )
