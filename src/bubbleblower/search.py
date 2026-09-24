"""Greedy search over bubble edits.

A candidate is accepted only when the global score rises by more than epsilon
and the instance count stays inside the hard limit.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from bubbleblower.classify import Posterior, classify_bubble
from bubbleblower.detect import Bubble, detect_bubbles
from bubbleblower.edits import Edit, pop_branch, split_instance
from bubbleblower.features import extract_features
from bubbleblower.graph import AssemblyGraph
from bubbleblower.score import Score, score_graph


@dataclass
class SearchStep:
    """One accepted or rejected comparison."""

    iteration: int
    edit: Edit | None
    score_before: float
    score_after: float
    accepted: bool


@dataclass
class SearchResult:
    """Best graph, history, and per-bubble posteriors from the initial graph."""

    graph: AssemblyGraph
    edits: list[Edit] = field(default_factory=list)
    steps: list[SearchStep] = field(default_factory=list)
    scores: list[Score] = field(default_factory=list)
    bubbles: list[Bubble] = field(default_factory=list)
    posteriors: list[Posterior] = field(default_factory=list)


def _candidate_pop(graph: AssemblyGraph, bubble: Bubble) -> AssemblyGraph | None:
    if len(bubble.branches) < 2:
        return None
    coverages = []
    for branch in bubble.branches:
        if branch.path:
            coverages.append(sum(graph.node_coverage[unitig_id] for unitig_id in branch.path) / len(branch.path))
        else:
            coverages.append(sum(graph.link_coverage[link_id] for link_id in branch.links) / len(branch.links))
    index = min(range(len(coverages)), key=lambda item: coverages[item])
    if coverages[index] <= 0:
        return None
    edited, _edit = pop_branch(graph, bubble, index)
    return edited


def _groups_for(bubble: Bubble, end: str) -> list[list[str]] | None:
    groups: list[list[str]] = []
    for branch in bubble.branches:
        link_id = branch.links[0] if end == "source" else branch.links[-1]
        groups.append([link_id])
    if len(groups) < 2:
        return None
    return groups


def propose_graphs(graph: AssemblyGraph, bubble: Bubble) -> list[AssemblyGraph]:
    """Pop the weakest branch, or split the source, or split the sink."""
    candidates = []
    popped = _candidate_pop(graph, bubble)
    if popped is not None:
        candidates.append(popped)
    for end in ("source", "sink"):
        groups = _groups_for(bubble, end)
        if groups is None:
            continue
        instance_id = bubble.source if end == "source" else bubble.sink
        try:
            edited, _edit = split_instance(graph, instance_id, groups)
        except ValueError:
            continue
        candidates.append(edited)
    return candidates


def greedy_search(
    graph: AssemblyGraph,
    *,
    max_iterations: int = 20,
    min_score_gain: float = 1e-3,
    max_instances: int = 500,
    patience: int = 1,
    mode: str = "coverage",
) -> SearchResult:
    """Accept the single best improving edit at each iteration."""
    current = graph.copy()
    result = SearchResult(graph=current)
    result.scores.append(score_graph(current))
    stalled = 0
    for iteration in range(max_iterations):
        bubbles = detect_bubbles(current)
        if iteration == 0:
            result.bubbles = bubbles
            result.posteriors = [
                classify_bubble(extract_features(current, bubble), mode=mode) for bubble in bubbles
            ]
        if not bubbles:
            break
        base = score_graph(current)
        best_graph: AssemblyGraph | None = None
        best_delta = min_score_gain
        for bubble in bubbles:
            for candidate in propose_graphs(current, bubble):
                if len(candidate.cdbg.unitigs) > max_instances:
                    continue
                delta = score_graph(candidate).total - base.total
                if delta > best_delta:
                    best_delta = delta
                    best_graph = candidate
        if best_graph is None:
            stalled += 1
            result.steps.append(
                SearchStep(iteration, None, base.total, base.total, False)
            )
            if stalled >= patience:
                break
            continue
        stalled = 0
        after = score_graph(best_graph)
        result.steps.append(SearchStep(iteration, None, base.total, after.total, True))
        result.scores.append(after)
        current = best_graph
        result.graph = current
    return result
