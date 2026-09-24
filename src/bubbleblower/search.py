"""Greedy search over bubble edits.

A candidate is accepted only when the global score rises by more than epsilon
and the instance count stays inside the hard limit.
"""

from __future__ import annotations

import math
import random
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


def _candidate_pop(graph: AssemblyGraph, bubble: Bubble) -> tuple[AssemblyGraph, Edit] | None:
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
    return pop_branch(graph, bubble, index)


def _groups_for(bubble: Bubble, end: str) -> list[list[str]] | None:
    groups: list[list[str]] = []
    for branch in bubble.branches:
        link_id = branch.links[0] if end == "source" else branch.links[-1]
        groups.append([link_id])
    if len(groups) < 2:
        return None
    return groups


def propose_edits(graph: AssemblyGraph, bubble: Bubble) -> list[tuple[AssemblyGraph, Edit]]:
    """Pop the weakest branch, or split the source, or split the sink."""
    candidates: list[tuple[AssemblyGraph, Edit]] = []
    popped = _candidate_pop(graph, bubble)
    if popped is not None:
        candidates.append(popped)
    for end in ("source", "sink"):
        groups = _groups_for(bubble, end)
        if groups is None:
            continue
        instance_id = bubble.source if end == "source" else bubble.sink
        try:
            candidates.append(split_instance(graph, instance_id, groups))
        except ValueError:
            continue
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
        best_edit: Edit | None = None
        best_delta = min_score_gain
        for bubble in bubbles:
            for candidate, edit in propose_edits(current, bubble):
                if len(candidate.cdbg.unitigs) > max_instances:
                    continue
                delta = score_graph(candidate).total - base.total
                if delta > best_delta:
                    best_delta = delta
                    best_graph = candidate
                    best_edit = edit
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
        result.steps.append(SearchStep(iteration, best_edit, base.total, after.total, True))
        if best_edit is not None:
            result.edits.append(best_edit)
        result.scores.append(after)
        current = best_graph
        result.graph = current
    return result


def beam_search(
    graph: AssemblyGraph,
    *,
    beam_width: int = 4,
    max_iterations: int = 8,
    min_score_gain: float = 1e-3,
    max_instances: int = 500,
    mode: str = "coverage",
) -> SearchResult:
    """Keep the top ``beam_width`` graphs at each iteration."""
    if beam_width < 1:
        raise ValueError("beam_width must be >= 1")
    start = graph.copy()
    beam: list[tuple[AssemblyGraph, list[Edit], list[Score]]] = [(start, [], [score_graph(start)])]
    result = SearchResult(graph=start, scores=[beam[0][2][0]])
    bubbles = detect_bubbles(start)
    result.bubbles = bubbles
    result.posteriors = [classify_bubble(extract_features(start, bubble), mode=mode) for bubble in bubbles]
    for iteration in range(max_iterations):
        pool: list[tuple[float, AssemblyGraph, list[Edit], list[Score], Edit | None, float]] = []
        for state, edits, scores in beam:
            base = scores[-1]
            pool.append((base.total, state, edits, scores, None, base.total))
            for bubble in detect_bubbles(state):
                for candidate, edit in propose_edits(state, bubble):
                    if len(candidate.cdbg.unitigs) > max_instances:
                        continue
                    scored = score_graph(candidate)
                    if scored.total - base.total <= min_score_gain:
                        continue
                    pool.append((scored.total, candidate, [*edits, edit], [*scores, scored], edit, base.total))
        pool.sort(key=lambda item: item[0], reverse=True)
        chosen = pool[:beam_width]
        best_total, best_graph, best_edits, best_scores, _edit, _before = chosen[0]
        improved = best_total > result.scores[-1].total + min_score_gain
        result.graph = best_graph
        result.edits = best_edits
        result.scores = best_scores
        result.steps.append(
            SearchStep(iteration, best_edits[-1] if best_edits else None, result.scores[0].total, best_total, improved)
        )
        beam = [(item[1], item[2], item[3]) for item in chosen]
        if not improved:
            break
    return result


def mcmc_search(
    graph: AssemblyGraph,
    *,
    n_steps: int = 20,
    seed: int = 1,
    max_instances: int = 500,
    mode: str = "coverage",
) -> SearchResult:
    """Metropolis chain. The returned graph is the best scoring state seen."""
    rng = random.Random(seed)
    current = graph.copy()
    current_score = score_graph(current)
    result = SearchResult(graph=current.copy(), scores=[current_score])
    result.bubbles = detect_bubbles(current)
    result.posteriors = [
        classify_bubble(extract_features(current, bubble), mode=mode) for bubble in result.bubbles
    ]
    best = current.copy()
    best_score = current_score.total
    best_edits: list[Edit] = []
    accepted_edits: list[Edit] = []
    for iteration in range(n_steps):
        bubbles = detect_bubbles(current)
        proposals = [item for bubble in bubbles for item in propose_edits(current, bubble)]
        proposals = [item for item in proposals if len(item[0].cdbg.unitigs) <= max_instances]
        if not proposals:
            result.steps.append(SearchStep(iteration, None, current_score.total, current_score.total, False))
            break
        candidate, edit = proposals[rng.randrange(len(proposals))]
        proposed = score_graph(candidate)
        delta = proposed.total - current_score.total
        accept = delta >= 0 or rng.random() < math.exp(max(delta, -700.0))
        result.steps.append(SearchStep(iteration, edit, current_score.total, proposed.total, accept))
        if not accept:
            continue
        current = candidate
        current_score = proposed
        accepted_edits.append(edit)
        result.scores.append(proposed)
        if proposed.total > best_score:
            best = current.copy()
            best_score = proposed.total
            best_edits = list(accepted_edits)
    result.graph = best
    result.edits = best_edits
    return result
