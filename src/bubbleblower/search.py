"""Search over bubble edits: greedy, beam, and a seeded Metropolis chain.

A candidate is accepted only when the global score rises by more than epsilon,
the instance count stays inside the hard limit, and runtime has not expired.
Classifier posteriors adjust pop candidates; they never edit the graph.
"""

from __future__ import annotations

import math
import random
import time
from dataclasses import dataclass, field

from bubbleblower.classify import Posterior, classify_bubble
from bubbleblower.detect import Bubble, detect_bubbles
from bubbleblower.edits import Edit, pop_branch, split_instance
from bubbleblower.features import extract_features
from bubbleblower.graph import AssemblyGraph
from bubbleblower.score import Score, score_graph

# Soft weight: classifier informs the optimiser, it does not pop.
_POSTERIOR_WEIGHT = 6.0


@dataclass
class SearchStep:
    """One accepted or rejected comparison."""

    iteration: int
    edit: Edit | None
    score_before: float
    score_after: float
    accepted: bool
    state_id: int = 0
    parent_state: int = 0


@dataclass
class SearchResult:
    """Best graph, history, and per-bubble posteriors."""

    graph: AssemblyGraph
    edits: list[Edit] = field(default_factory=list)
    steps: list[SearchStep] = field(default_factory=list)
    scores: list[Score] = field(default_factory=list)
    bubbles: list[Bubble] = field(default_factory=list)
    posteriors: list[Posterior] = field(default_factory=list)


def _branch_coverage(graph: AssemblyGraph, bubble: Bubble, index: int) -> float:
    branch = bubble.branches[index]
    if branch.path:
        return sum(graph.node_coverage[unitig_id] for unitig_id in branch.path) / len(branch.path)
    return sum(graph.link_coverage[link_id] for link_id in branch.links) / len(branch.links)


def _candidate_pop(graph: AssemblyGraph, bubble: Bubble) -> tuple[AssemblyGraph, Edit] | None:
    if len(bubble.branches) < 2:
        return None
    coverages = [_branch_coverage(graph, bubble, index) for index in range(len(bubble.branches))]
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
            edited, edit = split_instance(graph, instance_id, groups)
        except ValueError:
            continue
        edit.bubble_id = bubble.bubble_id
        candidates.append((edited, edit))
    return candidates


def _pop_adjust(edit: Edit, posterior: Posterior | None) -> float:
    if edit.edit_type != "pop" or posterior is None:
        return 0.0
    return _POSTERIOR_WEIGHT * (posterior.p_error - posterior.p_strain)


def _classify_map(graph: AssemblyGraph, bubbles: list[Bubble], mode: str) -> dict[tuple[str, str], Posterior]:
    return {
        (bubble.source, bubble.sink): classify_bubble(extract_features(graph, bubble), mode=mode)
        for bubble in bubbles
    }


def _stamp_parents(bubbles: list[Bubble], parent_id: str | None) -> list[Bubble]:
    if not parent_id:
        return bubbles
    stamped = []
    for bubble in bubbles:
        stamped.append(
            Bubble(
                bubble_id=bubble.bubble_id,
                source=bubble.source,
                sink=bubble.sink,
                branches=bubble.branches,
                parent_bubble_id=parent_id,
            )
        )
    return stamped


def _deadline(max_runtime_s: float | None, started: float) -> bool:
    if max_runtime_s is None:
        return False
    return time.monotonic() - started >= max_runtime_s


def greedy_search(
    graph: AssemblyGraph,
    *,
    max_iterations: int = 20,
    min_score_gain: float = 1e-3,
    max_instances: int = 500,
    patience: int = 1,
    mode: str = "coverage",
    max_runtime_s: float | None = 90.0,
) -> SearchResult:
    """Accept the single best improving edit at each iteration."""
    started = time.monotonic()
    current = graph.copy()
    result = SearchResult(graph=current)
    result.scores.append(score_graph(current))
    stalled = 0
    state_id = 0
    for iteration in range(max_iterations):
        if _deadline(max_runtime_s, started):
            break
        bubbles = detect_bubbles(current)
        posteriors = _classify_map(current, bubbles, mode)
        if iteration == 0:
            result.bubbles = bubbles
            result.posteriors = [posteriors[(bubble.source, bubble.sink)] for bubble in bubbles]
        if not bubbles:
            break
        base = result.scores[-1]
        best_graph: AssemblyGraph | None = None
        best_edit: Edit | None = None
        best_score: Score | None = None
        best_delta = min_score_gain
        for bubble in bubbles:
            posterior = posteriors.get((bubble.source, bubble.sink))
            for candidate, edit in propose_edits(current, bubble):
                if len(candidate.cdbg.unitigs) > max_instances:
                    continue
                scored = score_graph(candidate)
                raw = scored.total - base.total
                if raw <= min_score_gain:
                    continue
                delta = raw + _pop_adjust(edit, posterior)
                if delta > best_delta:
                    best_delta = delta
                    best_graph = candidate
                    best_edit = edit
                    best_score = scored
        parent = state_id
        if best_graph is None or best_score is None:
            stalled += 1
            result.steps.append(
                SearchStep(iteration, None, base.total, base.total, False, state_id, parent)
            )
            if stalled >= patience:
                break
            continue
        stalled = 0
        state_id += 1
        result.steps.append(
            SearchStep(iteration, best_edit, base.total, best_score.total, True, state_id, parent)
        )
        if best_edit is not None:
            result.edits.append(best_edit)
            result.bubbles = _stamp_parents(detect_bubbles(best_graph), best_edit.bubble_id)
        result.scores.append(best_score)
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
    max_runtime_s: float | None = 90.0,
) -> SearchResult:
    """Keep the top ``beam_width`` graphs at each iteration."""
    if beam_width < 1:
        raise ValueError("beam_width must be >= 1")
    started = time.monotonic()
    start = graph.copy()
    start_score = score_graph(start)
    next_state = 1
    beam: list[tuple[int, int, AssemblyGraph, list[Edit], list[Score]]] = [
        (0, 0, start, [], [start_score])
    ]
    result = SearchResult(graph=start, scores=[start_score])
    bubbles = detect_bubbles(start)
    result.bubbles = bubbles
    result.posteriors = list(_classify_map(start, bubbles, mode).values())
    for iteration in range(max_iterations):
        if _deadline(max_runtime_s, started):
            break
        pool: list[tuple[float, int, int, AssemblyGraph, list[Edit], list[Score], Edit | None]] = []
        for state_id, parent_id, state, edits, scores in beam:
            base = scores[-1]
            pool.append((base.total, state_id, parent_id, state, edits, scores, None))
            posteriors = _classify_map(state, detect_bubbles(state), mode)
            for bubble in detect_bubbles(state):
                posterior = posteriors.get((bubble.source, bubble.sink))
                for candidate, edit in propose_edits(state, bubble):
                    if len(candidate.cdbg.unitigs) > max_instances:
                        continue
                    scored = score_graph(candidate)
                    raw = scored.total - base.total
                    if raw <= min_score_gain:
                        continue
                    delta = raw + _pop_adjust(edit, posterior)
                    child = next_state
                    next_state += 1
                    pool.append(
                        (scored.total + _pop_adjust(edit, posterior), child, state_id, candidate, [*edits, edit], [*scores, scored], edit)
                    )
        pool.sort(key=lambda item: item[0], reverse=True)
        chosen = pool[:beam_width]
        best_adj, best_id, best_parent, best_graph, best_edits, best_scores, _edit = chosen[0]
        improved = best_scores[-1].total > result.scores[-1].total + min_score_gain
        result.graph = best_graph
        result.edits = best_edits
        result.scores = best_scores
        result.steps.append(
            SearchStep(
                iteration,
                best_edits[-1] if best_edits else None,
                result.scores[0].total,
                best_scores[-1].total,
                improved,
                best_id,
                best_parent,
            )
        )
        if best_edits:
            result.bubbles = _stamp_parents(detect_bubbles(best_graph), best_edits[-1].bubble_id)
        beam = [(item[1], item[2], item[3], item[4], item[5]) for item in chosen]
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
    max_runtime_s: float | None = 90.0,
) -> SearchResult:
    """Metropolis chain. The returned graph is the best scoring state seen."""
    started = time.monotonic()
    rng = random.Random(seed)
    current = graph.copy()
    current_score = score_graph(current)
    result = SearchResult(graph=current.copy(), scores=[current_score])
    result.bubbles = detect_bubbles(current)
    result.posteriors = list(_classify_map(current, result.bubbles, mode).values())
    best = current.copy()
    best_score = current_score.total
    best_edits: list[Edit] = []
    accepted_edits: list[Edit] = []
    state_id = 0
    for iteration in range(n_steps):
        if _deadline(max_runtime_s, started):
            break
        bubbles = detect_bubbles(current)
        posteriors = _classify_map(current, bubbles, mode)
        proposals = []
        for bubble in bubbles:
            for item in propose_edits(current, bubble):
                if len(item[0].cdbg.unitigs) <= max_instances:
                    proposals.append((item[0], item[1], posteriors.get((bubble.source, bubble.sink))))
        if not proposals:
            result.steps.append(
                SearchStep(iteration, None, current_score.total, current_score.total, False, state_id, state_id)
            )
            break
        candidate, edit, posterior = proposals[rng.randrange(len(proposals))]
        proposed = score_graph(candidate)
        raw = proposed.total - current_score.total
        delta = raw + _pop_adjust(edit, posterior)
        accept = raw >= 0 or rng.random() < math.exp(max(delta, -700.0))
        parent = state_id
        if accept:
            state_id += 1
        result.steps.append(
            SearchStep(iteration, edit, current_score.total, proposed.total, accept, state_id, parent)
        )
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
    if best_edits:
        result.bubbles = _stamp_parents(detect_bubbles(best), best_edits[-1].bubble_id)
    return result
