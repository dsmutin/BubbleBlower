"""Named resolution modes. Each mode returns a new graph.

``retain`` leaves the graph unchanged.
``colour_pop`` removes a branch only when its colours equal the other branch
and its coverage is a small fraction of that branch.
``greedy`` is the global score search.
``kmer_divergence`` and ``colour_topology`` are graph-only debubblers.
``read_colour_break`` keeps the graph and splits contigs where the read colour changes.
"""

from __future__ import annotations

from bubbleblower.classify import classify_bubble
from bubbleblower.debubblers import resolve_debubbler
from bubbleblower.detect import detect_bubbles
from bubbleblower.edits import pop_branch
from bubbleblower.features import extract_features
from bubbleblower.graph import AssemblyGraph
from bubbleblower.search import greedy_search


def colour_pop(graph: AssemblyGraph, *, max_ratio: float = 0.15) -> AssemblyGraph:
    """Pop same-colour low-coverage branches. Disjoint colours are kept."""
    current = graph.copy()
    changed = True
    while changed:
        changed = False
        for bubble in detect_bubbles(current):
            features = extract_features(current, bubble)
            if not features.colours_equal or features.coverage_ratio > max_ratio:
                continue
            coverages = [branch.mean_coverage for branch in features.branches]
            index = min(range(len(coverages)), key=lambda item: coverages[item])
            current, _edit = pop_branch(current, bubble, index)
            changed = True
            break
    return current


def resolve_mode(graph: AssemblyGraph, mode: str, **kwargs) -> AssemblyGraph:
    """Apply one named mode. Unknown names raise ``ValueError``."""
    if mode == "retain":
        return graph.copy()
    if mode == "colour_pop":
        return colour_pop(graph, max_ratio=kwargs.get("max_ratio", 0.15))
    if mode == "greedy":
        return greedy_search(
            graph,
            max_iterations=kwargs.get("max_iterations", 30),
            max_runtime_s=kwargs.get("max_runtime_s", 120.0),
            mode=kwargs.get("classifier", "coverage"),
        ).graph
    if mode in {"kmer_divergence", "colour_topology"}:
        return resolve_debubbler(graph, mode)
    raise ValueError(f"unknown mode: {mode}")


def decisions(graph: AssemblyGraph) -> list[dict[str, str | float]]:
    """Classify every current bubble. Does not edit the graph."""
    rows = []
    for bubble in detect_bubbles(graph):
        features = extract_features(graph, bubble)
        posterior = classify_bubble(features)
        rows.append(
            {
                "bubble_id": bubble.bubble_id,
                "source": bubble.source,
                "sink": bubble.sink,
                "decision": posterior.decision,
                "p_error": posterior.p_error,
                "coverage_ratio": features.coverage_ratio,
                "colours_equal": features.colours_equal,
            }
        )
    return rows
