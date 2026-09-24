"""Detect, classify, and greedily edit a coloured graph."""

from __future__ import annotations

from pathlib import Path

from metametro.formats.cdbg.io import load_cdbg

from bubbleblower.classify import classify_bubble
from bubbleblower.detect import detect_bubbles
from bubbleblower.features import extract_features
from bubbleblower.graph import AssemblyGraph, from_cdbg
from bubbleblower.search import SearchResult, greedy_search


def load_assembly(path: str | Path, node_coverage: dict[str, float] | None = None) -> AssemblyGraph:
    """Load a MetaMetro CDBG directory."""
    return from_cdbg(load_cdbg(path), node_coverage=node_coverage)


def resolve(graph: AssemblyGraph, **kwargs) -> SearchResult:
    """Run greedy debubbling. The input graph is copied inside the search."""
    return greedy_search(graph, **kwargs)


def summarize(graph: AssemblyGraph, *, mode: str = "coverage") -> dict:
    """Classify current bubbles without editing."""
    bubbles = detect_bubbles(graph)
    posteriors = [classify_bubble(extract_features(graph, bubble), mode=mode) for bubble in bubbles]
    return {
        "n_bubbles": len(bubbles),
        "bubbles": [
            {
                "bubble_id": posterior.bubble_id,
                "p_error": posterior.p_error,
                "p_strain": posterior.p_strain,
                "p_other": posterior.p_other,
                "decision": posterior.decision,
                "confidence": posterior.confidence,
            }
            for posterior in posteriors
        ],
    }
