"""Detect, classify, and greedily edit a coloured graph."""

from __future__ import annotations

from pathlib import Path

from metametro.formats.cdbg.io import load_cdbg

from bubbleblower.classify import classify_bubble
from bubbleblower.detect import detect_bubbles
from bubbleblower.features import extract_features
from bubbleblower.graph import AssemblyGraph, from_cdbg
from bubbleblower.search import SearchResult, greedy_search


def _read_coverage(path: Path) -> dict[str, float] | None:
    if not path.is_file():
        return None
    rows: dict[str, float] = {}
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"empty coverage table: {path}")
    for line in lines[1:]:
        ident, value = line.split("\t")
        rows[ident] = float(value)
    return rows


def load_assembly(path: str | Path, node_coverage: dict[str, float] | None = None) -> AssemblyGraph:
    """Load a MetaMetro CDBG directory.

    Optional ``node_coverage.tsv`` and ``link_coverage.tsv`` in that directory
    supply coverage. They are not inferred when absent.
    """
    root = Path(path)
    nodes = node_coverage if node_coverage is not None else _read_coverage(root / "node_coverage.tsv")
    links = _read_coverage(root / "link_coverage.tsv")
    return from_cdbg(load_cdbg(root), node_coverage=nodes, link_coverage=links)


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
