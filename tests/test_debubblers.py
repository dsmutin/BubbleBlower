"""Mandatory graph-only debubblers: error pop, strain retain, multi split."""

from __future__ import annotations

import pytest

from bubbleblower.debubblers import classify_debubble, resolve_debubbler
from bubbleblower.detect import detect_bubbles
from bubbleblower.fixtures import error_bubble, strain_bubble
from bubbleblower.graph import build_graph
from bubbleblower.modes import resolve_mode

pytestmark = pytest.mark.mandatory

_NAMES = ("kmer_divergence", "colour_topology")


def _three_branch():
    """Three equal-colour branches that share one source and one sink."""
    return build_graph(
        graph_id="three_branch",
        colors=[{"color_id": "0", "namespace": "taxon", "value": "taxon_1"}],
        nodes=[
            {"id": "S", "sequence": "ACGTACGT", "colors": [0], "coverage": 150.0},
            {"id": "A", "sequence": "AAAAAAAA", "colors": [0], "coverage": 50.0},
            {"id": "B", "sequence": "CCCCCCCC", "colors": [0], "coverage": 50.0},
            {"id": "C", "sequence": "GGGGGGGG", "colors": [0], "coverage": 50.0},
            {"id": "T", "sequence": "GGCCTTAA", "colors": [0], "coverage": 150.0},
        ],
        links=[
            {"id": "eSA", "source": "S", "target": "A", "colors": [0], "coverage": 50.0},
            {"id": "eAT", "source": "A", "target": "T", "colors": [0], "coverage": 50.0},
            {"id": "eSB", "source": "S", "target": "B", "colors": [0], "coverage": 50.0},
            {"id": "eBT", "source": "B", "target": "T", "colors": [0], "coverage": 50.0},
            {"id": "eSC", "source": "S", "target": "C", "colors": [0], "coverage": 50.0},
            {"id": "eCT", "source": "C", "target": "T", "colors": [0], "coverage": 50.0},
        ],
    )


def test_error_bubble_pops_weak_branch() -> None:
    """Both models remove the 1x branch and keep the 90x branch."""
    graph = error_bubble()
    for name in _NAMES:
        decision = classify_debubble(graph, detect_bubbles(graph)[0], name)
        assert decision.label == "error"
        assert decision.action == "pop"
        resolved = resolve_mode(graph, name)
        ids = {unitig.unitig_id for unitig in resolved.cdbg.unitigs}
        assert "E" not in ids
        assert "A" in ids
        resolved.validate()
    assert any(unitig.unitig_id == "E" for unitig in graph.cdbg.unitigs)


def test_strain_bubble_keeps_both_alleles() -> None:
    """30x and 90x disjoint alleles are variation. Neither model pops."""
    graph = strain_bubble()
    for name in _NAMES:
        decision = classify_debubble(graph, detect_bubbles(graph)[0], name)
        assert decision.label == "variation"
        assert decision.action == "retain"
        resolved = resolve_debubbler(graph, name)
        ids = {unitig.unitig_id for unitig in resolved.cdbg.unitigs}
        assert {"A", "B", "S", "T"} <= ids
        sequences = {unitig.sequence for unitig in resolved.cdbg.unitigs}
        assert "ATATCG" in sequences
        assert "CGCGTA" in sequences
        resolved.validate()


def test_three_branch_bubble_is_multi() -> None:
    """Three branches are multi. The edited graph still validates."""
    graph = _three_branch()
    bubble = detect_bubbles(graph)[0]
    assert len(bubble.branches) == 3
    for name in _NAMES:
        decision = classify_debubble(graph, bubble, name)
        assert decision.label == "multi"
        resolved = resolve_debubbler(graph, name)
        resolved.validate()
