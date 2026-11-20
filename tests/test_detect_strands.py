"""Mandatory: bubble detection reads links as arcs between oriented unitigs."""

from __future__ import annotations

import pytest

from bubbleblower.detect import detect_bubbles
from bubbleblower.fixtures import strain_bubble
from bubbleblower.graph import records_to_tocumg

pytestmark = pytest.mark.mandatory

NODES = [
    {"id": "S", "sequence": "ACGTACGT", "colors": [0, 1], "coverage": 120.0},
    {"id": "A", "sequence": "ATATCG", "colors": [0], "coverage": 30.0},
    {"id": "B", "sequence": "CGCGTA", "colors": [1], "coverage": 90.0},
    {"id": "T", "sequence": "GGCCTTAA", "colors": [0, 1], "coverage": 120.0},
]
COLORS = [
    {"color_id": "0", "namespace": "sample", "value": "s1"},
    {"color_id": "1", "namespace": "sample", "value": "s2"},
]
FORWARD = [("eSA", "S", "A"), ("eAT", "A", "T"), ("eSB", "S", "B"), ("eBT", "B", "T")]


def _graph(links: list[dict]):
    return records_to_tocumg(graph_id="strands", colors=COLORS, nodes=NODES, links=links)


def _name(graph, unitig_id: str) -> str:
    return "".join(graph.unitig(unitig_id).members)


def _link(link_id: str, source: str, target: str, orientation: str) -> dict:
    return {"id": link_id, "source": source, "target": target, "orientation": orientation, "colors": [0], "coverage": 1.0}


def test_both_strands_give_one_bubble_with_forward_signature() -> None:
    forward = [_link(i, s, t, "++") for i, s, t in FORWARD]
    reverse = [_link(f"r{i}", t, s, "--") for i, s, t in FORWARD]
    expected = detect_bubbles(_graph(forward))
    both = _graph(forward + reverse)
    found = detect_bubbles(both)
    assert len(expected) == 1
    assert len(found) == 1
    assert found[0].bubble_id == min(expected[0].bubble_id, detect_bubbles(_graph(reverse))[0].bubble_id)
    assert {_name(both, found[0].source), _name(both, found[0].sink)} == {"S", "T"}
    assert sorted(_name(both, branch.path[0]) for branch in found[0].branches) == ["A", "B"]


def test_reverse_strand_only_is_found_with_strands() -> None:
    reverse = [_link(f"r{i}", t, s, "--") for i, s, t in FORWARD]
    graph = _graph(reverse)
    found = detect_bubbles(graph)
    assert len(found) == 1
    bubble = found[0]
    assert (_name(graph, bubble.source), bubble.source_strand) == ("T", "-")
    assert (_name(graph, bubble.sink), bubble.sink_strand) == ("S", "-")
    assert all(branch.strands == ("-",) for branch in bubble.branches)


def test_forward_fixture_unchanged() -> None:
    graph = strain_bubble()
    (bubble,) = detect_bubbles(graph)
    assert (bubble.source_strand, bubble.sink_strand) == ("+", "+")
    assert all(branch.strands == ("+",) for branch in bubble.branches)


def test_ids_match_metametro_strand_bubbles() -> None:
    bubble_gold = pytest.importorskip("metametro.bench.bubble_gold")
    forward = [_link(i, s, t, "++") for i, s, t in FORWARD]
    reverse = [_link(f"r{i}", t, s, "--") for i, s, t in FORWARD]
    for links in (forward, reverse, forward + reverse):
        graph = _graph(links)
        ours = [bubble.bubble_id for bubble in detect_bubbles(graph)]
        theirs = [bubble.bubble_id for bubble in bubble_gold.strand_bubbles(graph.cdbg)]
        assert ours == theirs
