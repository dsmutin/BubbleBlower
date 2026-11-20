"""Mandatory: a pop removes each branch link with its reverse-complement twin (#1203)."""

from __future__ import annotations

import pytest

from bubbleblower.detect import detect_bubbles
from bubbleblower.edits import pop_branch, restore_pop
from bubbleblower.graph import records_to_tocumg

pytestmark = pytest.mark.mandatory

COLORS = [{"color_id": "0", "namespace": "sample", "value": "s1"}]
NODES = [
    {"id": "S", "sequence": "ACGTACGT", "colors": [0], "coverage": 120.0},
    {"id": "A", "sequence": "ATATCG", "colors": [0], "coverage": 30.0},
    {"id": "B", "sequence": "CGCGTA", "colors": [0], "coverage": 90.0},
    {"id": "T", "sequence": "GGCCTTAA", "colors": [0], "coverage": 120.0},
]
FORWARD = [("eSA", "S", "A"), ("eAT", "A", "T"), ("eSB", "S", "B"), ("eBT", "B", "T")]


def _link(link_id: str, source: str, target: str, orientation: str) -> dict:
    return {"id": link_id, "source": source, "target": target, "orientation": orientation, "colors": [0], "coverage": 1.0}


def _graph(links: list[dict]):
    return records_to_tocumg(graph_id="twins", colors=COLORS, nodes=NODES, links=links)


def _member(graph, unitig_id: str) -> str:
    return "".join(graph.unitig(unitig_id).members)


def _no_dangling(graph) -> None:
    unitigs = {unitig.unitig_id for unitig in graph.cdbg.unitigs}
    assert all(link.source in unitigs and link.target in unitigs for link in graph.cdbg.links)


@pytest.mark.parametrize(
    "reverse_orientation",
    ["--", "mixed"],
)
def test_pop_on_double_stranded_graph_leaves_no_dangling_link(reverse_orientation: str) -> None:
    forward = [_link(i, s, t, "++") for i, s, t in FORWARD]
    if reverse_orientation == "--":
        reverse = [_link(f"r{i}", t, s, "--") for i, s, t in FORWARD]
    else:
        # B stored on its minus strand: S+ -> B-, B- -> T+; twins B+ -> S-, T- -> B+.
        forward = [_link("eSA", "S", "A", "++"), _link("eAT", "A", "T", "++"), _link("eSB", "S", "B", "+-"), _link("eBT", "B", "T", "-+")]
        reverse = [_link("reSA", "A", "S", "--"), _link("reAT", "T", "A", "--"), _link("reSB", "B", "S", "+-"), _link("reBT", "T", "B", "-+")]
    graph = _graph(forward + reverse)
    (bubble,) = detect_bubbles(graph)
    for index, branch in enumerate(bubble.branches):
        popped, edit = pop_branch(graph, bubble, index)
        popped.validate()
        _no_dangling(popped)
        assert len(popped.cdbg.links) == len(graph.cdbg.links) - 2 * len(branch.links)
        assert set(popped.link_coverage) == {link.link_id for link in popped.cdbg.links}
        removed = {_member(graph, unitig_id) for unitig_id in branch.path}
        assert not any(_member(popped, unitig.unitig_id) in removed for unitig in popped.cdbg.unitigs)
        restored = restore_pop(popped, edit)
        assert sorted(link.link_id for link in restored.cdbg.links) == sorted(link.link_id for link in graph.cdbg.links)
        assert restored.link_coverage == graph.link_coverage


def test_single_strand_pop_removes_only_branch_links() -> None:
    graph = _graph([_link(i, s, t, "++") for i, s, t in FORWARD])
    (bubble,) = detect_bubbles(graph)
    popped, _edit = pop_branch(graph, bubble, 0)
    kept = {link.link_id for link in graph.cdbg.links} - set(bubble.branches[0].links)
    assert {link.link_id for link in popped.cdbg.links} == kept
