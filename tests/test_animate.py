"""Pinned-layout animation of iterative splits."""

from __future__ import annotations

from pathlib import Path

import pytest

from bubbleblower.animate import animate_states, colour_partition_states, fork_resolution_states, layout_states
from bubbleblower.fixtures import strain_bubble
from bubbleblower.graph import records_to_tocumg
from bubbleblower.search import greedy_search

pytestmark = pytest.mark.mandatory


def test_greedy_records_a_frame_after_each_accepted_edit() -> None:
    """The optional frame list starts at the input and grows once per edit."""
    frames: list = []
    result = greedy_search(strain_bubble(), max_iterations=5, frames=frames)
    assert len(frames) == len(result.edits) + 1
    assert frames[0].member_ids() == strain_bubble().member_ids()
    assert frames[-1].member_ids() == result.graph.member_ids()


def test_colour_partition_splits_the_strain_bubble(tmp_path: Path) -> None:
    """Different taxon colours on the two branches split the shared source."""
    frames, edits = colour_partition_states(strain_bubble(), namespace="taxon")
    assert edits
    assert edits[0].edit_type == "split"
    assert len(frames) == len(edits) + 1
    layouts = layout_states(frames, edits, seed=0)
    moved = set(edits[0].source_ids) | set(edits[0].target_ids)
    for node_id in layouts[0]:
        if node_id in layouts[1] and node_id not in moved and node_id not in _neighbours(frames[1], moved):
            assert layouts[1][node_id] == layouts[0][node_id]
    gif = animate_states(frames, edits, tmp_path / "strain.gif", namespace="taxon", seconds=5.0)
    assert gif.is_file()
    assert gif.stat().st_size > 0


def test_fork_resolution_splits_every_outgoing_fork() -> None:
    """Two separate out-degree-2 nodes are both split, and none remain."""
    graph = records_to_tocumg(
        graph_id="two_forks",
        colors=[
            {"color_id": "0", "namespace": "type", "value": "residential"},
            {"color_id": "1", "namespace": "type", "value": "cycleway"},
        ],
        nodes=[
            {"id": "S", "sequence": "AAAAAA", "colors": [0], "coverage": 1.0},
            {"id": "A", "sequence": "CCCCCC", "colors": [0], "coverage": 1.0},
            {"id": "B", "sequence": "GGGGGG", "colors": [1], "coverage": 1.0},
            {"id": "P", "sequence": "TTTTTT", "colors": [0], "coverage": 1.0},
            {"id": "C", "sequence": "ACACAC", "colors": [0], "coverage": 1.0},
            {"id": "D", "sequence": "GTGTGT", "colors": [1], "coverage": 1.0},
        ],
        links=[
            {"id": "eSA", "source": "S", "target": "A", "colors": [0], "coverage": 1.0},
            {"id": "eSB", "source": "S", "target": "B", "colors": [1], "coverage": 1.0},
            {"id": "ePC", "source": "P", "target": "C", "colors": [0], "coverage": 1.0},
            {"id": "ePD", "source": "P", "target": "D", "colors": [1], "coverage": 1.0},
        ],
    )
    frames, edits = fork_resolution_states(graph)
    assert len(edits) == 2
    assert len(frames) == 3
    from bubbleblower.animate import _incident_ids

    outgoing, incoming = _incident_ids(frames[-1])
    assert all(len(links) <= 1 for links in outgoing.values())
    assert all(len(links) <= 1 for links in incoming.values())


def test_fork_resolution_splits_an_incoming_junction() -> None:
    """A node with two incoming links and one outgoing link is split too."""
    graph = records_to_tocumg(
        graph_id="in_fork",
        colors=[{"color_id": "0", "namespace": "type", "value": "residential"}],
        nodes=[
            {"id": "A", "sequence": "AAAAAA", "colors": [0], "coverage": 1.0},
            {"id": "B", "sequence": "CCCCCC", "colors": [0], "coverage": 1.0},
            {"id": "S", "sequence": "GGGGGG", "colors": [0], "coverage": 1.0},
            {"id": "T", "sequence": "TTTTTT", "colors": [0], "coverage": 1.0},
        ],
        links=[
            {"id": "eAS", "source": "A", "target": "S", "colors": [0], "coverage": 1.0},
            {"id": "eBS", "source": "B", "target": "S", "colors": [0], "coverage": 1.0},
            {"id": "eST", "source": "S", "target": "T", "colors": [0], "coverage": 1.0},
        ],
    )
    frames, edits = fork_resolution_states(graph)
    assert len(edits) == 1
    assert len(frames) == 2
    from bubbleblower.animate import _incident_ids

    outgoing, incoming = _incident_ids(frames[-1])
    assert all(len(links) <= 1 for links in outgoing.values())
    assert all(len(links) <= 1 for links in incoming.values())


def test_relax_leaves_a_distant_node_fixed() -> None:
    """A node outside the free set keeps its coordinate."""
    from bubbleblower.animate import _relax

    positions = {"a": (0.0, 0.0), "b": (1.0, 0.0), "c": (3.0, 0.0)}
    relaxed = _relax(positions, [("a", "b"), ("b", "c")], {"a"})
    assert relaxed["c"] == positions["c"]
    assert relaxed["b"] == positions["b"]


def _neighbours(graph, seeds: set[str]) -> set[str]:
    found: set[str] = set()
    for link in graph.cdbg.links:
        if link.source in seeds or link.target in seeds:
            found.add(link.source)
            found.add(link.target)
    return found
