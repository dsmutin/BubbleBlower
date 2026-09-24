"""Contig walks respect link orientation and overlap."""

from __future__ import annotations

import pytest

from bubbleblower.contigs import contig_sequences
from bubbleblower.graph import build_graph

pytestmark = pytest.mark.mandatory


def _graph(links: list[dict]) -> object:
    return build_graph(
        graph_id="walk",
        colors=[{"color_id": "0", "namespace": "taxon", "value": "taxon_1"}],
        nodes=[
            {"id": "A", "sequence": "AAAA", "colors": [0], "coverage": 10.0},
            {"id": "B", "sequence": "CCCC", "colors": [0], "coverage": 10.0},
        ],
        links=links,
    )


def test_forward_link_consumes_overlap() -> None:
    """A ++ link appends the target minus the overlap."""
    graph = _graph(
        [{"id": "e", "source": "A", "target": "B", "orientation": "++", "colors": [0], "coverage": 10.0}]
    )
    graph.link("e").overlap = 1
    sequences = [sequence for _name, sequence in contig_sequences(graph)]
    assert sequences == ["AAAACCC"]


def test_reverse_target_is_complemented() -> None:
    """A +- link appends the reverse complement of the target."""
    graph = _graph(
        [{"id": "e", "source": "A", "target": "B", "orientation": "+-", "colors": [0], "coverage": 10.0}]
    )
    graph.link("e").overlap = 0
    sequences = [sequence for _name, sequence in contig_sequences(graph)]
    assert sequences == ["AAAAGGGG"]


def test_missing_overlap_is_refused() -> None:
    """A join without a stored overlap is an error."""
    graph = _graph(
        [{"id": "e", "source": "A", "target": "B", "orientation": "++", "colors": [0], "coverage": 10.0}]
    )
    with pytest.raises(ValueError, match="no overlap"):
        contig_sequences(graph)
