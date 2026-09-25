"""MetaMetro CDBG fixtures stay valid inputs."""

from __future__ import annotations

import pytest
from metametro.converters.cdbg_to_cfa import cdbg_to_cfa
from metametro.converters.cdbg_to_cgt import cdbg_to_cgt
from metametro.converters.cfa_to_cdbg import cfa_to_cdbg
from metametro.fixtures import mock_cdbg, mock_cfa

from bubbleblower.detect import detect_bubbles
from bubbleblower.graph import as_cgt, from_cdbg, from_cgt

pytestmark = pytest.mark.mandatory


def test_mock_cfa_round_trip_through_bubbleblower() -> None:
    """CFA to CDBG to CFA keeps sequences. Detection does not require coverage."""
    cfa = mock_cfa()
    cdbg = cfa_to_cdbg(cfa)
    graph = from_cdbg(cdbg)
    detect_bubbles(graph)
    restored = cdbg_to_cfa(graph.cdbg)
    assert restored.sequences == cfa.sequences


def test_cgt_is_a_view_of_the_tocumg() -> None:
    """A tensor does not become a second graph, and colours stay off its features."""
    cdbg = mock_cdbg()
    cgt = cdbg_to_cgt(cdbg)
    wrapped = from_cgt(cgt, cdbg)
    assert len(wrapped.cdbg.unitigs) == len(cdbg.unitigs)
    viewed = as_cgt(wrapped)
    assert viewed.node_features.shape[1] == 0
    assert viewed.node_colors.shape == cgt.node_colors.shape
    assert viewed.node_colors.tolist() == cgt.node_colors.tolist()
