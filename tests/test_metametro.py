"""MetaMetro CDBG fixtures stay valid inputs."""

from __future__ import annotations

import pytest
from metametro.converters.cdbg_to_cfa import cdbg_to_cfa
from metametro.converters.cfa_to_cdbg import cfa_to_cdbg
from metametro.fixtures import mock_cfa

from bubbleblower.detect import detect_bubbles
from bubbleblower.graph import from_cdbg

pytestmark = pytest.mark.mandatory


def test_mock_cfa_round_trip_through_bubbleblower() -> None:
    """CFA to CDBG to CFA keeps sequences. Detection does not require coverage."""
    cfa = mock_cfa()
    cdbg = cfa_to_cdbg(cfa)
    graph = from_cdbg(cdbg)
    detect_bubbles(graph)
    restored = cdbg_to_cfa(graph.cdbg)
    assert restored.sequences == cfa.sequences
