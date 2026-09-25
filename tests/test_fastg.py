"""MEGAHIT FASTG loader keeps forward nodes and finds a simple bubble."""

from __future__ import annotations

from pathlib import Path

import pytest

from bubbleblower.detect import detect_bubbles
from bubbleblower.fastg import load_fastg

pytestmark = pytest.mark.mandatory

# Assembly k=4. Branch sequences share that overlap with the source and the sink.
_FASTG = """\
>NODE_1_length_8_cov_10.0_ID_1:NODE_2_length_12_cov_5.0_ID_3,NODE_3_length_12_cov_5.0_ID_5;
AAAATTTT
>NODE_1_length_8_cov_10.0_ID_1';
AAAATTTT
>NODE_2_length_12_cov_5.0_ID_3:NODE_4_length_8_cov_10.0_ID_7;
TTTTACGTAAAA
>NODE_3_length_12_cov_5.0_ID_5:NODE_4_length_8_cov_10.0_ID_7;
TTTTGCGCAAAA
>NODE_4_length_8_cov_10.0_ID_7;
AAAACCCC
"""


def test_fastg_forward_bubble(tmp_path: Path) -> None:
    """Reverse-strand records are folded in and the forward bubble is detected."""
    path = tmp_path / "k4.fastg"
    path.write_text(_FASTG, encoding="utf-8")
    graph = load_fastg(path, k=4, graph_id="k4")
    assert len(graph.cdbg.unitigs) == 4
    assert graph.node_coverage[graph.unitig_id_for_member("n000001")] == 10.0
    bubbles = detect_bubbles(graph)
    assert len(bubbles) == 1
    assert len(bubbles[0].branches) == 2
