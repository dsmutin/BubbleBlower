"""MEGAHIT FASTG loader keeps forward nodes and finds a simple bubble."""

from __future__ import annotations

from pathlib import Path

import pytest

from bubbleblower.detect import detect_bubbles
from bubbleblower.fastg import load_fastg

pytestmark = pytest.mark.mandatory

_FASTG = """\
>NODE_1_length_8_cov_10.0_ID_1:NODE_2_length_8_cov_5.0_ID_3,NODE_3_length_8_cov_5.0_ID_5;
ACGTACGT
>NODE_1_length_8_cov_10.0_ID_1';
ACGTACGT
>NODE_2_length_8_cov_5.0_ID_3:NODE_4_length_8_cov_10.0_ID_7;
ATATATAT
>NODE_3_length_8_cov_5.0_ID_5:NODE_4_length_8_cov_10.0_ID_7;
CGCGCGCG
>NODE_4_length_8_cov_10.0_ID_7;
GGCCTTAA
"""


def test_fastg_forward_bubble(tmp_path: Path) -> None:
    """Reverse-strand records are dropped and the forward bubble is detected."""
    path = tmp_path / "k21.fastg"
    path.write_text(_FASTG, encoding="utf-8")
    graph = load_fastg(path, graph_id="k21")
    assert len(graph.cdbg.unitigs) == 4
    assert graph.node_coverage["NODE_1_length_8_cov_10.0_ID_1"] == 10.0
    bubbles = detect_bubbles(graph)
    assert len(bubbles) == 1
    assert len(bubbles[0].branches) == 2
