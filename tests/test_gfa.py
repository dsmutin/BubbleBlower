"""GFA load and colour-aware pop mode."""

from __future__ import annotations

from pathlib import Path

import pytest

from bubbleblower.colour_reads import colour_from_fastq
from bubbleblower.detect import detect_bubbles
from bubbleblower.gfa import load_gfa
from bubbleblower.modes import resolve_mode

pytestmark = pytest.mark.mandatory

_GFA = """S\tS\tACGTACGTACGTACGTACGTACGTACGTACGT\tDP:f:100
S\tA\tATATATATATATATATATATATATATATATAT\tDP:f:90
S\tE\tCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG\tDP:f:1
S\tT\tGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAA\tDP:f:91
L\tS\t+\tA\t+\t0M
L\tA\t+\tT\t+\t0M
L\tS\t+\tE\t+\t0M
L\tE\t+\tT\t+\t0M
"""


def _fastq(path: Path) -> None:
    lines = []
    for index in range(4):
        lines.extend(
            [
                f"@GCF_000000001_{index}",
                "ACGTACGTACGTACGTACGTACGTACGTACGTATATATATATATATATATATATATATATATATGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAA",
                "+",
                "I" * 96,
            ]
        )
    lines.extend(
        [
            "@GCF_000000001_err",
            "ACGTACGTACGTACGTACGTACGTACGTACGT" + "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG" + "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAA",
            "+",
            "I" * 96,
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_flye_lowercase_dp_is_coverage(tmp_path: Path) -> None:
    """Flye writes depth as ``dp:i:``, not ``DP``."""
    path = tmp_path / "flye.gfa"
    path.write_text("S\tedge_1\tACGTACGTACGTACGTACGTACGTACGTACGT\tdp:i:8\n", encoding="utf-8")
    graph, _ids = load_gfa(path, graph_id="flye")
    assert graph.node_coverage["edge_1"] == 8.0


def test_gfa_bubble_and_colour_pop(tmp_path: Path) -> None:
    """A same-colour 1x branch is popped; the 90x branch stays."""
    gfa = tmp_path / "graph.gfa"
    gfa.write_text(_GFA, encoding="utf-8")
    graph, _ids = load_gfa(gfa)
    assert len(detect_bubbles(graph)) == 1
    fastq = tmp_path / "reads.fastq"
    _fastq(fastq)
    coloured = colour_from_fastq(graph, [fastq], k=16, min_depth=1)
    assert coloured.cdbg.colors
    resolved = resolve_mode(coloured, "colour_pop", max_ratio=0.2)
    ids = {unitig.unitig_id for unitig in resolved.cdbg.unitigs}
    assert "E" not in ids
    assert "A" in ids
