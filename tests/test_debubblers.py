"""Mandatory graph-only debubblers: error pop, strain retain, multi split."""

from __future__ import annotations

import pytest

from bubbleblower.colour_break import break_read_colour_chimeras
from bubbleblower.debubblers import classify_debubble, compact_same_colour, resolve_debubbler
from bubbleblower.edits import merge_adjacent, revert
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


def test_read_colour_break_splits_a_chimera(tmp_path) -> None:
    """A contig painted by two accessions is emitted as two pieces."""
    left = "A" * 24
    right = "C" * 24
    chimera = left + right
    fastq = tmp_path / "reads.fastq"
    lines = []
    for index in range(3):
        lines.extend([f"@GCF_000000001_{index}", left + "T" * 8, "+", "I" * 32])
        lines.extend([f"@GCF_000000002_{index}", "G" * 8 + right, "+", "I" * 32])
    fastq.write_text("\n".join(lines) + "\n", encoding="utf-8")
    broken = break_read_colour_chimeras(
        [("chimera", chimera), ("plain", left + "T" * 8)],
        [fastq],
        k=8,
        min_run=4,
        min_piece=10,
    )
    pieces = [sequence for name, sequence in broken if name.startswith("chimera")]
    assert len(pieces) >= 2
    assert any(set(piece) == {"A"} for piece in pieces)
    assert any(set(piece) == {"C"} for piece in pieces)
    plain = [sequence for name, sequence in broken if name == "plain"]
    assert plain == [left + "T" * 8]


def test_read_colour_break_reads_badread_accession(tmp_path) -> None:
    """Badread puts the genome after the read uuid, not in the first token."""
    left = "A" * 24
    right = "C" * 24
    fastq = tmp_path / "ont.fastq"
    lines = []
    for index in range(3):
        lines.extend(
            [
                f"@uuid-{index} GCF_000000001|0,+strand,1-40 length=40",
                left + "T" * 8,
                "+",
                "I" * 32,
            ]
        )
        lines.extend(
            [
                f"@uuid-{index}b GCF_000000002|0,+strand,1-40 length=40",
                "G" * 8 + right,
                "+",
                "I" * 32,
            ]
        )
    fastq.write_text("\n".join(lines) + "\n", encoding="utf-8")
    broken = break_read_colour_chimeras(
        [("chimera", left + right)],
        [fastq],
        k=8,
        min_run=4,
        min_piece=10,
    )
    pieces = [sequence for _name, sequence in broken]
    assert len(pieces) >= 2
    assert any(set(piece) == {"A"} for piece in pieces)
    assert any(set(piece) == {"C"} for piece in pieces)


def test_adjacent_same_colour_nodes_merge() -> None:
    """A simple same-colour path collapses. Different colours stay apart."""
    graph = build_graph(
        graph_id="linear",
        colors=[
            {"color_id": "0", "namespace": "taxon", "value": "taxon_1"},
            {"color_id": "1", "namespace": "taxon", "value": "taxon_2"},
        ],
        nodes=[
            {"id": "A", "sequence": "ACGTACGT", "colors": [0], "coverage": 10.0},
            {"id": "B", "sequence": "TTGGTTGG", "colors": [0], "coverage": 10.0},
            {"id": "C", "sequence": "GGCCAAGG", "colors": [1], "coverage": 10.0},
        ],
        links=[
            {"id": "eAB", "source": "A", "target": "B", "colors": [0], "coverage": 10.0},
            {"id": "eBC", "source": "B", "target": "C", "colors": [0, 1], "coverage": 10.0},
        ],
    )
    for link in graph.cdbg.links:
        link.overlap = 0
    merged, edit = merge_adjacent(graph, "eAB")
    assert {unitig.unitig_id for unitig in graph.cdbg.unitigs} == {"A", "B", "C"}
    assert len(merged.cdbg.unitigs) == 2
    sequences = {unitig.sequence for unitig in merged.cdbg.unitigs}
    assert "ACGTACGTTTGGTTGG" in sequences
    restored = revert(merged, edit)
    assert {unitig.unitig_id for unitig in restored.cdbg.unitigs} == {"A", "B", "C"}
    scratch = graph.copy()
    inplace, _inplace_edit = merge_adjacent(scratch, "eAB", copy_graph=False)
    assert inplace is scratch
    assert any(unitig.sequence == "ACGTACGTTTGGTTGG" for unitig in inplace.cdbg.unitigs)
    compacted = compact_same_colour(graph)
    assert any(unitig.sequence == "ACGTACGTTTGGTTGG" for unitig in compacted.cdbg.unitigs)
    assert any(unitig.sequence == "GGCCAAGG" for unitig in compacted.cdbg.unitigs)
