"""Mandatory bubble contracts: detection, edits, classifier, search."""

from __future__ import annotations

import pytest

from bubbleblower.classify import classify_bubble
from bubbleblower.detect import detect_bubbles
from bubbleblower.edits import duplicate_instance, merge_instances, pop_branch, revert, split_instance
from bubbleblower.evaluate import classification_metrics
from bubbleblower.features import extract_features
from bubbleblower.fixtures import error_bubble, nested_bubbles, shared_duplicate_node, strain_bubble
from bubbleblower.generate import generate_bubble_benchmark
from bubbleblower.graph import semantic_signature, sequence_hash
from bubbleblower.search import greedy_search

pytestmark = pytest.mark.mandatory


def test_strain_bubble_is_retained() -> None:
    graph = strain_bubble()
    bubbles = detect_bubbles(graph)
    assert len(bubbles) == 1
    features = extract_features(graph, bubbles[0])
    posterior = classify_bubble(features)
    assert posterior.decision == "strain"
    assert features.colours_disjoint is True
    assert {node.unitig_id for node in graph.cdbg.unitigs} >= {"S", "T"}


def test_error_bubble_is_classified_error() -> None:
    graph = error_bubble()
    bubbles = detect_bubbles(graph)
    assert len(bubbles) == 1
    posterior = classify_bubble(extract_features(graph, bubbles[0]))
    assert posterior.decision == "error"
    assert posterior.p_error > posterior.p_strain


def test_duplicate_keeps_sequence_and_changes_instance() -> None:
    graph = shared_duplicate_node()
    edited, edit = duplicate_instance(graph, "X")
    original = edited.unitig("X")
    created = edited.unitig(edit.target_ids[1])
    assert original.sequence == created.sequence
    assert sequence_hash(original.sequence) == sequence_hash(created.sequence)
    assert original.unitig_id != created.unitig_id
    assert semantic_signature(graph) != semantic_signature(edited)


def test_duplicate_merge_round_trip() -> None:
    graph = shared_duplicate_node()
    edited, edit = duplicate_instance(graph, "X")
    restored, _merge = merge_instances(edited, list(edit.target_ids), "X")
    assert semantic_signature(restored) == semantic_signature(graph)


def test_pop_restore_round_trip() -> None:
    graph = error_bubble()
    bubble = detect_bubbles(graph)[0]
    edited, edit = pop_branch(graph, bubble, 1)
    assert semantic_signature(edited) != semantic_signature(graph)
    restored = revert(edited, edit)
    assert semantic_signature(restored) == semantic_signature(graph)


def test_split_restore_round_trip() -> None:
    graph = strain_bubble()
    bubble = detect_bubbles(graph)[0]
    groups = [[branch.links[0]] for branch in bubble.branches]
    edited, edit = split_instance(graph, bubble.source, groups)
    assert len(edited.cdbg.unitigs) == len(graph.cdbg.unitigs) + 1
    restored = revert(edited, edit)
    assert semantic_signature(restored) == semantic_signature(graph)


def test_fifty_bubble_classifier() -> None:
    graph, truth = generate_bubble_benchmark(seed=42)
    bubbles = detect_bubbles(graph)
    assert len(bubbles) == 50
    by_source = {bubble.source: bubble for bubble in bubbles}
    rows = []
    for row in truth:
        bubble = by_source[row["source_id"]]
        rows.append((row["type"], classify_bubble(extract_features(graph, bubble))))
    metrics = classification_metrics(rows)
    assert metrics["auroc_error"] >= 0.9
    assert metrics["f1_error"] >= 0.75


def test_rare_strain_needs_more_than_coverage() -> None:
    from bubbleblower.graph import build_graph

    graph = build_graph(
        graph_id="rare",
        colors=[
            {"color_id": "0", "namespace": "taxon", "value": "dominant"},
            {"color_id": "1", "namespace": "taxon", "value": "rare"},
        ],
        nodes=[
            {"id": "S", "sequence": "ACGTACGT", "colors": [0, 1], "coverage": 105.0},
            {"id": "A", "sequence": "ATATATAT", "colors": [0], "coverage": 100.0},
            {"id": "R", "sequence": "CGCGCGCG", "colors": [1], "coverage": 5.0},
            {"id": "T", "sequence": "GGCCTTAA", "colors": [0, 1], "coverage": 105.0},
        ],
        links=[
            {"id": "e1", "source": "S", "target": "A", "colors": [0], "coverage": 100.0},
            {"id": "e2", "source": "A", "target": "T", "colors": [0], "coverage": 100.0},
            {"id": "e3", "source": "S", "target": "R", "colors": [1], "coverage": 5.0},
            {"id": "e4", "source": "R", "target": "T", "colors": [1], "coverage": 5.0},
        ],
    )
    bubble = detect_bubbles(graph)[0]
    coverage = classify_bubble(extract_features(graph, bubble), mode="coverage")
    multimodal = classify_bubble(
        extract_features(graph, bubble, linkage=True),
        mode="multimodal",
    )
    assert coverage.decision != "strain" or coverage.p_strain < multimodal.p_strain
    assert multimodal.decision == "strain"


def test_equal_coverage_is_uncertain_without_linkage() -> None:
    from bubbleblower.graph import build_graph

    graph = build_graph(
        graph_id="equal",
        colors=[
            {"color_id": "0", "namespace": "taxon", "value": "taxon_1"},
            {"color_id": "1", "namespace": "taxon", "value": "taxon_2"},
        ],
        nodes=[
            {"id": "S", "sequence": "ACGTACGT", "colors": [0, 1], "coverage": 102.0},
            {"id": "A", "sequence": "ATATATAT", "colors": [0], "coverage": 50.0},
            {"id": "B", "sequence": "CGCGCGCG", "colors": [1], "coverage": 52.0},
            {"id": "T", "sequence": "GGCCTTAA", "colors": [0, 1], "coverage": 102.0},
        ],
        links=[
            {"id": "e1", "source": "S", "target": "A", "colors": [0], "coverage": 50.0},
            {"id": "e2", "source": "A", "target": "T", "colors": [0], "coverage": 50.0},
            {"id": "e3", "source": "S", "target": "B", "colors": [1], "coverage": 52.0},
            {"id": "e4", "source": "B", "target": "T", "colors": [1], "coverage": 52.0},
        ],
    )
    bubble = detect_bubbles(graph)[0]
    coverage = classify_bubble(extract_features(graph, bubble))
    linked = classify_bubble(extract_features(graph, bubble, linkage=True), mode="multimodal")
    assert coverage.decision == "other"
    assert linked.decision == "strain"


def test_misleading_coverage_uses_kmer_evidence() -> None:
    from bubbleblower.graph import build_graph

    graph = build_graph(
        graph_id="misleading",
        colors=[
            {"color_id": "0", "namespace": "taxon", "value": "taxon_1"},
            {"color_id": "1", "namespace": "taxon", "value": "taxon_2"},
        ],
        nodes=[
            {"id": "S", "sequence": "ACGTACGT", "colors": [0], "coverage": 18.0},
            {"id": "A", "sequence": "ATATATAT", "colors": [0], "coverage": 10.0},
            {"id": "E", "sequence": "CGCGCGCG", "colors": [1], "coverage": 8.0},
            {"id": "T", "sequence": "GGCCTTAA", "colors": [0], "coverage": 18.0},
        ],
        links=[
            {"id": "e1", "source": "S", "target": "A", "colors": [0], "coverage": 10.0},
            {"id": "e2", "source": "A", "target": "T", "colors": [0], "coverage": 10.0},
            {"id": "e3", "source": "S", "target": "E", "colors": [1], "coverage": 8.0},
            {"id": "e4", "source": "E", "target": "T", "colors": [1], "coverage": 8.0},
        ],
    )
    bubble = detect_bubbles(graph)[0]
    coverage = classify_bubble(extract_features(graph, bubble))
    full = classify_bubble(extract_features(graph, bubble, low_kmer=True), mode="multimodal")
    assert coverage.decision != "error"
    assert full.decision == "error"


def test_greedy_improves_score_and_pops_error() -> None:
    error = error_bubble()
    strain = strain_bubble()
    before_error = greedy_search(error, max_iterations=5)
    assert before_error.scores[-1].total > before_error.scores[0].total
    assert all(unitig.unitig_id != "E" for unitig in before_error.graph.cdbg.unitigs)
    resolved = greedy_search(strain, max_iterations=5)
    sequences = {unitig.sequence for unitig in resolved.graph.cdbg.unitigs}
    assert "ATATCG" in sequences
    assert "CGCGTA" in sequences
    assert resolved.scores[-1].total >= resolved.scores[0].total


def test_nested_bubbles_are_both_found() -> None:
    bubbles = detect_bubbles(nested_bubbles())
    assert len(bubbles) >= 2
