"""Search variants and contracted output tables."""

from __future__ import annotations

from pathlib import Path

import pytest

from bubbleblower.fixtures import error_bubble
from bubbleblower.pipeline import load_assembly
from bubbleblower.report import write_result
from bubbleblower.search import beam_search, greedy_search, mcmc_search

pytestmark = pytest.mark.mandatory


def test_beam_and_mcmc_do_not_lose_score() -> None:
    graph = error_bubble()
    greedy = greedy_search(graph, max_iterations=5)
    beam = beam_search(graph, beam_width=4, max_iterations=5)
    chain = mcmc_search(graph, n_steps=15, seed=1)
    assert greedy.scores[-1].total > greedy.scores[0].total
    assert beam.scores[-1].total >= greedy.scores[0].total
    assert chain.scores[-1].total >= chain.scores[0].total - 1e-6
    assert greedy.edits
    assert greedy.edits[0].edit_type == "pop"


def test_result_tables_round_trip_coverage(tmp_path: Path) -> None:
    result = greedy_search(error_bubble(), max_iterations=4)
    write_result(result, tmp_path)
    assert (tmp_path / "bubble_results.tsv").is_file()
    assert (tmp_path / "edit_history.tsv").is_file()
    assert (tmp_path / "state_scores.tsv").is_file()
    loaded = load_assembly(tmp_path / "resolved")
    assert loaded.node_coverage.keys() == result.graph.node_coverage.keys()
    for unitig_id, coverage in result.graph.node_coverage.items():
        assert loaded.node_coverage[unitig_id] == pytest.approx(coverage)
