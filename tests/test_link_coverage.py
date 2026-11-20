"""Mandatory: links without coverage get the ``split`` share of the source (#1202)."""

from __future__ import annotations

import pytest

from bubbleblower.graph import from_cdbg, records_to_tocumg, split_link_coverage
from bubbleblower.score import score_graph

pytestmark = pytest.mark.mandatory

COLORS = [{"color_id": "0", "namespace": "sample", "value": "s1"}]
NODES = [
    {"id": "S", "sequence": "ACGTACGT", "colors": [0], "coverage": 120.0},
    {"id": "A", "sequence": "ATATCG", "colors": [0], "coverage": 30.0},
    {"id": "B", "sequence": "CGCGTA", "colors": [0], "coverage": 90.0},
    {"id": "T", "sequence": "GGCCTTAA", "colors": [0], "coverage": 120.0},
]
FORWARD = [("eSA", "S", "A"), ("eAT", "A", "T"), ("eSB", "S", "B"), ("eBT", "B", "T")]


def _link(link_id: str, source: str, target: str, orientation: str) -> dict:
    return {"id": link_id, "source": source, "target": target, "orientation": orientation, "colors": [0], "coverage": 1.0}


def _without_link_coverage(links: list[dict]):
    built = records_to_tocumg(graph_id="split", colors=COLORS, nodes=NODES, links=links)
    return from_cdbg(built.cdbg, node_coverage=dict(built.node_coverage))


def _uid(graph, name: str) -> str:
    return graph.unitig_id_for_member(name)


def test_graph_without_link_coverage_scores_after_split() -> None:
    forward = [_link(i, s, t, "++") for i, s, t in FORWARD]
    reverse = [_link(f"r{i}", t, s, "--") for i, s, t in FORWARD]
    graph = _without_link_coverage(forward + reverse)
    assert graph.link_coverage == {}
    with pytest.raises(KeyError):
        score_graph(graph)
    graph.link_coverage = split_link_coverage(graph.cdbg, graph.node_coverage, graph.link_coverage)
    graph.validate()
    score_graph(graph)


def test_fork_sums_to_source_coverage_on_each_end() -> None:
    forward = [_link(i, s, t, "++") for i, s, t in FORWARD]
    reverse = [_link(f"r{i}", t, s, "--") for i, s, t in FORWARD]
    graph = _without_link_coverage(forward + reverse)
    coverage = split_link_coverage(graph.cdbg, graph.node_coverage, {})
    source, sink = _uid(graph, "S"), _uid(graph, "T")
    for unitig, strand in ((source, "+"), (sink, "-")):
        fork = [link for link in graph.cdbg.links if link.source == unitig and link.orientation[0] == strand]
        assert len(fork) == 2
        assert sum(coverage[link.link_id] for link in fork) == pytest.approx(graph.node_coverage[unitig])
    by_target = {link.target: coverage[link.link_id] for link in graph.cdbg.links if link.source == source}
    assert by_target[_uid(graph, "A")] == pytest.approx(120.0 * 30.0 / 120.0)
    assert by_target[_uid(graph, "B")] == pytest.approx(120.0 * 90.0 / 120.0)


def test_existing_link_coverage_is_kept() -> None:
    graph = _without_link_coverage([_link(i, s, t, "++") for i, s, t in FORWARD])
    first = graph.cdbg.links[0].link_id
    coverage = split_link_coverage(graph.cdbg, graph.node_coverage, {first: 7.0})
    assert coverage[first] == 7.0
    assert set(coverage) == {link.link_id for link in graph.cdbg.links}


def test_bench_without_link_coverage_loads_and_scores(tmp_path, monkeypatch) -> None:
    pytest.importorskip("metametro.bench.build")
    from bubbleblower.bench_input import BENCH_ROOT_ENV, load_bench_graph

    monkeypatch.delenv(BENCH_ROOT_ENV, raising=False)
    _root, graph = load_bench_graph("bubble_reads_2", root=tmp_path)
    assert set(graph.link_coverage) == {link.link_id for link in graph.cdbg.links}
    assert "split" in graph.coverage_source
    graph.validate()
    score_graph(graph)
