"""``ensure_bench`` and ``load_bench_graph`` honour a bench root (#1165) and reject a missing namespace (#1169)."""

from __future__ import annotations

import pytest

pytest.importorskip("metametro.bench.build")

from metametro.errors import ContractError  # noqa: E402

from bubbleblower.bench_input import BENCH_ROOT_ENV, ensure_bench, load_bench_graph  # noqa: E402

BENCH = "mock_bubble"


def test_root_argument_builds_under_root(tmp_path, monkeypatch):
    monkeypatch.delenv(BENCH_ROOT_ENV, raising=False)
    root, graph = load_bench_graph(BENCH, root=tmp_path)
    assert root.is_relative_to(tmp_path / "data" / "bench" / BENCH)
    assert (root / "cfa" / "metadata.yaml").is_file()
    assert graph.cdbg.unitigs
    assert ensure_bench(BENCH, root=tmp_path) == root


def test_env_root_is_used_without_argument(tmp_path, monkeypatch):
    monkeypatch.setenv(BENCH_ROOT_ENV, str(tmp_path))
    root = ensure_bench(BENCH)
    assert root.is_relative_to(tmp_path / "data" / "bench" / BENCH)


def test_argument_wins_over_env(tmp_path, monkeypatch):
    monkeypatch.setenv(BENCH_ROOT_ENV, str(tmp_path / "env"))
    root = ensure_bench(BENCH, root=tmp_path / "arg")
    assert root.is_relative_to(tmp_path / "arg")
    assert not (tmp_path / "env").exists()


def test_missing_namespace_raises(tmp_path, monkeypatch):
    monkeypatch.delenv(BENCH_ROOT_ENV, raising=False)
    with pytest.raises(ContractError, match="load_bench_gold"):
        load_bench_graph(BENCH, namespaces=("taxon",), root=tmp_path)
