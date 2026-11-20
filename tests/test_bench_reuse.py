"""Mandatory: ``ensure_bench`` rejects a stale reused build (#1209)."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.mandatory

pytest.importorskip("metametro.bench.build")
canonical = pytest.importorskip("metametro.bench").__dict__.get("canonical_identity")

from metametro.errors import ContractError  # noqa: E402
from metametro.tables import read_yaml, write_yaml  # noqa: E402

from bubbleblower.bench_input import BENCH_ROOT_ENV, ensure_bench, load_bench_graph  # noqa: E402

BENCH = "mock_bubble"
STALE = "0" * 64


@pytest.fixture
def built(tmp_path, monkeypatch):
    monkeypatch.delenv(BENCH_ROOT_ENV, raising=False)
    return ensure_bench(BENCH, root=tmp_path)


def _needs_canonical() -> None:
    if canonical is None or canonical(BENCH) is None:
        pytest.skip("this MetaMetro has no canonical identity for mock_bubble")


def test_fresh_default_build_is_reused(tmp_path, built) -> None:
    assert ensure_bench(BENCH, root=tmp_path) == built
    load_bench_graph(BENCH, root=tmp_path)


def test_stale_default_build_raises(tmp_path, built) -> None:
    _needs_canonical()
    (built / "identity.sha256").write_text(f"{STALE}  identity\n", encoding="utf-8")
    with pytest.raises(ContractError, match="rerun metametro benchbuild mock_bubble"):
        ensure_bench(BENCH, root=tmp_path)
    with pytest.raises(ContractError, match="canonical"):
        load_bench_graph(BENCH, root=tmp_path)


def test_manifest_without_gold_namespaces_raises(tmp_path, built) -> None:
    manifest = read_yaml(built / "manifest.yaml")
    manifest.pop("gold_namespaces")
    write_yaml(built / "manifest.yaml", manifest)
    with pytest.raises(ContractError, match="gold_namespaces"):
        ensure_bench(BENCH, root=tmp_path)


def test_build_with_explicit_colourings_is_accepted(tmp_path, monkeypatch) -> None:
    _needs_canonical()
    from metametro.bench.build import build
    from metametro.bench.paths import default_outdir
    from metametro.bench.registry import resolve

    monkeypatch.delenv(BENCH_ROOT_ENV, raising=False)
    destination = default_outdir(resolve(BENCH), root=tmp_path)
    result = build(BENCH, outdir=destination, colourings=("composition_kmeans",))
    assert result.identity != canonical(BENCH)
    assert ensure_bench(BENCH, root=tmp_path) == destination


def test_bench_without_canonical_identity_is_accepted(tmp_path, built, monkeypatch) -> None:
    import metametro.bench

    monkeypatch.setattr(metametro.bench, "canonical_identity", lambda _name: None, raising=False)
    (built / "identity.sha256").write_text(f"{STALE}  identity\n", encoding="utf-8")
    assert ensure_bench(BENCH, root=tmp_path) == built


@pytest.mark.parametrize(
    ("applied", "default"),
    [
        (["as_built", "composition_kmeans", "decaying"], True),
        (["as_built", "composition_kmeans", "kraken2", "decaying"], True),
        (["composition_kmeans"], False),
        (["as_built", "decaying", "composition_kmeans"], False),
        (["as_built", "read_accession"], False),
        ([], False),
    ],
)
def test_default_colourings_rule(applied: list[str], default: bool) -> None:
    from bubbleblower.bench_input import _default_colourings

    assert _default_colourings(applied) is default
