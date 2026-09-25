"""Optional integration against an existing totally coloured graph."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

pytestmark = pytest.mark.optional


def _tocumg() -> Path | None:
    """Return a held-out ToCUMG directory, or None when this machine has not built one."""
    env = os.environ.get("BUBBLEBLOWER_TOCUMG", "")
    if env:
        return Path(env)
    try:
        from metametro.bench.paths import default_outdir
        from metametro.bench.registry import resolve
    except ImportError:
        return None
    path = default_outdir(resolve("bacteria_species_20_heldout")) / "work" / "reprofile" / "tocumg"
    return path if path.is_dir() else None


def test_heldout_cdbg_loads() -> None:
    """The existing ToCUMG loads. This export is not a bubble graph.

    Both ``cfa/edges.tsv`` and ``cdbg/links.tsv`` contain 24 adjacencies on
    thousands of nodes (``graph_type: knn``). Detection finishing with zero
    bubbles is the observed topology, not a missing importer.
    """
    root = _tocumg()
    if root is None or not (root / "cdbg").is_dir():
        pytest.skip("held-out ToCUMG is not built; set BUBBLEBLOWER_TOCUMG or run benchbuild")
    from bubbleblower.detect import detect_bubbles
    from bubbleblower.pipeline import load_assembly

    graph = load_assembly(root / "cdbg")
    bubbles = detect_bubbles(graph)
    assert graph.cdbg.unitigs
    assert graph.cdbg.metadata.get("graph_type") == "knn"
    assert len(graph.cdbg.links) < 100
    assert isinstance(bubbles, list)
    cfa_edges = root / "cfa" / "edges.tsv"
    if cfa_edges.is_file():
        n_cfa = max(0, len(cfa_edges.read_text(encoding="utf-8").splitlines()) - 1)
        assert n_cfa == len(graph.cdbg.links)
