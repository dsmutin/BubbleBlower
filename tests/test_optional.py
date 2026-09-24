"""Optional integration against an existing totally coloured graph."""

from __future__ import annotations

from pathlib import Path

import pytest

pytestmark = pytest.mark.optional

_TOCUMG = Path(
    "/mnt/tank/scratch/partition-metagenomics/smuteam/vaegbin_improved/examples/heldout_genera/work/reprofile/tocumg"
)
_HELDOUT = _TOCUMG / "cdbg"


@pytest.mark.skipif(not _HELDOUT.is_dir(), reason="heldout_genera CDBG is not on this machine")
def test_heldout_cdbg_loads() -> None:
    """The existing ToCUMG loads. This export is not a bubble graph.

    Both ``cfa/edges.tsv`` and ``cdbg/links.tsv`` contain 24 adjacencies on
    thousands of nodes (``graph_type: knn``). Detection finishing with zero
    bubbles is the observed topology, not a missing importer.
    """
    from bubbleblower.detect import detect_bubbles
    from bubbleblower.pipeline import load_assembly

    graph = load_assembly(_HELDOUT)
    bubbles = detect_bubbles(graph)
    assert graph.cdbg.unitigs
    assert graph.cdbg.metadata.get("graph_type") == "knn"
    assert len(graph.cdbg.links) < 100
    assert isinstance(bubbles, list)
    cfa_edges = _TOCUMG / "cfa" / "edges.tsv"
    if cfa_edges.is_file():
        n_cfa = max(0, len(cfa_edges.read_text(encoding="utf-8").splitlines()) - 1)
        assert n_cfa == len(graph.cdbg.links)
