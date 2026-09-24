"""Optional integration against an existing totally coloured graph."""

from __future__ import annotations

from pathlib import Path

import pytest

pytestmark = pytest.mark.optional

_HELDOUT = Path(
    "/mnt/tank/scratch/partition-metagenomics/smuteam/vaegbin_improved/examples/heldout_genera/work/reprofile/tocumg/cdbg"
)


@pytest.mark.skipif(not _HELDOUT.is_dir(), reason="heldout_genera CDBG is not on this machine")
def test_heldout_cdbg_loads() -> None:
    """The existing ToCUMG loads and bubble detection finishes."""
    from bubbleblower.detect import detect_bubbles
    from bubbleblower.pipeline import load_assembly

    graph = load_assembly(_HELDOUT)
    bubbles = detect_bubbles(graph)
    assert graph.cdbg.unitigs
    assert isinstance(bubbles, list)
