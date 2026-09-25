"""Load a MEGAHIT FASTG as a MetaMetro ToCUMG.

``fastg_to_cfa`` keeps forward sequences, stores strand in the edge
orientation, and compacts chains that overlap by the assembly ``k``.
Coverage is the ``cov_`` field already in the header.
"""

from __future__ import annotations

from pathlib import Path

from metametro.contracts.assembly import fastg_to_cfa

from bubbleblower.graph import AssemblyGraph, adopt_cfa


def load_fastg(path: str | Path, *, k: int, graph_id: str = "fastg") -> AssemblyGraph:
    """Read a ``contig2fastg`` file. ``k`` is the assembly k MEGAHIT used."""
    if k < 1:
        raise ValueError("k must be positive")
    cfa = fastg_to_cfa(Path(path), k=k, graph_id=graph_id)
    graph = adopt_cfa(cfa)
    graph.coverage_source = "megahit_fastg_cov"
    return graph
