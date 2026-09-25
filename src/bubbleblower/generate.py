"""Deterministic graph-only bubble benchmark.

Ground truth is written beside the graph and is not read by the resolver.
The graph itself is built by MetaMetro ``bubble_strain_3_n50``.
"""

from __future__ import annotations

from pathlib import Path

from bubbleblower.graph import AssemblyGraph, adopt_cfa


def generate_bubble_benchmark(
    *,
    n_strains: int = 3,
    n_bubbles: int = 50,
    n_error_bubbles: int = 25,
    abundance: tuple[float, ...] = (60.0, 40.0, 20.0),
    error_rate: float = 0.01,
    seed: int = 42,
) -> tuple[AssemblyGraph, list[dict[str, str]]]:
    """Build ``n_bubbles`` disjoint bubbles and a ground-truth table.

    Sequences and the retain/pop labels come from MetaMetro. This function
    adopts the compacted graph and rewrites truth ids to unitig ids.
    """
    from metametro.bench.data.universal.bubbles import synthetic_bubbles

    cfa, truth = synthetic_bubbles(
        n_strains=n_strains,
        n_bubbles=n_bubbles,
        n_error_bubbles=n_error_bubbles,
        abundance=abundance,
        error_rate=error_rate,
        seed=seed,
    )
    graph = adopt_cfa(cfa)
    graph.coverage_source = "simulated_poisson"
    unitig_of = {
        member: unitig.unitig_id
        for unitig in graph.cdbg.unitigs
        for member in unitig.members
    }
    for row in truth:
        row["source_id"] = unitig_of[row["source_id"]]
        row["sink_id"] = unitig_of[row["sink_id"]]
        row["branch_ids"] = ",".join(unitig_of[branch] for branch in row["branch_ids"].split(","))
    return graph, truth


def write_ground_truth(rows: list[dict[str, str]], path: Path) -> None:
    """Write ``bubbles.tsv``. This file is for the evaluator only."""
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ["bubble_id", "type", "strain_ids", "source_id", "sink_id", "branch_ids", "expected_action"]
    lines = ["\t".join(header)]
    for row in rows:
        lines.append("\t".join(row[column] for column in header))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
