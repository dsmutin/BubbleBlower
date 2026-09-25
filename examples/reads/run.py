#!/usr/bin/env python3
"""Read-level example. The graph is the MetaMetro ToCUMG of the simulated reads.

Reads, the de Bruijn CFA, and its colours come from MetaMetro. BubbleBlower
adopts that compacted graph and does not build a second allele graph.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
_METAMETRO = ROOT.parent / "metametro" / "src"
if _METAMETRO.is_dir():
    sys.path.insert(0, str(_METAMETRO))

from bubbleblower.bench_input import load_bench_graph  # noqa: E402
from bubbleblower.detect import detect_bubbles  # noqa: E402


def run() -> int:
    """Load the MetaMetro two-strain read graph and count bubbles on that ToCUMG."""
    from metametro.formats.cfa.io import load_cfa

    root, graph = load_bench_graph("bubble_reads_2")
    cfa = load_cfa(root / "cfa")
    truth_lines = (root / "ground_truth" / "read_to_genome.tsv").read_text(encoding="utf-8").splitlines()
    n_reads = sum(1 for line in truth_lines if line.strip()) - 1
    bubbles = detect_bubbles(graph)
    payload = {
        "n_reads": n_reads,
        "n_cfa_nodes": len(cfa.nodes),
        "n_unitigs": len(graph.cdbg.unitigs),
        "n_bubbles": len(bubbles),
        "n_colours": len(cfa.colors or []),
        "coverage_source": graph.coverage_source,
        "graph_source": "metametro.benchbuild bubble_reads_2",
    }
    out = Path(__file__).resolve().parent / "data"
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if payload["n_reads"] != 40 or payload["n_colours"] != 2 or payload["n_unitigs"] < 1:
        return 1
    if graph.cdbg.metadata.get("contract") != "cfa_to_cdbg":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
