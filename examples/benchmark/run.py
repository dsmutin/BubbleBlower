#!/usr/bin/env python3
"""Run the 50-bubble benchmark. Primary quality metric is AMBER F1."""

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
from bubbleblower.classify import classify_bubble  # noqa: E402
from bubbleblower.detect import detect_bubbles  # noqa: E402
from bubbleblower.evaluate import classification_metrics, resolution_metrics  # noqa: E402
from bubbleblower.features import extract_features  # noqa: E402
from bubbleblower.report import write_result  # noqa: E402
from bubbleblower.search import greedy_search  # noqa: E402

# Primary metric. Classifier and resolver must both clear this floor.
AMBER_F1_MIN = 0.85


def run() -> int:
    """Score the MetaMetro 50-bubble build. Primary quality metric is AMBER F1."""
    root, graph = load_bench_graph("bubble_strain_3_n50")
    from metametro.tables import read_tsv

    _header, truth = read_tsv(root / "ground_truth" / "bubbles.tsv")
    unitig_of = {
        member: unitig.unitig_id
        for unitig in graph.cdbg.unitigs
        for member in unitig.members
    }
    for row in truth:
        row["source_id"] = unitig_of[row["source_id"]]
        row["sink_id"] = unitig_of[row["sink_id"]]
        row["branch_ids"] = ",".join(unitig_of[branch] for branch in row["branch_ids"].split(","))
    out = Path(__file__).resolve().parent / "data"
    by_source = {bubble.source: bubble for bubble in detect_bubbles(graph)}
    rows = []
    table = ["bubble_id\ttype\tdecision\tp_error\tp_strain"]
    for row in truth:
        posterior = classify_bubble(extract_features(graph, by_source[row["source_id"]]))
        rows.append((row["type"], posterior))
        table.append(
            f"{row['bubble_id']}\t{row['type']}\t{posterior.decision}\t{posterior.p_error:.4f}\t{posterior.p_strain:.4f}"
        )
    classified = classification_metrics(rows)
    resolved = greedy_search(graph, max_iterations=60, patience=1, max_runtime_s=60.0)
    write_result(resolved, out / "resolved_run")
    resolution = resolution_metrics(
        resolved.graph,
        truth,
        score_before=resolved.scores[0],
        score_after=resolved.scores[-1],
    )
    metrics = {
        "amber_f1": resolution["amber_f1"],
        "classification": classified,
        "resolution": resolution,
    }
    (out / "bubble_results.tsv").write_text("\n".join(table) + "\n", encoding="utf-8")
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metrics, indent=2))
    if classified["amber_f1"] < AMBER_F1_MIN or classified["auroc_error"] < 0.9:
        return 1
    if resolution["amber_f1"] < AMBER_F1_MIN or resolution["false_pops"] != 0:
        return 1
    if resolution["delta_score"] < 0:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
