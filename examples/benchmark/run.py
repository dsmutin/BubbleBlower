#!/usr/bin/env python3
"""Classify the 50-bubble graph-only benchmark and require a usable error F1."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
_METAMETRO = ROOT.parent / "metametro" / "src"
if _METAMETRO.is_dir():
    sys.path.insert(0, str(_METAMETRO))

from bubbleblower.classify import classify_bubble  # noqa: E402
from bubbleblower.detect import detect_bubbles  # noqa: E402
from bubbleblower.evaluate import classification_metrics  # noqa: E402
from bubbleblower.features import extract_features  # noqa: E402
from bubbleblower.generate import generate_bubble_benchmark, write_ground_truth  # noqa: E402


def run() -> int:
    """Write ground truth and metrics. Fail when error F1 is below 0.75."""
    graph, truth = generate_bubble_benchmark(seed=42)
    out = Path(__file__).resolve().parent / "data"
    write_ground_truth(truth, out / "ground_truth" / "bubbles.tsv")
    by_source = {bubble.source: bubble for bubble in detect_bubbles(graph)}
    rows = []
    table = ["bubble_id\ttype\tdecision\tp_error\tp_strain"]
    for row in truth:
        posterior = classify_bubble(extract_features(graph, by_source[row["source_id"]]))
        rows.append((row["type"], posterior))
        table.append(
            f"{row['bubble_id']}\t{row['type']}\t{posterior.decision}\t{posterior.p_error:.4f}\t{posterior.p_strain:.4f}"
        )
    metrics = classification_metrics(rows)
    (out / "bubble_results.tsv").write_text("\n".join(table) + "\n", encoding="utf-8")
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metrics, indent=2))
    if metrics["f1_error"] < 0.75 or metrics["auroc_error"] < 0.9:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
