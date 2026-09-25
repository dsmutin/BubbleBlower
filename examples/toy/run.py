#!/usr/bin/env python3
"""Run bubbleblower baseline on the toy example and check the contract."""

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
from bubbleblower.pipeline import summarize  # noqa: E402


def run() -> int:
    """Classify the MetaMetro two-taxon bubble and require a strain decision."""
    _root, graph = load_bench_graph("bubble_strain_2", namespaces=("taxon",))
    summary = summarize(graph)
    decision = summary["bubbles"][0]["decision"] if summary["bubbles"] else None
    payload = {
        "status": "resolved",
        "ok": decision == "strain",
        "decision": decision,
        "n_bubbles": summary["n_bubbles"],
        "bench": "bubble_strain_2",
        "input_path": str(_root),
    }
    out = Path(__file__).resolve().parent / "data" / "toy_result.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if payload["ok"] is not True or payload["n_bubbles"] != 1:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
