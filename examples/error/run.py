"""Error-bubble example. The low-coverage branch is removed and the score rises."""

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
from bubbleblower.report import write_result  # noqa: E402
from bubbleblower.search import greedy_search  # noqa: E402


def run() -> int:
    """Resolve the sequencing-error bubble from MetaMetro and require the error node to be gone."""
    _root, graph = load_bench_graph("bubble_error_1", namespaces=("taxon",))
    result = greedy_search(graph, max_iterations=5)
    out = Path(__file__).resolve().parent / "data"
    write_result(result, out)
    members = result.graph.member_ids()
    payload = {
        "score_before": result.scores[0].total,
        "score_after": result.scores[-1].total,
        "removed_error_branch": "E" not in members,
        "kept_true_branch": "A" in members,
    }
    (out / "summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if not payload["removed_error_branch"] or not payload["kept_true_branch"]:
        return 1
    if payload["score_after"] <= payload["score_before"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
