"""Baseline implementations that hold contracts until real logic lands."""

from __future__ import annotations


def run_pipeline(input_path: str | None = None) -> dict:
    """Classify the built-in strain bubble, or a CDBG directory.

    The return value keeps ``status``, ``ok``, and ``input_path``.
    """
    from bubbleblower.fixtures import strain_bubble
    from bubbleblower.pipeline import load_assembly, summarize

    if input_path:
        graph = load_assembly(input_path)
        summary = summarize(graph)
        summary["status"] = "resolved"
        summary["ok"] = True
        summary["input_path"] = input_path
        return summary
    summary = summarize(strain_bubble())
    decision = summary["bubbles"][0]["decision"] if summary["bubbles"] else None
    return {
        "status": "resolved",
        "ok": decision == "strain",
        "input_path": input_path,
        "decision": decision,
        "n_bubbles": summary["n_bubbles"],
    }
