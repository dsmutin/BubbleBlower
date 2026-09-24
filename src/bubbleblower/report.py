"""Write the contracted BubbleBlower tables and a resolved CDBG."""

from __future__ import annotations

from pathlib import Path

from metametro.formats.cdbg.io import dump_cdbg

from bubbleblower.search import SearchResult


def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_result(result: SearchResult, directory: str | Path) -> None:
    """Write ``resolved/``, ``bubble_results.tsv``, ``edit_history.tsv``, and ``state_scores.tsv``."""
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    resolved = root / "resolved"
    dump_cdbg(result.graph.cdbg, resolved)
    _write_tsv(
        resolved / "node_coverage.tsv",
        ["unitig_id", "coverage"],
        [
            [unitig_id, f"{coverage:.6g}"]
            for unitig_id, coverage in sorted(result.graph.node_coverage.items())
        ],
    )
    _write_tsv(
        resolved / "link_coverage.tsv",
        ["link_id", "coverage"],
        [
            [link_id, f"{coverage:.6g}"]
            for link_id, coverage in sorted(result.graph.link_coverage.items())
        ],
    )
    bubble_rows = []
    for index, posterior in enumerate(result.posteriors):
        edit = result.edits[0].edit_type if result.edits else ""
        bubble_rows.append(
            [
                posterior.bubble_id,
                "0",
                f"{posterior.p_error:.6g}",
                f"{posterior.p_strain:.6g}",
                posterior.decision,
                edit if index == 0 else "",
                f"{posterior.confidence:.6g}",
            ]
        )
    _write_tsv(
        root / "bubble_results.tsv",
        ["bubble_id", "iteration", "P_error", "P_strain", "decision", "selected_edit", "confidence"],
        bubble_rows,
    )
    history = []
    for step_index, step in enumerate(result.steps):
        edit = step.edit
        history.append(
            [
                str(step.iteration),
                str(step_index + 1),
                str(step_index),
                "" if edit is None else edit.edit_id,
                "" if edit is None else edit.edit_type,
                "" if edit is None else ",".join(edit.source_ids),
                "" if edit is None else ",".join(edit.target_ids),
                f"{step.score_before:.6g}",
                f"{step.score_after:.6g}",
                "yes" if step.accepted else "no",
            ]
        )
    _write_tsv(
        root / "edit_history.tsv",
        [
            "iteration",
            "state_id",
            "parent_state",
            "edit_id",
            "edit_type",
            "source_ids",
            "target_ids",
            "score_before",
            "score_after",
            "accepted",
        ],
        history,
    )
    score_rows = []
    for index, score in enumerate(result.scores):
        score_rows.append(
            [
                str(index),
                str(index),
                f"{score.total:.6g}",
                f"{score.coverage:.6g}",
                f"{score.flow:.6g}",
                f"{-score.topology:.6g}",
                f"{-score.complexity:.6g}",
                str(score.n_bubbles),
                str(score.n_instances),
            ]
        )
    _write_tsv(
        root / "state_scores.tsv",
        [
            "iteration",
            "state_id",
            "score",
            "coverage_score",
            "likelihood_score",
            "topology_penalty",
            "complexity_penalty",
            "n_bubbles",
            "n_instances",
        ],
        score_rows,
    )
