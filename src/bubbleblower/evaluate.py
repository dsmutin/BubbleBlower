"""Classifier and resolution metrics. Ground truth is an evaluator input."""

from __future__ import annotations

from bubbleblower.classify import Posterior
from bubbleblower.graph import AssemblyGraph
from bubbleblower.score import Score


def _safe_div(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return 0.0
    return numerator / denominator


def amber_f1(precision: float, recall: float) -> float:
    """AMBER-style F1: harmonic mean of purity (precision) and completeness (recall)."""
    return _safe_div(2 * precision * recall, precision + recall)


def classification_metrics(rows: list[tuple[str, Posterior]]) -> dict[str, float]:
    """Precision, recall, AMBER F1, AUROC, and AUPRC for the error class.

    Each row is ``(truth_type, posterior)`` with truth ``error`` or ``strain``.
    ``amber_f1`` is the primary quality metric (same value as ``f1_error``).
    """
    tp = fp = fn = tn = 0
    for truth, posterior in rows:
        predicted_error = posterior.decision == "error"
        actual_error = truth == "error"
        if predicted_error and actual_error:
            tp += 1
        elif predicted_error and not actual_error:
            fp += 1
        elif not predicted_error and actual_error:
            fn += 1
        else:
            tn += 1
    precision = _safe_div(tp, tp + fp)
    recall = _safe_div(tp, tp + fn)
    f1 = amber_f1(precision, recall)
    strain_precision = _safe_div(tn, tn + fn)
    strain_recall = _safe_div(tn, tn + fp)
    return {
        "amber_f1": f1,
        "precision_error": precision,
        "recall_error": recall,
        "f1_error": f1,
        "precision_strain": strain_precision,
        "recall_strain": strain_recall,
        "f1_strain": amber_f1(strain_precision, strain_recall),
        "auroc_error": _auroc(rows),
        "auprc_error": _auprc(rows),
        "tp": float(tp),
        "fp": float(fp),
        "fn": float(fn),
        "tn": float(tn),
    }


def resolution_metrics(
    graph: AssemblyGraph,
    truth: list[dict[str, str]],
    *,
    score_before: Score | None = None,
    score_after: Score | None = None,
) -> dict[str, float]:
    """Graph-resolution counts and AMBER F1.

    A true pop is an error bubble whose last listed branch is gone. A false
    pop is a strain bubble that lost any branch. ``amber_f1`` is the primary
    quality metric: purity of pops times completeness of error removal.
    """
    ids = {unitig.unitig_id for unitig in graph.cdbg.unitigs}
    true_pops = 0
    false_pops = 0
    missed_errors = 0
    retained_strains = 0
    for row in truth:
        branches = row["branch_ids"].split(",")
        gone = [branch for branch in branches if branch not in ids]
        if row["type"] == "error":
            error_branch = branches[-1]
            if error_branch not in ids:
                true_pops += 1
            else:
                missed_errors += 1
        else:
            if gone:
                false_pops += 1
            else:
                retained_strains += 1
    precision = _safe_div(true_pops, true_pops + false_pops)
    recall = _safe_div(true_pops, true_pops + missed_errors)
    f1 = amber_f1(precision, recall)
    delta = 0.0
    if score_before is not None and score_after is not None:
        delta = score_after.total - score_before.total
    return {
        "amber_f1": f1,
        "precision": precision,
        "recall": recall,
        "true_pops": float(true_pops),
        "false_pops": float(false_pops),
        "missed_errors": float(missed_errors),
        "retained_strains": float(retained_strains),
        "delta_score": delta,
    }


def _auroc(rows: list[tuple[str, Posterior]]) -> float:
    positives = [posterior.p_error for truth, posterior in rows if truth == "error"]
    negatives = [posterior.p_error for truth, posterior in rows if truth != "error"]
    if not positives or not negatives:
        return 0.0
    wins = 0.0
    for positive in positives:
        for negative in negatives:
            if positive > negative:
                wins += 1.0
            elif positive == negative:
                wins += 0.5
    return wins / (len(positives) * len(negatives))


def _auprc(rows: list[tuple[str, Posterior]]) -> float:
    """Average precision of ``P_error`` ranked high-to-low."""
    ranked = sorted(rows, key=lambda item: item[1].p_error, reverse=True)
    n_pos = sum(1 for truth, _posterior in ranked if truth == "error")
    if n_pos == 0:
        return 0.0
    seen_pos = 0
    area = 0.0
    for index, (truth, _posterior) in enumerate(ranked, start=1):
        if truth != "error":
            continue
        seen_pos += 1
        area += seen_pos / index
    return area / n_pos
