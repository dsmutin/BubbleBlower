"""Classifier metrics. Ground truth is an evaluator input, not a graph input."""

from __future__ import annotations

from bubbleblower.classify import Posterior


def _safe_div(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return 0.0
    return numerator / denominator


def classification_metrics(rows: list[tuple[str, Posterior]]) -> dict[str, float]:
    """Precision, recall, and F1 for the error class, plus AUROC of ``P_error``.

    Each row is ``(truth_type, posterior)`` with truth ``error`` or ``strain``.
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
    f1 = _safe_div(2 * precision * recall, precision + recall)
    return {
        "precision_error": precision,
        "recall_error": recall,
        "f1_error": f1,
        "auroc_error": _auroc(rows),
        "tp": float(tp),
        "fp": float(fp),
        "fn": float(fn),
        "tn": float(tn),
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
