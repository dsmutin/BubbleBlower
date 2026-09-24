"""Coverage Bayesian bubble classifier.

The classifier returns probabilities. It does not edit the graph.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from bubbleblower.features import BubbleFeatures

MODEL_VERSION = "coverage-poisson-1"
_ERROR_EPS = (0.005, 0.01, 0.02, 0.03, 0.05)
_STRAIN_RHO = (0.15, 0.25, 0.33, 0.5, 0.67, 0.85, 1.0)


@dataclass(frozen=True)
class Posterior:
    """``P(error)``, ``P(strain)``, ``P(other)`` for one bubble."""

    bubble_id: str
    p_error: float
    p_strain: float
    p_other: float
    decision: str
    confidence: float
    model_version: str = MODEL_VERSION


def _log_poisson(count: float, lam: float) -> float:
    k = int(round(max(count, 0.0)))
    rate = max(lam, 1e-9)
    return k * math.log(rate) - rate - math.lgamma(k + 1)


def _logsumexp(values: list[float]) -> float:
    top = max(values)
    return top + math.log(sum(math.exp(value - top) for value in values))


def _error_ll(minor: float, major: float) -> float:
    terms = [_log_poisson(minor, eps * major) for eps in _ERROR_EPS]
    return _logsumexp(terms) - math.log(len(terms))


def _strain_ll(minor: float, major: float) -> float:
    terms = [
        _log_poisson(minor, rho * major) + _log_poisson(major, major)
        for rho in _STRAIN_RHO
    ]
    return _logsumexp(terms) - math.log(len(terms))


def classify_bubble(features: BubbleFeatures, *, mode: str = "coverage") -> Posterior:
    """Posterior over error, strain, and other.

    ``mode='coverage'`` uses branch coverage only.
    ``mode='multimodal'`` also uses colour equality, long-read linkage, and
    k-mer support when those fields are present.
    """
    if mode not in {"coverage", "multimodal"}:
        raise ValueError(f"unknown classifier mode: {mode}")
    coverages = [branch.mean_coverage for branch in features.branches]
    major = max(coverages) if coverages else 0.0
    minor = min(coverages) if coverages else 0.0
    ll_error = _error_ll(minor, major)
    ll_strain = _strain_ll(minor, major)
    ratio = features.coverage_ratio
    ll_other = min(ll_error, ll_strain) - 3.0
    if ratio >= 0.85:
        ll_other = max(ll_error, ll_strain) + 0.25
    if mode == "multimodal":
        if features.colours_disjoint and features.colour_consistent:
            ll_strain += 6.0
        if features.colours_equal:
            ll_error += 3.0
        if not features.colour_consistent:
            ll_error += 8.0
        if features.linkage is True:
            ll_strain += 6.0
            ll_other -= 5.0
        if features.low_kmer is True:
            ll_error += 8.0
    weights = [ll_error, ll_strain, ll_other]
    top = max(weights)
    exp = [math.exp(value - top) for value in weights]
    total = sum(exp)
    probs = [value / total for value in exp]
    labels = ("error", "strain", "other")
    decision = labels[max(range(3), key=lambda index: probs[index])]
    return Posterior(
        bubble_id=features.bubble_id,
        p_error=probs[0],
        p_strain=probs[1],
        p_other=probs[2],
        decision=decision,
        confidence=max(probs),
    )
