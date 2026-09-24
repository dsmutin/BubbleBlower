"""Graph-only debubblers. They never read ground-truth labels.

Each simple bubble is labelled ``error``, ``variation``, or ``multi``.
The edit is then pop the weak branch, retain both branches, or split the
shared source and sink.

``kmer_divergence`` uses pairwise branch identity, length ratio, and coverage
ratio. ``colour_topology`` uses colour entropy, colour equality, colour
disjointness, colour consistency, source and sink degree, and branch count.
"""

from __future__ import annotations

from dataclasses import dataclass

from metametro.errors import ContractError

from bubbleblower.detect import Bubble, detect_bubbles
from bubbleblower.edits import merge_adjacent, pop_branch, split_instance
from bubbleblower.features import BubbleFeatures, extract_features
from bubbleblower.graph import AssemblyGraph

# Pairwise identity at or above this, with unbalanced coverage, is an error.
_HIGH_IDENTITY = 0.85
# Identity below this is moderate sequence divergence.
_MODERATE_IDENTITY = 0.95
# Coverage ratio at or above this is similar enough to keep a strain allele.
_SIMILAR_COVERAGE = 0.25
# Coverage ratio at or below this is unbalanced when identity is high.
_UNBALANCED_COVERAGE = 0.5
# Below the strain-ratio overlap (~0.4). A 1x branch is an error.
_ERROR_COVERAGE_FLOOR = 0.05
# Equal colours at or below this ratio are an error for colour_topology.
_VERY_SMALL_COVERAGE = 0.15
# Every branch must reach this depth before a split is considered.
_HIGH_BRANCH_COVERAGE = 10.0

DEBUBBLERS = ("kmer_divergence", "colour_topology")


@dataclass(frozen=True)
class DebubbleDecision:
    """Label and edit for one bubble. The label is not a ground-truth class."""

    label: str
    action: str
    branch_index: int | None = None


def _branch_sequence(graph: AssemblyGraph, path: tuple[str, ...]) -> str:
    return "".join(graph.unitig(unitig_id).sequence for unitig_id in path)


def _identity(left: str, right: str) -> float:
    """Pairwise identity in ``[0, 1]``. Empty parallel links are identical."""
    if not left and not right:
        return 1.0
    if not left or not right:
        return 0.0
    width = max(len(left), len(right))
    hamming = sum(base == other for base, other in zip(left, right)) / width
    k = 4
    if len(left) < k or len(right) < k:
        return hamming
    left_kmers = {left[index : index + k] for index in range(len(left) - k + 1)}
    right_kmers = {right[index : index + k] for index in range(len(right) - k + 1)}
    union = left_kmers | right_kmers
    jaccard = len(left_kmers & right_kmers) / len(union) if union else 1.0
    return (hamming + jaccard) / 2.0


def _length_ratio(left: int, right: int) -> float:
    if left == 0 and right == 0:
        return 1.0
    if left == 0 or right == 0:
        return 0.0
    return min(left, right) / max(left, right)


def _pair_stats(graph: AssemblyGraph, features: BubbleFeatures) -> tuple[float, float]:
    """Minimum pairwise identity and minimum length ratio across branches."""
    sequences = [_branch_sequence(graph, branch.path) for branch in features.branches]
    identities: list[float] = []
    ratios: list[float] = []
    for index, left in enumerate(sequences):
        for right in sequences[index + 1 :]:
            identities.append(_identity(left, right))
            ratios.append(_length_ratio(len(left), len(right)))
    if not identities:
        return 1.0, 1.0
    return min(identities), min(ratios)


def _weak_index(features: BubbleFeatures) -> int:
    coverages = [branch.mean_coverage for branch in features.branches]
    return min(range(len(coverages)), key=lambda index: coverages[index])


def _dangling_links(graph: AssemblyGraph) -> int:
    """Links whose source or target unitig is missing."""
    unitig_ids = {unitig.unitig_id for unitig in graph.cdbg.unitigs}
    return sum(
        1
        for link in graph.cdbg.links
        if link.source not in unitig_ids or link.target not in unitig_ids
    )


def _groups(bubble: Bubble, end: str) -> list[list[str]]:
    groups: list[list[str]] = []
    for branch in bubble.branches:
        link_id = branch.links[0] if end == "source" else branch.links[-1]
        groups.append([link_id])
    return groups


def _split_end(graph: AssemblyGraph, instance_id: str, groups: list[list[str]]) -> AssemblyGraph | None:
    """Split one node. Return None when the partition is rejected."""
    if len(groups) < 2:
        return None
    before = _dangling_links(graph)
    try:
        edited, _edit = split_instance(graph, instance_id, groups)
    except (ValueError, ContractError):
        return None
    if _dangling_links(edited) > before:
        return None
    return edited


def _split_source_and_sink(graph: AssemblyGraph, bubble: Bubble) -> AssemblyGraph | None:
    """Split both ends. If either partition is rejected, keep the input graph."""
    updated = _split_end(graph, bubble.source, _groups(bubble, "source"))
    if updated is None:
        return None
    updated = _split_end(updated, bubble.sink, _groups(bubble, "sink"))
    if updated is None:
        return None
    return updated


def _classify_kmer(graph: AssemblyGraph, features: BubbleFeatures) -> DebubbleDecision:
    """High identity plus unbalanced coverage pops. Three branches split."""
    if features.n_branches >= 3:
        return DebubbleDecision("multi", "split")
    identity, length_ratio = _pair_stats(graph, features)
    ratio = features.coverage_ratio
    weak = _weak_index(features)
    high_identity_error = (
        identity >= _HIGH_IDENTITY
        and length_ratio >= 0.5
        and ratio <= _UNBALANCED_COVERAGE
    )
    floor_error = ratio <= _ERROR_COVERAGE_FLOOR and length_ratio >= 0.5
    if high_identity_error or floor_error:
        return DebubbleDecision("error", "pop", weak)
    divergent = identity < _MODERATE_IDENTITY or length_ratio < 0.9
    if ratio >= _SIMILAR_COVERAGE and divergent:
        return DebubbleDecision("variation", "retain")
    return DebubbleDecision("variation", "retain")


def _high_coverage(features: BubbleFeatures) -> bool:
    if features.n_branches < 2:
        return False
    return all(branch.mean_coverage >= _HIGH_BRANCH_COVERAGE for branch in features.branches)


def _classify_colour(features: BubbleFeatures) -> DebubbleDecision:
    """Disjoint consistent colours are kept. Equal colours may pop or split."""
    if features.colours_disjoint and features.colour_consistent and features.colour_entropy > 0:
        return DebubbleDecision("variation", "retain")
    if features.colours_equal and features.coverage_ratio <= _VERY_SMALL_COVERAGE:
        return DebubbleDecision("error", "pop", _weak_index(features))
    if features.colours_equal and _high_coverage(features):
        label = "multi" if features.n_branches >= 3 else "variation"
        source_ok = features.source_degree >= 2
        sink_ok = features.sink_degree >= 2
        if source_ok or sink_ok:
            return DebubbleDecision(label, "split")
        return DebubbleDecision(label, "retain")
    return DebubbleDecision("variation", "retain")


def classify_debubble(graph: AssemblyGraph, bubble: Bubble, name: str) -> DebubbleDecision:
    """Classify one bubble. ``name`` is ``kmer_divergence`` or ``colour_topology``."""
    if name not in DEBUBBLERS:
        raise ValueError(f"unknown debubbler: {name}")
    features = extract_features(graph, bubble)
    if name == "kmer_divergence":
        return _classify_kmer(graph, features)
    return _classify_colour(features)


def _apply_split(graph: AssemblyGraph, bubble: Bubble, name: str, features: BubbleFeatures) -> AssemblyGraph | None:
    if name == "kmer_divergence":
        return _split_source_and_sink(graph, bubble)
    updated = graph
    changed = False
    if features.source_degree >= 2:
        source = _split_end(updated, bubble.source, _groups(bubble, "source"))
        if source is not None and _dangling_links(source) <= _dangling_links(graph):
            updated = source
            changed = True
    if features.sink_degree >= 2:
        sink = _split_end(updated, bubble.sink, _groups(bubble, "sink"))
        if sink is not None and _dangling_links(sink) <= _dangling_links(graph):
            updated = sink
            changed = True
    if not changed:
        return None
    return updated


def compact_same_colour(graph: AssemblyGraph) -> AssemblyGraph:
    """Merge simple same-colour links until none remain.

    A link is simple when its source has no other outgoing link and its
    target has no other incoming link. Only forward-forward links are merged.
    Different colours are left apart.
    """
    current = graph.copy()
    for _ in range(len(current.cdbg.links) + 1):
        chosen = None
        outgoing: dict[str, int] = {}
        incoming: dict[str, int] = {}
        for link in current.cdbg.links:
            outgoing[link.source] = outgoing.get(link.source, 0) + 1
            incoming[link.target] = incoming.get(link.target, 0) + 1
        for link in current.cdbg.links:
            if outgoing.get(link.source, 0) != 1 or incoming.get(link.target, 0) != 1:
                continue
            if link.overlap is None or link.orientation != "++":
                continue
            left = current.unitig(link.source).color_ids
            right = current.unitig(link.target).color_ids
            if left and list(left) == list(right):
                chosen = link.link_id
                break
        if chosen is None:
            break
        current, _edit = merge_adjacent(current, chosen, copy_graph=False)
    return current


def resolve_debubbler(graph: AssemblyGraph, name: str) -> AssemblyGraph:
    """Return a new graph edited by the named debubbler. The input is copied."""
    if name not in DEBUBBLERS:
        raise ValueError(f"unknown debubbler: {name}")
    current = graph.copy()
    skipped: set[str] = set()
    for _ in range(max(8, len(current.cdbg.unitigs) * 4)):
        acted = False
        for bubble in detect_bubbles(current):
            if bubble.bubble_id in skipped:
                continue
            decision = classify_debubble(current, bubble, name)
            if decision.action == "retain":
                skipped.add(bubble.bubble_id)
                continue
            if decision.action == "pop" and decision.branch_index is not None:
                current, _edit = pop_branch(current, bubble, decision.branch_index)
                acted = True
                break
            if decision.action == "split":
                features = extract_features(current, bubble)
                edited = _apply_split(current, bubble, name, features)
                if edited is None:
                    skipped.add(bubble.bubble_id)
                    continue
                current = edited
                acted = True
                break
            skipped.add(bubble.bubble_id)
        if not acted:
            break
    return compact_same_colour(current)
