"""Reversible graph edits: pop, duplicate, split, merge.

Each edit returns a new graph. The input graph is not mutated.
Inverse edits restore sequence, topology, colours, coverage, and provenance.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field

from metametro.formats.cdbg.model import Link, NodeMap, Unitig

from bubbleblower.detect import Bubble
from bubbleblower.graph import AssemblyGraph, sequence_hash


@dataclass
class Edit:
    """One applied operation and the payload that inverts it."""

    edit_id: str
    edit_type: str
    source_ids: tuple[str, ...]
    target_ids: tuple[str, ...]
    reason: str
    payload: dict = field(default_factory=dict)


def _remember(graph: AssemblyGraph, edit_type: str) -> str:
    graph._seq += 1
    return f"{edit_type}-{graph._seq:06d}"


def _add_unitig(
    graph: AssemblyGraph,
    *,
    sequence: str,
    color_ids: list[int],
    coverage: float,
    parent: str | None,
    unitig_id: str | None = None,
    cfa_id: str | None = None,
) -> str:
    unitig_id = unitig_id or graph.alloc_unitig_id()
    cfa_id = cfa_id or graph.alloc_cfa_id()
    graph.cdbg.unitigs.append(
        Unitig(unitig_id=unitig_id, sequence=sequence, members=[cfa_id], color_ids=list(color_ids))
    )
    graph.cdbg.mapping.append(
        NodeMap(
            cfa_node_id=cfa_id,
            unitig_id=unitig_id,
            ordinal=0,
            length=len(sequence),
            color_ids=list(color_ids),
        )
    )
    graph.node_coverage[unitig_id] = coverage
    graph.lineage[unitig_id] = {
        "sequence_hash": sequence_hash(sequence),
        "parent_instance": parent,
    }
    return unitig_id


def _drop_unitig(graph: AssemblyGraph, unitig_id: str) -> None:
    graph.cdbg.unitigs = [unitig for unitig in graph.cdbg.unitigs if unitig.unitig_id != unitig_id]
    graph.cdbg.mapping = [row for row in graph.cdbg.mapping if row.unitig_id != unitig_id]
    graph.node_coverage.pop(unitig_id, None)
    graph.lineage.pop(unitig_id, None)


def pop_branch(graph: AssemblyGraph, bubble: Bubble, branch_index: int) -> tuple[AssemblyGraph, Edit]:
    """Remove one bubble branch. Source and sink coverage stay as observed."""
    if branch_index < 0 or branch_index >= len(bubble.branches):
        raise IndexError(branch_index)
    updated = graph.copy()
    branch = bubble.branches[branch_index]
    removed_ids = set(branch.path)
    link_ids = set(branch.links)
    payload = {
        "unitigs": [copy.deepcopy(updated.unitig(unitig_id)) for unitig_id in branch.path],
        "links": [copy.deepcopy(updated.link(link_id)) for link_id in branch.links],
        "mapping": [
            copy.deepcopy(row) for row in updated.cdbg.mapping if row.unitig_id in removed_ids
        ],
        "node_coverage": {unitig_id: updated.node_coverage[unitig_id] for unitig_id in branch.path},
        "link_coverage": {link_id: updated.link_coverage[link_id] for link_id in branch.links},
        "lineage": {unitig_id: dict(updated.lineage[unitig_id]) for unitig_id in branch.path},
    }
    updated.cdbg.unitigs = [unitig for unitig in updated.cdbg.unitigs if unitig.unitig_id not in removed_ids]
    updated.cdbg.links = [link for link in updated.cdbg.links if link.link_id not in link_ids]
    updated.cdbg.mapping = [row for row in updated.cdbg.mapping if row.unitig_id not in removed_ids]
    for unitig_id in branch.path:
        updated.node_coverage.pop(unitig_id, None)
        updated.lineage.pop(unitig_id, None)
    for link_id in branch.links:
        updated.link_coverage.pop(link_id, None)
    edit = Edit(
        edit_id=_remember(updated, "pop"),
        edit_type="pop",
        source_ids=tuple(branch.path),
        target_ids=(),
        reason=f"pop branch {branch_index} of {bubble.bubble_id}",
        payload=payload,
    )
    updated.validate()
    return updated, edit


def restore_pop(graph: AssemblyGraph, edit: Edit) -> AssemblyGraph:
    """Inverse of ``pop_branch``."""
    if edit.edit_type != "pop":
        raise ValueError(edit.edit_type)
    updated = graph.copy()
    updated.cdbg.unitigs.extend(copy.deepcopy(edit.payload["unitigs"]))
    updated.cdbg.links.extend(copy.deepcopy(edit.payload["links"]))
    updated.cdbg.mapping.extend(copy.deepcopy(edit.payload["mapping"]))
    updated.node_coverage.update(edit.payload["node_coverage"])
    updated.link_coverage.update(edit.payload["link_coverage"])
    updated.lineage.update({key: dict(value) for key, value in edit.payload["lineage"].items()})
    updated.validate()
    return updated


def duplicate_instance(graph: AssemblyGraph, instance_id: str) -> tuple[AssemblyGraph, Edit]:
    """Second instance of the same sequence. Coverage is split in half.

    Links stay on the original instance. Sequence hashes match; instance ids do not.
    """
    updated = graph.copy()
    unitig = updated.unitig(instance_id)
    coverage = updated.node_coverage[instance_id]
    left = coverage / 2.0
    right = coverage - left
    updated.node_coverage[instance_id] = left
    new_id = _add_unitig(
        updated,
        sequence=unitig.sequence,
        color_ids=list(unitig.color_ids),
        coverage=right,
        parent=instance_id,
    )
    edit = Edit(
        edit_id=_remember(updated, "duplicate"),
        edit_type="duplicate",
        source_ids=(instance_id,),
        target_ids=(instance_id, new_id),
        reason=f"duplicate {instance_id}",
        payload={"original": instance_id, "created": new_id, "coverage": coverage},
    )
    updated.validate()
    return updated, edit


def merge_instances(graph: AssemblyGraph, instance_ids: list[str], target_id: str | None = None) -> tuple[AssemblyGraph, Edit]:
    """Merge instances that share one sequence. Coverage is summed. Colours are united.

    Parallel links that land on the same endpoints after the merge are collapsed.
    """
    if len(instance_ids) < 2:
        raise ValueError("merge needs two instances")
    updated = graph.copy()
    ids = list(dict.fromkeys(instance_ids))
    sequences = {updated.unitig(unitig_id).sequence for unitig_id in ids}
    if len(sequences) != 1:
        raise ValueError("merge requires identical sequences")
    target = target_id or ids[0]
    if target not in ids:
        raise ValueError("target must be one of the merged instances")
    others = [unitig_id for unitig_id in ids if unitig_id != target]
    snapshot = {
        "unitigs": [copy.deepcopy(updated.unitig(unitig_id)) for unitig_id in ids],
        "mapping": [copy.deepcopy(row) for row in updated.cdbg.mapping if row.unitig_id in set(ids)],
        "links": [
            copy.deepcopy(link)
            for link in updated.cdbg.links
            if link.source in ids or link.target in ids
        ],
        "node_coverage": {unitig_id: updated.node_coverage[unitig_id] for unitig_id in ids},
        "link_coverage": {
            link.link_id: updated.link_coverage[link.link_id]
            for link in updated.cdbg.links
            if link.source in ids or link.target in ids
        },
        "lineage": {unitig_id: dict(updated.lineage[unitig_id]) for unitig_id in ids},
    }
    colours: set[int] = set()
    coverage = 0.0
    for unitig_id in ids:
        colours.update(updated.unitig(unitig_id).color_ids)
        coverage += updated.node_coverage[unitig_id]
    target_unitig = updated.unitig(target)
    target_unitig.color_ids = sorted(colours)
    for row in updated.cdbg.mapping:
        if row.unitig_id == target:
            row.color_ids = sorted(colours)
    updated.node_coverage[target] = coverage
    remap = {unitig_id: target for unitig_id in others}
    for link in updated.cdbg.links:
        if link.source in remap:
            link.source = remap[link.source]
        if link.target in remap:
            link.target = remap[link.target]
    _collapse_parallel(updated)
    for unitig_id in others:
        _drop_unitig(updated, unitig_id)
    edit = Edit(
        edit_id=_remember(updated, "merge"),
        edit_type="merge",
        source_ids=tuple(ids),
        target_ids=(target,),
        reason=f"merge {','.join(ids)}",
        payload=snapshot,
    )
    updated.validate()
    return updated, edit


def _collapse_parallel(graph: AssemblyGraph) -> None:
    """Sum coverage of links that share endpoints, orientation, and colours."""
    grouped: dict[tuple, list[Link]] = {}
    for link in graph.cdbg.links:
        key = (link.source, link.target, link.orientation, tuple(link.color_ids))
        grouped.setdefault(key, []).append(link)
    kept: list[Link] = []
    coverage: dict[str, float] = {}
    for links in grouped.values():
        head = links[0]
        total = sum(graph.link_coverage[link.link_id] for link in links)
        kept.append(head)
        coverage[head.link_id] = total
        for link in links[1:]:
            graph.link_coverage.pop(link.link_id, None)
    graph.cdbg.links = kept
    graph.link_coverage = coverage


def split_instance(graph: AssemblyGraph, instance_id: str, groups: list[list[str]]) -> tuple[AssemblyGraph, Edit]:
    """Replace one instance with one instance per link group.

    Each group lists link ids that move to that instance. Incident links
    left out of every group are copied onto each new instance, with coverage
    split in proportion to the grouped link coverage.
    """
    if len(groups) < 2:
        raise ValueError("split needs two groups")
    updated = graph.copy()
    unitig = updated.unitig(instance_id)
    named = [link_id for group in groups for link_id in group]
    if len(named) != len(set(named)):
        raise ValueError("a link cannot sit in two split groups")
    incident = [
        link.link_id
        for link in updated.cdbg.links
        if link.source == instance_id or link.target == instance_id
    ]
    unknown = [link_id for link_id in named if link_id not in incident]
    if unknown:
        raise ValueError(f"links are not incident to {instance_id}: {unknown}")
    leftover = [link_id for link_id in incident if link_id not in set(named)]
    weights = []
    for group in groups:
        weights.append(sum(updated.link_coverage[link_id] for link_id in group) or 1.0)
    weight_total = sum(weights)
    snapshot = {
        "unitig": copy.deepcopy(unitig),
        "mapping": [copy.deepcopy(row) for row in updated.cdbg.mapping if row.unitig_id == instance_id],
        "links": [copy.deepcopy(updated.link(link_id)) for link_id in incident],
        "node_coverage": updated.node_coverage[instance_id],
        "link_coverage": {link_id: updated.link_coverage[link_id] for link_id in incident},
        "lineage": dict(updated.lineage[instance_id]),
        "created": [],
    }
    parent_coverage = updated.node_coverage[instance_id]
    new_ids: list[str] = []
    for weight in weights:
        new_ids.append(
            _add_unitig(
                updated,
                sequence=unitig.sequence,
                color_ids=list(unitig.color_ids),
                coverage=parent_coverage * weight / weight_total,
                parent=instance_id,
            )
        )
    snapshot["created"] = list(new_ids)
    for new_id, group in zip(new_ids, groups):
        for link_id in group:
            link = updated.link(link_id)
            if link.source == instance_id:
                link.source = new_id
            if link.target == instance_id:
                link.target = new_id
    for link_id in leftover:
        original = updated.link(link_id)
        total = updated.link_coverage[link_id]
        for index, new_id in enumerate(new_ids):
            if index == 0:
                clone = original
                clone_id = link_id
            else:
                clone_id = updated.alloc_link_id()
                clone = Link(
                    link_id=clone_id,
                    source=original.source,
                    target=original.target,
                    orientation=original.orientation,
                    color_ids=list(original.color_ids),
                    overlap=original.overlap,
                )
                updated.cdbg.links.append(clone)
            if clone.source == instance_id:
                clone.source = new_id
            if clone.target == instance_id:
                clone.target = new_id
            updated.link_coverage[clone_id] = total * weights[index] / weight_total
        if link_id not in {link.link_id for link in updated.cdbg.links}:
            pass
    _drop_unitig(updated, instance_id)
    edit = Edit(
        edit_id=_remember(updated, "split"),
        edit_type="split",
        source_ids=(instance_id,),
        target_ids=tuple(new_ids),
        reason=f"split {instance_id}",
        payload=snapshot,
    )
    updated.validate()
    return updated, edit


def restore_split(graph: AssemblyGraph, edit: Edit) -> AssemblyGraph:
    """Inverse of ``split_instance``. Restores the original instance id."""
    if edit.edit_type != "split":
        raise ValueError(edit.edit_type)
    updated = graph.copy()
    created = set(edit.payload["created"])
    touched = {
        link.link_id
        for link in updated.cdbg.links
        if link.source in created or link.target in created
    }
    updated.cdbg.links = [link for link in updated.cdbg.links if link.link_id not in touched]
    for link_id in touched:
        updated.link_coverage.pop(link_id, None)
    for unitig_id in created:
        _drop_unitig(updated, unitig_id)
    unitig = copy.deepcopy(edit.payload["unitig"])
    updated.cdbg.unitigs.append(unitig)
    updated.cdbg.mapping.extend(copy.deepcopy(edit.payload["mapping"]))
    updated.cdbg.links.extend(copy.deepcopy(edit.payload["links"]))
    updated.node_coverage[unitig.unitig_id] = edit.payload["node_coverage"]
    updated.link_coverage.update(edit.payload["link_coverage"])
    updated.lineage[unitig.unitig_id] = dict(edit.payload["lineage"])
    updated.validate()
    return updated


def revert(graph: AssemblyGraph, edit: Edit) -> AssemblyGraph:
    """Apply the inverse edit."""
    if edit.edit_type == "pop":
        return restore_pop(graph, edit)
    if edit.edit_type == "duplicate":
        merged, _merge_edit = merge_instances(graph, list(edit.target_ids), edit.payload["original"])
        return merged
    if edit.edit_type == "split":
        return restore_split(graph, edit)
    if edit.edit_type == "merge":
        return _restore_merge(graph, edit)
    raise ValueError(edit.edit_type)


def _restore_merge(graph: AssemblyGraph, edit: Edit) -> AssemblyGraph:
    updated = graph.copy()
    payload = edit.payload
    target = edit.target_ids[0]
    _drop_unitig(updated, target)
    incident = {
        link.link_id
        for link in updated.cdbg.links
        if link.source == target or link.target == target
    }
    updated.cdbg.links = [link for link in updated.cdbg.links if link.link_id not in incident]
    for link_id in incident:
        updated.link_coverage.pop(link_id, None)
    updated.cdbg.unitigs.extend(copy.deepcopy(payload["unitigs"]))
    updated.cdbg.mapping.extend(copy.deepcopy(payload["mapping"]))
    updated.cdbg.links.extend(copy.deepcopy(payload["links"]))
    updated.node_coverage.update(payload["node_coverage"])
    updated.link_coverage.update(payload["link_coverage"])
    updated.lineage.update({key: dict(value) for key, value in payload["lineage"].items()})
    updated.validate()
    return updated
