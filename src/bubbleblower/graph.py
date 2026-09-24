"""MetaMetro CDBG plus coverage and instance lineage.

Unitig ids are graph-instance ids. Sequence hashes are never used as ids.
Coverage lives beside the CDBG: MetaMetro unitigs do not store it.
"""

from __future__ import annotations

import copy
import hashlib
import re
from dataclasses import dataclass, field

from metametro.formats.cdbg.model import SCHEMA_VERSION, Cdbg, Link, NodeMap, Unitig
from metametro.formats.cdbg.validator import validate_cdbg

_COV = re.compile(r"_cov_([0-9]+(?:\.[0-9]+)?)")


def sequence_hash(sequence: str) -> str:
    """Return the SHA-256 hex digest of a DNA sequence."""
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


@dataclass
class AssemblyGraph:
    """One editable totally coloured graph state."""

    cdbg: Cdbg
    node_coverage: dict[str, float] = field(default_factory=dict)
    link_coverage: dict[str, float] = field(default_factory=dict)
    lineage: dict[str, dict[str, str | None]] = field(default_factory=dict)
    _seq: int = 0

    def copy(self) -> AssemblyGraph:
        """Return a deep copy. The caller keeps the previous state."""
        return copy.deepcopy(self)

    def unitig(self, unitig_id: str) -> Unitig:
        for unitig in self.cdbg.unitigs:
            if unitig.unitig_id == unitig_id:
                return unitig
        raise KeyError(unitig_id)

    def link(self, link_id: str) -> Link:
        for link in self.cdbg.links:
            if link.link_id == link_id:
                return link
        raise KeyError(link_id)

    def alloc_unitig_id(self) -> str:
        used = {unitig.unitig_id for unitig in self.cdbg.unitigs}
        while True:
            self._seq += 1
            unitig_id = f"bb{self._seq:06d}"
            if unitig_id not in used:
                return unitig_id

    def alloc_link_id(self) -> str:
        used = {link.link_id for link in self.cdbg.links}
        while True:
            self._seq += 1
            link_id = f"be{self._seq:06d}"
            if link_id not in used:
                return link_id

    def alloc_cfa_id(self) -> str:
        used = {row.cfa_node_id for row in self.cdbg.mapping}
        while True:
            self._seq += 1
            cfa_id = f"bn{self._seq:06d}"
            if cfa_id not in used:
                return cfa_id

    def validate(self) -> None:
        """Reject dangling ids, missing coverage, and broken CDBG invariants."""
        validate_cdbg(self.cdbg)
        unitig_ids = {unitig.unitig_id for unitig in self.cdbg.unitigs}
        link_ids = {link.link_id for link in self.cdbg.links}
        if set(self.node_coverage) != unitig_ids:
            raise ValueError("node coverage keys do not match unitig ids")
        if set(self.link_coverage) != link_ids:
            raise ValueError("link coverage keys do not match link ids")
        if set(self.lineage) != unitig_ids:
            raise ValueError("lineage keys do not match unitig ids")


def _coverage_from_name(*names: str) -> float | None:
    for name in names:
        match = _COV.search(name)
        if match:
            return float(match.group(1))
    return None


def from_cdbg(
    cdbg: Cdbg,
    node_coverage: dict[str, float] | None = None,
    link_coverage: dict[str, float] | None = None,
) -> AssemblyGraph:
    """Wrap a CDBG. Missing coverage stays absent until the caller sets it.

    A ``_cov_`` token in a unitig id or its CFA member id is read when the
    caller does not pass ``node_coverage``. Link coverage is never invented.
    """
    graph = AssemblyGraph(cdbg=copy.deepcopy(cdbg))
    supplied_nodes = node_coverage is not None
    supplied_links = link_coverage is not None
    graph.node_coverage = dict(node_coverage or {})
    graph.link_coverage = dict(link_coverage or {})
    for unitig in graph.cdbg.unitigs:
        graph.lineage[unitig.unitig_id] = {
            "sequence_hash": sequence_hash(unitig.sequence),
            "parent_instance": None,
        }
        if unitig.unitig_id not in graph.node_coverage and not supplied_nodes:
            parsed = _coverage_from_name(unitig.unitig_id, *unitig.members)
            if parsed is not None:
                graph.node_coverage[unitig.unitig_id] = parsed
    if not supplied_links:
        for link in graph.cdbg.links:
            parsed = _coverage_from_name(link.link_id)
            if parsed is not None:
                graph.link_coverage[link.link_id] = parsed
    graph._seq = len(graph.cdbg.unitigs) + len(graph.cdbg.links) + len(graph.cdbg.mapping)
    return graph


def build_graph(
    *,
    graph_id: str,
    nodes: list[dict],
    links: list[dict],
    colors: list[dict[str, str]],
) -> AssemblyGraph:
    """Build an identity-compaction CDBG from plain node and link records.

    Each node dict has ``id``, ``sequence``, ``colors``, and ``coverage``.
    Each link dict has ``id``, ``source``, ``target``, ``colors``, and ``coverage``.
    """
    unitigs: list[Unitig] = []
    mapping: list[NodeMap] = []
    node_coverage: dict[str, float] = {}
    for node in nodes:
        unitig_id = str(node["id"])
        sequence = str(node["sequence"])
        color_ids = [int(color) for color in node["colors"]]
        cfa_id = str(node.get("cfa_id", unitig_id))
        unitigs.append(
            Unitig(
                unitig_id=unitig_id,
                sequence=sequence,
                members=[cfa_id],
                color_ids=color_ids,
            )
        )
        mapping.append(
            NodeMap(
                cfa_node_id=cfa_id,
                unitig_id=unitig_id,
                ordinal=0,
                length=len(sequence),
                color_ids=list(color_ids),
            )
        )
        node_coverage[unitig_id] = float(node["coverage"])
    link_rows: list[Link] = []
    link_coverage: dict[str, float] = {}
    for link in links:
        link_id = str(link["id"])
        color_ids = [int(color) for color in link["colors"]]
        link_rows.append(
            Link(
                link_id=link_id,
                source=str(link["source"]),
                target=str(link["target"]),
                orientation=link.get("orientation", "++"),
                color_ids=color_ids,
            )
        )
        link_coverage[link_id] = float(link["coverage"])
    cdbg = Cdbg(
        metadata={
            "schema_version": SCHEMA_VERSION,
            "graph_id": graph_id,
            "graph_type": "repeat",
            "contract": "bubbleblower",
            "contract_version": "1.0",
            "compaction": "identity",
        },
        k=None,
        unitigs=unitigs,
        links=link_rows,
        mapping=mapping,
        colors=colors,
    )
    graph = from_cdbg(cdbg, node_coverage, link_coverage)
    graph.validate()
    return graph


def semantic_signature(graph: AssemblyGraph) -> tuple:
    """Identity-independent view used for round-trip checks.

    Nodes are matched by sequence, colours, coverage, and CFA provenance.
    Links are matched by the provenance of their endpoints, not by id order.
    """
    prov = {unitig.unitig_id: tuple(unitig.members) for unitig in graph.cdbg.unitigs}
    nodes = tuple(
        sorted(
            (
                unitig.sequence,
                tuple(unitig.color_ids),
                round(graph.node_coverage[unitig.unitig_id], 6),
                tuple(unitig.members),
                graph.lineage[unitig.unitig_id]["sequence_hash"],
            )
            for unitig in graph.cdbg.unitigs
        )
    )
    links = tuple(
        sorted(
            (
                prov[link.source],
                prov[link.target],
                tuple(link.color_ids),
                round(graph.link_coverage[link.link_id], 6),
                link.orientation,
            )
            for link in graph.cdbg.links
        )
    )
    return nodes, links
