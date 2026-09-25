"""MetaMetro CDBG plus coverage and instance lineage.

Unitig ids are graph-instance ids. Sequence hashes are never used as ids.
Coverage lives beside the CDBG: MetaMetro unitigs do not store it.
"""

from __future__ import annotations

import copy
import hashlib
import re
from dataclasses import dataclass, field

from metametro.converters.cdbg_to_cgt import cdbg_to_cgt
from metametro.converters.cfa_to_cdbg import cfa_to_cdbg
from metametro.formats.cdbg.model import Cdbg, Unitig
from metametro.formats.cdbg.validator import validate_cdbg
from metametro.formats.cfa.model import CfaGraph
from metametro.formats.cgt.model import Cgt
from metametro.identity import assert_cgt_matches_cdbg

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
    coverage_source: str = "caller"
    _seq: int = 0

    def copy(self) -> AssemblyGraph:
        """Return a deep copy. The caller keeps the previous state."""
        return copy.deepcopy(self)

    def unitig(self, unitig_id: str) -> Unitig:
        for unitig in self.cdbg.unitigs:
            if unitig.unitig_id == unitig_id:
                return unitig
        raise KeyError(unitig_id)

    def unitig_id_for_member(self, cfa_node_id: str) -> str:
        """Return the ToCUMG unitig that contains one CFA node id."""
        hits = [unitig.unitig_id for unitig in self.cdbg.unitigs if cfa_node_id in unitig.members]
        if len(hits) != 1:
            raise KeyError(cfa_node_id)
        return hits[0]

    def member_ids(self) -> set[str]:
        """CFA node ids stored on this ToCUMG. Unitig ids are MetaMetro's."""
        return {member for unitig in self.cdbg.unitigs for member in unitig.members}

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
    if supplied_nodes:
        graph.coverage_source = "caller"
    elif any(unitig.unitig_id in graph.node_coverage for unitig in graph.cdbg.unitigs):
        graph.coverage_source = "unitig_name_cov_token"
    else:
        graph.coverage_source = "absent"
    return graph


def _column(rows: list[dict[str, str]], name: str) -> dict[str, str] | None:
    if not rows or name not in rows[0]:
        return None
    key = "node_id" if "node_id" in rows[0] else "edge_id"
    return {row[key]: row[name] for row in rows}


def adopt_cfa(cfa: CfaGraph) -> AssemblyGraph:
    """Wrap the ToCUMG that MetaMetro compacts from ``cfa``.

    Unitig ids come from ``cfa_to_cdbg``. This function does not invent nodes
    or links. A ``coverage`` column, when the CFA declares one, is copied onto
    unitigs (mean of member nodes) and onto links. A link with no coverage
    column takes the source unitig's coverage.
    """
    cdbg = cfa_to_cdbg(cfa)
    node_column = _column(cfa.nodes, "coverage")
    edge_column = _column(cfa.edges, "coverage")
    node_coverage: dict[str, float] = {}
    if node_column is not None:
        for unitig in cdbg.unitigs:
            values = [float(node_column[member]) for member in unitig.members]
            node_coverage[unitig.unitig_id] = sum(values) / len(values)
    link_coverage: dict[str, float] = {}
    if edge_column is not None:
        for link in cdbg.links:
            if link.link_id in edge_column:
                link_coverage[link.link_id] = float(edge_column[link.link_id])
    if node_coverage and len(link_coverage) != len(cdbg.links):
        for link in cdbg.links:
            link_coverage.setdefault(link.link_id, node_coverage[link.source])
    graph = from_cdbg(
        cdbg,
        node_coverage or None,
        link_coverage or None,
    )
    if node_coverage:
        graph.coverage_source = "cfa_coverage_column"
    graph.validate()
    return graph


def from_cgt(cgt: Cgt, cdbg: Cdbg) -> AssemblyGraph:
    """Bind a graph tensor to the ToCUMG it was built from.

    The tensor is checked against that ToCUMG and is not turned into a second
    graph. Topology, sequences, and colours stay on the ToCUMG.
    """
    assert_cgt_matches_cdbg(cgt, cdbg)
    return from_cdbg(cdbg)


def as_cgt(graph: AssemblyGraph) -> Cgt:
    """Tensor view of an existing ToCUMG. Colours stay off the feature matrix."""
    return cdbg_to_cgt(graph.cdbg)


def _color_set(colors: list[int]) -> str:
    return ",".join(str(color) for color in colors)


def records_to_tocumg(
    *,
    graph_id: str,
    nodes: list[dict],
    links: list[dict],
    colors: list[dict[str, str]],
) -> AssemblyGraph:
    """Hand node and link records to MetaMetro and adopt the compacted ToCUMG.

    Records are a CFA. ``cfa_to_cdbg`` assigns unitig ids. There is no overlap
    column, so compaction stays one unitig per CFA node.
    """
    if not nodes:
        raise ValueError("a ToCUMG needs at least one CFA node")
    sequences = {str(node["id"]): str(node["sequence"]) for node in nodes}
    node_rows = [
        {
            "node_id": str(node["id"]),
            "coverage": str(float(node["coverage"])),
            "color_set": _color_set([int(color) for color in node["colors"]]),
        }
        for node in nodes
    ]
    edge_rows = [
        {
            "edge_id": str(link["id"]),
            "source": str(link["source"]),
            "target": str(link["target"]),
            "orientation": str(link.get("orientation", "++")),
            "coverage": str(float(link["coverage"])),
            "color_set": _color_set([int(color) for color in link["colors"]]),
        }
        for link in links
    ]
    cfa = CfaGraph(
        metadata={
            "schema_version": "1.0",
            "graph_id": graph_id,
            "graph_type": "repeat",
            "contract": "metagenome_to_graph",
            "contract_version": "1.0",
            "features": {
                "node": {"coverage": "float", "color_set": "color_set"},
                "edge": {
                    "orientation": "orientation",
                    "coverage": "float",
                    "color_set": "color_set",
                },
            },
        },
        sequences=sequences,
        nodes=node_rows,
        edges=edge_rows,
        colors=colors,
        node_header=["node_id", "coverage", "color_set"],
        edge_header=["edge_id", "source", "target", "orientation", "coverage", "color_set"],
    )
    return adopt_cfa(cfa)


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
