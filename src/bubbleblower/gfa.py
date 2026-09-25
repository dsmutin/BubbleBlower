"""Load an assembly GFA as a MetaMetro ToCUMG.

Topology comes from ``gfa_to_cfa`` and ``cfa_to_cdbg``. Coverage is read from
``DP``, Flye's ``dp``, or ``KC`` tags and stored beside that ToCUMG. Missing
tags stay absent and are not imputed.
"""

from __future__ import annotations

from pathlib import Path

from metametro.contracts.assembly import gfa_to_cfa
from metametro.converters.cfa_to_cdbg import cfa_to_cdbg

from bubbleblower.graph import AssemblyGraph, from_cdbg


def _tags(fields: list[str]) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for field in fields:
        if ":" not in field:
            continue
        key, _type, value = field.split(":", 2)
        parsed[key] = value
    return parsed


def _coverage(tags: dict[str, str], length: int) -> float | None:
    for key in ("DP", "dp"):
        if key in tags:
            return float(tags[key])
    if "KC" in tags and length > 0:
        return float(tags["KC"]) / length
    return None


def _segment_coverage(path: Path) -> dict[str, float]:
    """Read depth tags. A segment without one is an error."""
    coverage: dict[str, float] = {}
    missing: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw or raw.startswith("#"):
            continue
        fields = raw.split("\t")
        if fields[0] != "S" or len(fields) < 3:
            continue
        name = fields[1]
        sequence = fields[2]
        value = _coverage(_tags(fields[3:]), len(sequence))
        if value is None:
            missing.append(name)
        else:
            coverage[name] = value
    if missing:
        raise ValueError(
            f"{len(missing)} GFA segments have no DP or KC coverage; refusing to impute"
        )
    return coverage


def load_gfa(path: str | Path, *, graph_id: str = "gfa") -> tuple[AssemblyGraph, dict[str, str]]:
    """Return the MetaMetro ToCUMG and the CFA-segment to unitig map.

    ``gfa_to_cfa`` builds the graph. BubbleBlower only attaches coverage that
    the GFA already stored.
    """
    gfa_path = Path(path)
    coverage = _segment_coverage(gfa_path)
    cfa = gfa_to_cfa(gfa_path, graph_id=graph_id, graph_type="repeat")
    cdbg = cfa_to_cdbg(cfa)
    node_coverage: dict[str, float] = {}
    for unitig in cdbg.unitigs:
        values = [coverage[member] for member in unitig.members]
        node_coverage[unitig.unitig_id] = sum(values) / len(values)
    link_coverage = {link.link_id: node_coverage[link.source] for link in cdbg.links}
    graph = from_cdbg(cdbg, node_coverage, link_coverage)
    graph.coverage_source = "gfa_dp_or_kc; links without RC use source DP"
    graph.validate()
    return graph, {row.cfa_node_id: row.unitig_id for row in graph.cdbg.mapping}
