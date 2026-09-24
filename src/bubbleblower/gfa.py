"""Load an assembly GFA into a MetaMetro CDBG.

Segments become unitigs. Links become edges. Coverage is read from ``DP``,
Flye's ``dp``, or ``KC`` tags when the GFA stores them. Missing tags stay absent.
"""

from __future__ import annotations

from pathlib import Path

from bubbleblower.graph import AssemblyGraph, build_graph


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


def load_gfa(path: str | Path, *, graph_id: str = "gfa") -> tuple[AssemblyGraph, dict[str, str]]:
    """Return a graph and the raw segment-id map.

    Sequences that contain bases outside ACGTN are skipped, and links that
    touch a skipped segment are skipped. Coverage is not invented.
    """
    segments: dict[str, str] = {}
    node_cov: dict[str, float] = {}
    missing_node_cov: list[str] = []
    links: list[dict] = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        if not raw or raw.startswith("#"):
            continue
        fields = raw.split("\t")
        kind = fields[0]
        if kind == "S" and len(fields) >= 3:
            name = fields[1]
            sequence = fields[2].upper()
            if any(base not in "ACGTN" for base in sequence) or sequence == "":
                continue
            segments[name] = sequence
            coverage = _coverage(_tags(fields[3:]), len(sequence))
            if coverage is None:
                missing_node_cov.append(name)
            else:
                node_cov[name] = coverage
        elif kind == "L" and len(fields) >= 6:
            tags = _tags(fields[6:])
            coverage = float(tags["RC"]) if "RC" in tags else None
            overlap_token = "".join(ch for ch in fields[5] if ch.isdigit())
            overlap = int(overlap_token) if overlap_token else None
            links.append(
                {
                    "id": f"L{len(links) + 1:07d}",
                    "source": fields[1],
                    "target": fields[3],
                    "orientation": f"{fields[2]}{fields[4]}",
                    "colors": [],
                    "coverage": coverage,
                    "overlap": overlap,
                }
            )
    if missing_node_cov:
        raise ValueError(
            f"{len(missing_node_cov)} GFA segments have no DP or KC coverage; refusing to impute"
        )
    kept_links = []
    for link in links:
        if link["source"] not in segments or link["target"] not in segments:
            continue
        if link["coverage"] is None:
            link["coverage"] = node_cov[link["source"]]
        kept_links.append(link)
    nodes = [
        {"id": name, "sequence": sequence, "colors": [], "coverage": node_cov[name]}
        for name, sequence in segments.items()
    ]
    # Empty colour dictionaries are not valid rows; colours are attached later.
    graph = build_graph(graph_id=graph_id, nodes=nodes, links=kept_links, colors=[])
    overlaps = {link["id"]: link.get("overlap") for link in kept_links}
    for link in graph.cdbg.links:
        link.overlap = overlaps.get(link.link_id)
    graph.coverage_source = "gfa_dp_or_kc; links without RC use source DP"
    graph.cdbg.metadata["graph_type"] = "assembly"
    graph.cdbg.metadata["gfa"] = Path(path).name
    return graph, {name: name for name in segments}
