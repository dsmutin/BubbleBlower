"""Load a MEGAHIT FASTG into a MetaMetro CDBG.

Only forward records are kept. A neighbour whose name ends with ``'`` is the
reverse strand and is stored with orientation ``+-``. Coverage is the
``cov_`` field MEGAHIT writes into the header. Records without that field
are refused.
"""

from __future__ import annotations

import re
from pathlib import Path

from bubbleblower.graph import AssemblyGraph, build_graph

_COV = re.compile(r"_cov_([0-9]+(?:\.[0-9]+)?)")


def _header_parts(header: str) -> tuple[str, list[str]]:
    body = header[1:].strip().rstrip(";")
    if ":" not in body:
        return body, []
    name, rest = body.split(":", 1)
    neighbours = [item for item in rest.split(",") if item]
    return name, neighbours


def load_fastg(path: str | Path, *, graph_id: str = "fastg") -> AssemblyGraph:
    """Read forward unitigs and their overlap edges from a FASTG file."""
    nodes: dict[str, tuple[str, float]] = {}
    edges: dict[str, list[str]] = {}
    name: str | None = None
    chunks: list[str] = []
    pending: list[str] = []

    def flush() -> None:
        if name is None or name.endswith("'"):
            return
        sequence = "".join(chunks).upper()
        if not sequence or any(base not in "ACGTN" for base in sequence):
            return
        match = _COV.search(name)
        if match is None:
            raise ValueError(f"FASTG record {name} has no cov_ field; refusing to impute coverage")
        nodes[name] = (sequence, float(match.group(1)))
        edges[name] = list(pending)

    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            name, pending = _header_parts(line)
            chunks = []
            continue
        chunks.append(line)
    flush()
    if not nodes:
        raise ValueError(f"no forward FASTG records in {path}")
    links = []
    seen: set[tuple[str, str, str]] = set()
    for source, neighbours in edges.items():
        if source not in nodes:
            continue
        for neighbour in neighbours:
            reverse = neighbour.endswith("'")
            target = neighbour[:-1] if reverse else neighbour
            if target not in nodes or target == source:
                continue
            orientation = "+-" if reverse else "++"
            key = (source, target, orientation)
            if key in seen:
                continue
            seen.add(key)
            links.append(
                {
                    "id": f"L{len(links) + 1:07d}",
                    "source": source,
                    "target": target,
                    "orientation": orientation,
                    "colors": [],
                    "coverage": nodes[source][1],
                    "overlap": None,
                }
            )
    graph = build_graph(
        graph_id=graph_id,
        nodes=[
            {"id": node_id, "sequence": sequence, "colors": [], "coverage": coverage}
            for node_id, (sequence, coverage) in nodes.items()
        ],
        links=links,
        colors=[],
    )
    graph.coverage_source = "megahit_fastg_cov"
    graph.cdbg.metadata["graph_type"] = "assembly"
    graph.cdbg.metadata["fastg"] = Path(path).name
    return graph
