"""Load a generalised example from a MetaMetro benchbuild directory.

Mandatory tests keep their own fixtures. These helpers are for the example
runners, which must use the bench build rather than a second copy of the graph.
"""

from __future__ import annotations

from pathlib import Path

from bubbleblower.graph import AssemblyGraph, from_cdbg


def ensure_bench(name: str) -> Path:
    """Return the bench directory, building an in-process graph when it is absent."""
    try:
        from metametro.bench.build import build
        from metametro.bench.paths import default_outdir
        from metametro.bench.registry import resolve
    except ImportError as exc:
        raise SystemExit(
            "metametro.bench is not importable. Install MetaMetro or set PYTHONPATH to its src."
        ) from exc
    spec = resolve(name)
    destination = default_outdir(spec)
    if (destination / "cfa" / "metadata.yaml").is_file() and (destination / "cdbg" / "metadata.yaml").is_file():
        return destination
    if destination.exists() and any(destination.iterdir()):
        raise SystemExit(
            f"bench {spec.name} at {destination} is incomplete. Remove that directory and rerun "
            f"metametro benchbuild {spec.name}"
        )
    build(spec.name, outdir=destination)
    return destination


def load_bench_graph(name: str) -> tuple[Path, AssemblyGraph]:
    """Build or reuse ``name`` and copy CFA coverage onto the compacted graph."""
    from metametro.formats.cdbg.io import load_cdbg
    from metametro.formats.cfa.io import load_cfa

    root = ensure_bench(name)
    cfa = load_cfa(root / "cfa")
    cdbg = load_cdbg(root / "cdbg")
    node_column = {row["node_id"]: float(row["coverage"]) for row in cfa.nodes if "coverage" in row}
    edge_column = {row["edge_id"]: float(row["coverage"]) for row in cfa.edges if "coverage" in row}
    node_coverage: dict[str, float] = {}
    for unitig in cdbg.unitigs:
        values = [node_column[member] for member in unitig.members if member in node_column]
        if values:
            node_coverage[unitig.unitig_id] = sum(values) / len(values)
    link_coverage = {
        link.link_id: edge_column[link.link_id] for link in cdbg.links if link.link_id in edge_column
    }
    return root, from_cdbg(
        cdbg,
        node_coverage=node_coverage or None,
        link_coverage=link_coverage or None,
    )
