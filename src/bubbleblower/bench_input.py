"""Load a generalised example from a MetaMetro benchbuild directory.

Mandatory tests keep their own fixtures. These helpers are for the example
runners, which must use the bench build rather than a second copy of the graph.
"""

from __future__ import annotations

import os
from pathlib import Path

from bubbleblower.graph import AssemblyGraph, from_cdbg


BENCH_ROOT_ENV = "BUBBLEBLOWER_BENCH_ROOT"


def _bench_root(root: Path | str | None) -> Path | None:
    """``root``, else ``$BUBBLEBLOWER_BENCH_ROOT``, else ``None`` (MetaMetro repository)."""
    if root is not None:
        return Path(root)
    value = os.environ.get(BENCH_ROOT_ENV)
    return Path(value) if value else None


def ensure_bench(name: str, *, root: Path | str | None = None) -> Path:
    """Return the bench directory, building an in-process graph when it is absent.

    The directory is ``{root}/data/bench/{name}/...`` as in MetaMetro
    ``default_outdir``. ``root`` falls back to ``$BUBBLEBLOWER_BENCH_ROOT`` and
    then to the MetaMetro repository.
    """
    try:
        from metametro.bench.build import build
        from metametro.bench.paths import default_outdir
        from metametro.bench.registry import resolve
    except ImportError as exc:
        raise SystemExit(
            "metametro.bench is not importable. Install MetaMetro or set PYTHONPATH to its src."
        ) from exc
    spec = resolve(name)
    destination = default_outdir(spec, root=_bench_root(root))
    if (destination / "cfa" / "metadata.yaml").is_file() and (destination / "cdbg" / "metadata.yaml").is_file():
        return destination
    if destination.exists() and any(destination.iterdir()):
        raise SystemExit(
            f"bench {spec.name} at {destination} is incomplete. Remove that directory and rerun "
            f"metametro benchbuild {spec.name}"
        )
    build(spec.name, outdir=destination)
    return destination


def load_bench_graph(
    name: str, *, namespaces: tuple[str, ...] | None = None, root: Path | str | None = None
) -> tuple[Path, AssemblyGraph]:
    """Build or reuse ``name`` and copy CFA coverage onto the compacted graph.

    ``namespaces`` keeps those colour layers from the MetaMetro ToCUMG.
    ``root`` is passed to :func:`ensure_bench`.
    """
    from metametro.contracts.colour_filter import filter_colours
    from metametro.formats.cdbg.io import load_cdbg
    from metametro.formats.cfa.io import load_cfa

    root = ensure_bench(name, root=root)
    cfa = load_cfa(root / "cfa")
    cdbg = load_cdbg(root / "cdbg")
    if namespaces:
        cfa = filter_colours(cfa, namespaces=namespaces)
        cdbg = filter_colours(cdbg, namespaces=namespaces)
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
