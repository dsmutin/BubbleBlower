"""Load a generalised example from a MetaMetro benchbuild directory.

Mandatory tests keep their own fixtures. These helpers are for the example
runners, which must use the bench build rather than a second copy of the graph.
"""

from __future__ import annotations

import os
from pathlib import Path

from bubbleblower.graph import AssemblyGraph, from_cdbg, split_link_coverage


BENCH_ROOT_ENV = "BUBBLEBLOWER_BENCH_ROOT"


def _bench_root(root: Path | str | None) -> Path | None:
    """``root``, else ``$BUBBLEBLOWER_BENCH_ROOT``, else ``None`` (MetaMetro repository)."""
    if root is not None:
        return Path(root)
    value = os.environ.get(BENCH_ROOT_ENV)
    return Path(value) if value else None


def _default_colourings(applied: list[str]) -> bool:
    """Whether manifest ``colourings`` look like a default build (``colourings=None``).

    MetaMetro does not record the request in the manifest. A default build runs
    every automatic colouring that can run, in registration order, starting
    with ``as_built``. An explicit selection either omits ``as_built``, names a
    non-automatic colouring, or breaks that order.
    """
    from metametro.bench.colourings import listed_colourings

    auto = [colouring for colouring, _namespace, is_auto in listed_colourings() if is_auto]
    if not applied or applied[0] != "as_built" or any(item not in auto for item in applied):
        return False
    positions = [auto.index(item) for item in applied]
    return positions == sorted(positions)


def check_reused_build(name: str, destination: Path) -> None:
    """Reject a bench directory that a current MetaMetro would not have written.

    A manifest without ``gold_namespaces`` predates gold ground truth outside
    the colour layers (MetaMetro 0.23.0). When ``canonical_identity(name)`` is
    known and the build used the default colourings, its ``identity.sha256``
    (or the manifest ``identity``) must equal it. A build with an explicit
    colouring selection, or a bench without a canonical identity (roxel, spb,
    communities), is not compared. Raises ``ContractError`` with a rebuild hint.
    """
    from metametro.errors import ContractError
    from metametro.tables import read_yaml

    hint = f"remove {destination} and rerun metametro benchbuild {name}"
    manifest_path = destination / "manifest.yaml"
    manifest = read_yaml(manifest_path) if manifest_path.is_file() else {}
    if "gold_namespaces" not in manifest:
        raise ContractError([f"bench {name} at {destination} has no gold_namespaces in manifest.yaml (stale build); {hint}"])
    try:
        from metametro.bench import canonical_identity
    except ImportError:
        return
    expected = canonical_identity(name)
    if expected is None or not _default_colourings([str(item) for item in manifest.get("colourings") or []]):
        return
    identity_path = destination / "identity.sha256"
    if identity_path.is_file():
        found = (identity_path.read_text(encoding="utf-8").split() or [""])[0]
    else:
        found = str(manifest.get("identity", ""))
    if found != expected:
        raise ContractError(
            [f"bench {name} at {destination} has identity {found or 'none'}, not the canonical {expected} (stale build); {hint}"]
        )


def ensure_bench(name: str, *, root: Path | str | None = None) -> Path:
    """Return the bench directory, building an in-process graph when it is absent.

    The directory is ``{root}/data/bench/{name}/...`` as in MetaMetro
    ``default_outdir``. ``root`` falls back to ``$BUBBLEBLOWER_BENCH_ROOT`` and
    then to the MetaMetro repository. A directory that is already built is
    checked with :func:`check_reused_build` before it is reused.
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
        check_reused_build(spec.name, destination)
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

    Links without a CFA coverage (MEGAHIT builds) get one from
    :func:`bubbleblower.graph.split_link_coverage` when every unitig has coverage.

    ``namespaces`` keeps those colour layers from the MetaMetro ToCUMG. ``None``
    keeps every layer; an empty tuple drops every layer (a colourless baseline).
    A requested namespace that the build does not have raises ``ContractError``.
    ``taxon`` is gold ground truth since MetaMetro 0.20.0; read it with
    ``metametro.bench.load_bench_gold``. ``root`` is passed to :func:`ensure_bench`.
    """
    from metametro.contracts.colour_filter import filter_colours
    from metametro.formats.cdbg.io import load_cdbg
    from metametro.formats.cfa.io import load_cfa

    root = ensure_bench(name, root=root)
    cfa = load_cfa(root / "cfa")
    cdbg = load_cdbg(root / "cdbg")
    if namespaces is not None and len(namespaces) == 0:
        cfa = filter_colours(cfa, color_ids=[])
        cdbg = filter_colours(cdbg, color_ids=[])
    elif namespaces:
        present = {row["namespace"] for row in cfa.colors}
        missing = [name for name in namespaces if name not in present]
        if missing:
            from metametro.errors import ContractError

            hint = "; taxon is gold ground truth, read it with load_bench_gold" if "taxon" in missing else ""
            raise ContractError(
                [f"bench {name} has no colour namespace {', '.join(missing)} (has {', '.join(sorted(present))}){hint}"]
            )
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
    missing = len(cdbg.links) - len(link_coverage)
    if node_coverage and missing and len(node_coverage) == len(cdbg.unitigs):
        link_coverage = split_link_coverage(cdbg, node_coverage, link_coverage)
    graph = from_cdbg(
        cdbg,
        node_coverage=node_coverage or None,
        link_coverage=link_coverage or None,
    )
    if node_coverage and missing and len(node_coverage) == len(cdbg.unitigs):
        graph.coverage_source = f"cfa_coverage_column; {missing} links without coverage split source coverage by target coverage"
    return root, graph
