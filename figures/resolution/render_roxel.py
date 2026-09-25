"""Animate fork resolution of the MetaMetro roxel street graph.

A fork is a unitig with two or more outgoing links. The simple-bubble
detector only keeps pairs of paths that meet again, and roxel has a handful
of those. This script splits every fork, one per iteration, until none remain.

The CFA has no coverage column. Every unitig and every link is given
coverage 1 so ``split_instance`` can run. That constant is not a measured
depth.

Usage::

    python figures/resolution/render_roxel.py --cfa /path/to/roxel/cfa --out figures/resolution/roxel.mp4
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from metametro.converters.cfa_to_cdbg import cfa_to_cdbg  # noqa: E402
from metametro.formats.cfa.io import load_cfa  # noqa: E402

from bubbleblower.animate import animate_states, fork_resolution_states  # noqa: E402
from bubbleblower.graph import from_cdbg  # noqa: E402


def _with_unit_coverage(graph):
    """Copy coverage 1 onto every unitig and link. The caller asked for this assumption."""
    missing_nodes = [unitig.unitig_id for unitig in graph.cdbg.unitigs if unitig.unitig_id not in graph.node_coverage]
    missing_links = [link.link_id for link in graph.cdbg.links if link.link_id not in graph.link_coverage]
    if not missing_nodes and not missing_links:
        return graph
    if graph.node_coverage or graph.link_coverage:
        raise ValueError("roxel coverage is partial; refusing to fill the gaps")
    for unitig_id in missing_nodes:
        graph.node_coverage[unitig_id] = 1.0
    for link_id in missing_links:
        graph.link_coverage[link_id] = 1.0
    graph.coverage_source = "unit_assumption"
    return graph


def main(argv: list[str] | None = None) -> int:
    """Load a roxel CFA, split every outgoing fork, and write the video."""
    parser = argparse.ArgumentParser(description="Animate roxel fork resolution.")
    parser.add_argument("--cfa", type=Path, required=True, help="MetaMetro roxel CFA directory")
    parser.add_argument("--out", type=Path, required=True, help="GIF or MP4 path")
    parser.add_argument("--seconds", type=float, default=8.0, help="animation length in seconds")
    args = parser.parse_args(argv)
    if not args.cfa.is_dir():
        print(f"missing CFA directory: {args.cfa}", file=sys.stderr)
        return 2
    graph = _with_unit_coverage(from_cdbg(cfa_to_cdbg(load_cfa(args.cfa))))
    frames, edits = fork_resolution_states(graph)
    print(f"edits {len(edits)} frames {len(frames)}", flush=True)
    animate_states(frames, edits, args.out, namespace="type", seconds=args.seconds, seed=0)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
