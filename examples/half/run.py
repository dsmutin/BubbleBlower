#!/usr/bin/env python3
"""Colour the pre-pop metaSPAdes graph and compare contig metrics.

BubbleBlower walks the graph with bulge removal left off. The baseline is
metaSPAdes with its own bulge remover. Both are scored with minimap2 against
the three source genomes. Missing files stop the run.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
_METAMETRO = ROOT.parent / "metametro" / "src"
if _METAMETRO.is_dir():
    sys.path.insert(0, str(_METAMETRO))

from bubbleblower.assembly_metrics import quast_like  # noqa: E402
from bubbleblower.colour_reads import colour_from_fastq  # noqa: E402
from bubbleblower.contigs import contig_sequences, write_fasta  # noqa: E402
from bubbleblower.detect import detect_bubbles  # noqa: E402
from bubbleblower.gfa import load_gfa  # noqa: E402
from bubbleblower.modes import resolve_mode  # noqa: E402


def run() -> int:
    """Write per-mode metrics. Exit 2 when an input path is missing."""
    config = yaml.safe_load((Path(__file__).resolve().parent / "config.yaml").read_text(encoding="utf-8"))
    paths = {
        name: (ROOT / config[name]).resolve() if not str(config[name]).startswith("/") else Path(config[name])
        for name in ("reads_r1", "reads_r2", "graph", "baseline_contigs", "references", "minimap2")
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        print("missing inputs:", ", ".join(missing), file=sys.stderr)
        return 2
    graph, _ids = load_gfa(paths["graph"], graph_id="close_strains")
    coloured = colour_from_fastq(
        graph,
        [paths["reads_r1"], paths["reads_r2"]],
        k=int(config["k"]),
        min_depth=int(config["min_depth"]),
    )
    out = Path(__file__).resolve().parent / "data"
    out.mkdir(parents=True, exist_ok=True)
    summary = {
        "n_unitigs": len(coloured.cdbg.unitigs),
        "n_links": len(coloured.cdbg.links),
        "n_bubbles": len(detect_bubbles(coloured)),
        "n_colours": len(coloured.cdbg.colors or []),
        "modes": {},
    }
    for mode in ("retain", "colour_pop"):
        resolved = resolve_mode(coloured, mode)
        fasta = out / f"{mode}.fasta"
        write_fasta(contig_sequences(resolved), fasta)
        summary["modes"][mode] = quast_like(
            fasta,
            paths["references"],
            minimap2=str(paths["minimap2"]),
            min_align=int(config["min_align"]),
        )
        summary["modes"][mode]["n_bubbles"] = len(detect_bubbles(resolved))
    summary["metaspades"] = quast_like(
        paths["baseline_contigs"],
        paths["references"],
        minimap2=str(paths["minimap2"]),
        min_align=int(config["min_align"]),
    )
    (out / "metrics.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if summary["n_bubbles"] < 1:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
