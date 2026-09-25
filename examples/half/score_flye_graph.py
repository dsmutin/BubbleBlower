#!/usr/bin/env python3
"""Run both graph debubblers on the half_strains metaFlye GFA.

The baseline is Flye's polished contig FASTA. Coverage comes from Flye's
``dp`` tag. Read colours are not attached.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from bubbleblower.assembly_metrics import quast_like  # noqa: E402
from bubbleblower.contigs import contig_sequences, write_fasta  # noqa: E402
from bubbleblower.debubblers import DEBUBBLERS, classify_debubble, resolve_debubbler  # noqa: E402
from bubbleblower.detect import detect_bubbles  # noqa: E402
from bubbleblower.gfa import load_gfa  # noqa: E402

from bench_paths import minimap2, work_dir  # noqa: E402

MINIMAP = minimap2()
WORK = ROOT / "examples" / "half" / "work"


def wins_against(mode: dict, baseline: dict) -> list[str]:
    """Main scores the mode wins. N50 is not a main score."""
    won = []
    if mode["misassemblies"] < baseline["misassemblies"]:
        won.append("misassemblies")
    if mode["mismatches_per_100kbp"] < baseline["mismatches_per_100kbp"]:
        won.append("mismatches_per_100kbp")
    if mode["duplication"] < baseline["duplication"]:
        won.append("duplication")
    if mode["genome_fraction"] > baseline["genome_fraction"]:
        won.append("genome_fraction")
    return won


def main() -> int:
    """Score both debubblers against the polished Flye contigs."""
    gfa = WORK / "flye" / "assembly_graph.gfa"
    assembly = WORK / "flye" / "assembly.fasta"
    ref = WORK / "megahit_k21" / "half_strains.references.fna"
    for path in (gfa, assembly, ref):
        if not path.is_file():
            raise SystemExit(f"missing {path}")
    graph, _ids = load_gfa(gfa, graph_id="half_strains_flye")
    bubbles = detect_bubbles(graph)
    print("unitigs", len(graph.cdbg.unitigs), "links", len(graph.cdbg.links), "bubbles", len(bubbles), flush=True)
    row: dict = {"n_bubbles": len(bubbles), "modes": {}}
    for name in DEBUBBLERS:
        counts: Counter[str] = Counter()
        for bubble in bubbles:
            decision = classify_debubble(graph, bubble, name)
            counts[f"{decision.label}:{decision.action}"] += 1
        print(name, dict(counts), flush=True)
        n_edits = sum(count for key, count in counts.items() if not key.endswith(":retain"))
        resolved = resolve_debubbler(graph, name, max_edits=n_edits)
        fasta = ROOT / "examples" / "half" / "data" / f"half_strains_flye_{name}.fasta"
        write_fasta(contig_sequences(resolved), fasta)
        metrics = quast_like(fasta, ref, minimap2=MINIMAP, min_align=200)
        row["modes"][name] = {"decisions": dict(counts), "metrics": metrics}
    baseline = quast_like(assembly, ref, minimap2=MINIMAP, min_align=200)
    row["metaflye"] = baseline
    for mode in row["modes"].values():
        mode["wins"] = wins_against(mode["metrics"], baseline)
    dest = ROOT / "examples" / "half" / "data" / "flye_debubble_metrics.json"
    dest.write_text(json.dumps({"half_strains": row}, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"half_strains": row}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
