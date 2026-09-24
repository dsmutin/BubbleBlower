#!/usr/bin/env python3
"""Run both graph debubblers on local metaSPAdes keep-graphs and score them.

The baseline is metaSPAdes with bulge removal. The input graph is the
matching run with bulge removal disabled. Contigs follow link orientation.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from bubbleblower.assembly_metrics import quast_like  # noqa: E402
from bubbleblower.colour_reads import colour_from_fastq  # noqa: E402
from bubbleblower.contigs import contig_sequences, write_fasta  # noqa: E402
from bubbleblower.debubblers import DEBUBBLERS, classify_debubble, resolve_debubbler  # noqa: E402
from bubbleblower.detect import detect_bubbles  # noqa: E402
from bubbleblower.gfa import load_gfa  # noqa: E402

MINIMAP = "/mnt/tank/scratch/dsmutin/partition-metagenomics/envs/vaegbin_env/bin/minimap2"
WORK = ROOT / "examples" / "half" / "work"
EX = Path("/mnt/tank/scratch/partition-metagenomics/smuteam/vaegbin_improved/examples")


def decision_counts(graph, name: str) -> dict[str, int]:
    """Count label:action pairs. The graph is not edited."""
    counts: Counter[str] = Counter()
    for bubble in detect_bubbles(graph):
        decision = classify_debubble(graph, bubble, name)
        counts[f"{decision.label}:{decision.action}"] += 1
    return dict(counts)


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


def score_one(name: str, graph_path: Path, reads: list[Path], pop: Path, ref: Path, out_dir: Path) -> dict:
    """Colour one keep-graph, resolve both debubblers, and score the contigs."""
    for path in (graph_path, *reads, pop, ref):
        if not path.is_file():
            raise SystemExit(f"missing {path}")
    print("loading", name, flush=True)
    graph, _ids = load_gfa(graph_path, graph_id=name)
    print("colouring", name, "unitigs", len(graph.cdbg.unitigs), flush=True)
    graph = colour_from_fastq(graph, reads, k=21, min_depth=2)
    row: dict = {"n_bubbles": len(detect_bubbles(graph)), "modes": {}}
    for debubbler in DEBUBBLERS:
        print("classifying", name, debubbler, flush=True)
        counts = decision_counts(graph, debubbler)
        print("decisions", debubbler, counts, flush=True)
        resolved = resolve_debubbler(graph, debubbler)
        fasta = out_dir / f"{name}_{debubbler}.fasta"
        write_fasta(contig_sequences(resolved), fasta)
        print("scoring", name, debubbler, flush=True)
        metrics = quast_like(fasta, ref, minimap2=MINIMAP, min_align=200)
        row["modes"][debubbler] = {"decisions": counts, "metrics": metrics}
    print("scoring pop", name, flush=True)
    baseline = quast_like(pop, ref, minimap2=MINIMAP, min_align=200)
    row["metaspades_pop"] = baseline
    for mode in row["modes"].values():
        mode["wins"] = wins_against(mode["metrics"], baseline)
    return row


def main() -> int:
    """Score the local keep-graphs named on the command line."""
    available = {
        "close_k33": (
            WORK / "close" / "spades_keep" / "assembly_graph_after_simplification.gfa",
            [WORK / "close" / "R1.fastq", WORK / "close" / "R2.fastq"],
            WORK / "close" / "spades_pop" / "contigs.fasta",
            WORK / "close" / "references.fna",
        ),
        "half_strains_k55": (
            WORK / "metaspades_nobulge" / "assembly_graph_after_simplification.gfa",
            [
                EX / "half_strains" / "work" / "iss" / "initial" / "sample_full_R1.fastq",
                EX / "half_strains" / "work" / "iss" / "initial" / "sample_full_R2.fastq",
            ],
            WORK / "metaspades" / "contigs.fasta",
            WORK / "megahit_k21" / "half_strains.references.fna",
        ),
    }
    names = sys.argv[1:] or list(available)
    out_dir = ROOT / "examples" / "half" / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / "spades_debubble_metrics.json"
    payload: dict = {}
    if dest.is_file():
        payload = json.loads(dest.read_text(encoding="utf-8"))
    for name in names:
        if name not in available:
            raise SystemExit(f"unknown dataset {name}")
        payload[name] = score_one(name, *available[name], out_dir)
        dest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({name: payload[name]}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
