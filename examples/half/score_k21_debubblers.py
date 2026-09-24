#!/usr/bin/env python3
"""Run both graph debubblers on MEGAHIT k21 FASTG and score the unitigs.

The FASTG is ``contig2fastg 21`` of ``k21.contigs.fa``. Forward links whose
sequences share a (k-1) suffix/prefix get that overlap, so same-colour
compaction can join them. Reverse-strand links keep a missing overlap and
are not concatenated. Unitig sequences are scored against the final MEGAHIT
contigs. Read colours come from the example's Illumina FASTQ.
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
from bubbleblower.debubblers import DEBUBBLERS, classify_debubble, resolve_debubbler  # noqa: E402
from bubbleblower.detect import detect_bubbles  # noqa: E402
from bubbleblower.fastg import load_fastg  # noqa: E402

EX = Path("/mnt/tank/scratch/partition-metagenomics/smuteam/vaegbin_improved/examples")
MINIMAP = "/mnt/tank/scratch/dsmutin/partition-metagenomics/envs/vaegbin_env/bin/minimap2"
K = 21


def assign_forward_overlaps(graph, k: int = K) -> int:
    """Store a (k-1) overlap on forward links that actually share it."""
    overlap = k - 1
    sequences = {unitig.unitig_id: unitig.sequence for unitig in graph.cdbg.unitigs}
    matched = 0
    for link in graph.cdbg.links:
        if link.orientation != "++":
            continue
        left = sequences[link.source]
        right = sequences[link.target]
        if len(left) >= overlap and len(right) >= overlap and left[-overlap:] == right[:overlap]:
            link.overlap = overlap
            matched += 1
    return matched


def decision_counts(graph, name: str) -> dict[str, int]:
    """Count label:action pairs. The graph is not edited."""
    counts: Counter[str] = Counter()
    for bubble in detect_bubbles(graph):
        decision = classify_debubble(graph, bubble, name)
        counts[f"{decision.label}:{decision.action}"] += 1
    return dict(counts)


def write_unitigs(graph, path: Path) -> None:
    """Write one FASTA record per unitig. No extra joins are invented."""
    lines = []
    for unitig in graph.cdbg.unitigs:
        lines.append(f">{unitig.unitig_id}")
        lines.append(unitig.sequence)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


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
    """Score each dataset named on the command line."""
    names = sys.argv[1:] or ["high100"]
    out_dir = ROOT / "examples" / "half" / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    payload: dict = {}
    dest = out_dir / "megahit_k21_debubble_metrics.json"
    if dest.is_file():
        payload = json.loads(dest.read_text(encoding="utf-8"))
    for name in names:
        fastg = ROOT / "examples" / "half" / "work" / "megahit_k21" / f"{name}.k21.fastg"
        final = EX / name / "work" / "megahit" / "final.contigs.fa"
        reads = [
            EX / name / "work" / "iss" / "initial" / "sample_full_R1.fastq",
            EX / name / "work" / "iss" / "initial" / "sample_full_R2.fastq",
        ]
        ref = ROOT / "examples" / "half" / "work" / "megahit_k21" / f"{name}.references.fna"
        for path in (fastg, final, *reads, ref):
            if not path.is_file():
                raise SystemExit(f"missing {path}")
        print("loading", name, flush=True)
        graph = load_fastg(fastg, graph_id=name)
        matched = assign_forward_overlaps(graph)
        print("overlaps", matched, "links", len(graph.cdbg.links), flush=True)
        print("colouring", name, flush=True)
        graph = colour_from_fastq(graph, reads, k=K, min_depth=2)
        row: dict = {"forward_overlaps": matched, "modes": {}}
        for debubbler in DEBUBBLERS:
            print("classifying", name, debubbler, flush=True)
            counts = decision_counts(graph, debubbler)
            print("decisions", debubbler, counts, flush=True)
            print("resolving", name, debubbler, flush=True)
            resolved = resolve_debubbler(graph, debubbler)
            fasta = out_dir / f"{name}_{debubbler}.fasta"
            write_unitigs(resolved, fasta)
            print("scoring", name, debubbler, "unitigs", len(resolved.cdbg.unitigs), flush=True)
            metrics = quast_like(fasta, ref, minimap2=MINIMAP, min_align=200)
            row["modes"][debubbler] = {"decisions": counts, "metrics": metrics, "n_unitigs": len(resolved.cdbg.unitigs)}
        print("scoring final", name, flush=True)
        baseline = quast_like(final, ref, minimap2=MINIMAP, min_align=200)
        row["megahit_final"] = baseline
        for debubbler, mode in row["modes"].items():
            mode["wins"] = wins_against(mode["metrics"], baseline)
        payload[name] = row
        dest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({name: row}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
