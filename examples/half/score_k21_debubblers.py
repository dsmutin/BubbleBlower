#!/usr/bin/env python3
"""Run both graph debubblers on MEGAHIT k21 FASTG and score the unitigs.

The FASTG is ``contig2fastg 21`` of ``k21.contigs.fa``. MetaMetro compacts
it into a ToCUMG and keeps the assembly overlap. Unitig sequences are scored
against the final MEGAHIT contigs. Read colours come from the example's
Illumina FASTQ.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from bubbleblower.assembly_metrics import quast_like  # noqa: E402
from bubbleblower.colour_reads import genome_id  # noqa: E402
from bubbleblower.debubblers import DEBUBBLERS, classify_debubble, resolve_debubbler  # noqa: E402
from bubbleblower.detect import detect_bubbles  # noqa: E402
from bubbleblower.fastg import load_fastg  # noqa: E402

from bench_paths import minimap2, work_dir  # noqa: E402

MINIMAP = minimap2()
K = 21


def colour_bubble_unitigs(graph, fastq_paths: list[Path], *, k: int = K, min_depth: int = 2) -> int:
    """Paint only unitigs that sit in a simple bubble. Returns how many were painted."""
    wanted: set[str] = set()
    for bubble in detect_bubbles(graph):
        wanted.add(bubble.source)
        wanted.add(bubble.sink)
        for branch in bubble.branches:
            wanted.update(branch.path)
    sequences = {unitig.unitig_id: unitig.sequence for unitig in graph.cdbg.unitigs if unitig.unitig_id in wanted}
    index: dict[str, list[str]] = {}
    for unitig_id, sequence in sequences.items():
        if len(sequence) < k:
            continue
        for index_at in range(len(sequence) - k + 1):
            index.setdefault(sequence[index_at : index_at + k], []).append(unitig_id)
    hits: dict[str, Counter[str]] = {unitig_id: Counter() for unitig_id in sequences}
    for path in fastq_paths:
        lines = path.read_text(encoding="utf-8").splitlines()
        for offset in range(0, len(lines), 4):
            genome = genome_id(lines[offset][1:].split()[0])
            sequence = lines[offset + 1].strip().upper()
            seen: set[str] = set()
            if len(sequence) < k:
                continue
            for index_at in range(len(sequence) - k + 1):
                for unitig_id in index.get(sequence[index_at : index_at + k], ()):
                    if unitig_id in seen:
                        continue
                    seen.add(unitig_id)
                    hits[unitig_id][genome] += 1
    genomes = sorted({genome for counts in hits.values() for genome in counts})
    if not genomes:
        raise SystemExit("no read k-mers hit bubble unitigs")
    palette = {genome: colour_index for colour_index, genome in enumerate(genomes)}
    graph.cdbg.colors = [
        {"color_id": str(palette[genome]), "namespace": "genome", "value": genome} for genome in genomes
    ]
    by_id = {unitig.unitig_id: unitig for unitig in graph.cdbg.unitigs}
    painted = 0
    for unitig_id, counts in hits.items():
        kept = [palette[genome] for genome, depth in counts.items() if depth >= min_depth]
        if not kept:
            continue
        by_id[unitig_id].color_ids = sorted(kept)
        painted += 1
    for row in graph.cdbg.mapping:
        if row.unitig_id in by_id and by_id[row.unitig_id].color_ids:
            row.color_ids = list(by_id[row.unitig_id].color_ids)
    graph.validate()
    return painted


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
        final = work_dir(name) / "megahit" / "final.contigs.fa"
        reads = [
            work_dir(name) / "iss" / "initial" / "sample_full_R1.fastq",
            work_dir(name) / "iss" / "initial" / "sample_full_R2.fastq",
        ]
        ref = ROOT / "examples" / "half" / "work" / "megahit_k21" / f"{name}.references.fna"
        for path in (fastg, final, *reads, ref):
            if not path.is_file():
                raise SystemExit(f"missing {path}")
        print("loading", name, flush=True)
        graph = load_fastg(fastg, k=K, graph_id=name)
        matched = sum(1 for link in graph.cdbg.links if link.overlap is not None)
        print("overlaps", matched, "links", len(graph.cdbg.links), flush=True)
        print("colouring bubbles", name, flush=True)
        painted = colour_bubble_unitigs(graph, reads)
        print("painted", painted, flush=True)
        row: dict = {"forward_overlaps": matched, "painted_unitigs": painted, "modes": {}}
        for debubbler in DEBUBBLERS:
            print("classifying", name, debubbler, flush=True)
            counts = decision_counts(graph, debubbler)
            print("decisions", debubbler, counts, flush=True)
            print("resolving", name, debubbler, flush=True)
            n_edits = sum(count for key, count in counts.items() if not key.endswith(":retain"))
            resolved = resolve_debubbler(graph, debubbler, max_edits=n_edits)
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
