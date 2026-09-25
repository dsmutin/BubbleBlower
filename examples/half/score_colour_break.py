#!/usr/bin/env python3
"""Split contigs where the read colour changes and score them against metaSPAdes."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from bubbleblower.assembly_metrics import quast_like  # noqa: E402
from bubbleblower.colour_break import break_read_colour_chimeras  # noqa: E402

from bench_paths import minimap2, work_dir  # noqa: E402

MINIMAP = minimap2()
ISS = work_dir("half_strains") / "iss" / "initial"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Read FASTA records."""
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def score_pair(name: str, assembly: Path, reads: list[Path], reference: Path, dest: Path) -> dict:
    """Break colour chimeras in ``assembly`` and compare with the untouched file."""
    print("breaking", name, flush=True)
    broken = break_read_colour_chimeras(read_fasta(assembly), reads, k=21, min_run=40, min_piece=80)
    text = "".join(f">{ident}\n{sequence}\n" for ident, sequence in broken)
    dest.write_text(text, encoding="utf-8")
    print("scoring", name, flush=True)
    baseline = quast_like(assembly, reference, minimap2=MINIMAP, min_align=200)
    mode = quast_like(dest, reference, minimap2=MINIMAP, min_align=200)
    wins = []
    if mode["misassemblies"] < baseline["misassemblies"]:
        wins.append("misassemblies")
    if mode["mismatches_per_100kbp"] < baseline["mismatches_per_100kbp"]:
        wins.append("mismatches_per_100kbp")
    if mode["duplication"] < baseline["duplication"]:
        wins.append("duplication")
    if mode["genome_fraction"] > baseline["genome_fraction"]:
        wins.append("genome_fraction")
    return {"baseline": baseline, "read_colour_break": mode, "wins": wins, "n_pieces": len(broken)}


def main() -> int:
    """Score the close-strain walk and the full half_strains metaSPAdes contigs."""
    work = ROOT / "examples" / "half"
    out = work / "data"
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "close_k33": score_pair(
            "close_k33",
            out / "colour_topology.fasta",
            [work / "work" / "close" / "R1.fastq", work / "work" / "close" / "R2.fastq"],
            work / "work" / "close" / "references.fna",
            out / "close_read_colour_break.fasta",
        ),
        "half_strains_k55": score_pair(
            "half_strains_k55",
            work / "work" / "metaspades" / "contigs.fasta",
            [ISS / "sample_full_R1.fastq", ISS / "sample_full_R2.fastq"],
            work / "work" / "megahit_k21" / "half_strains.references.fna",
            out / "half_strains_read_colour_break.fasta",
        ),
    }
    dest = out / "colour_break_metrics.json"
    dest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
