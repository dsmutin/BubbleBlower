#!/usr/bin/env python3
"""Apply read-colour breaks to MEGAHIT k21 contigs and score them.

The baseline is the final MEGAHIT contig set, after later k-mer iterations
and MEGAHIT's own later bubble removal. The mode input is
``intermediate_contigs/k21.contigs.fa``, not the final FASTA.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from bubbleblower.assembly_metrics import quast_like  # noqa: E402
from bubbleblower.colour_break import break_read_colour_chimeras  # noqa: E402

EX = Path("/mnt/tank/scratch/partition-metagenomics/smuteam/vaegbin_improved/examples")
MINIMAP = "/mnt/tank/scratch/dsmutin/partition-metagenomics/envs/vaegbin_env/bin/minimap2"


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


def main() -> int:
    """Score each dataset named on the command line."""
    names = sys.argv[1:] or ["half_strains"]
    out_dir = ROOT / "examples" / "half" / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    payload: dict = {}
    dest = out_dir / "megahit_k21_colour_break_metrics.json"
    if dest.is_file():
        payload = json.loads(dest.read_text(encoding="utf-8"))
    for name in names:
        k21 = EX / name / "work" / "megahit" / "intermediate_contigs" / "k21.contigs.fa"
        final = EX / name / "work" / "megahit" / "final.contigs.fa"
        reads = [
            EX / name / "work" / "iss" / "initial" / "sample_full_R1.fastq",
            EX / name / "work" / "iss" / "initial" / "sample_full_R2.fastq",
        ]
        ref = ROOT / "examples" / "half" / "work" / "megahit_k21" / f"{name}.references.fna"
        for path in (k21, final, *reads, ref):
            if not path.is_file():
                raise SystemExit(f"missing {path}")
        before = read_fasta(k21)
        print("breaking", name, len(before), flush=True)
        broken = break_read_colour_chimeras(before, reads, k=21, min_run=40, min_piece=80)
        fasta = out_dir / f"{name}_k21_colour_break.fasta"
        fasta.write_text("".join(f">{ident}\n{seq}\n" for ident, seq in broken), encoding="utf-8")
        print("scoring", name, "pieces", len(broken), flush=True)
        baseline = quast_like(final, ref, minimap2=MINIMAP, min_align=200)
        mode = quast_like(fasta, ref, minimap2=MINIMAP, min_align=200)
        wins = []
        if mode["misassemblies"] < baseline["misassemblies"]:
            wins.append("misassemblies")
        if mode["mismatches_per_100kbp"] < baseline["mismatches_per_100kbp"]:
            wins.append("mismatches_per_100kbp")
        if mode["duplication"] < baseline["duplication"]:
            wins.append("duplication")
        if mode["genome_fraction"] > baseline["genome_fraction"]:
            wins.append("genome_fraction")
        payload[name] = {
            "megahit_final": baseline,
            "k21_read_colour_break": mode,
            "n_k21_contigs": len(before),
            "n_pieces": len(broken),
            "wins": wins,
        }
        dest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({name: payload[name]}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
