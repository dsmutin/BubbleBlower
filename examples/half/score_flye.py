#!/usr/bin/env python3
"""Score metaFlye contigs and a read-colour break of those contigs.

The reads are the badread 8× ONT simulation of the ten half_strains genomes.
The baseline is Flye's own contig FASTA. Missing files stop the run.
"""

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
WORK = ROOT / "examples" / "half" / "work"


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
    """Score the half_strains metaFlye run."""
    assembly = WORK / "flye" / "assembly.fasta"
    reads = WORK / "ont" / "reads.fastq"
    ref = WORK / "megahit_k21" / "half_strains.references.fna"
    for path in (assembly, reads, ref):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"missing {path}")
    records = read_fasta(assembly)
    print("breaking", len(records), "flye contigs", flush=True)
    broken = break_read_colour_chimeras(records, [reads], k=21, min_run=40, min_piece=80)
    fasta = ROOT / "examples" / "half" / "data" / "half_strains_flye_colour_break.fasta"
    fasta.write_text("".join(f">{ident}\n{seq}\n" for ident, seq in broken), encoding="utf-8")
    print("scoring", len(broken), "pieces", flush=True)
    baseline = quast_like(assembly, ref, minimap2=MINIMAP, min_align=200)
    mode = quast_like(fasta, ref, minimap2=MINIMAP, min_align=200)
    payload = {
        "half_strains": {
            "metaflye": baseline,
            "read_colour_break": mode,
            "n_contigs": len(records),
            "n_pieces": len(broken),
            "wins": wins_against(mode, baseline),
        }
    }
    dest = ROOT / "examples" / "half" / "data" / "flye_colour_break_metrics.json"
    dest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
