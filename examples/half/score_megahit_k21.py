#!/usr/bin/env python3
"""Score MEGAHIT k21 contigs against the final MEGAHIT contigs.

References are the simulated genomes for that example. The k21 FASTA is the
unitig set after MEGAHIT's own k21 bubble removal. The final FASTA is the
baseline after every later k and every later bubble removal.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from bubbleblower.assembly_metrics import quast_like  # noqa: E402

from bench_paths import minimap2, work_dir  # noqa: E402

MINIMAP = minimap2()
DATASETS = ("half_strains", "low75", "low75half", "high100", "heldout_genera")


def write_references(dataset: str, dest: Path) -> None:
    """Concatenate simulated genomes. Skip files that are not FASTA."""
    sim = work_dir(dataset) / "sim"
    if not sim.is_dir():
        raise FileNotFoundError(sim)
    chunks: list[str] = []
    for path in sorted(sim.glob("*.fna")):
        text = path.read_text(encoding="utf-8", errors="replace")
        if not text.startswith(">"):
            continue
        chunks.append(text if text.endswith("\n") else text + "\n")
    if not chunks:
        raise FileNotFoundError(f"no genomes in {sim}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("".join(chunks), encoding="utf-8")


def score_dataset(dataset: str) -> dict:
    """Return metric dicts for k21 contigs and final contigs."""
    work = work_dir(dataset) / "megahit"
    k21 = work / "intermediate_contigs" / "k21.contigs.fa"
    final = work / "final.contigs.fa"
    if not k21.is_file() or not final.is_file():
        raise FileNotFoundError(dataset)
    ref = ROOT / "examples" / "half" / "work" / "megahit_k21" / f"{dataset}.references.fna"
    if not ref.is_file():
        write_references(dataset, ref)
    return {
        "megahit_k21": quast_like(k21, ref, minimap2=MINIMAP, min_align=200),
        "megahit_final": quast_like(final, ref, minimap2=MINIMAP, min_align=200),
    }


def main() -> int:
    """Score every dataset named on the command line, or half_strains."""
    names = sys.argv[1:] or ["half_strains"]
    out = ROOT / "examples" / "half" / "data"
    out.mkdir(parents=True, exist_ok=True)
    payload = {}
    for name in names:
        print("scoring", name, flush=True)
        payload[name] = score_dataset(name)
        print(json.dumps({name: payload[name]}, indent=2), flush=True)
    dest = out / "megahit_k21_metrics.json"
    existing: dict = {}
    if dest.is_file():
        existing = json.loads(dest.read_text(encoding="utf-8"))
    existing.update(payload)
    dest.write_text(json.dumps(existing, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
