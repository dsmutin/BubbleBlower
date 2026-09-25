#!/usr/bin/env python3
"""Score local metaSPAdes contigs that are already on disk.

Bulge-on contigs are the baseline. Bulge-off contigs are the unreduced graph.
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
WORK = ROOT / "examples" / "half" / "work"


def main() -> int:
    """Write metrics for the two local metaSPAdes runs."""
    pairs = {
        "half_strains_k55": (
            WORK / "metaspades" / "contigs.fasta",
            WORK / "metaspades_nobulge" / "contigs.fasta",
            WORK / "megahit_k21" / "half_strains.references.fna",
        ),
        "close_k33": (
            WORK / "close" / "spades_pop" / "contigs.fasta",
            WORK / "close" / "spades_keep" / "contigs.fasta",
            WORK / "close" / "references.fna",
        ),
    }
    payload = {}
    for name, (pop, keep, ref) in pairs.items():
        for path in (pop, keep, ref):
            if not path.is_file():
                raise SystemExit(f"missing {path}")
        print("scoring", name, flush=True)
        payload[name] = {
            "metaspades_pop": quast_like(pop, ref, minimap2=MINIMAP, min_align=200),
            "metaspades_keep": quast_like(keep, ref, minimap2=MINIMAP, min_align=200),
        }
        print(json.dumps({name: payload[name]}, indent=2), flush=True)
    dest = ROOT / "examples" / "half" / "data" / "local_spades_metrics.json"
    dest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
