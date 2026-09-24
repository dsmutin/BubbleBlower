#!/usr/bin/env python3
"""Read-level example using MetaMetro's mock metagenome contract.

Reads and a de Bruijn CFA are written with MetaMetro. BubbleBlower sees a
four-node allele graph of the same two strains, coloured by those reads.
A compacted de Bruijn graph of short homopolymer-rich sequences is not a
simple bubble; that topology is stored but not used as the resolver input.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
_METAMETRO = ROOT.parent / "metametro" / "src"
if _METAMETRO.is_dir():
    sys.path.insert(0, str(_METAMETRO))

from metametro.contracts.assembly import dbg_from_sequences, read_fastq, simulate_metagenome  # noqa: E402
from metametro.contracts.colouring import colour_by_reads  # noqa: E402
from metametro.formats.cfa.io import dump_cfa  # noqa: E402

from bubbleblower.detect import detect_bubbles  # noqa: E402
from bubbleblower.graph import build_graph  # noqa: E402


def _hits(reads: list[tuple[str, str]], motif: str) -> list[int]:
    colours: set[int] = set()
    for sample, sequence in reads:
        if motif in sequence:
            colours.add(0 if sample == "sample_A" else 1)
    return sorted(colours)


def run() -> int:
    """Simulate two strains, write MetaMetro reads/CFA, colour a simple bubble."""
    out = Path(__file__).resolve().parent / "data"
    left = "ACGTAGCTTG"
    allele_a = "CATGCA"
    allele_b = "GATCGA"
    right = "TGCCTAAGGC"
    genomes = {
        "strain_A": left + allele_a + right,
        "strain_B": left + allele_b + right,
    }
    genome_dir = out / "genomes"
    genome_dir.mkdir(parents=True, exist_ok=True)
    for name, sequence in genomes.items():
        (genome_dir / f"{name}.fna").write_text(f">{name}\n{sequence}\n", encoding="utf-8")
    paths = simulate_metagenome(
        genomes,
        [("sample_A", "strain_A", 20), ("sample_B", "strain_B", 20)],
        out / "reads",
        read_length=12,
        seed=1,
    )
    fastq = read_fastq(paths["reads"])
    sample_reads = []
    coloured_reads = []
    for read_id, sequence in fastq:
        sample = "sample_A" if read_id.startswith("sample_A") else "sample_B"
        sample_reads.append((sample, sequence))
        coloured_reads.append((read_id, sample, sequence))
    cfa = dbg_from_sequences(fastq, k=5, graph_id="read_level")
    cfa = colour_by_reads(cfa, coloured_reads, ["sample_A", "sample_B"], min_edge_kmer_density=1)
    dump_cfa(cfa, out / "cfa")
    graph = build_graph(
        graph_id="read_level_alleles",
        colors=[
            {"color_id": "0", "namespace": "sample", "value": "sample_A"},
            {"color_id": "1", "namespace": "sample", "value": "sample_B"},
        ],
        nodes=[
            {"id": "S", "sequence": left, "colors": _hits(sample_reads, left[:8]), "coverage": 40.0},
            {"id": "A", "sequence": allele_a, "colors": _hits(sample_reads, allele_a), "coverage": 20.0},
            {"id": "B", "sequence": allele_b, "colors": _hits(sample_reads, allele_b), "coverage": 20.0},
            {"id": "T", "sequence": right, "colors": _hits(sample_reads, right[:8]), "coverage": 40.0},
        ],
        links=[
            {"id": "eSA", "source": "S", "target": "A", "colors": [0], "coverage": 20.0},
            {"id": "eAT", "source": "A", "target": "T", "colors": [0], "coverage": 20.0},
            {"id": "eSB", "source": "S", "target": "B", "colors": [1], "coverage": 20.0},
            {"id": "eBT", "source": "B", "target": "T", "colors": [1], "coverage": 20.0},
        ],
    )
    graph.coverage_source = "read_counts"
    bubbles = detect_bubbles(graph)
    payload = {
        "n_reads": len(fastq),
        "n_cfa_nodes": len(cfa.nodes),
        "n_allele_unitigs": len(graph.cdbg.unitigs),
        "n_bubbles": len(bubbles),
        "n_colours": len(cfa.colors or []),
        "coverage_source": graph.coverage_source,
    }
    (out / "summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if payload["n_reads"] != 40 or payload["n_bubbles"] != 1 or payload["n_colours"] != 2:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
