"""Read-level example using MetaMetro's deterministic read and de Bruijn builders."""

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
from metametro.converters.cfa_to_cdbg import cfa_to_cdbg  # noqa: E402
from metametro.formats.cfa.io import dump_cfa  # noqa: E402

from bubbleblower.detect import detect_bubbles  # noqa: E402
from bubbleblower.graph import from_cdbg  # noqa: E402


def run() -> int:
    """Simulate two strains, colour the de Bruijn graph, and detect bubbles."""
    out = Path(__file__).resolve().parent / "data"
    shared = "ACGTACGTACGTACGT"
    genomes = {
        "strain_A": "AAAACCCC" + shared + "GGGGTTTT",
        "strain_B": "TTTTGGGG" + shared + "CCCCAAAA",
    }
    genome_dir = out / "genomes"
    genome_dir.mkdir(parents=True, exist_ok=True)
    for name, sequence in genomes.items():
        (genome_dir / f"{name}.fna").write_text(f">{name}\n{sequence}\n", encoding="utf-8")
    paths = simulate_metagenome(
        genomes,
        [("sample_A", "strain_A", 6), ("sample_B", "strain_B", 6)],
        out / "reads",
        read_length=16,
        seed=1,
    )
    fastq = read_fastq(paths["reads"])
    coloured_reads = []
    for read_id, sequence in fastq:
        sample = "sample_A" if read_id.startswith("sample_A") else "sample_B"
        coloured_reads.append((read_id, sample, sequence))
    cfa = dbg_from_sequences(fastq, k=5, graph_id="read_level")
    cfa = colour_by_reads(cfa, coloured_reads, ["sample_A", "sample_B"])
    dump_cfa(cfa, out / "cfa")
    cdbg = cfa_to_cdbg(cfa)
    coverage = {}
    for unitig in cdbg.unitigs:
        coverage[unitig.unitig_id] = float(sum(sequence.count(unitig.sequence) for _read_id, sequence in fastq))
    graph = from_cdbg(cdbg, node_coverage=coverage)
    bubbles = detect_bubbles(graph)
    payload = {
        "n_reads": len(fastq),
        "n_unitigs": len(cdbg.unitigs),
        "n_bubbles": len(bubbles),
        "n_colours": len(cfa.colors or []),
    }
    (out / "summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if payload["n_reads"] != 12 or payload["n_unitigs"] < 1 or payload["n_colours"] != 2:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
