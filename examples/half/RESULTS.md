# Close-strain graph

Three half_strains genomes were assembled with metaSPAdes 4.2.0 (`-k 21,33`):
Escherichia coli `GCF_000952955`, Shigella sonnei `GCF_001518855`, and
Escherichia albertii `GCF_001549955`. Reads are the matching InSilicoSeq
records from `examples/half_strains` (27,491 pairs).

The ten-strain k=55 graph has 4 simple bubbles (140 links on 20,839 unitigs).
Turning the bulge remover off there raised that to 10. The three-strain graph
is the one with a real bubble set: 27 bubbles, 341 links, 7,866 unitigs when
bulge removal is disabled, against 68 links when metaSPAdes pops bulges.

## What the bubbles are

Branch sequences were aligned to the three genomes (`minimap2 -k 11`).
16 bubbles have branches that hit different genomes. 9 hit the same genome.
2 are parallel links with empty internal paths. Coverage ratio and colour
equality do not separate the two classes: both branches usually carry the
same read colour, and ratios run from about 0.4 to 1.0 in both classes.

## Metrics

Contigs are minimap2 `asm5` alignments, minimum 200 bp. Higher is better for
N50 and genome fraction. Lower is better for misassemblies, mismatches per
100 kbp, and duplication.

| assembly | N50 | genome fraction | misassemblies | mismatches / 100 kbp | duplication |
| --- | ---: | ---: | ---: | ---: | ---: |
| metaSPAdes contigs (bulges popped) | 289 | 0.1234 | 9 | 599.8 | 0.1236 |
| BubbleBlower `retain` (walk the unreduced graph) | 286 | 0.1226 | 9 | 212.8 | 0.1228 |

`retain` cuts mismatches by about 3x and does not add misassemblies. It does
not beat metaSPAdes on N50 or genome fraction. `colour_pop` does not remove
these bubbles, because the branches share a colour and the coverage ratio is
above 0.15.

metaFlye is not in this table. `half_strains` has Illumina reads only, so
there is no ONT graph to resolve.

Numbers are recomputed by `examples/half/run.py`.
