# Local metaSPAdes: keeping every bubble

These assemblies were already on disk. Bulge-on contigs are the metaSPAdes baseline. Bulge-off contigs keep the bubbles the remover would have deleted. minimap2 `asm5`, alignments of at least 200 bp.

| dataset | assembly | N50 | genome fraction | misassemblies | mismatches / 100 kbp | duplication |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| half_strains k55 | bulge on | 358 | 0.1457 | 18 | 1017.5 | 0.1462 |
| half_strains k55 | bulge off | 358 | 0.1457 | 19 | 1021.0 | 0.1462 |
| close strains k33 | bulge on | 289 | 0.1234 | 9 | 599.8 | 0.1236 |
| close strains k33 | bulge off | 287 | 0.1228 | 11 | 606.6 | 0.1230 |

Leaving every bubble in the metaSPAdes graph does not cut misassemblies. It adds one on the ten-strain graph and two on the three-strain graph, and mismatches rise slightly. This is the opposite of the MEGAHIT k21 result, where later iterations create long two-genome contigs. A mode that keeps every branch is the wrong edit for these metaSPAdes graphs.

## Read-colour break

`read_colour_break` splits a contig where the accession painted on its k-mers changes. N50 falls. Three of the four main scores improve against the metaSPAdes bulge-on contigs.

| dataset | assembly | N50 | genome fraction | misassemblies | mismatches / 100 kbp | duplication |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| half_strains k55 | bulge on | 358 | 0.1457 | 18 | 1017.5 | 0.1462 |
| half_strains k55 | colour break | 337 | 0.1375 | 6 | 977.3 | 0.1380 |
| close strains k33 | bulge on | 289 | 0.1234 | 9 | 599.8 | 0.1236 |
| close strains k33 | colour break of the unreduced walk | 244 | 0.1028 | 4 | 158.4 | 0.1029 |

On `half_strains` the break is applied to the bulge-on contigs. On the close strains the break is applied to the walk of the unreduced graph, whose own mismatch rate was already 212.8 rather than 599.8. In both cases misassemblies, mismatches, and duplication improve, and genome fraction drops. Numbers are in `examples/half/data/colour_break_metrics.json`.

## Graph debubblers on the close-strain keep graph

Contigs are walked with link orientation and scored against bulge-on contigs.

| dataset | bubbles | kmer_divergence | colour_topology | misassemblies | mismatches / 100 kbp |
| --- | ---: | --- | --- | ---: | ---: |
| close strains k33 | 27 | 26 retain, 1 pop | 27 retain | 10 (bulge-on 9) | 213 (bulge-on 600) |
| half_strains k55 | 10 | 8 retain, 2 pop | 10 retain | 19 (bulge-on 18) | 131 (bulge-on 1017) |
| low75 k21,33 | 27 | 25 retain, 1 pop, 1 split | 27 retain | 0 (bulge-on 2) | 91 (bulge-on 502) |

On the close-strain and ten-strain graphs both modes win mismatches and duplication and lose misassemblies and genome fraction. That is two of the four main scores. On `low75` both win misassemblies (0 versus 2), mismatches (91 versus 502), and duplication. Genome fraction is slightly lower (0.0145 versus 0.0146). That is three of the four main scores. N50 stays 245 versus 246. Numbers are in `examples/half/data/spades_debubble_metrics.json`.
