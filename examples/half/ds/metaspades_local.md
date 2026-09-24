# Local metaSPAdes: keeping every bubble

These assemblies were already on disk. Bulge-on contigs are the metaSPAdes baseline. Bulge-off contigs keep the bubbles the remover would have deleted. minimap2 `asm5`, alignments of at least 200 bp.

| dataset | assembly | N50 | genome fraction | misassemblies | mismatches / 100 kbp | duplication |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| half_strains k55 | bulge on | 358 | 0.1457 | 18 | 1017.5 | 0.1462 |
| half_strains k55 | bulge off | 358 | 0.1457 | 19 | 1021.0 | 0.1462 |
| close strains k33 | bulge on | 289 | 0.1234 | 9 | 599.8 | 0.1236 |
| close strains k33 | bulge off | 287 | 0.1228 | 11 | 606.6 | 0.1230 |

Leaving every bubble in the metaSPAdes graph does not cut misassemblies. It adds one on the ten-strain graph and two on the three-strain graph, and mismatches rise slightly. This is the opposite of the MEGAHIT k21 result, where later iterations create long two-genome contigs. A mode that keeps every branch is the wrong edit for these metaSPAdes graphs. The useful edit is to split a contig where the read colour changes, not to retain the bubble.
