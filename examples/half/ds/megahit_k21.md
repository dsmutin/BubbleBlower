# Initial MEGAHIT graphs

The final k141 FASTG in these examples has almost no edges. The bubble-rich
objects are the k=21 stage, before later k-mer iterations collapse them.

`intermediate_contigs/k21.bubble_seq.fa` is the set of sequences MEGAHIT
removed as bubbles. `k21.contigs.fa` converted with `megahit_toolkit
contig2fastg 21` is the unitig graph that remains after that removal.
Forward-only counts below use `bubbleblower.fastg.load_fastg`.

## Sequences MEGAHIT removed at k=21

| dataset | bubble sequences |
| --- | ---: |
| half_strains | 50045 |
| heldout_genera | 38843 |
| low75half | 21690 |
| low75 | 16875 |
| high100 | 10291 |

`k141.bubble_seq.fa` is empty in every example. The final FASTG is not the
training graph.

## Unitig graph after that removal

| dataset | forward unitigs | links | simple bubbles | kmer_divergence retain / pop / split |
| --- | ---: | ---: | ---: | --- |
| high100 | 149792 | 12767 | 100 | 98 / 2 / 0 |
| half_strains | 93070 | 41742 | 384 | 370 / 14 / 0 |
| low75 | 122002 | 13669 | 135 | 134 / 1 / 0 |
| low75half | 108039 | 32218 | 144 | 138 / 6 / 0 |
| heldout_genera | 140306 | 67156 | 525 | 512 / 12 / 1 |

`kmer_divergence` labels almost every remaining simple bubble as variation and keeps it. `heldout_genera` is the richest graph (525 bubbles). FASTG files are under `examples/half/work/megahit_k21/`. `half100half` has no MEGAHIT run.

Grouping bubble sequences by a shared (k-1)-prefix and (k-1)-suffix does
not recover one bubble per removed sequence: almost every sequence is its
own group. The FASTG, not that grouping, is the graph to debubble.

## k21 contigs versus MEGAHIT final contigs

Scored with minimap2 `asm5`, alignments of at least 200 bp, against `work/sim/*.fna`. Lower is better for misassemblies, mismatches per 100 kbp, and duplication. Higher is better for genome fraction. N50 is secondary and is expected to fall.

| dataset | assembly | N50 | genome fraction | misassemblies | mismatches / 100 kbp | duplication |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| half_strains | k21 contigs | 187 | 0.1312 | 9 | 116.8 | 0.1316 |
| half_strains | final contigs | 460 | 0.0990 | 24 | 202.8 | 0.0993 |
| low75 | k21 contigs | 126 | 0.0139 | 0 | 82.6 | 0.0140 |
| low75 | final contigs | 465 | 0.0066 | 1 | 139.2 | 0.0066 |
| low75half | k21 contigs | 126 | 0.0333 | 0 | 94.0 | 0.0334 |
| low75half | final contigs | 429 | 0.0180 | 2 | 232.7 | 0.0182 |
| heldout_genera | k21 contigs | 126 | 0.0456 | 11 | 148.0 | 0.0457 |
| heldout_genera | final contigs | 390 | 0.0276 | 45 | 318.9 | 0.0276 |
| high100 | k21 contigs | 126 | 0.0020 | 0 | 70.0 | 0.0020 |
| high100 | final contigs | 428 | 0.0006 | 2 | 103.6 | 0.0006 |

On all five finished sets the k21 contigs win misassemblies, mismatches, and genome fraction. Duplication is higher by about the same amount as the genome-fraction gain, so the extra alignment is new reference coverage rather than repeated mapping of the same bases. N50 is lower, as expected. That is three of the four main scores. `half100half` has no MEGAHIT output.

## Where the final assembly adds misassemblies

A misassembly here is one contig with two alignments of at least 200 bp on different reference records. On `half_strains` the final graph has 24 such contigs and k21 has 9. On `heldout_genera` the counts are 45 and 11. The extra final contigs are longer: the longest `half_strains` chimera is `k141_1194` at 1334 bp, and the longest `heldout_genera` chimera is `k141_299` at 2893 bp. The longest k21 chimeras in those two sets are 837 bp and 622 bp. Later k-mer iterations are joining sequences from two genomes into one contig. Keeping the k21 contigs avoids those joins. The two debubblers, run on the k21 FASTG, mostly retain the remaining simple bubbles (`heldout_genera`: 512 retain, 12 pop, 1 split), which is the same direction as keeping the k21 sequences.
