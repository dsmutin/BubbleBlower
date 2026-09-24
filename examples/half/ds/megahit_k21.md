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
