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

## Colour break on the final contigs

`read_colour_break` applied to the final FASTA, using the example's own Illumina reads, cuts the long two-genome contigs without going back to k21. On every finished set it wins misassemblies, mismatches, and duplication. Genome fraction falls. N50 falls a little. `half100half` has no reads and no final contigs.

| dataset | assembly | N50 | genome fraction | misassemblies | mismatches / 100 kbp | duplication |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| half_strains | final | 460 | 0.0990 | 24 | 202.8 | 0.0993 |
| half_strains | colour break | 436 | 0.0924 | 9 | 167.7 | 0.0927 |
| heldout_genera | final | 390 | 0.0276 | 45 | 318.9 | 0.0276 |
| heldout_genera | colour break | 343 | 0.0220 | 12 | 193.2 | 0.0220 |
| low75 | final | 465 | 0.0066 | 1 | 139.2 | 0.0066 |
| low75 | colour break | 462 | 0.0066 | 0 | 137.4 | 0.0066 |
| low75half | final | 429 | 0.0180 | 2 | 232.7 | 0.0182 |
| low75half | colour break | 426 | 0.0180 | 1 | 230.8 | 0.0181 |
| high100 | final | 428 | 0.0006 | 2 | 103.6 | 0.0006 |
| high100 | colour break | 421 | 0.0006 | 0 | 100.0 | 0.0006 |

## Colour break on the k21 contigs

`read_colour_break` on `k21.contigs.fa`, with the same Illumina reads, is a mode on the initial MEGAHIT sequences rather than on the final FASTA. Pieces shorter than 80 bp are dropped. Against the final contigs it wins misassemblies, mismatches, and genome fraction on every finished set. Duplication rises with that extra coverage. N50 falls to the k21 length.

| dataset | assembly | N50 | genome fraction | misassemblies | mismatches / 100 kbp | duplication |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| half_strains | k21 colour break | 176 | 0.1243 | 3 | 100.4 | 0.1248 |
| heldout_genera | k21 colour break | 126 | 0.0405 | 2 | 113.1 | 0.0405 |
| low75 | k21 colour break | 126 | 0.0139 | 0 | 82.5 | 0.0139 |
| low75half | k21 colour break | 126 | 0.0332 | 0 | 94.1 | 0.0333 |
| high100 | k21 colour break | 126 | 0.0020 | 0 | 69.6 | 0.0020 |

The k21 contigs already contained two-genome joins. On `half_strains` the break cuts those from 9 to 3, and mismatches from 116.8 to 100.4, while genome fraction stays above the final assembly (0.1243 versus 0.0990; uncut k21 was 0.1312). On `heldout_genera` the same cut is 11 to 2 misassemblies and 148.0 to 113.1 mismatches per 100 kbp, with genome fraction 0.0405 versus 0.0276 for the final contigs (uncut k21 was 0.0456). `low75`, `low75half`, and `high100` already had no k21 misassembly, and the break does not add one. Their scores stay on the k21 side of the comparison.
