# Bubble census

Simple-bubble counts are for the full ISS read set of each example, assembled with metaSPAdes 4.2.0 (`--meta -k 21,33 -t 8 -m 8`) on SLURM partition `main`. Bulge-off jobs pass `--configs-dir examples/half/work/spades_configs`. Counts are read from `assembly_graph_after_simplification.gfa` with `bubbleblower.gfa.load_gfa` and `bubbleblower.detect.detect_bubbles`. No 3-genome subset is used.

The first submission (`892928`–`892935`, 32G each) was cancelled: together with the running 510G `cami3_illumina` job it exceeded the per-user memory limit. The live jobs are `892936`–`892943` at 8G. They are `PENDING` with reason `QOSMaxMemoryPerUser` until that large job releases memory. Link and bubble columns stay empty until those jobs finish. Nothing below is a stand-in graph.

| dataset | reads | genomes | links bulge-on | links bulge-off | n_bubbles | SLURM pop | SLURM keep |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| low75 | 100000 | 75 | pending | pending | pending | 892943 | 892942 |
| low75half | 100000 | 38 | pending | pending | pending | 892941 | 892940 |
| high100 | 100000 | 100 | pending | pending | pending | 892939 | 892938 |
| half100half | absent | absent | | | | | |
| half_strains | 100000 | 10 | pending | pending | pending | 892937 | 892936 |

`half100half` has no `work/sim` genomes and no FASTQ. It was not assembled.

`high100` has no `sample_full_R*.fastq`. The pair that is present, and size-stable at 100000 reads per file, is `examples/high100/work/iss/initial/.iss_full/pool_full_R1.fastq` and `pool_full_R2.fastq`. Those files are the inputs for jobs 892939 and 892938. The InSilicoSeq Snakemake process was still running when the jobs were submitted.

Finished graphs, once the jobs complete, will be:

- `examples/half/work/bench/low75_pop/assembly_graph_after_simplification.gfa`
- `examples/half/work/bench/low75_keep/assembly_graph_after_simplification.gfa`
- `examples/half/work/bench/low75half_pop/assembly_graph_after_simplification.gfa`
- `examples/half/work/bench/low75half_keep/assembly_graph_after_simplification.gfa`
- `examples/half/work/bench/high100_pop/assembly_graph_after_simplification.gfa`
- `examples/half/work/bench/high100_keep/assembly_graph_after_simplification.gfa`
- `examples/half/work/bench/half_strains_pop/assembly_graph_after_simplification.gfa`
- `examples/half/work/bench/half_strains_keep/assembly_graph_after_simplification.gfa`

Which of these full read sets is bubble-rich is not known yet. The keep graphs do not exist, so no dataset can be chosen for a full downstream run from this census.
