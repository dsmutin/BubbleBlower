# bubbleblower features

Check a box only after mandatory tests pass. This is the implementation contract.

## Graph source

- [x] Edit a MetaMetro ToCUMG (`cfa_to_cdbg` / `load_cdbg`). Do not assemble unitigs in BubbleBlower (VERSION 0.9.0)
- [x] A CGT is a view of that ToCUMG (`from_cgt`, `as_cgt`). It is not a second graph
- [x] GFA and MEGAHIT FASTG enter through MetaMetro `gfa_to_cfa` and `fastg_to_cfa`
- [x] Read colouring is MetaMetro `colour_by_read_accessions`; examples select ToCUMG namespaces

## Edits and launch criteria

- [x] Animate each resolution state on a pinned Fruchterman–Reingold layout (VERSION 0.10.0)
- [ ] On each iteration, either merge two nodes that sit in one bubble or next to each other, unioning their colours, or split one node
- [ ] A split partitions that node's colour: a binary 0/1 mask, or a gradual weight in `[0, 1]`, chosen by the method and by whether the colour channel is a mask or a score
- [ ] A launch criterion chooses merge versus split. Bayesian criteria that exist today are posteriors only: `coverage` (Poisson, `coverage-poisson-1`) and `multimodal` (those coverage terms plus colour equality, long-read linkage, and k-mer support). They adjust a pop candidate. They do not choose merge versus split
- [ ] An ML launch criterion on a CGT (NumPy GCN or PyG GCN from MetaMetro). Not implemented. `kmer_divergence` and `colour_topology` are fixed thresholds, not a trained model

## Benchmarks

- [ ] Debubbling with those two edits is compared with the assembler's own bubble removal
- [ ] The comparison reports AMBER F1, contig F1, misassemblies, mismatches, N50, and runtime
- [ ] That comparison is reproducible from pinned inputs, a seed, and recorded commands
- [ ] Benchmark colours use the taxid assigned to each read before the metagenome was generated. See CONTRIBUTING.md

## Already checked, and still true

- [x] CLI resolves the strain bubble (`status=resolved`)
- [x] Detect, classify, reversible edits, greedy search
- [x] Beam search and MCMC
- [x] Read-level example adopts the MetaMetro ToCUMG (`examples/reads/run.py`). That graph has no simple bubble
- [x] Contract tables: `bubble_results.tsv`, `edit_history.tsv`, `state_scores.tsv`
- [x] Colour-aware score
- [x] AMBER F1 on the synthetic 50-bubble graph (seed 42). Those colours are fixture labels, not pre-generation taxids
- [x] Graph-only debubblers `kmer_divergence` and `colour_topology`
- [x] Adjacent same-colour merge (equal colours; not a union of different colours)
- [x] In-place adjacent merge
- [x] Orientation-aware contig walk
- [ ] Beat metaSPAdes, MEGAHIT k21 bubble removal, and metaFlye on the majority of main scores

## Not in the code

These are absent. Do not treat the lines above as done.

- The iteration is not "merge two nodes and union their colours" or "split one node and partition its colour". `merge_adjacent` requires equal colours. `merge_instances` unions colours only for copies of one sequence. `split_instance` copies the full colour set onto every piece.
- There is no binary 0/1 colour mask and no gradual colour weight in `[0, 1]`.
- Nothing chooses merge versus split. `coverage` and `multimodal` are Bayesian posteriors that reweight a pop. `kmer_divergence` and `colour_topology` are fixed thresholds.
- No trained launch criterion. The MetaMetro NumPy GCN and the PyTorch Geometric GCN are not called.
- No comparison of those two edits with assembler bubble removal on AMBER F1, contig F1, misassemblies, mismatches, N50, and runtime, and no pinned reproducible run of that comparison.
- The 50-bubble AMBER F1 run uses CFA fixture labels. It does not use the taxid assigned to each read before metagenome generation.
- Accepted edits still allocate unitig records on the ToCUMG. They do not go through MetaMetro `EditProposal`.
- `examples/reads/run.py` adopts the MetaMetro ToCUMG. That graph has no simple bubble.
