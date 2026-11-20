# bubbleblower features

Check a box only after mandatory tests pass. This is the implementation contract.

## Graph source

- [x] Edit a MetaMetro ToCUMG (`cfa_to_cdbg` / `load_cdbg`). Do not assemble unitigs in BubbleBlower (VERSION 0.9.0)
- [x] A CGT is a view of that ToCUMG (`from_cgt`, `as_cgt`). It is not a second graph
- [x] GFA and MEGAHIT FASTG enter through MetaMetro `gfa_to_cfa` and `fastg_to_cfa`
- [x] Read colouring is MetaMetro `colour_by_read_accessions`; examples select ToCUMG namespaces

## Edits and launch criteria

- [x] Animate each resolution state on a pinned Fruchterman–Reingold layout (VERSION 0.10.2)
- [ ] On each iteration, either merge two nodes that sit in one bubble or next to each other, unioning their colours, or split one node
- [ ] A split partitions that node's colour: a binary 0/1 mask, or a gradual weight in `[0, 1]`, chosen by the method and by whether the colour channel is a mask or a score
- [ ] A launch criterion chooses merge versus split. Bayesian criteria that exist today are posteriors only: `coverage` (Poisson, `coverage-poisson-1`) and `multimodal` (those coverage terms plus colour equality, long-read linkage, and k-mer support). They adjust a pop candidate. They do not choose merge versus split
- [ ] An ML launch criterion on a CGT (NumPy GCN or PyG GCN from MetaMetro). Not implemented. `kmer_divergence` and `colour_topology` are fixed thresholds, not a trained model
- [ ] Train debubblers on MetaMetro benches and check them on initial MEGAHIT graphs (k21 FASTG plus bubble_seq.fa) from half_strains, low75, low75half, high100, heldout_genera. Uses the ML launch criterion above
- [ ] Cap debubbler edits at the first-pass decision count in the library. Today the cap lives only in `examples/half/score_*`; the library default is `max(8, 3 × bubbles)` (`debubblers.py`)

## Environment and CLI

- [x] Install MetaMetro 0.17.1 in the GitHub conda env (VERSION 0.10.3)
- [x] Score toy, error, and benchmark benches on the taxon layer (VERSION 0.10.4)
- [x] CLI `--search` defaults to `greedy`; `none` is a deprecated alias (VERSION 0.10.5)
- [x] Restore the MEGAHIT training and edit-cap items, and move environment items out of edits (VERSION 0.10.6)
- [x] Remove the `baseline` module; the no-graph CLI path is `pipeline.run_demo` (VERSION 0.10.7)
- [x] `state_scores.tsv` reports `n_forks`, detected `n_bubbles`, and `colour_score` separately (VERSION 0.11.0)
- [x] Example outputs under `examples/{toy,error,reads,benchmark}/data` are untracked and ignored; `examples/half/data` metrics JSONs stay until the #231 manifest (VERSION 0.11.1)
- [x] `ensure_bench` and `load_bench_graph` take `root=` or `$BUBBLEBLOWER_BENCH_ROOT` for the bench directory (#1165, VERSION 0.12.0)
- [x] Pin MetaMetro 0.20.1 (19a3edb) in environment.yml, CI, and the examples rule (#1169, VERSION 0.12.1)
- [x] Pin MetaMetro 0.20.2 (97fd453) (#1169, VERSION 0.12.2)
- [x] Pin MetaMetro 0.21.0 (4ed0762) (VERSION 0.13.1)
- [x] `examples/benchmark/run.py --namespace` (#1179, VERSION 0.14.0)
- [x] Examples rule loads `composition`, not gold `taxon` (#1182, VERSION 0.14.1)
- [x] Pin MetaMetro 0.22.1 (fdd0519) (VERSION 0.14.2)
- [x] Strand-aware `detect_bubbles`: oriented unitigs, one bubble per strand pair; phage_10 ids match gold 94/94 (#1190, VERSION 0.15.0)
- [x] Colourless baseline: `--namespace none` / `namespaces=()` drops every colour layer (#1172 stage 1, VERSION 0.15.1)
- [x] `.cursor/rules/` mirrors `agents/rules/` byte for byte; mandatory test `test_rules_mirror.py` (#1191, VERSION 0.15.2)
- [x] Pin MetaMetro 0.24.0 (f268a53); the #1190 parity test runs on the pin (#1206, VERSION 0.15.3)
- [x] Links without coverage (MEGAHIT builds) get the `split` share of source coverage in `load_bench_graph`; a fork sums to the source (#1202, VERSION 0.15.4)
- [x] `pop_branch` removes each branch link with its reverse-complement twin; no dangling link on double-stranded graphs (#1203, VERSION 0.15.5)
- [x] Strand fixture with source and sink on different strands; parity test checks ids and both ends against MetaMetro `strand_bubbles` (#1190 review, VERSION 0.15.6)
- [x] `tests/test_bench_input.py` is mandatory; the colourless test also checks link, internal-edge, and mapping colours (#1172 stage 1 review, VERSION 0.15.7)
- [x] Examples read the `composition` colouring, never gold `taxon`; `load_bench_graph` raises on a missing namespace; benchmark gates AMBER F1, AUROC, and delta_score, not false_pops (#1169, VERSION 0.13.0)

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
