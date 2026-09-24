# bubbleblower features

Scaffold checklist. Check a box only after mandatory tests pass.

- [x] CLI resolves the strain bubble (`status=resolved`)
- [x] Detect, classify, reversible edits, greedy search (VERSION 0.1.0)
- [x] Beam search and MCMC
- [x] Read-level example via MetaMetro (`examples/reads/run.py`)
- [x] Contract tables: `bubble_results.tsv`, `edit_history.tsv`, `state_scores.tsv`
- [x] Colour-aware score; no strain false pops (VERSION 0.3.0)
- [x] AMBER F1 is the primary quality metric
- [x] Close-strain metaSPAdes graph coloured and scored (VERSION 0.4.0)
- [x] Graph-only debubblers kmer_divergence and colour_topology (VERSION 0.5.0)
- [x] Adjacent same-colour merge, inverse of a linear split (VERSION 0.6.0)
- [x] In-place adjacent merge so compaction can run on a large graph (VERSION 0.6.1)
- [x] Orientation-aware contig walk (VERSION 0.6.2)
- [x] Do not re-split the same bubble signature (VERSION 0.6.3)
- [ ] Beat metaSPAdes, MEGAHIT k21 bubble removal, and metaFlye on the majority of main scores
- [ ] Train debubblers on initial MEGAHIT graphs (k21 FASTG plus bubble_seq.fa) from half_strains, low75, low75half, high100, heldout_genera
