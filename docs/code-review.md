# Code review: BubbleBlower 0.2.0

Read-only review against the BubbleBlower specification (contracts BB-1 to
BB-13, test matrix in section 34, evaluation metrics in section 35, and the
acceptance criterion in section 36). Commits reviewed: `ed9d243`, `eaad827`.

No `method-decision.md` exists in the repository. Decisions were taken from
the specification alone.

## Baseline

Stop criterion from the task: all tests pass; code is organised and
reproducible; every example runs and returns the expected output.

Measured on 2026-09-24 with `PYTHONPATH=src:../metametro/src`:

| Check | Result |
|-------|--------|
| `pytest -m mandatory` | 20 passed |
| `pytest` (all) | 21 passed |
| `examples/toy/run.py` | `decision=strain`, exit 0 |
| `examples/error/run.py` | branch `E` popped, `A` kept, score −7.08 → −4.90 |
| `examples/benchmark/run.py` | F1(error) 0.958, AUROC 0.995, exit 0 |
| `examples/reads/run.py` | 12 reads, 11 unitigs, 2 colours, 0 bubbles |
| `heldout_genera` CDBG | loads, validates, 0 bubbles (24 links on 6439 unitigs) |

Classifier misses on the benchmark: `B049` (23× vs 10.7×, hard band) and
`B050` (49× vs 7.7×, hard band). Both are the two "hard" error bubbles the
specification asks for, so the coverage-only model fails exactly where the
spec predicts it should.

## Findings

### Critical

**C1. The global score pops true strain branches on the 50-bubble graph.**

- Files: `src/bubbleblower/score.py`, `src/bubbleblower/search.py`
- Evidence: `greedy_search(benchmark, max_iterations=60)` accepts 37 pops.
  25 remove error branches, **12 remove strain branches** (B003, B005,
  B008, B009, B012, B013, B016, B021–B025). Every false pop removes the
  lowest-abundance strain (strain_C at 20×, or strain_B when C is absent).
  Decomposition for B003 (64× vs 21×): Δcoverage +3.60, Δflow −3.59,
  Δtopology +0.05, Δcomplexity +0.19, Δtotal **+0.25**.
- Cause: `coverage_score` rewards removing any low-coverage instance from
  the mixture. The only counter-force is `flow_penalty`, and it is divided
  by `|coverage| + 1`, so a 20× branch lost from a 85× source costs very
  little. The bubble-count term (`−0.05` per bubble) and the complexity
  refund (`+0.15` per removed instance) both push towards popping.
- Expected (spec section 12 and 36): `S(G_{t+1}) > S(G_t)` must reflect
  evidence fit. Removing a 21× branch with its own colour is not a better
  explanation of the data.
- Recommended fix: make the likelihood term dominant and colour-aware.
  Score the popped branch's coverage as unexplained mass (Poisson or
  NegBin log-likelihood of observing 21× under an error model), and add a
  colour term: removing the only instance of a colour on a source→sink
  path should cost, because that colour's reads are then unexplained.
  Re-run the audit below and require `false_pops == 0` on the benchmark
  before this finding is closed.

**C2. The optimiser is not tested on the benchmark.**

- Files: `tests/test_bubbles.py`, `examples/benchmark/run.py`
- Evidence: `test_greedy_improves_score_and_pops_error` uses one error
  bubble and one strain bubble. `examples/benchmark/run.py` only runs the
  classifier. Spec section 35 asks for `number of false pops`, `number of
  missed errors`, and `ΔS` on the benchmark; none are computed anywhere.
  C1 was invisible to the suite for that reason.
- Recommended fix: add a mandatory test that runs greedy on
  `generate_bubble_benchmark(seed=42)` (or a 10-bubble slice for speed,
  spec row "iterative 5–10 bubbles") and asserts zero false pops and at
  least the very-low and medium error bands popped. Write
  `false_pops`, `missed_errors`, `delta_score` into
  `examples/benchmark/data/metrics.json`.

### Major

**M1. Greedy on the 50-bubble graph takes ~60 s; `max_runtime` is missing.**

- Files: `src/bubbleblower/search.py`, `src/bubbleblower/score.py`
- Evidence: 38 iterations in 61.7 s. Each iteration re-scores every
  candidate with `score_graph`, which calls `detect_bubbles` and a full
  k-means on all instances. MCMC with 40 steps took 74.7 s.
- Spec section 16 lists `max_runtime` as a hard limit. It is not
  implemented in any search.
- Recommended fix: add `max_runtime_s` to all three searches; cache the
  base score; score the candidate incrementally (only the touched bubble's
  neighbourhood changes flow and bubble count).

**M2. `bubble_results.tsv` does not record the edit per bubble.**

- File: `src/bubbleblower/report.py` lines 40–53
- Evidence: `selected_edit` is the type of the *first* accepted edit,
  written only on row 0. Every other row is blank. `iteration` is always
  `0`. Posteriors are computed once on the initial graph and never
  refreshed after edits.
- Spec section 17: each row needs `iteration`, `selected_edit`, and
  `confidence` for that bubble.
- Recommended fix: keep `bubble_id` on `Edit` (spec section 8.1
  `bubble_id:` field), map accepted edits back to bubbles by signature,
  and write one row per (bubble, iteration) when the bubble is
  re-detected after an edit.

**M3. Beam search's `state_id`/`parent_state` are not tracked.**

- File: `src/bubbleblower/report.py` lines 59–75
- Evidence: `state_id = step_index + 1`, `parent_state = step_index`.
  For beam search the parent of a kept state is not the previous row.
- Recommended fix: give `SearchStep` real `state_id` and `parent_state`
  fields set by each search.

**M4. `heldout_genera` is not a debubbling target as loaded.**

- Evidence: graph_type `knn`, 6439 unitigs, **24 links**, 0 bubbles.
  Coverage comes from the `_cov_` token in node names; link coverage
  comes from the same regex on `link_id`, which happens to match.
- The optional test `test_heldout_cdbg_loads` therefore only proves the
  importer works. Spec section 33 expects a graph "ready for debubbling".
  Either the CFA under `tocumg/cfa/` has more edges, or the export dropped
  them. This must be checked with the graph's owner; do not fabricate
  edges.
- Recommended fix: report the gap to the user. If the CFA carries the
  full edge set, load from `cfa/` via `cfa_to_cdbg` instead of `cdbg/`.

**M5. `environment.yml` pins nothing and pulls MetaMetro from `main`.**

- File: `environment.yml`
- Evidence: `python=3.12`, `numpy`, `pyyaml`, `pytest` unpinned;
  `pip: git+https://github.com/dsmutin/MetaMetro.git` has no ref.
  MetaMetro's own env pins `numpy=2.5.3`, `pyyaml=6.0.3`,
  `pytest=9.1.1`, `python=3.12.14`.
- Rule `reproducibility`: pin versions. CI will drift with MetaMetro
  `main`.
- Recommended fix: copy MetaMetro's pins and add `@<commit>` to the git
  URL. Also note that tests run locally against `../metametro/src` via
  `conftest.py`, while CI would use the pip install; the two can diverge.

### Minor

**m1. `generate.py::math_exp` is a pointless wrapper with an inline import.**
Replace with a module-level `import math`.

**m2. Split candidates are never accepted, so `duplicate`/`split`/`merge`
are exercised only by round-trip tests.** On the benchmark all 37 accepted
edits are pops. The "two paths through a shared node" case from spec
section 3 has a fixture (`shared_duplicate_node`) but no search test that
should end in a duplicate. Add a fixture where the correct answer is a split
and assert the optimiser finds it.

**m3. `split_instance` lines 319–320 contain a dead `if ...: pass`.**

**m4. `Bubble.parent_bubble_id` is declared but never set.** Spec 5 asks
for lineage of bubbles created by edits.

**m5. `Edit` has no `bubble_id` field** (spec 8.1). `reason` carries it as
free text only.

**m6. Coverage on the real graph is parsed from node names with a regex.**
This is documented in `from_cdbg`, but the coverage source is not written to
provenance (spec section 6 "coverage source"). Record it in
`AssemblyGraph` or the result metadata.

**m7. `examples/reads/run.py` reports `n_bubbles=0`** and does not assert
on it. The example demonstrates the MetaMetro pipeline, not debubbling.
Either choose two genomes that differ inside a shared flank (a real
bubble) or state in the README that this example is import-only.

**m8. `tests/test_search.py::test_beam_and_mcmc_do_not_lose_score`** asserts
`beam.scores[-1] >= greedy.scores[0]`, which is always true because beam
starts from the same graph and only accepts improvements. It does not test
that beam does anything.

### Suggestions

- `classification_metrics` computes AUROC only. Spec 35 also asks for
  AUPRC and per-class metrics for `strain`.
- `_means` k-means with 12 fixed iterations is fine for the MVP; if it stays,
  seed the initialisation explicitly in the docstring (it is deterministic
  via sorting, which is good).
- `examples/*/data/` outputs are tracked in git. That is acceptable for a
  toy, but the read-level FASTQ and CFA are regenerated on each run and
  will produce noisy diffs.

## What matches the specification

- CDBG is the only graph model; `graph_type` is never consulted (spec 37).
- Sequence hash and instance id are separated in `AssemblyGraph.lineage`;
  `duplicate` produces equal hashes with different ids (spec 3).
- `pop`, `duplicate`, `split`, `merge` are all reversible and
  `semantic_signature` compares by sequence, colours, coverage, and CFA
  provenance rather than id order (spec 18). Round-trip tests pass.
- The classifier returns probabilities and never edits (spec 7).
- Fixtures: strain, error (same colour on both branches), nested, shared
  node, rare-strain, equal-coverage, misleading-coverage (spec 20–28).
- Benchmark composition matches spec 23–25: 10 two-strain, 10
  strain-specific, 5 three-strain; error bands 10/8/5/2; Poisson coverage.
- Ground truth is written to a separate file and never read by the
  resolver (spec 29).
- Hard limits `max_iterations` and `max_instances` exist; `patience` exists
  for greedy.

## Conclusion

The classifier meets the benchmark thresholds (F1 0.96, AUROC 0.995), and
the CLI, examples, and tests run green. The stop criterion "all examples
return the expected output" is met only for the examples as written.

The acceptance criterion in spec 36 is **not met**: on the same 50-bubble
graph the greedy optimiser increases `S(G)` by removing 12 of 25 strain
branches. The score is currently a variance-reduction objective with a weak
flow term, which is the failure mode the specification warns about in
section 12. Findings C1 and C2 should be fixed before this version is
considered a working MVP. This review cannot rule out defects not covered by
the checks above.

## Audit script used

```python
from bubbleblower.generate import generate_bubble_benchmark
from bubbleblower.search import greedy_search

graph, truth = generate_bubble_benchmark(seed=42)
result = greedy_search(graph, max_iterations=60, patience=1)
ids = {u.unitig_id for u in result.graph.cdbg.unitigs}
false_pops = sum(
    1 for r in truth if r["type"] == "strain"
    and any(b not in ids for b in r["branch_ids"].split(","))
)
missed = sum(
    1 for r in truth if r["type"] == "error"
    and r["branch_ids"].split(",")[1] in ids
)
print(false_pops, missed, result.scores[0].total, result.scores[-1].total)
# observed: 12 0 -860.87 -722.50
```
