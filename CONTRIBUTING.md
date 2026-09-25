# Contributing to bubbleblower

All development follows this guide. Project rules require it (`follow-contributing`).

English is required for every public function, class, module, CLI flag, data file, and user-facing document. Other languages or missing documentation are not allowed.

## Architecture

```
CLI (bubbleblower.cli)
  → baseline pipeline (bubbleblower.baseline.run_pipeline)
      → JSON result {status, ok, input_path}

tests/          mandatory vs optional pytest
examples/toy/   end-to-end run of the current (baseline) tool
cite/           BibTeX for integrated third-party tools
agents/         portable rules and skills (any IDE)
```

The wiki is the [GitHub wiki](https://github.com/dsmutin/BubbleBlower/wiki), not a directory in this repository. Keep the documented return keys until you change the contract on the [Contracts](https://github.com/dsmutin/BubbleBlower/wiki/Contracts) page and the tests together.

## Testing architecture

| Kind | Marker | Command | When |
|------|--------|---------|------|
| Required | `mandatory` | `pytest -m mandatory` | every commit; GitHub Action `required-tests` on every push and pull request |
| Optional | `optional` | `pytest` (all) | release or workflow_dispatch; Action `full-tests` |
| Examples | — | `examples/toy`, `examples/error`, `examples/reads`, `examples/benchmark` | full CI; toy also after CLI changes |
| Local only | — | `examples/half/run.py` | needs reads, a graph, and minimap2; not in CI |

Do not mark a contract test `optional`. Optional tests are slow, extra, or nice-to-have.

After **any new feature**, run the **mandatory** suite (and toy if the CLI changed) before you stop.

## Feature checklist (`todo.md`)

Track work in `todo.md` (checkboxes). One line per feature or fix. Check it off only when mandatory tests pass. This is a **feature list**, not a `/do` analysis graph.

## Versioning

Edit **only** `VERSION`. Everything else reads it.

Starting value: `0.0.1`.

| Change | Bump |
|--------|------|
| New feature | **minor** (`0.0.1` → `0.1.0`) |
| Fix or update of an existing feature | **patch** (`0.1.0` → `0.1.1`) |
| Release | **major** (`0.1.1` → `1.0.0`) |

## Install (conda only)

```bash
conda env create -f environment.yml
conda activate bubbleblower
```

## GitHub

Never `git push` unless the human explicitly asks. CI runs on GitHub after they push.

## Benchmarks

Do not add a new benchmark or pinned community in this repository. Add it in MetaMetro (`metametro benchbuild`) and open a pull request there. This package reads that build.

Do not hard-code a machine path. Use the MetaMetro bench directory, a CLI argument, or an environment variable (`METAMETRO_SRC`, `MINIMAP2`). If it is missing, stop.

Do not mock a benchmark graph, a taxonomy label, or a metric. A tiny graph that exists only inside a test file stays in that test.

Do not copy an evaluation target into graph features, colours, or any file a model reads as input. Ground truth stays in `ground_truth/` and is used only by the scorer.

## Colourings

Do not add a new colouring method in this repository. Add it in MetaMetro (`metametro.bench.colourings`) and open a pull request there.

Do not mock a colouring. Paint graphs with MetaMetro colourings, or load the ToCUMG namespaces the bench already wrote.

Check that a selected colouring exists on the current MetaMetro build (`metametro benchbuild --list-colourings` and the bench `manifest.yaml`). If the namespace is missing, stop.

Use MetaMetro to select layers:

```python
from bubbleblower.bench_input import load_bench_graph

root, graph = load_bench_graph("bubble_strain_2", namespaces=("taxon", "composition"))
```

## Taxon identifiers

Two checks apply to every benchmark and to every graph tensor a model trains on.

**Pre-generation taxid.** A benchmark colour on a read is the genome id assigned before the metagenome was simulated. That id is the simulator manifest or MetaMetro `ground_truth/read_to_genome.tsv`. It is not a mock label written onto nodes, not a Kraken call, and not a model prediction. If that table is missing, stop. Do not invent a substitute.

**No taxid leakage into the GCN.** The true taxid is for evaluation only. It must not be copied into the feature matrix (`X_node`, `X_edge`) or into the colour channels of the tensor (`C_node`, `C_edge`, and any colour-weight channel). A colour the debubbler or the GCN is allowed to see is a sample colour or a prediction, not that ground-truth taxid.

## Citations

Add a `.bib` entry in `cite/` only for tools this package actually integrates. Do not invent papers.
