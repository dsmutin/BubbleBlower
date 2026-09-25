# BubbleBlower🫧 <img src="BubbleBlower.png" align="right" width="150" alt="BubbleBlower logo">

[![version](https://img.shields.io/badge/dynamic/yaml?url=https%3A%2F%2Fraw.githubusercontent.com%2Fdsmutin%2FBubbleBlower%2Fmaster%2FVERSION&query=%24&label=version&color=blue)](VERSION)
[![required tests](https://img.shields.io/github/actions/workflow/status/dsmutin/BubbleBlower/required-tests.yml?branch=master&label=required%20tests)](https://github.com/dsmutin/BubbleBlower/actions/workflows/required-tests.yml)
[![full tests](https://img.shields.io/github/actions/workflow/status/dsmutin/BubbleBlower/full-tests.yml?branch=master&label=full%20tests)](https://github.com/dsmutin/BubbleBlower/actions/workflows/full-tests.yml)
[![license](https://img.shields.io/github/license/dsmutin/BubbleBlower)](LICENSE)

Iterative debubbling of a MetaMetro ToCUMG (totally coloured universal metagenomic graph).

The graph is an input. BubbleBlower adopts a ToCUMG from `cfa_to_cdbg` or `load_cdbg`, or binds a coloured graph tensor to the ToCUMG it was built from. It does not assemble unitigs, and a tensor is not turned into another graph. GFA and MEGAHIT FASTG are loaded with MetaMetro's `gfa_to_cfa` and `fastg_to_cfa`.

Strain variation and sequencing error both appear as bubbles: two or more paths between the same source and sink. A colour-aware score keeps an edit only when that score rises.

Concepts and contracts: [wiki](https://github.com/dsmutin/BubbleBlower/wiki). Version: [`VERSION`](VERSION).

## How it works

1. **Load** a ToCUMG directory (`--cdbg`), or the built-in strain bubble when no graph is given. That bubble is a CFA compacted by MetaMetro.
2. **Detect** bubbles between a shared source and sink.
3. **Classify** each bubble. Two Bayesian posteriors exist: `coverage` (Poisson coverage, model `coverage-poisson-1`) and `multimodal` (coverage plus colour equality, long-read linkage, and k-mer support). They return probabilities. They do not edit the graph and they do not choose merge versus split.
4. **Search** reversible edits — pop a branch, split a shared node, or merge adjacent unitigs of the same colour — and keep an edit only when the global score rises. Search is greedy, beam, or MCMC.
5. **Write** the resolved ToCUMG and three tables: `bubble_results.tsv`, `edit_history.tsv`, `state_scores.tsv`.

Graph-only modes `kmer_divergence` and `colour_topology` are fixed thresholds on sequence and colour. They label a simple bubble as error, variation, or multi, then pop, retain, or split. They are not a trained model.

## What the implementation is expected to do

Each iteration either merges two nodes (inside a bubble or adjacent) and unions their colours, or splits one node and partitions its colour. A binary channel becomes a 0/1 mask. A gradual channel keeps a weight in `[0, 1]`. A launch criterion chooses the edit. The Bayesian posteriors above are not that switch. An ML criterion on the MetaMetro GCN (NumPy or PyTorch Geometric) is not implemented.

Debubbling with those edits is to be compared with the assembler's own bubble removal on AMBER F1, contig F1, misassemblies, mismatches, N50, and runtime, from pinned inputs. That comparison is not in the tree. The synthetic benchmark (`examples/benchmark/run.py`, seed 42) scores AMBER F1 on 50 fixture bubbles. Its colours are labels on the CFA, not taxids assigned before read simulation. See [CONTRIBUTING.md](CONTRIBUTING.md).

## Install

Conda is the only supported install:

```bash
conda env create -f environment.yml
conda activate bubbleblower
```

`environment.yml` sets `PYTHONPATH=src`.

## Usage

```bash
bubbleblower --version
bubbleblower
bubbleblower --cdbg path/to/cdbg --out-dir out --search greedy
bubbleblower -o result.json
```

With no graph, `bubbleblower` classifies the built-in strain bubble and prints JSON (`status`, `ok`, `decision`). `--search` is `greedy` (default), `beam`, or `mcmc`. Optional coverage tables `node_coverage.tsv` and `link_coverage.tsv` in the CDBG directory are read when present; coverage is not inferred when they are absent.

```bash
python examples/toy/run.py
python examples/benchmark/run.py
python examples/error/run.py
python examples/reads/run.py
python examples/half/run.py
```

## Tests

Required tests run on every push and every pull request.

```bash
pytest -m mandatory    # every commit
pytest                 # mandatory and optional
```

Optional tests and the examples run in `full-tests.yml` on a published release or a manual dispatch.

## License

MIT. See [CONTRIBUTING.md](CONTRIBUTING.md).
