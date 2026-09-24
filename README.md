# BubbleBlower🫧 <img src="BubbleBlower.png" align="right" width="150" alt="BubbleBlower logo">

[![version](https://img.shields.io/badge/dynamic/yaml?url=https%3A%2F%2Fraw.githubusercontent.com%2Fdsmutin%2FBubbleBlower%2Fmaster%2FVERSION&query=%24&label=version&color=blue)](VERSION)
[![required tests](https://img.shields.io/github/actions/workflow/status/dsmutin/BubbleBlower/required-tests.yml?branch=master&label=required%20tests)](https://github.com/dsmutin/BubbleBlower/actions/workflows/required-tests.yml)
[![full tests](https://img.shields.io/github/actions/workflow/status/dsmutin/BubbleBlower/full-tests.yml?branch=master&label=full%20tests)](https://github.com/dsmutin/BubbleBlower/actions/workflows/full-tests.yml)
[![license](https://img.shields.io/github/license/dsmutin/BubbleBlower)](LICENSE)

Iterative debubbling of totally coloured assembly graphs.

Strain variation and sequencing error both appear as bubbles: two or more paths between the same source and sink. BubbleBlower tells those cases apart on a graph whose nodes and links already carry taxon colours, then removes an error branch only when that edit raises a colour-aware score. True strain branches stay in the graph.

Concepts and contracts: [wiki](https://github.com/dsmutin/BubbleBlower/wiki). Version: [`VERSION`](VERSION).

## How it works

1. **Load** a MetaMetro compacted de Bruijn graph (`--cdbg`), or the built-in two-taxon strain bubble when no graph is given.
2. **Detect** bubbles between a shared source and sink.
3. **Classify** each bubble. The coverage model returns probabilities of sequencing error, strain variation, and other. Classification does not edit the graph.
4. **Search** reversible edits — pop a branch, split a shared node, or merge adjacent unitigs of the same colour — and keep an edit only when the global score rises. Search is greedy, beam, or MCMC.
5. **Write** the resolved graph and three tables: `bubble_results.tsv`, `edit_history.tsv`, `state_scores.tsv`.

Graph-only modes `kmer_divergence` and `colour_topology` label a simple bubble as error, variation, or multi from sequence and colour features, then pop, retain, or split. They do not read ground-truth labels.

The primary quality metric is **AMBER F1**: the harmonic mean of purity and completeness. The graph-only benchmark (`examples/benchmark/run.py`) builds 50 bubbles (25 strain, 25 sequencing error, 3 strains). A false pop of a strain branch fails that run.

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
