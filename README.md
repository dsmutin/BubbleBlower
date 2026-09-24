# bubbleblower

[![version](https://img.shields.io/badge/dynamic/yaml?url=https%3A%2F%2Fraw.githubusercontent.com%2Fdsmutin%2Fdivide-et-impera%2Fmain%2FVERSION&query=%24&label=version&color=blue)](VERSION)
[![required tests](https://img.shields.io/github/actions/workflow/status/dsmutin/divide-et-impera/required-tests.yml?branch=main&label=required%20tests)](https://github.com/dsmutin/divide-et-impera/actions/workflows/required-tests.yml)
[![full tests](https://img.shields.io/github/actions/workflow/status/dsmutin/divide-et-impera/full-tests.yml?branch=main&label=full%20tests)](https://github.com/dsmutin/divide-et-impera/actions/workflows/full-tests.yml)
[![warning](https://img.shields.io/badge/warning-in%20development-yellow)](https://shields.io/badges/static-badge)

Iterative debubbling of totally coloured assembly graphs

**Warning: in development.** Interfaces may change. See `VERSION` (single source of truth).

## Install

Conda is the only supported install:

```bash
conda env create -f environment.yml
conda activate bubbleblower
```

`environment.yml` sets `PYTHONPATH=src`. Do not publish a pip-first install path.

## Usage

```bash
bubbleblower --version
bubbleblower
python examples/toy/run.py
python examples/benchmark/run.py
python examples/error/run.py
python examples/reads/run.py
```

`bubbleblower` classifies the built-in two-taxon strain bubble and prints JSON.
Pass `-o` to write that JSON to a file. A CDBG directory can be passed later through `bubbleblower.pipeline.load_assembly`.
The benchmark writes `examples/benchmark/data/` and gates on **AMBER F1**
(harmonic mean of pop purity and error-removal completeness) for both
classification and greedy resolution. The graph has 50 bubbles (25 strain,
25 sequencing error, 3 strains). False pops of strain branches fail the run.

## Tests

```bash
pytest -m mandatory    # every commit
pytest                 # mandatory + optional (release / manual CI)
```

## License

MIT. See [CONTRIBUTING.md](CONTRIBUTING.md).
