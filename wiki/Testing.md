# Testing

## Test data

| Path | Role |
|------|------|
| `examples/toy/data/` | strain-bubble JSON from `run.py` |
| `examples/benchmark/data/` | 50-bubble ground truth; AMBER F1 for classify and resolve |
| `examples/error/data/` | resolved error bubble, score rises, branch `E` removed |
| `examples/reads/data/` | MetaMetro reads, coloured CFA, bubble count |
| `tests/` | contract tests; graphs are built in code |

Add real fixtures under `tests/data/` or `examples/toy/data/` and document them here. Do not invent datasets.

## Integrative testing

1. Mandatory unit/contract tests: `pytest -m mandatory`
2. CLI write path: `tests/test_integration.py`
3. Toy: `python examples/toy/run.py` (strain bubble, `decision=strain`)
4. Benchmark: `python examples/benchmark/run.py` (error F1 and AUROC)
5. Full suite (optional + toy + any vignettes): GitHub Action `full-tests.yml` on **release** or **workflow_dispatch**

Required CI on every push/PR: `required-tests.yml` (mandatory only), conda from `environment.yml`.
