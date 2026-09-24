# Contracts

Scaffold contracts from `/start-dev`. Add tool-specific rows when the product spec is known. Change a row only together with tests.

## Implementation table

| Contract | Implementation |
|----------|----------------|
| Single version source (`VERSION`) | baseline (`bubbleblower.__version__` reads `VERSION`) |
| CLI `--version` and JSON run | baseline (`bubbleblower.cli.main`) |
| Pipeline result keys `status`, `ok`, `input_path` | baseline (`bubbleblower.baseline.run_pipeline`) |
| Conda-only install (`environment.yml`) | baseline (env file + `PYTHONPATH=src`) |
| Mandatory pytest on every commit | baseline (`pytest -m mandatory`) |
| Optional pytest on release / manual | baseline (`pytest -m optional` placeholder `pass`) |
| Toy example runs the tool | `examples/toy/run.py` classifies the strain bubble |
| Bubble detection | `bubbleblower.detect.detect_bubbles` |
| Coverage classifier does not edit | `bubbleblower.classify.classify_bubble` |
| Reversible pop, duplicate, split, merge | `bubbleblower.edits` |
| Global score and greedy search | `bubbleblower.score`, `bubbleblower.search.greedy_search` |
| 50-bubble graph-only benchmark | `bubbleblower.generate.generate_bubble_benchmark` |
| Beam search and MCMC | `bubbleblower.search.beam_search`, `mcmc_search` |
| Output tables and resolved CDBG | `bubbleblower.report.write_result` |
| Read-level graph via MetaMetro | `examples/reads/run.py` |
| English docs on public APIs | baseline (module/function docstrings) |
| No `git push` unless the human asks | baseline (rule `no-push`) |
