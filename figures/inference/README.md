# Inference figures

Publication figures for checked-in BubbleBlower examples. Built with `cnsplots` 0.6.0 by `plot_inference.py`. The script reads existing tables and does not rerun inference.

| File | What it shows | Source |
|------|----------------|--------|
| `benchmark_metrics` | AMBER F1, precision, and recall for classification and resolution; strain precision, strain recall, strain F1, AUROC, and AUPRC; confusion counts; resolution counts; global score change | `examples/benchmark/data/metrics.json` |
| `benchmark_classification` | Decision confusion and ROC of rounded `p_error` for the error class | `examples/benchmark/data/bubble_results.tsv` |
| `benchmark_search` | Global score at each accepted greedy iteration | `examples/benchmark/data/resolved_run/state_scores.tsv` |
| `error_score` | Global score before and after resolving the error bubble | `examples/error/data/summary.json` |
| `spades_methods` | Genome fraction, misassemblies, mismatches per 100 kbp, duplication ratio, and N50 for metaSPAdes modes and baselines | `examples/half/data/metrics.json`, `spades_debubble_metrics.json`, `colour_break_metrics.json`, `local_spades_metrics.json` |
| `flye_methods` | The same scores for metaFlye modes and the metaFlye baseline | `examples/half/data/flye_debubble_metrics.json`, `flye_colour_break_metrics.json` |
| `megahit_methods` | The same scores for MEGAHIT k21 modes, colour breaks, the k21 contig set, and the final MEGAHIT baseline | `examples/half/data/megahit_k21_debubble_metrics.json`, `megahit_k21_metrics.json`, `megahit_colour_break_metrics.json`, `megahit_k21_colour_break_metrics.json` |

Classification `f1_error` in the benchmark file is the same value as classification AMBER F1, so it is the AMBER F1 bar. The ROC uses the four-decimal `p_error` column. `auroc_error` on the metrics figure is the unrounded value stored in `examples/benchmark/data/metrics.json`.

AMBER F1 is stored only for that 50-bubble benchmark. The assembly examples are compared on the four main scores (genome fraction, misassemblies, mismatches per 100 kbp, duplication) plus N50. Each metric has its own axis. A baseline whose scores are identical to a named method on the same dataset is not drawn a second time. `retain` and `colour_pop` on `close_k33` are both drawn. The toy, reads, and error examples do not store these scores.

```bash
python figures/inference/plot_inference.py
```
