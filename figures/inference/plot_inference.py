"""Publication figures for BubbleBlower inference on the checked-in examples.

Reads existing example tables and metric JSON. Does not rerun inference and
does not impute missing values. Outputs SVG and PDF under ``figures/inference/``.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pandas as pd

import cnsplots as cns

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
EXAMPLES = ROOT / "examples"

BENCHMARK_BUBBLES = EXAMPLES / "benchmark" / "data" / "bubble_results.tsv"
BENCHMARK_SCORES = EXAMPLES / "benchmark" / "data" / "resolved_run" / "state_scores.tsv"
BENCHMARK_METRICS = EXAMPLES / "benchmark" / "data" / "metrics.json"
ERROR_SUMMARY = EXAMPLES / "error" / "data" / "summary.json"
HALF = EXAMPLES / "half" / "data"
CLOSE_MODES = HALF / "metrics.json"
SPADES_DEBURBLE = HALF / "spades_debubble_metrics.json"
SPADES_COLOUR_BREAK = HALF / "colour_break_metrics.json"
SPADES_LOCAL = HALF / "local_spades_metrics.json"
FLYE_DEBURBLE = HALF / "flye_debubble_metrics.json"
FLYE_COLOUR_BREAK = HALF / "flye_colour_break_metrics.json"
MEGAHIT_DEBURBLE = HALF / "megahit_k21_debubble_metrics.json"
MEGAHIT_K21 = HALF / "megahit_k21_metrics.json"
MEGAHIT_COLOUR_BREAK = HALF / "megahit_colour_break_metrics.json"
MEGAHIT_K21_COLOUR_BREAK = HALF / "megahit_k21_colour_break_metrics.json"

REQUIRED = (
    BENCHMARK_BUBBLES,
    BENCHMARK_SCORES,
    BENCHMARK_METRICS,
    ERROR_SUMMARY,
    CLOSE_MODES,
    SPADES_DEBURBLE,
    SPADES_COLOUR_BREAK,
    SPADES_LOCAL,
    FLYE_DEBURBLE,
    FLYE_COLOUR_BREAK,
    MEGAHIT_DEBURBLE,
    MEGAHIT_K21,
    MEGAHIT_COLOUR_BREAK,
    MEGAHIT_K21_COLOUR_BREAK,
)

# Main scores used to decide whether a mode beats its assembler baseline.
# N50 is reported with them. It is not one of those four win criteria.
_COMPARED_METRICS = (
    ("genome_fraction", "Genome fraction"),
    ("misassemblies", "Misassemblies"),
    ("mismatches_per_100kbp", "Mismatches per 100 kbp"),
    ("duplication", "Duplication ratio"),
    ("n50", "N50 (bp)"),
)
_METRIC_KEYS = tuple(key for key, _label in _COMPARED_METRICS)


def _require_inputs() -> None:
    """Stop if a figure input is missing or empty."""
    missing = [str(path.relative_to(ROOT)) for path in REQUIRED if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("missing or empty inference inputs:\n" + "\n".join(missing))


def _load_json(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not payload:
        raise ValueError(f"expected a non-empty object: {path}")
    return payload


def _save(stem: str) -> None:
    cns.savefig(OUT / f"{stem}.svg")
    cns.savefig(OUT / f"{stem}.pdf")


def _benchmark_classification() -> None:
    """Confusion matrix and ROC for the 50-bubble graph."""
    bubbles = pd.read_csv(BENCHMARK_BUBBLES, sep="\t")
    expected = {"bubble_id", "type", "decision", "p_error", "p_strain"}
    if set(bubbles.columns) != expected:
        raise ValueError(f"unexpected bubble_results columns: {list(bubbles.columns)}")
    if bubbles[["type", "decision", "p_error"]].isna().any().any():
        raise ValueError("bubble_results has missing type, decision, or p_error")

    roc = bubbles.copy()
    if set(roc["type"]) != {"error", "strain"}:
        raise ValueError(f"ROC expects truth types error and strain, found {sorted(roc['type'].unique())}")
    roc["error"] = (roc["type"] == "error").astype(int)
    roc = roc.rename(columns={"p_error": "P(error)"})

    mp = cns.multipanel(max_width=520)
    ax_confusion = mp.panel("A", width=210, height=190, margin_right=20, margin_bottom=36)
    cns.confusionplot(
        bubbles,
        x="decision",
        y="type",
        x_order=["error", "strain", "other"],
        y_order=["error", "strain"],
        add_pvalue=False,
        ax=ax_confusion,
    )
    ax_confusion.set(xlabel="Predicted decision", ylabel="True type")

    ax_roc = mp.panel("B", width=200, height=190, margin_right=16, margin_bottom=28)
    cns.rocplot(roc, true_label_col="error", pred_prob_cols="P(error)", ax=ax_roc)
    ax_roc.set(xlabel="False positive rate", ylabel="True positive rate")
    roc_legend = ax_roc.get_legend()
    if roc_legend is not None:
        roc_legend.set_loc("lower right")
    _save("benchmark_classification")


def _rate_frame(metrics: dict, keys: list[tuple[str, str]], task: str) -> pd.DataFrame:
    """One row per named rate in a metrics object."""
    missing = [key for key, _label in keys if key not in metrics]
    if missing:
        raise KeyError(f"{task} metrics missing {missing}")
    return pd.DataFrame(
        {
            "metric": [label for _key, label in keys],
            "value": [float(metrics[key]) for key, _label in keys],
            "task": task,
        }
    )


def _benchmark_metrics() -> None:
    """Every rate and count stored for the 50-bubble benchmark."""
    metrics = _load_json(BENCHMARK_METRICS)
    classification = metrics["classification"]
    resolution = metrics["resolution"]
    if "amber_f1" not in metrics:
        raise KeyError("benchmark metrics missing amber_f1")

    shared = _rate_frame(
        classification,
        [("amber_f1", "AMBER F1"), ("precision_error", "Precision"), ("recall_error", "Recall")],
        "Classification",
    )
    shared = pd.concat(
        [
            shared,
            _rate_frame(
                resolution,
                [("amber_f1", "AMBER F1"), ("precision", "Precision"), ("recall", "Recall")],
                "Resolution",
            ),
        ],
        ignore_index=True,
    )
    extra = _rate_frame(
        classification,
        [
            ("precision_strain", "Strain precision"),
            ("recall_strain", "Strain recall"),
            ("f1_strain", "Strain F1"),
            ("auroc_error", "AUROC"),
            ("auprc_error", "AUPRC"),
        ],
        "Classification",
    )
    class_counts = _rate_frame(
        classification,
        [("tp", "True positives"), ("fp", "False positives"), ("fn", "False negatives"), ("tn", "True negatives")],
        "Classification",
    )
    resolve_counts = _rate_frame(
        resolution,
        [
            ("true_pops", "True pops"),
            ("false_pops", "False pops"),
            ("missed_errors", "Missed errors"),
            ("retained_strains", "Retained strains"),
        ],
        "Resolution",
    )
    if "delta_score" not in resolution:
        raise KeyError("resolution metrics missing delta_score")
    delta = pd.DataFrame({"metric": ["Score change"], "value": [float(resolution["delta_score"])]})

    mp = cns.multipanel(max_width=760)
    ax_shared = mp.panel("A", width=220, height=190, margin_right=120, margin_bottom=40)
    cns.barplot(
        shared,
        x="metric",
        y="value",
        hue="task",
        order=["AMBER F1", "Precision", "Recall"],
        hue_order=["Classification", "Resolution"],
        ax=ax_shared,
    )
    ax_shared.set(xlabel="Metric", ylabel="Value")
    ax_shared.set_ylim(0, 1.08)
    cns.take_legend_out(ax=ax_shared, title="Task")

    ax_extra = mp.panel("B", width=250, height=190, margin_right=16, margin_bottom=52)
    cns.lollipopplot(
        extra,
        x="metric",
        y="value",
        order=["Strain precision", "Strain recall", "Strain F1", "AUROC", "AUPRC"],
        ax=ax_extra,
    )
    ax_extra.set(xlabel="Classification metric", ylabel="Value")
    ax_extra.set_ylim(0, 1.08)

    ax_class = mp.panel("C", width=220, height=180, margin_right=16, margin_bottom=48)
    cns.lollipopplot(
        class_counts,
        x="metric",
        y="value",
        order=["True positives", "False positives", "False negatives", "True negatives"],
        ax=ax_class,
    )
    ax_class.set(xlabel="Classification count", ylabel="Bubbles")

    ax_resolve = mp.panel("D", width=230, height=180, margin_right=16, margin_bottom=48)
    cns.lollipopplot(
        resolve_counts,
        x="metric",
        y="value",
        order=["True pops", "False pops", "Missed errors", "Retained strains"],
        ax=ax_resolve,
    )
    ax_resolve.set(xlabel="Resolution count", ylabel="Bubbles")

    ax_delta = mp.panel("E", width=120, height=180, margin_right=16, margin_bottom=40)
    cns.barplot(delta, x="metric", y="value", ax=ax_delta)
    ax_delta.set(xlabel="Resolution", ylabel="Global score change")
    _save("benchmark_metrics")


def _benchmark_search() -> None:
    """Global score accepted during greedy resolution of the benchmark graph."""
    scores = pd.read_csv(BENCHMARK_SCORES, sep="\t")
    for column in ("iteration", "score"):
        if column not in scores.columns:
            raise ValueError(f"state_scores missing {column}")
    if scores.empty:
        raise ValueError("state_scores is empty")
    cns.figure(width=240, height=180)
    ax = cns.lineplot(scores, x="iteration", y="score", errorbar=None, marker="o")
    ax.set(xlabel="Accepted iteration", ylabel="Global score")
    _save("benchmark_search")


def _error_score() -> None:
    """Score of the single error-bubble example before and after the edit."""
    summary = _load_json(ERROR_SUMMARY)
    for key in ("score_before", "score_after"):
        if key not in summary:
            raise KeyError(f"error summary missing {key}")
    frame = pd.DataFrame(
        {
            "example": ["Error bubble", "Error bubble"],
            "stage": ["Before", "After"],
            "score": [float(summary["score_before"]), float(summary["score_after"])],
            "pair": ["error", "error"],
        }
    )
    cns.figure(width=180, height=180)
    ax = cns.slopeplot(
        frame,
        x="example",
        y="score",
        hue="stage",
        pair="pair",
        hue_order=["Before", "After"],
    )
    ax.set(xlabel="Example", ylabel="Global score")
    cns.take_legend_out(ax=ax, title="Stage")
    _save("error_score")


def _metric_signature(metrics: dict, source: str) -> tuple[float, ...]:
    """Require the compared metric keys and return their values."""
    missing = [key for key in _METRIC_KEYS if key not in metrics]
    if missing:
        raise KeyError(f"{source} lacks {missing}")
    return tuple(float(metrics[key]) for key in _METRIC_KEYS)


class _MethodTable:
    """One score vector per graph, dataset, and method."""

    def __init__(self) -> None:
        self._values: dict[tuple[str, str, str], tuple[float, ...]] = {}

    def add(self, graph: str, dataset: str, method: str, metrics: dict, source: str) -> None:
        """Store a method. The same method must carry the same scores."""
        signature = _metric_signature(metrics, f"{source} {dataset} {method}")
        key = (graph, dataset, method)
        previous = self._values.get(key)
        if previous is not None and previous != signature:
            raise ValueError(f"conflicting scores for {key} from {source}")
        self._values[key] = signature

    def add_unless_duplicate(self, graph: str, dataset: str, method: str, metrics: dict, source: str) -> None:
        """Skip a baseline whose scores already belong to a named method on that dataset."""
        signature = _metric_signature(metrics, f"{source} {dataset} {method}")
        for (stored_graph, stored_dataset, _stored_method), stored in self._values.items():
            if stored_graph == graph and stored_dataset == dataset and stored == signature:
                return
        self.add(graph, dataset, method, metrics, source)

    def frame(self, graph: str) -> pd.DataFrame:
        """Long table for one assembler family."""
        rows: list[dict[str, object]] = []
        for (stored_graph, dataset, method), signature in self._values.items():
            if stored_graph != graph:
                continue
            for (key, label), value in zip(_COMPARED_METRICS, signature, strict=True):
                rows.append(
                    {
                        "dataset": dataset,
                        "method": method,
                        "metric": label,
                        "metric_key": key,
                        "value": value,
                    }
                )
        frame = pd.DataFrame(rows)
        if frame.empty:
            raise ValueError(f"no methods for {graph}")
        return frame


def _add_mode_file(table: _MethodTable, path: Path, graph: str, baseline_key: str, baseline_label: str) -> None:
    """Add debubbler modes and the named assembler baseline from one JSON file."""
    payload = _load_json(path)
    for dataset, block in payload.items():
        modes = block.get("modes")
        baseline = block.get(baseline_key)
        if not isinstance(modes, dict) or not isinstance(baseline, dict):
            raise ValueError(f"{path.name} dataset {dataset} lacks modes or {baseline_key}")
        for mode, mode_block in modes.items():
            metrics = mode_block.get("metrics") if isinstance(mode_block, dict) else None
            if not isinstance(metrics, dict):
                raise ValueError(f"{path.name} {dataset} {mode} lacks metrics")
            table.add(graph, dataset, mode, metrics, path.name)
        table.add(graph, dataset, baseline_label, baseline, path.name)


def _add_pair_file(table: _MethodTable, path: Path, graph: str, labels: dict[str, str]) -> None:
    """Add every named score block in a dataset object."""
    payload = _load_json(path)
    for dataset, block in payload.items():
        if not isinstance(block, dict):
            raise ValueError(f"{path.name} dataset {dataset} is not an object")
        found = False
        for key, label in labels.items():
            metrics = block.get(key)
            if isinstance(metrics, dict) and _METRIC_KEYS[0] in metrics:
                table.add(graph, dataset, label, metrics, path.name)
                found = True
        if not found:
            raise ValueError(f"{path.name} dataset {dataset} has none of {list(labels)}")


def _comparison_table() -> _MethodTable:
    """Baselines and BubbleBlower modes for every scored assembly example."""
    table = _MethodTable()
    _add_mode_file(table, SPADES_DEBURBLE, "metaSPAdes", "metaspades_pop", "metaSPAdes pop")
    _add_mode_file(table, FLYE_DEBURBLE, "metaFlye", "metaflye", "metaFlye")
    _add_mode_file(table, MEGAHIT_DEBURBLE, "MEGAHIT", "megahit_final", "MEGAHIT final")
    close = _load_json(CLOSE_MODES)
    close_pop = _load_json(SPADES_DEBURBLE)["close_k33"]["metaspades_pop"]
    if _metric_signature(close["metaspades"], "metrics.json metaspades") != _metric_signature(
        close_pop, "spades close_k33 metaSPAdes pop"
    ):
        raise ValueError("metrics.json metaSPAdes scores do not match close_k33 metaSPAdes pop")
    for mode, metrics in close["modes"].items():
        table.add("metaSPAdes", "close_k33", mode, metrics, "metrics.json")
    _add_pair_file(
        table,
        SPADES_LOCAL,
        "metaSPAdes",
        {"metaspades_pop": "metaSPAdes pop", "metaspades_keep": "metaSPAdes keep"},
    )
    _add_pair_file(table, MEGAHIT_K21, "MEGAHIT", {"megahit_k21": "MEGAHIT k21", "megahit_final": "MEGAHIT final"})
    for path, graph in (
        (SPADES_COLOUR_BREAK, "metaSPAdes"),
        (FLYE_COLOUR_BREAK, "metaFlye"),
        (MEGAHIT_COLOUR_BREAK, "MEGAHIT"),
        (MEGAHIT_K21_COLOUR_BREAK, "MEGAHIT"),
    ):
        payload = _load_json(path)
        for dataset, block in payload.items():
            for key, metrics in block.items():
                if not isinstance(metrics, dict) or "genome_fraction" not in metrics:
                    continue
                label = {
                    "read_colour_break": "read_colour_break",
                    "k21_read_colour_break": "k21_read_colour_break",
                    "metaflye": "metaFlye",
                    "megahit_final": "MEGAHIT final",
                }.get(key, key)
                if key == "baseline":
                    table.add_unless_duplicate(graph, dataset, "baseline", metrics, path.name)
                else:
                    table.add(graph, dataset, label, metrics, path.name)
    return table


def _plot_method_family(frame: pd.DataFrame, graph: str, stem: str, width: int) -> None:
    """One panel per compared metric. Methods and baselines share the axis."""
    methods = [method for method in _METHOD_ORDER if method in set(frame["method"])]
    extra = sorted(set(frame["method"]) - set(methods))
    if extra:
        raise ValueError(f"unordered methods for {graph}: {extra}")
    datasets = list(dict.fromkeys(frame["dataset"]))
    mp = cns.multipanel(max_width=width + 180)
    first = True
    for _key, ylabel in _COMPARED_METRICS:
        subset = frame[frame["metric"] == ylabel]
        margin_right = 160 if first else 16
        ax = mp.panel(width=width, height=180, margin_right=margin_right, margin_bottom=48, margin_top=8)
        cns.barplot(
            subset,
            x="dataset",
            y="value",
            hue="method",
            order=datasets,
            hue_order=methods,
            ax=ax,
        )
        ax.set(xlabel=graph, ylabel=ylabel)
        if first:
            cns.take_legend_out(ax=ax, title="Method")
            first = False
        elif ax.get_legend() is not None:
            ax.get_legend().remove()
    _save(stem)


_METHOD_ORDER = (
    "retain",
    "colour_pop",
    "kmer_divergence",
    "colour_topology",
    "read_colour_break",
    "k21_read_colour_break",
    "MEGAHIT k21",
    "metaSPAdes keep",
    "baseline",
    "metaSPAdes pop",
    "metaFlye",
    "MEGAHIT final",
)


def _method_comparisons() -> None:
    """Compare every stored mode with the assembler baseline on each example."""
    table = _comparison_table()
    _plot_method_family(table.frame("metaSPAdes"), "metaSPAdes", "spades_methods", 340)
    _plot_method_family(table.frame("metaFlye"), "metaFlye", "flye_methods", 220)
    _plot_method_family(table.frame("MEGAHIT"), "MEGAHIT", "megahit_methods", 380)


def main() -> None:
    """Write the inference figures next to this script."""
    _require_inputs()
    cns.settings.font_sans_serif = ["DejaVu Sans"]
    _benchmark_classification()
    _benchmark_metrics()
    _benchmark_search()
    _error_score()
    _method_comparisons()


if __name__ == "__main__":
    main()
