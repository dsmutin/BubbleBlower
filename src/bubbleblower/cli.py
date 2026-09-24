"""Command-line entry for bubbleblower."""

from __future__ import annotations

import argparse
import json

from bubbleblower import __version__
from bubbleblower.baseline import run_pipeline
from bubbleblower.fixtures import strain_bubble
from bubbleblower.pipeline import load_assembly
from bubbleblower.report import write_result
from bubbleblower.search import beam_search, greedy_search, mcmc_search


def _resolve_cli(cdbg: str | None, out_dir: str | None, search: str) -> dict:
    """Resolve one graph and optionally write the contracted tables."""
    graph = load_assembly(cdbg) if cdbg else strain_bubble()
    if search == "beam":
        found = beam_search(graph)
    elif search == "mcmc":
        found = mcmc_search(graph)
    else:
        found = greedy_search(graph)
    if out_dir:
        write_result(found, out_dir)
    before = found.scores[0].total
    after = found.scores[-1].total
    return {
        "status": "resolved",
        "ok": after + 1e-9 >= before,
        "input_path": cdbg,
        "search": "greedy" if search == "none" else search,
        "score_before": before,
        "score_after": after,
        "n_edits": len(found.edits),
        "n_bubbles": len(found.bubbles),
        "out_dir": out_dir,
    }


def main(argv: list[str] | None = None) -> int:
    """Parse CLI arguments and run the baseline pipeline."""
    parser = argparse.ArgumentParser(prog="bubbleblower", description="Iterative debubbling of totally coloured assembly graphs")
    parser.add_argument("--version", action="store_true", help="print version and exit")
    parser.add_argument("-o", "--output", default="-", help="JSON output path or - for stdout")
    parser.add_argument("--cdbg", default=None, help="MetaMetro CDBG directory to resolve")
    parser.add_argument("--out-dir", default=None, help="directory for resolved graph and TSV tables")
    parser.add_argument(
        "--search",
        choices=("none", "greedy", "beam", "mcmc"),
        default="none",
        help="search applied when --cdbg or --out-dir is set",
    )
    args = parser.parse_args(argv)
    if args.version:
        print(__version__)
        return 0
    if args.cdbg or args.out_dir:
        result = _resolve_cli(args.cdbg, args.out_dir, args.search)
    else:
        result = run_pipeline()
    text = json.dumps(result, indent=2)
    if args.output in {"", "-"}:
        print(text)
    else:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
