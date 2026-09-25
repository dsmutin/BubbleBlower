"""Locate a MetaMetro bench directory and minimap2 without a machine path."""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path


def _metametro_src() -> None:
    try:
        import metametro.bench  # noqa: F401
    except ImportError:
        sibling = Path(__file__).resolve().parents[3] / "metametro" / "src"
        env = os.environ.get("METAMETRO_SRC", "")
        for candidate in (Path(env) if env else None, sibling):
            if candidate is not None and candidate.is_dir():
                sys.path.insert(0, str(candidate))
                return
        raise SystemExit("metametro is not installed, METAMETRO_SRC is unset, and ../metametro/src is missing")


def work_dir(legacy_name: str) -> Path:
    """Return the ``work/`` directory of a MetaMetro community build."""
    _metametro_src()
    from metametro.bench.paths import default_outdir
    from metametro.bench.registry import resolve

    return default_outdir(resolve(legacy_name)) / "work"


def minimap2() -> str:
    """Return the minimap2 executable from ``MINIMAP2`` or ``PATH``."""
    found = os.environ.get("MINIMAP2") or shutil.which("minimap2")
    if not found:
        raise SystemExit("minimap2 is not on PATH and MINIMAP2 is unset")
    return found
