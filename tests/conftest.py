"""Shared pytest fixtures."""

from __future__ import annotations

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_METAMETRO = ROOT.parent / "metametro" / "src"
if _METAMETRO.is_dir():
    import sys

    sys.path.insert(0, str(_METAMETRO))


@pytest.fixture
def repo_root() -> Path:
    """Return the repository root (directory that contains VERSION)."""
    return ROOT
