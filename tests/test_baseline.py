"""Mandatory: baseline pipeline contract."""

from __future__ import annotations

import pytest

from bubbleblower.baseline import run_pipeline

pytestmark = pytest.mark.mandatory


def test_run_pipeline_keys() -> None:
    """Baseline result has status, ok, and input_path."""
    result = run_pipeline()
    assert set(result) >= {"status", "ok", "input_path"}
    assert result["status"] == "resolved"
    assert result["ok"] is True
    assert result["decision"] == "strain"
