"""Mandatory: built-in strain bubble contract."""

from __future__ import annotations

import pytest

from bubbleblower.pipeline import run_demo

pytestmark = pytest.mark.mandatory


def test_run_demo_keys() -> None:
    """Demo result has status, ok, and input_path."""
    result = run_demo()
    assert set(result) >= {"status", "ok", "input_path"}
    assert result["status"] == "resolved"
    assert result["ok"] is True
    assert result["decision"] == "strain"
