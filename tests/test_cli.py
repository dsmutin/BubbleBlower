"""Mandatory: CLI entry works."""

from __future__ import annotations

import json

import pytest

from bubbleblower.cli import main

pytestmark = pytest.mark.mandatory


def test_cli_version(capsys: pytest.CaptureFixture[str]) -> None:
    """--version prints the package version."""
    assert main(["--version"]) == 0
    out = capsys.readouterr().out.strip()
    assert out


def test_cli_baseline_stdout(capsys: pytest.CaptureFixture[str]) -> None:
    """Default CLI run prints JSON for the built-in strain bubble."""
    assert main([]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["status"] == "resolved"
    assert payload["ok"] is True
    assert payload["decision"] == "strain"


def test_cli_search_defaults_to_greedy(tmp_path, capsys: pytest.CaptureFixture[str]) -> None:
    """Without --search the CLI runs and reports greedy search."""
    assert main(["--out-dir", str(tmp_path / "out")]) == 0
    assert json.loads(capsys.readouterr().out)["search"] == "greedy"


def test_cli_search_none_is_greedy_alias(tmp_path, capsys: pytest.CaptureFixture[str]) -> None:
    """The deprecated value none still runs greedy search."""
    assert main(["--out-dir", str(tmp_path / "out"), "--search", "none"]) == 0
    assert json.loads(capsys.readouterr().out)["search"] == "greedy"
