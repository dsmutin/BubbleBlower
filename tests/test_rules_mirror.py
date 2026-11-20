"""Mandatory: `.cursor/rules/` is a byte-identical mirror of `agents/rules/`."""

from __future__ import annotations

from pathlib import Path

import pytest

pytestmark = pytest.mark.mandatory


def _rule_names(directory: Path) -> set[str]:
    return {path.name for path in directory.glob("*.mdc")}


def test_cursor_rules_mirror_has_same_files(repo_root: Path) -> None:
    """Every rule in `agents/rules/` has a mirror, and the mirror has no extra rule."""
    source = _rule_names(repo_root / "agents" / "rules")
    mirror = _rule_names(repo_root / ".cursor" / "rules")
    assert source, "agents/rules/ has no .mdc rules"
    assert source == mirror, (
        f"missing in .cursor/rules: {sorted(source - mirror)}; "
        f"extra in .cursor/rules: {sorted(mirror - source)}"
    )


@pytest.mark.parametrize(
    "name",
    sorted(_rule_names(Path(__file__).resolve().parents[1] / "agents" / "rules")),
)
def test_cursor_rule_is_byte_identical(repo_root: Path, name: str) -> None:
    """`.cursor/rules/<name>` must equal `agents/rules/<name>` byte for byte."""
    source = (repo_root / "agents" / "rules" / name).read_bytes()
    mirror_path = repo_root / ".cursor" / "rules" / name
    assert mirror_path.is_file(), f"{mirror_path} is missing; copy agents/rules/{name}"
    assert mirror_path.read_bytes() == source, (
        f".cursor/rules/{name} differs from agents/rules/{name}; "
        f"copy the source of truth: cp agents/rules/{name} .cursor/rules/{name}"
    )
