"""Small coloured graphs with known bubble labels."""

from __future__ import annotations

from bubbleblower.graph import AssemblyGraph, build_graph

TAXA = [
    {"color_id": "0", "namespace": "taxon", "value": "taxon_1"},
    {"color_id": "1", "namespace": "taxon", "value": "taxon_2"},
    {"color_id": "2", "namespace": "taxon", "value": "taxon_3"},
]


def _colors(names: list[str]) -> list[dict[str, str]]:
    wanted = set(names)
    return [row for row in TAXA if row["value"] in wanted]


def strain_bubble() -> AssemblyGraph:
    """Shared source and sink, one branch per taxon. Ground truth: retain."""
    return build_graph(
        graph_id="strain_bubble",
        colors=_colors(["taxon_1", "taxon_2"]),
        nodes=[
            {"id": "S", "sequence": "ACGTACGT", "colors": [0, 1], "coverage": 120.0},
            {"id": "A", "sequence": "ATATCG", "colors": [0], "coverage": 30.0},
            {"id": "B", "sequence": "CGCGTA", "colors": [1], "coverage": 90.0},
            {"id": "T", "sequence": "GGCCTTAA", "colors": [0, 1], "coverage": 120.0},
        ],
        links=[
            {"id": "eSA", "source": "S", "target": "A", "colors": [0], "coverage": 30.0},
            {"id": "eAT", "source": "A", "target": "T", "colors": [0], "coverage": 30.0},
            {"id": "eSB", "source": "S", "target": "B", "colors": [1], "coverage": 90.0},
            {"id": "eBT", "source": "B", "target": "T", "colors": [1], "coverage": 90.0},
        ],
    )


def error_bubble() -> AssemblyGraph:
    """Same colour on the true branch and the error branch. Ground truth: pop."""
    return build_graph(
        graph_id="error_bubble",
        colors=_colors(["taxon_2"]),
        nodes=[
            {"id": "S", "sequence": "ACGTACGT", "colors": [1], "coverage": 91.0},
            {"id": "A", "sequence": "ATATCGAT", "colors": [1], "coverage": 90.0},
            {"id": "E", "sequence": "CGCGTACG", "colors": [1], "coverage": 1.0},
            {"id": "T", "sequence": "GGCCTTAA", "colors": [1], "coverage": 91.0},
        ],
        links=[
            {"id": "eSA", "source": "S", "target": "A", "colors": [1], "coverage": 90.0},
            {"id": "eAT", "source": "A", "target": "T", "colors": [1], "coverage": 90.0},
            {"id": "eSE", "source": "S", "target": "E", "colors": [1], "coverage": 1.0},
            {"id": "eET", "source": "E", "target": "T", "colors": [1], "coverage": 1.0},
        ],
    )


def nested_bubbles() -> AssemblyGraph:
    """One outer simple bubble and one inner bubble on a side path."""
    return build_graph(
        graph_id="nested_bubbles",
        colors=_colors(["taxon_1", "taxon_2"]),
        nodes=[
            {"id": "S", "sequence": "AAAAAA", "colors": [0, 1], "coverage": 40.0},
            {"id": "A", "sequence": "CCCCCC", "colors": [0], "coverage": 20.0},
            {"id": "U", "sequence": "GGGGGG", "colors": [1], "coverage": 20.0},
            {"id": "P", "sequence": "TTTTTT", "colors": [1], "coverage": 10.0},
            {"id": "Q", "sequence": "ACACAC", "colors": [1], "coverage": 10.0},
            {"id": "V", "sequence": "GTGTGT", "colors": [1], "coverage": 20.0},
            {"id": "B", "sequence": "CACACA", "colors": [1], "coverage": 20.0},
            {"id": "T", "sequence": "ATATAT", "colors": [0, 1], "coverage": 40.0},
        ],
        links=[
            {"id": "eSA", "source": "S", "target": "A", "colors": [0], "coverage": 20.0},
            {"id": "eAT", "source": "A", "target": "T", "colors": [0], "coverage": 20.0},
            {"id": "eSB", "source": "S", "target": "B", "colors": [1], "coverage": 20.0},
            {"id": "eBT", "source": "B", "target": "T", "colors": [1], "coverage": 20.0},
            {"id": "eTU", "source": "T", "target": "U", "colors": [1], "coverage": 20.0},
            {"id": "eUP", "source": "U", "target": "P", "colors": [1], "coverage": 10.0},
            {"id": "eUQ", "source": "U", "target": "Q", "colors": [1], "coverage": 10.0},
            {"id": "ePV", "source": "P", "target": "V", "colors": [1], "coverage": 10.0},
            {"id": "eQV", "source": "Q", "target": "V", "colors": [1], "coverage": 10.0},
            {"id": "eVT", "source": "V", "target": "T", "colors": [1], "coverage": 20.0},
        ],
    )


def shared_duplicate_node() -> AssemblyGraph:
    """One shared sequence instance sitting on two coloured paths."""
    return build_graph(
        graph_id="shared_node",
        colors=_colors(["taxon_1", "taxon_2"]),
        nodes=[
            {"id": "S", "sequence": "ACGTACGT", "colors": [0, 1], "coverage": 100.0},
            {"id": "X", "sequence": "ATGCCATG", "colors": [0, 1], "coverage": 100.0},
            {"id": "T", "sequence": "GGCCTTAA", "colors": [0, 1], "coverage": 100.0},
        ],
        links=[
            {"id": "eSX", "source": "S", "target": "X", "colors": [0, 1], "coverage": 100.0},
            {"id": "eXT", "source": "X", "target": "T", "colors": [0, 1], "coverage": 100.0},
        ],
    )
