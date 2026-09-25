"""Small coloured graphs with known bubble labels."""

from __future__ import annotations

from bubbleblower.graph import AssemblyGraph, adopt_cfa


def strain_bubble() -> AssemblyGraph:
    """Shared source and sink, one branch per taxon. Ground truth: retain."""
    from metametro.bench.data.universal.bubbles import strain_bubble as cfa

    return adopt_cfa(cfa())


def error_bubble() -> AssemblyGraph:
    """Same colour on the true branch and the error branch. Ground truth: pop."""
    from metametro.bench.data.universal.bubbles import error_bubble as cfa

    return adopt_cfa(cfa())


def nested_bubbles() -> AssemblyGraph:
    """One outer simple bubble and one inner bubble on a side path."""
    from metametro.bench.data.universal.bubbles import nested_bubbles as cfa

    return adopt_cfa(cfa())


mock_bubble_strain = strain_bubble
mock_bubble_error = error_bubble
mock_bubble_nested = nested_bubbles


def shared_duplicate_node() -> AssemblyGraph:
    """One shared sequence instance sitting on two coloured paths."""
    from metametro.bench.data.universal.bubbles import shared_node as cfa

    return adopt_cfa(cfa())


mock_bubble_shared = shared_duplicate_node
mock_bubble_duplicate = shared_duplicate_node
