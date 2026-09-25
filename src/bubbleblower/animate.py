"""Animate iterative graph edits on a pinned Fruchterman–Reingold layout.

The first frame is placed with MetaMetro's Fruchterman–Reingold layout.
Later frames keep that drawing. A split or a merge moves only the new nodes
and their graph neighbours; every other node stays where it was.
"""

from __future__ import annotations

import math
import shutil
from pathlib import Path

import numpy as np

from bubbleblower.detect import detect_bubbles
from bubbleblower.edits import Edit, split_instance
from bubbleblower.graph import AssemblyGraph

# RColorBrewer Set1. Other is last and grey; Unclassified is near-black.
_SET1 = (
    "#E41A1C",
    "#377EB8",
    "#4DAF4A",
    "#984EA3",
    "#FF7F00",
    "#FFFF33",
    "#A65628",
    "#F781BF",
    "#999999",
)
_OTHER = "#CCCCCC"
_UNCLASSIFIED = "#333333"


def fork_resolution_states(
    graph: AssemblyGraph,
    *,
    max_edits: int | None = None,
) -> tuple[list[AssemblyGraph], list[Edit]]:
    """Split until every node has at most one incoming and one outgoing link.

    One node is split per iteration. A node with several outgoing links is
    split first: each outgoing link is its own group, and incoming links
    stay on the first group so they are not copied. A node that only has
    several incoming links is split the same way, and its single outgoing
    link stays on the first group. The loop stops when no such node remains,
    or when ``max_edits`` splits have been accepted. Coverage must already
    be present. This function does not fill it in.
    """
    _require_coverage(graph)
    current = graph.copy()
    frames = [current.copy()]
    edits: list[Edit] = []
    skipped: set[str] = set()
    limit = max_edits if max_edits is not None else max(len(current.cdbg.unitigs) * 4, 1)
    for _ in range(limit):
        outgoing, incoming = _incident_ids(current)
        forks = []
        for node_id in outgoing:
            if node_id in skipped:
                continue
            out_count = len(outgoing[node_id])
            in_count = len(incoming[node_id])
            if out_count >= 2 or in_count >= 2:
                forks.append((-max(out_count, in_count), -out_count, node_id))
        if not forks:
            break
        forks.sort()
        node_id = forks[0][2]
        out_links = outgoing[node_id]
        in_links = [link_id for link_id in incoming[node_id] if link_id not in out_links]
        if len(out_links) >= 2:
            groups = [[link_id] for link_id in out_links]
            if in_links:
                groups[0] = groups[0] + in_links
        else:
            groups = [[link_id] for link_id in in_links]
            if out_links:
                groups[0] = groups[0] + out_links
        try:
            updated, edit = split_instance(current, node_id, groups)
        except (KeyError, ValueError):
            skipped.add(node_id)
            continue
        current = updated
        edits.append(edit)
        frames.append(current.copy())
    return frames, edits


def colour_partition_states(
    graph: AssemblyGraph,
    *,
    namespace: str,
    max_edits: int = 12,
) -> tuple[list[AssemblyGraph], list[Edit]]:
    """Split one bubble source per iteration when branch colours differ.

    Each branch contributes its first link. Those links are the split groups.
    A bubble whose branches carry the same values in ``namespace`` is left
    in place. Coverage must already be present on every unitig and link;
    this function does not fill it in.

    The returned list starts with a copy of ``graph`` and gains one copy
    after each accepted split.
    """
    _require_coverage(graph)
    current = graph.copy()
    frames = [current.copy()]
    edits: list[Edit] = []
    skipped: set[str] = set()
    for _ in range(max_edits):
        chosen = None
        for bubble in detect_bubbles(current):
            if bubble.bubble_id in skipped:
                continue
            if _branch_labels(current, bubble, namespace) is None:
                skipped.add(bubble.bubble_id)
                continue
            chosen = bubble
            break
        if chosen is None:
            break
        groups = [[branch.links[0]] for branch in chosen.branches if branch.links]
        if len(groups) < 2:
            skipped.add(chosen.bubble_id)
            continue
        try:
            updated, edit = split_instance(current, chosen.source, groups)
        except (KeyError, ValueError):
            skipped.add(chosen.bubble_id)
            continue
        edit.bubble_id = chosen.bubble_id
        current = updated
        edits.append(edit)
        frames.append(current.copy())
    return frames, edits


def layout_states(
    frames: list[AssemblyGraph],
    edits: list[Edit],
    *,
    seed: int = 0,
) -> list[dict[str, tuple[float, float]]]:
    """Place each state. Distant nodes keep the coordinates of the first layout.

    ``frames`` has one more graph than ``edits``: the graph before any edit,
    then the graph after each edit.
    """
    if not frames:
        raise ValueError("layout needs at least one graph")
    if len(edits) != len(frames) - 1:
        raise ValueError("each edit must sit between two graphs")
    positions = _initial_layout(frames[0], seed=seed)
    placed = [dict(positions)]
    for _before, after, edit in zip(frames[:-1], frames[1:], edits, strict=True):
        positions = _place_edit(_before, after, edit, positions)
        placed.append(dict(positions))
    return placed


def animate_states(
    frames: list[AssemblyGraph],
    edits: list[Edit],
    path: str | Path,
    *,
    namespace: str,
    seconds: float = 5.0,
    seed: int = 0,
    hold_seconds: float = 0.0,
) -> Path:
    """Write a GIF or MP4 of about ``seconds``, then hold the last frame.

    ``.gif`` uses Pillow. ``.mp4`` uses ffmpeg, an optional dependency.
    ``hold_seconds`` repeats the settled final layout so the end state stays
    on screen. It may be zero.

    Edge colour is one value of ``namespace`` (Set1). A link with several
    values of that namespace is Other. A link with none is Unclassified.
    Node colour is degree on a YlGnBu scale fixed from the first graph.
    """
    if seconds <= 0:
        raise ValueError("seconds must be positive")
    if hold_seconds < 0:
        raise ValueError("hold_seconds must be non-negative")
    output = Path(path)
    if output.suffix.lower() not in {".gif", ".mp4"}:
        raise ValueError("animation path must be a .gif or a .mp4 file")
    if output.suffix.lower() == ".mp4" and shutil.which("ffmpeg") is None:
        raise RuntimeError("ffmpeg is not on PATH; it is an optional dependency")
    if not frames:
        raise ValueError("animation needs at least one graph")
    output.parent.mkdir(parents=True, exist_ok=True)
    layouts = layout_states(frames, edits, seed=seed)
    palette = _namespace_palette(frames[0], namespace)
    degree_limit = max(_degrees(frames[0]).values(), default=1.0)
    pictures = _tween(frames, edits, layouts, steps=_steps_per_edit(len(edits)))
    _write_animation(
        pictures,
        palette,
        namespace,
        degree_limit,
        output,
        seconds=seconds,
        hold_seconds=hold_seconds,
    )
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"animation was not written: {output}")
    return output


def _incident_ids(graph: AssemblyGraph) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """Outgoing and incoming link ids. A self-loop is outgoing only."""
    outgoing = {unitig.unitig_id: [] for unitig in graph.cdbg.unitigs}
    incoming = {unitig.unitig_id: [] for unitig in graph.cdbg.unitigs}
    for link in graph.cdbg.links:
        outgoing[link.source].append(link.link_id)
        if link.target != link.source:
            incoming[link.target].append(link.link_id)
    for links in outgoing.values():
        links.sort()
    for links in incoming.values():
        links.sort()
    return outgoing, incoming


def _require_coverage(graph: AssemblyGraph) -> None:
    unitig_ids = {unitig.unitig_id for unitig in graph.cdbg.unitigs}
    link_ids = {link.link_id for link in graph.cdbg.links}
    if set(graph.node_coverage) != unitig_ids or set(graph.link_coverage) != link_ids:
        raise ValueError("colour partition needs coverage on every unitig and every link")


def _branch_labels(graph: AssemblyGraph, bubble, namespace: str) -> set[frozenset[str]] | None:
    """Return one label set per branch, or None when they are not all different."""
    labels = []
    for branch in bubble.branches:
        if not branch.links:
            return None
        link = graph.link(branch.links[0])
        labels.append(frozenset(_namespace_values(graph, link.color_ids, namespace)))
    if len(labels) < 2 or len(set(labels)) < 2:
        return None
    return set(labels)


def _namespace_values(graph: AssemblyGraph, color_ids: list[int], namespace: str) -> list[str]:
    table = {
        int(row["color_id"]): row["value"]
        for row in graph.cdbg.colors or []
        if row.get("namespace") == namespace
    }
    if not table and graph.cdbg.colors:
        known = {row.get("namespace") for row in graph.cdbg.colors}
        if namespace not in known:
            raise ValueError(f"colour namespace {namespace!r} is not on the graph")
    return sorted(table[color_id] for color_id in color_ids if color_id in table)


def _edges(graph: AssemblyGraph) -> list[tuple[str, str]]:
    return [(link.source, link.target) for link in graph.cdbg.links]


def _ids(graph: AssemblyGraph) -> list[str]:
    return [unitig.unitig_id for unitig in graph.cdbg.unitigs]


def _initial_layout(graph: AssemblyGraph, *, seed: int) -> dict[str, tuple[float, float]]:
    from metametro.viz.cfa_colouring import spring_positions

    class _Edges:
        def __init__(self, node_ids: list[str], edges: list[tuple[str, str]]) -> None:
            self._ids = node_ids
            self.edges = [{"source": source, "target": target} for source, target in edges]

        def node_ids(self) -> list[str]:
            return list(self._ids)

    return {
        node_id: (float(x), float(y))
        for node_id, (x, y) in spring_positions(
            _Edges(_ids(graph), _edges(graph)),
            seed=seed,
            iterations=40,
            spread=1.0,
        ).items()
    }


def _place_edit(
    _before: AssemblyGraph,
    after: AssemblyGraph,
    edit: Edit,
    positions: dict[str, tuple[float, float]],
) -> dict[str, tuple[float, float]]:
    """Start new nodes on their parents, then relax those nodes and their neighbours."""
    after_ids = set(_ids(after))
    placed = {node_id: positions[node_id] for node_id in after_ids if node_id in positions}
    parents = [node_id for node_id in edit.source_ids if node_id in positions]
    origin = _centroid(parents, positions) if parents else (0.0, 0.0)
    born = [node_id for node_id in after_ids if node_id not in placed]
    born_set = set(born)
    for index, node_id in enumerate(born):
        placed[node_id] = _birth_offset(after, node_id, origin, placed, born_set, index, len(born))
    free = set(born)
    free.update(_neighbours(after, set(edit.source_ids) | set(edit.target_ids) | set(born)))
    free &= set(placed)
    return _relax(placed, _edges(after), free)


def _birth_offset(
    after: AssemblyGraph,
    node_id: str,
    origin: tuple[float, float],
    placed: dict[str, tuple[float, float]],
    born: set[str],
    index: int,
    n_born: int,
) -> tuple[float, float]:
    """Place a new node part-way toward the neighbour that is unique to it."""
    anchors = []
    for link in after.cdbg.links:
        other = link.target if link.source == node_id else link.source if link.target == node_id else None
        if other is None or other == node_id or other in born or other not in placed:
            continue
        anchors.append(placed[other])
    if not anchors:
        angle = 2.0 * math.pi * index / max(n_born, 1)
        return (origin[0] + 0.15 * math.cos(angle), origin[1] + 0.15 * math.sin(angle))
    ax = sum(point[0] for point in anchors) / len(anchors)
    ay = sum(point[1] for point in anchors) / len(anchors)
    dx = ax - origin[0]
    dy = ay - origin[1]
    length = math.hypot(dx, dy) or 1.0
    nudge = 0.05 * (index - (n_born - 1) / 2.0)
    return (
        origin[0] + 0.45 * dx + nudge * (-dy / length),
        origin[1] + 0.45 * dy + nudge * (dx / length),
    )


def _centroid(node_ids: list[str], positions: dict[str, tuple[float, float]]) -> tuple[float, float]:
    xs = [positions[node_id][0] for node_id in node_ids]
    ys = [positions[node_id][1] for node_id in node_ids]
    return sum(xs) / len(xs), sum(ys) / len(ys)


def _neighbours(graph: AssemblyGraph, seeds: set[str]) -> set[str]:
    found: set[str] = set()
    for link in graph.cdbg.links:
        if link.source in seeds:
            found.add(link.target)
        if link.target in seeds:
            found.add(link.source)
    return found


def _relax(
    positions: dict[str, tuple[float, float]],
    edges: list[tuple[str, str]],
    free: set[str],
    *,
    iterations: int = 20,
) -> dict[str, tuple[float, float]]:
    """Fruchterman–Reingold on ``free`` nodes. Every other coordinate is fixed."""
    order = list(positions)
    if not order or not free:
        return dict(positions)
    index = {node_id: slot for slot, node_id in enumerate(order)}
    pos = np.asarray([positions[node_id] for node_id in order], dtype=float)
    movable = np.asarray([index[node_id] for node_id in free if node_id in index], dtype=int)
    if movable.size == 0:
        return dict(positions)
    count = len(order)
    ideal = math.sqrt(1.0 / count)
    pairs = [(index[source], index[target]) for source, target in edges if source in index and target in index]
    temperature0 = 0.08
    for step in range(iterations):
        temperature = temperature0 * (1.0 - step / iterations)
        for slot in movable:
            delta = pos[slot] - pos
            distance = np.maximum(np.linalg.norm(delta, axis=1), 1e-6)
            repulsion = ((ideal * ideal) / distance)[:, None] * (delta / distance[:, None])
            repulsion[slot] = 0.0
            force = repulsion.sum(axis=0)
            for source, target in pairs:
                other = target if source == slot else source if target == slot else None
                if other is None:
                    continue
                edge = pos[slot] - pos[other]
                length = max(float(np.linalg.norm(edge)), 1e-6)
                force -= (length / ideal) * (edge / length)
            length = max(float(np.linalg.norm(force)), 1e-12)
            pos[slot] = pos[slot] + force / length * min(length, temperature)
    return {node_id: (float(pos[slot, 0]), float(pos[slot, 1])) for slot, node_id in enumerate(order)}


def _steps_per_edit(n_edits: int) -> int:
    if n_edits <= 0:
        return 1
    if n_edits > 40:
        return 2
    return max(4, min(12, 150 // n_edits))


def _tween(
    frames: list[AssemblyGraph],
    edits: list[Edit],
    layouts: list[dict[str, tuple[float, float]]],
    *,
    steps: int,
) -> list[tuple[AssemblyGraph, dict[str, tuple[float, float]]]]:
    if not edits:
        return [(frames[0], layouts[0])]
    pictures: list[tuple[AssemblyGraph, dict[str, tuple[float, float]]]] = [(frames[0], layouts[0])]
    for before, after, edit, start, end in zip(
        frames[:-1], frames[1:], edits, layouts[:-1], layouts[1:], strict=True
    ):
        origin = _birth_positions(before, after, edit, start)
        for step in range(steps):
            weight = 1.0 if steps == 1 else step / (steps - 1)
            pictures.append((after, _lerp(origin, end, weight)))
    return pictures


def _birth_positions(
    _before: AssemblyGraph,
    after: AssemblyGraph,
    edit: Edit,
    start: dict[str, tuple[float, float]],
) -> dict[str, tuple[float, float]]:
    """New nodes begin on the nodes they replace. Survivors begin where they were."""
    after_ids = _ids(after)
    parents = [node_id for node_id in edit.source_ids if node_id in start]
    origin = _centroid(parents, start) if parents else (0.0, 0.0)
    born = {
        node_id: origin
        for node_id in after_ids
        if node_id not in start
    }
    placed: dict[str, tuple[float, float]] = {}
    for node_id in after_ids:
        if node_id in born:
            placed[node_id] = born[node_id]
        elif node_id in start:
            placed[node_id] = start[node_id]
    return placed


def _lerp(
    start: dict[str, tuple[float, float]],
    end: dict[str, tuple[float, float]],
    weight: float,
) -> dict[str, tuple[float, float]]:
    mixed: dict[str, tuple[float, float]] = {}
    for node_id, destination in end.items():
        source = start.get(node_id, destination)
        mixed[node_id] = (
            source[0] + (destination[0] - source[0]) * weight,
            source[1] + (destination[1] - source[1]) * weight,
        )
    return mixed


def _namespace_palette(graph: AssemblyGraph, namespace: str) -> dict[str, str]:
    values = {
        row["value"]
        for row in graph.cdbg.colors or []
        if row.get("namespace") == namespace
    }
    if not values:
        raise ValueError(f"colour namespace {namespace!r} is not on the graph")
    ordered = sorted(values)
    colours = {label: _SET1[index % len(_SET1)] for index, label in enumerate(ordered)}
    colours["Other"] = _OTHER
    colours["Unclassified"] = _UNCLASSIFIED
    return colours


def _edge_label(graph: AssemblyGraph, color_ids: list[int], namespace: str) -> str:
    values = _namespace_values(graph, color_ids, namespace)
    if not values:
        return "Unclassified"
    if len(values) > 1:
        return "Other"
    return values[0]


def _degrees(graph: AssemblyGraph) -> dict[str, float]:
    degree = {unitig.unitig_id: 0.0 for unitig in graph.cdbg.unitigs}
    for link in graph.cdbg.links:
        degree[link.source] = degree.get(link.source, 0.0) + 1.0
        degree[link.target] = degree.get(link.target, 0.0) + 1.0
    return degree


def _write_animation(
    pictures: list[tuple[AssemblyGraph, dict[str, tuple[float, float]]]],
    palette: dict[str, str],
    namespace: str,
    degree_limit: float,
    path: Path,
    *,
    seconds: float,
    hold_seconds: float,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FFMpegWriter, PillowWriter
    from matplotlib.cm import ScalarMappable
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D

    figure, axis = plt.subplots(figsize=(7.2, 5.4))
    figure.subplots_adjust(left=0.08, right=0.76, bottom=0.1, top=0.96)
    norm = Normalize(vmin=0.0, vmax=max(degree_limit, 1.0))
    cmap = plt.get_cmap("YlGnBu")
    bar_axis = figure.add_axes((0.78, 0.55, 0.02, 0.36))
    bar = figure.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=bar_axis)
    bar.set_label("Degree (count)", fontsize=8)
    bar.ax.tick_params(labelsize=8)
    legend_labels = [label for label in palette if label not in {"Other", "Unclassified"}]
    legend_labels.append("Other")
    legend_labels.append("Unclassified")
    figure.legend(
        handles=[
            Line2D([0], [0], color=palette[label], linewidth=1.5, label=label) for label in legend_labels
        ],
        title=namespace,
        loc="center left",
        bbox_to_anchor=(0.78, 0.28),
        frameon=False,
        fontsize=7,
        title_fontsize=8,
    )
    fps = max(round(len(pictures) / seconds), 1)
    hold_frames = round(hold_seconds * fps)
    if path.suffix.lower() == ".mp4":
        writer = FFMpegWriter(fps=fps, codec="libx264")
    else:
        writer = PillowWriter(fps=fps)
    drawn = pictures + [pictures[-1]] * hold_frames
    with writer.saving(figure, str(path), dpi=120):
        for graph, positions in drawn:
            axis.clear()
            order = [node_id for node_id in _ids(graph) if node_id in positions]
            coordinates = np.asarray([positions[node_id] for node_id in order], dtype=float)
            index = {node_id: slot for slot, node_id in enumerate(order)}
            segments = []
            colours = []
            for link in graph.cdbg.links:
                if link.source not in index or link.target not in index:
                    continue
                segments.append([coordinates[index[link.source]], coordinates[index[link.target]]])
                colours.append(palette[_edge_label(graph, link.color_ids, namespace)])
            if segments:
                axis.add_collection(LineCollection(segments, colors=colours, linewidths=0.6, zorder=1))
            degree = _degrees(graph)
            axis.scatter(
                coordinates[:, 0],
                coordinates[:, 1],
                s=12 if len(order) > 200 else 28,
                c=[cmap(norm(degree.get(node_id, 0.0))) for node_id in order],
                linewidths=0,
                zorder=3,
            )
            axis.set_xlabel("Layout x", fontsize=9)
            axis.set_ylabel("Layout y", fontsize=9)
            axis.tick_params(labelsize=8)
            axis.set_aspect("equal")
            axis.autoscale_view()
            axis.margins(0.05)
            writer.grab_frame()
    plt.close(figure)
