"""
planetgen/galaxy/window.py

The Galaxy Map's data window (MAP.151, docs/design/drilldown-region-sizes.md
sections 1 and 5.3): `g` cells across, centred as near as the lattice
allows on a chosen cell of one level of the block ladder (`galaxy.drill`).
The frame is a round region of the lattice (`galaxy.region.Region`, the
shape the radial fills use) and the data is whole cells, so every window
is built from the stage files the cache already holds -- the stage of each
parent block the window touches (at most 8 for 9 cells across) -- cropped to
the window's cells.

The lattice is the sector grid scaled by the cell size, so a window needs no
edge length: a block `m` sectors a side has its ring centreline at
`(ring + 1/2) * m` edges, its wedge centre at `(wedge + 1/2)` of
`drill_wedge_count(m, ring)` turns, and its slab midplane at `slab * m` edges.
"""
import math

from planetgen.galaxy.drill import DRILL_LEVELS, DrillBlock, drill_parent, drill_wedge_count, format_drill_key
from planetgen.galaxy.region import Region

WINDOW_WIDTHS = (3, 5, 7, 9, 11)
"""tuple: Cells across a window, odd so the chosen cell sits in the middle
(docs/design/drilldown-region-sizes.md section 5.5)."""

DEFAULT_WINDOW = 9
"""int: Nine across reproduces today's stage of 729 children at most."""


def parse_window_key(key):
    """
    A cell key `m.ring.wedge.slab` of any level of the ladder, a sector
    (`1.ring.slot.layer`) included, checked like `drill.parse_drill_key`.

    Raises:
        ValueError: On anything else.
    """
    parts = str(key or "").split(".")
    if len(parts) != 4:
        raise ValueError(f"not a block key: {key!r}")
    try:
        m, ring, wedge, slab = (int(p) for p in parts)
    except ValueError:
        raise ValueError(f"not a block key: {key!r}") from None
    if m not in DRILL_LEVELS:
        raise ValueError(f"block size must be one of {DRILL_LEVELS}, got {m}")
    if ring < 0 or not 0 <= wedge < drill_wedge_count(m, ring):
        raise ValueError(f"no such block: {key!r}")
    return DrillBlock(m, ring, wedge, slab)


def check_width(g):
    """`g` as a whole number of cells across, one of `WINDOW_WIDTHS`.

    Raises:
        ValueError: On any other.
    """
    try:
        width = int(g)
    except (TypeError, ValueError):
        raise ValueError(f"cells across must be one of {WINDOW_WIDTHS}, got {g!r}") from None
    if width not in WINDOW_WIDTHS:
        raise ValueError(f"cells across must be one of {WINDOW_WIDTHS}, got {g!r}")
    return width


def window_region(centre, g=DEFAULT_WINDOW):
    """The window's frame: a `Region` of the lattice of `centre`'s level,
    in units of one sector edge."""
    m, ring, wedge, slab = centre
    wedges = drill_wedge_count(m, ring)
    theta = (wedge + 0.5) * 2 * math.pi / wedges
    radius = (ring + 0.5) * m
    position = (radius * math.cos(theta), radius * math.sin(theta), float(slab * m))
    return Region(position, g * m / 2.0, (g - 1) // 2, float(m), slot_count=lambda r: drill_wedge_count(m, r))


def window_cells(centre, g=DEFAULT_WINDOW):
    """
    The cells of the window, `[DrillBlock, ...]`, a slab at a time, rings out
    from the axis: every cell of `centre`'s level whose centre lies within
    `g/2` cells of `centre`'s across the plane and `(g-1)/2` slabs of its
    layer. `centre` itself is always among them.
    """
    m = centre.m
    return [DrillBlock(m, ring, wedge, slab) for ring, slab, wedge in window_region(centre, g).addresses()]


def window_parents(cells):
    """The stages a window is built from: `[None or block, ...]` -- `None`
    for the galaxy's own (the level-243 blocks), else the parent block --
    without repeats, in order of first use."""
    seen = {}
    for cell in cells:
        seen.setdefault(drill_parent(cell), None)
    return list(seen)


def crop_stages(centre, g, stages):
    """
    The window of `centre` and `g` cut from the stages of its parents.

    Args:
        centre (DrillBlock): The middle cell.
        g (int): Cells across.
        stages (dict): `{key: stage}` for each of `window_parents`, `key` the
            parent's `format_drill_key`, or `"galaxy"`; each a
            `queryDb.galaxy_stage` result.

    Returns:
        dict: `at` (the centre's key), `g`, `child_m` (the cells' size),
            `cells` (every cell: `[[ring, wedge, slab], ...]`), `children`
            (the cells that hold a generated sector, as a stage lists
            them) and `sectors` (a window of sectors lists the generated
            ones with their names; `None` for bigger cells).
    """
    cells = window_cells(centre, g)
    wanted = {(c.ring, c.wedge, c.slab) for c in cells}
    children = []
    sectors = [] if centre.m == 1 else None
    for stage in stages.values():
        children += [c for c in stage.get("children") or () if (c["ring"], c["wedge"], c["slab"]) in wanted]
        if sectors is not None:
            sectors += [s for s in stage.get("sectors") or () if (s["ring"], s["slot"], s["layer"]) in wanted]
    children.sort(key=lambda c: (c["slab"], c["ring"], c["wedge"]))
    if sectors is not None:
        sectors.sort(key=lambda s: (s["layer"], s["ring"], s["slot"]))
    return {
        "at": format_drill_key(centre), "g": g, "child_m": centre.m,
        "cells": [[c.ring, c.wedge, c.slab] for c in cells], "children": children, "sectors": sectors,
    }
