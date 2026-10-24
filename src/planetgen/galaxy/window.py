"""
planetgen/galaxy/window.py

The Galaxy Map's data window (MAP.151, docs/design/interface/drilldown-region-sizes.md
sections 1 and 5.3; on shells, docs/design/galaxy/spherical-block-ladder.md
section 4): `g` cells across, centred on a chosen cell of one level of the
block ladder (`galaxy.drill`). The window is a ball: every cell of the
centre's level whose centre lies within `g * m / 2` edges of the centre cell's
centre. The data is whole cells, so a window of blocks is read from the stored
block totals and a window of sectors is built from the stage files the cache
already holds -- the stage of each level-3 block the window touches -- cropped
to the window's cells.

The cells are found in closed form per shell: the ball meets a sphere in a cap
of angular radius `psi`, so the rows or bands to look at are those within
`psi` of the centre's colatitude, and in each the wedges are those within a
longitude half-width. No walk of the ball, and the cost is a handful of
arithmetic steps per (shell, band) the ball reaches.
"""
import math
from bisect import bisect_left, bisect_right

from planetgen.galaxy import spherical
from planetgen.galaxy.drill import (
    STAGE_LEVELS, DrillBlock, block_radius, cell_colatitude, drill_band_count, drill_wedge_count, format_drill_key,
    stage_parent,
)

WINDOW_WIDTHS = (3, 5, 7, 9, 11)
"""tuple: Cells across a window, odd so the chosen cell sits in the middle
(docs/design/interface/drilldown-region-sizes.md section 5.5)."""

DEFAULT_WINDOW = 9
"""int: Nine across: at most a few hundred cells."""

TWO_PI = 2 * math.pi
_TOLERANCE_PC = 1e-9


def parse_window_key(key):
    """
    A cell key `m.shell.band.wedge` of any level of the ladder, a sector
    (`1.shell.band.wedge`) included, checked like `drill.parse_drill_key`.

    Raises:
        ValueError: On anything else.
    """
    parts = str(key or "").split(".")
    if len(parts) != 4:
        raise ValueError(f"not a block key: {key!r}")
    try:
        m, shell, band, wedge = (int(p) for p in parts)
    except ValueError:
        raise ValueError(f"not a block key: {key!r}") from None
    if m not in STAGE_LEVELS:
        raise ValueError(f"block size must be one of {STAGE_LEVELS}, got {m}")
    if shell < 0:
        raise ValueError(f"no such block: {key!r}")
    if m == 1:
        low, high = spherical.band_range(shell)
        if not low <= band <= high or not 0 <= wedge < spherical.wedge_count(shell, band):
            raise ValueError(f"no such block: {key!r}")
    elif abs(band) > drill_band_count(m, shell) // 2 or not 0 <= wedge < drill_wedge_count(m, shell, band):
        raise ValueError(f"no such block: {key!r}")
    return DrillBlock(m, shell, band, wedge)


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


def _block_colatitude(m, shell, band):
    rows = drill_band_count(m, shell)
    return ((rows - 1) // 2 - band + 0.5) * math.pi / rows


def cell_position_pc(cell, edge_pc):
    """The `(x, y, z)` of a cell's centre, parsecs: a sector's own centre, a
    block's at its block shell's centreline radius, its row's centre
    colatitude and its wedge's centre longitude."""
    m, shell, band, wedge = cell
    if m == 1:
        return spherical.sector_position_pc(shell, band, wedge, edge_pc)
    radius = block_radius(m, shell) * edge_pc
    colat = _block_colatitude(m, shell, band)
    lon = (wedge + 0.5) * TWO_PI / drill_wedge_count(m, shell, band)
    plane = radius * math.sin(colat)
    return (plane * math.cos(lon), plane * math.sin(lon), radius * math.cos(colat))


def _levels_shells(m, low_pc, high_pc, edge_pc):
    """The shell ids of level `m` (block shells, or sector shells for 1)
    whose centreline radius lies in `[low_pc, high_pc]`, with that radius."""
    if m == 1:
        first = max(0, math.floor(low_pc / edge_pc))
        last = math.ceil(high_pc / edge_pc)
        return [(k, k * edge_pc) for k in range(first, last + 1) if low_pc <= k * edge_pc <= high_pc]
    first = max(0, math.floor((low_pc / edge_pc - (m - 1) / 2) / m))
    last = math.ceil((high_pc / edge_pc - (m - 1) / 2) / m)
    out = []
    for s in range(first, last + 1):
        radius = block_radius(m, s) * edge_pc
        if low_pc <= radius <= high_pc:
            out.append((s, radius))
    return out


def _bands(m, shell, low_colat, high_colat):
    """`(band, centre colatitude)` of the rows (blocks) or bands (sectors) of
    one shell whose centre colatitude lies in `[low_colat, high_colat]`,
    north to south."""
    if m == 1:
        if shell == 0:
            return [(0, math.pi / 2)]
        low, high = spherical.band_range(shell)
        bands = list(range(high, low - 1, -1))
        centres = [cell_colatitude(shell, b) for b in bands]
    else:
        rows = drill_band_count(m, shell)
        bands = list(range(rows // 2, -(rows // 2) - 1, -1))
        centres = [_block_colatitude(m, shell, b) for b in bands]
    start = bisect_left(centres, low_colat)
    stop = bisect_right(centres, high_colat)
    return [(bands[i], centres[i]) for i in range(start, stop)]


def _wedge_candidates(count, lon0, half_width):
    """The wedges of a band of `count` whose centre longitude may lie within
    `half_width` of `lon0` (a superset; the caller tests the distance)."""
    if half_width >= math.pi or count == 1:
        return range(count)
    step = TWO_PI / count
    first = math.floor((lon0 - half_width) / step - 0.5)
    last = math.ceil((lon0 + half_width) / step - 0.5)
    return sorted({w % count for w in range(first, last + 1)})


def cells_in_ball(m, centre_pc, reach_pc, edge_pc):
    """
    The cells of level `m` (a block size, or 1 for sectors) whose centre lies
    within `reach_pc` of `centre_pc`, as `DrillBlock`s ordered by shell, then
    band north to south, then wedge.
    """
    x0, y0, z0 = centre_pc
    rho0 = math.sqrt(x0 * x0 + y0 * y0 + z0 * z0)
    theta0 = math.acos(max(-1.0, min(1.0, z0 / rho0))) if rho0 else 0.0
    lon0 = math.atan2(y0, x0)
    cells = []
    for shell, radius in _levels_shells(m, rho0 - reach_pc - _TOLERANCE_PC, rho0 + reach_pc + _TOLERANCE_PC, edge_pc):
        if rho0 == 0 or radius == 0:
            psi, cos_need = math.pi, -1.0
        else:
            cos_need = (radius * radius + rho0 * rho0 - reach_pc * reach_pc) / (2 * radius * rho0)
            psi = math.pi if cos_need <= -1 else (0.0 if cos_need >= 1 else math.acos(cos_need))
        # Block rows are wide: look one row's width beyond the cap for rows whose centre may still be inside.
        slack = math.pi / drill_band_count(m, shell) if m > 1 else math.pi / max(1, spherical.band_range(shell)[1] * 2 + 1)
        for band, colat in _bands(m, shell, theta0 - psi - slack, theta0 + psi + slack):
            count = drill_wedge_count(m, shell, band) if m > 1 else spherical.wedge_count(shell, band)
            denominator = math.sin(colat) * math.sin(theta0)
            if abs(denominator) < 1e-9 or radius == 0:
                wedges = range(count)
            else:
                cos_lon = (cos_need - math.cos(colat) * math.cos(theta0)) / denominator
                half = math.pi if cos_lon <= -1 else (0.0 if cos_lon >= 1 else math.acos(cos_lon))
                step = TWO_PI / count
                wedges = _wedge_candidates(count, lon0, half + step)
            for wedge in wedges:
                cell = DrillBlock(m, shell, band, wedge)
                x, y, z = cell_position_pc(cell, edge_pc)
                if math.sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2) <= reach_pc + _TOLERANCE_PC:
                    cells.append(cell)
    return cells


def window_cells(centre, edge_pc, g=DEFAULT_WINDOW):
    """
    The cells of the window, `[DrillBlock, ...]`: every cell of `centre`'s
    level whose centre lies within `g * m / 2` edges of `centre`'s. `centre`
    itself is always among them.
    """
    centre = DrillBlock(*centre)
    reach = g * centre.m * edge_pc / 2.0
    cells = cells_in_ball(centre.m, cell_position_pc(centre, edge_pc), reach, edge_pc)
    if centre not in cells:
        cells.append(centre)
        cells.sort(key=lambda c: (c.shell, -c.band, c.wedge))
    return cells


def window_parents(cells):
    """The stages a window is built from: `[None or block, ...]` -- `None`
    for the galaxy's own (the level-243 blocks), else the parent block --
    without repeats, in order of first use."""
    seen = {}
    for cell in cells:
        seen.setdefault(stage_parent(cell), None)
    return list(seen)


def crop_stages(centre, g, stages, edge_pc):
    """
    The window of `centre` and `g` cut from the stages of its parents.

    Args:
        centre (DrillBlock): The middle cell.
        g (int): Cells across.
        stages (dict): `{key: stage}` for each of `window_parents`, `key` the
            parent's `format_drill_key`, or `"galaxy"`; each a
            `queryDb.galaxy_stage` result.
        edge_pc (float): The sector edge in parsecs.

    Returns:
        dict: `at` (the centre's key), `g`, `child_m` (the cells' size),
            `cells` (every cell: `[[shell, band, wedge], ...]`), `children`
            (the cells that hold a generated sector, as a stage lists them)
            and `sectors` (a window of sectors lists the generated ones with
            their names; `None` for bigger cells).
    """
    cells = window_cells(centre, edge_pc, g)
    wanted = {(c.shell, c.band, c.wedge) for c in cells}
    children = []
    sectors = [] if centre.m == 1 else None
    for stage in stages.values():
        children += [c for c in stage.get("children") or () if (c["shell"], c["band"], c["wedge"]) in wanted]
        if sectors is not None:
            sectors += [s for s in stage.get("sectors") or () if (s["shell"], s["band"], s["wedge"]) in wanted]
    children.sort(key=lambda c: (c["shell"], -c["band"], c["wedge"]))
    if sectors is not None:
        sectors.sort(key=lambda s: (s["shell"], -s["band"], s["wedge"]))
    return {
        "at": format_drill_key(centre), "g": g, "child_m": centre.m,
        "cells": [[c.shell, c.band, c.wedge] for c in cells], "children": children, "sectors": sectors,
    }
