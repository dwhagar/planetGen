"""MAP.151, MAP.205: the Galaxy Map's data window on the shell grid -- the cells of one level whose centre lies in a
ball, found in closed form, and the window cut from the parent stages."""
import math

import pytest

from planetgen.galaxy import spherical
from planetgen.galaxy.drill import DRILL_LEVELS, DrillBlock, drill_band_count, drill_parent, drill_wedge_count
from planetgen.galaxy.window import (
    DEFAULT_WINDOW, WINDOW_WIDTHS, cell_position_pc, cells_in_ball, check_width, crop_stages, parse_window_key,
    window_cells, window_parents,
)

EDGE = 4.0


def _brute_force(m, centre, reach, shells):
    """Every cell of level `m` in `shells` (sector shells for m=1, block shells otherwise) whose centre lies within
    `reach` of `centre`, by walking them all."""
    found = []
    for shell in shells:
        if m == 1:
            cells = [DrillBlock(1, *c) for c in spherical.shell_cells(shell)]
        else:
            rows = drill_band_count(m, shell)
            cells = [DrillBlock(m, shell, b, w) for b in range(rows // 2, -(rows // 2) - 1, -1)
                     for w in range(drill_wedge_count(m, shell, b))]
        for cell in cells:
            if math.dist(cell_position_pc(cell, EDGE), centre) <= reach + 1e-9:
                found.append(cell)
    return found


@pytest.mark.parametrize("m, shell, band, wedge, g", [
    (1, 12, 0, 3, 9),     # sectors on the plane
    (1, 12, 4, 20, 5),    # sectors off the plane
    (1, 9, 6, 0, 7),      # near a pole: the caps are in reach
    (1, 3, 0, 0, 9),      # the core's neighbourhood
    (3, 20, 2, 5, 9),
    (3, 4, 0, 0, 11),     # block shell 1, around the origin's blocks
    (9, 8, -2, 2, 9),
    (27, 3, 1, 0, 7),
])
def test_the_closed_form_ball_equals_a_walk_of_every_cell(m, shell, band, wedge, g):
    centre = DrillBlock(m, shell, band, wedge)
    reach = g * m * EDGE / 2
    point = cell_position_pc(centre, EDGE)
    radius = math.sqrt(sum(v * v for v in point))
    unit = 1 if m == 1 else m
    shells = range(max(0, int((radius - reach) / EDGE / unit) - 2), int((radius + reach) / EDGE / unit) + 4)
    expected = _brute_force(m, point, reach, shells)
    got = cells_in_ball(m, point, reach, EDGE)
    assert sorted(got) == sorted(expected)
    assert got == sorted(got, key=lambda c: (c.shell, -c.band, c.wedge)), "in shell, band and wedge order"


@pytest.mark.parametrize("m", DRILL_LEVELS)
def test_a_window_holds_its_centre_and_is_about_a_ball_of_cells(m):
    shell = 40 if m < 243 else 3
    if m == 1:
        centre = DrillBlock(1, shell, 5, 7)
    else:
        centre = DrillBlock(m, shell, 3, 2)
        assert 0 <= centre.wedge < drill_wedge_count(m, centre.shell, centre.band)
    cells = window_cells(centre, EDGE, DEFAULT_WINDOW)
    assert centre in cells and len(set(cells)) == len(cells)
    # A ball of radius g/2 cells holds about (4/3) pi (g/2)^3 = 382 cells of unit volume.
    assert 200 <= len(cells) <= 700, len(cells)


def test_the_window_is_centred_on_the_centre_cell():
    centre = DrillBlock(27, 20, 10, 3)
    cells = window_cells(centre, EDGE, 9)
    mean = [sum(c) / len(cells) for c in zip(*(cell_position_pc(c, EDGE) for c in cells))]
    assert math.dist(mean, cell_position_pc(centre, EDGE)) < 0.5 * 27 * EDGE


def test_widths_and_keys_are_checked():
    assert WINDOW_WIDTHS == (3, 5, 7, 9, 11) and DEFAULT_WINDOW == 9
    assert check_width("5") == 5
    for bad in (4, 13, "x", None):
        with pytest.raises(ValueError):
            check_width(bad)
    assert parse_window_key("81.12.-3.40") == DrillBlock(81, 12, -3, 40)
    assert parse_window_key("1.12.2.7") == DrillBlock(1, 12, 2, 7)
    for bad in ("", "81.12.-3", "729.0.0.0", "81.-1.0.0", "81.12.999.0", "1.12.99.0", "1.12.0.9999"):
        with pytest.raises(ValueError):
            parse_window_key(bad)


def test_a_window_of_sectors_is_cut_from_the_stages_of_its_level_three_parents():
    centre = DrillBlock(1, 30, 3, 11)
    cells = window_cells(centre, EDGE, 5)
    parents = window_parents(cells)
    assert None not in parents and all(p.m == 3 for p in parents)
    assert len(parents) == len(set(parents))
    stages = {}
    for parent in parents:
        stages[f"{parent.m}.{parent.shell}.{parent.band}.{parent.wedge}"] = {
            "children": [],
            "sectors": [{"shell": s.shell, "band": s.band, "wedge": s.wedge, "name": "x"} for c, in [(0,)]
                        for s in cells if drill_parent(s) == parent],
        }
    window = crop_stages(centre, 5, stages, EDGE)
    assert window["at"] == "1.30.3.11" and window["child_m"] == 1 and window["g"] == 5
    assert len(window["cells"]) == len(cells)
    assert {(s["shell"], s["band"], s["wedge"]) for s in window["sectors"]} == {(c.shell, c.band, c.wedge) for c in cells}
    assert window["sectors"] == sorted(window["sectors"], key=lambda s: (s["shell"], -s["band"], s["wedge"]))
