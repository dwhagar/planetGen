"""MAP.151: the Galaxy Map's data window -- cells of one level cut from the parent stages."""
import math

import pytest

from planetgen.galaxy.drill import DRILL_LEVELS, DrillBlock, drill_chain_of, drill_wedge_count, format_drill_key
from planetgen.galaxy.window import (
    DEFAULT_WINDOW, WINDOW_WIDTHS, check_width, crop_stages, parse_window_key, window_cells, window_parents,
)

CENTRES = [
    DrillBlock(243, 0, 0, 0), DrillBlock(243, 7, 14, 0), DrillBlock(27, 0, 0, 0), DrillBlock(27, 40, 11, -3),
    DrillBlock(3, 5, 2, 1), DrillBlock(3, 300, 100, -2), DrillBlock(1, 0, 1, 0), DrillBlock(1, 25, 40, 2),
]


def _centre_of(cell):
    m = cell.m
    theta = (cell.wedge + 0.5) * 2 * math.pi / drill_wedge_count(m, cell.ring)
    radius = (cell.ring + 0.5) * m
    return radius * math.cos(theta), radius * math.sin(theta), cell.slab * m


def _brute(centre, g):
    """Every cell of the level whose centre is within g/2 cells across and (g-1)/2 slabs of the centre's."""
    m = centre.m
    cx, cy, _ = _centre_of(centre)
    half = g * m / 2.0
    found = set()
    top_ring = int(math.hypot(cx, cy) / m + g + 2)
    for ring in range(0, top_ring + 1):
        for wedge in range(drill_wedge_count(m, ring)):
            cell = DrillBlock(m, ring, wedge, centre.slab)
            x, y, _z = _centre_of(cell)
            if math.hypot(x - cx, y - cy) <= half + 1e-9:
                for slab in range(centre.slab - (g - 1) // 2, centre.slab + (g - 1) // 2 + 1):
                    found.add(DrillBlock(m, ring, wedge, slab))
    return found


@pytest.mark.parametrize("centre", CENTRES)
@pytest.mark.parametrize("g", WINDOW_WIDTHS)
def test_a_window_is_the_cells_within_reach_of_the_middle_one(centre, g):
    cells = window_cells(centre, g)
    assert len(cells) == len(set(cells))
    assert set(cells) == _brute(centre, g)
    assert centre in cells
    assert all(isinstance(c, DrillBlock) and c.m == centre.m for c in cells)


def test_a_window_is_about_g_cubed_away_from_the_axis():
    cells = window_cells(DrillBlock(27, 40, 11, 0), 9)
    assert 0.7 * 729 <= len(cells) <= 1.05 * 729
    assert {c.slab for c in cells} == set(range(-4, 5))


def test_a_window_over_the_axis_wraps_through_every_wedge():
    cells = window_cells(DrillBlock(27, 0, 0, 0), 5)
    assert {c.ring for c in cells} >= {0, 1, 2}
    assert all(0 <= c.wedge < drill_wedge_count(27, c.ring) for c in cells)


def test_a_window_touches_few_parents():
    for centre in CENTRES:
        parents = window_parents(window_cells(centre, 9))
        # nine cells across touch two parents an axis when the parent is nine cells wide, four when three
        assert 1 <= len(parents) <= (64 if centre.m == 1 else 18 if centre.ring >= 5 else 60)
        if centre.m > 1 and centre.ring >= 5:
            assert len(parents) <= 12
        if centre.m == 243:
            assert parents == [None]
        else:
            assert None not in parents and all(p.m == centre.m * (3 if centre.m == 1 else 9) for p in parents)


def test_keys_and_widths_are_checked():
    assert parse_window_key("27.3.5.-2") == DrillBlock(27, 3, 5, -2)
    assert parse_window_key("1.4.7.0") == DrillBlock(1, 4, 7, 0)
    for bad in ("", None, "27.3.5", "81.0.0.0", "27.3.999.0", "27.-1.0.0", "a.b.c.d"):
        with pytest.raises(ValueError):
            parse_window_key(bad)
    assert check_width("7") == 7 and check_width(DEFAULT_WINDOW) == DEFAULT_WINDOW
    for bad in (4, 13, "x", None, 0):
        with pytest.raises(ValueError):
            check_width(bad)


def _stage_of(parent, sectors):
    """A stage like `queryDb.galaxy_stage`'s, worked out from `sectors` (ring, layer, slot)."""
    level = 0 if parent is None else DRILL_LEVELS.index(parent.m) + 1
    children = {}
    for ring, layer, slot in sectors:
        chain = drill_chain_of(ring, layer, slot)
        if parent is not None and chain[level - 1] != parent:
            continue
        children.setdefault(chain[level], []).append((ring, layer, slot))
    cells = [{"ring": b.ring, "wedge": b.wedge, "slab": b.slab, "generated": len(v), "stats": {}}
             for b, v in children.items()]
    sector_rows = None
    if level == 3:
        sector_rows = [{"ring": r, "layer": ly, "slot": sl, "id": 100 * r + sl, "name": f"S{r}-{sl}"} for v in children.values()
                       for r, ly, sl in v]
    return {"child_m": DRILL_LEVELS[level], "children": cells, "sectors": sector_rows}


SECTORS = [(0, 0, 1), (0, 0, 2), (1, 0, 3), (2, 1, 5), (5, 0, 20), (30, -1, 100), (30, 0, 101), (60, 2, 3), (200, 0, 40)]


@pytest.mark.parametrize("centre", [DrillBlock(243, 0, 0, 0), DrillBlock(27, 0, 0, 0), DrillBlock(27, 1, 0, 0),
                                    DrillBlock(3, 0, 0, 0), DrillBlock(3, 5, 4, 0), DrillBlock(1, 0, 1, 0)])
def test_cropped_stages_equal_a_direct_count_of_the_window(centre):
    g = 7
    cells = window_cells(centre, g)
    stages = {(format_drill_key(p) if p else "galaxy"): _stage_of(p, SECTORS) for p in window_parents(cells)}
    window = crop_stages(centre, g, stages)
    want = {}
    for ring, layer, slot in SECTORS:
        cell = drill_chain_of(ring, layer, slot)[DRILL_LEVELS.index(centre.m)]
        if cell in set(cells):
            want[(cell.ring, cell.wedge, cell.slab)] = want.get((cell.ring, cell.wedge, cell.slab), 0) + 1
    assert {(c["ring"], c["wedge"], c["slab"]): c["generated"] for c in window["children"]} == want
    assert window["cells"] == [[c.ring, c.wedge, c.slab] for c in cells]
    assert window["at"] == format_drill_key(centre) and window["g"] == g and window["child_m"] == centre.m
    keys = [(c["slab"], c["ring"], c["wedge"]) for c in window["children"]]
    assert keys == sorted(keys)
    if centre.m == 1:
        assert {(s["ring"], s["layer"], s["slot"]) for s in window["sectors"]} == {
            (r, ly, sl) for r, ly, sl in SECTORS if DrillBlock(1, r, sl, ly) in set(cells)}
    else:
        assert window["sectors"] is None


# --- The page's copy agrees --------------------------------------------------

def test_the_page_lists_the_same_window_cells():
    from tests.test_galaxydrill import NODE, _run
    if NODE is None:
        pytest.skip("node is not installed")
    centres = CENTRES + [DrillBlock(243, 15, 3, 1), DrillBlock(27, 120, 300, 5), DrillBlock(3, 1000, 1500, -4),
                         DrillBlock(1, 1500, 200, 0), DrillBlock(243, 0, 1, 0), DrillBlock(1, 1, 2, 0)]
    got = _run(f"""
const centres = {[list(c) for c in centres]!r};
const out = [];
for (const [m, ring, wedge, slab] of centres) {{
  for (const g of {list(WINDOW_WIDTHS)!r}) {{
    out.push(P.windowCells({{m, ring, wedge, slab}}, g).map(b => [b.m, b.ring, b.wedge, b.slab]));
  }}
}}
out.push(P.windowCells({{m: 27, ring: 40, wedge: 11, slab: 0}}).length, P.DEFAULT_WINDOW);
console.log(JSON.stringify(out));
""")
    want = [[list(c) for c in window_cells(centre, g)] for centre in centres for g in WINDOW_WIDTHS]
    assert got[:-2] == want
    assert got[-2:] == [len(window_cells(DrillBlock(27, 40, 11, 0))), DEFAULT_WINDOW]
