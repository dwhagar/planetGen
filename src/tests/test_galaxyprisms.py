"""
Tests for `html/static/galaxyprisms.js`, the Galaxy Map's grid and prisms, run under node (skipped where node isn't
installed). It is the browser's twin of `planetgen.galaxy.spherical` (the shell grid), `drill` (the block ladder),
`window` (the data window) and `skeleton` (the outline), so each is checked against the Python on samples that
reach the caps, the core, the hand-cut shell 1 and the outer shells: grid layouts, addresses, designations, block
chains and bounds, windows, the outline and the sectors a block can hold (MAP.205, MAP.206). Its density must match
`planetgen.galaxy.density.relative_density` exactly and its geometry must be closed and wound outward. Also
`html/static/galaxyblocks.js`, which packs a stage's cells for the GPU: every cell once, arrays that agree with
each other, and colors and opacities that follow what a cell holds.
"""

import json
import math
import os
import random
import shutil
import subprocess

import pytest

from planetgen.galaxy import spherical
from planetgen.galaxy.density import build_galaxy_shape, relative_density, shape_with_terms
from planetgen.galaxy.drill import (
    DRILL_LEVELS, DrillBlock, block_addresses, drill_band_count, drill_block_bounds, drill_chain_of, drill_parent,
    drill_wedge_count, format_drill_key,
)
from planetgen.galaxy.skeleton import GalaxyBounds, build_shell_extents, expected_system_count_at_density_1
from planetgen.galaxy.window import window_cells
from planetgen.physics.units import pc_to_ly

NODE = shutil.which("node")
STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "html", "static")
MODULE = os.path.join(STATIC, "galaxyprisms.js")
BLOCKS_MODULE = os.path.join(STATIC, "galaxyblocks.js")

pytestmark = pytest.mark.skipif(NODE is None, reason="node is not installed")

SHAPE = build_galaxy_shape(2800.0, 350.0, 200.0, 1.0, 2, math.radians(15.0), 0.4)
EDGE_PC = 4.0


def _run(script, shape_used=None):
    """Runs `script` as an ES module with the prisms module imported as `P`, the block scene module as `B` and
    the test shape as `shape`; returns its JSON output. The source goes in on standard input: a big list of
    cases would overrun the command line."""
    source = (
        f"import * as P from {json.dumps('file://' + os.path.abspath(MODULE))};\n"
        f"import * as B from {json.dumps('file://' + os.path.abspath(BLOCKS_MODULE))};\n"
        f"const shape = {json.dumps(shape_with_terms(shape_used or SHAPE))};\n"
        f"{script}\n"
    )
    result = subprocess.run([NODE, "--input-type=module"], input=source, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-3000:]
    return json.loads(result.stdout)


def _shape_with_threshold():
    """The test shape with the galaxy's own sector threshold, as the page embeds it
    (`web/maps/galaxymap3d._density_shape`)."""
    return {**shape_with_terms(SHAPE), "sector_min_density": 1.0 / expected_system_count_at_density_1(pc_to_ly(EDGE_PC))}


def test_density_matches_the_python_model():
    points = [(8000.0, 100.0, 50.0), (0.0, 0.0, 0.0), (300.0, -200.0, 10.0), (-12000.0, 4000.0, -600.0),
              (5.0, 5.0, 900.0), (20000.0, -3.0, 0.0)]
    got = _run(f"console.log(JSON.stringify({json.dumps(points)}.map(p => P.relativeDensity(...p, shape))));")
    for point, value in zip(points, got):
        assert value == pytest.approx(relative_density(point, SHAPE), rel=1e-12)


def _sample_sectors(count, seed):
    """Sectors across the grid: the core, the hand-cut shell 1, small shells, the outer ones, each at its caps,
    its equator and a few other bands."""
    rng = random.Random(seed)
    found = []
    for _ in range(count):
        shell = rng.choice([0, 1, 2, 3, rng.randint(0, 60), rng.randint(60, 3764), 3764, 3765])
        low, high = spherical.band_range(shell)
        band = rng.choice([low, high, 0, rng.randint(low, high), rng.randint(low, high)])
        found.append((shell, band, rng.randrange(spherical.wedge_count(shell, band))))
    return found


def test_the_shell_grid_matches_python():
    """Wedge counts, band ranges, caps and every cell's bounds, for shells from the core out to the rim."""
    sectors = _sample_sectors(1500, 3)
    got = _run(f"""
const cells = {json.dumps(sectors)};
console.log(JSON.stringify(cells.map(([s, b, w]) => {{
  const x = P.sectorBounds(s, b, w, {EDGE_PC});
  return [P.wedgeCount(s, b), P.bandRange(s), P.isCap(s, b), [x.r0, x.r1, x.c0, x.c1, x.t0, x.t1],
    P.sectorPositionPc(s, b, w, {EDGE_PC}), P.bandCentrePc(s, b, {EDGE_PC})];
}})));
""")
    for (shell, band, wedge), (count, band_range, cap, bounds, position, centre) in zip(sectors, got):
        assert count == spherical.wedge_count(shell, band), (shell, band)
        assert band_range == list(spherical.band_range(shell))
        assert cap == spherical.is_cap(shell, band)
        inner, outer = spherical.shell_bounds_pc(shell, EDGE_PC)
        north, south = spherical.band_colatitude_bounds(shell, band)
        start, end = spherical.wedge_longitude_bounds(shell, band, wedge)
        assert bounds == pytest.approx([inner, outer, north, south, start, end], abs=1e-9, rel=1e-12), (shell, band, wedge)
        assert position == pytest.approx(spherical.sector_position_pc(shell, band, wedge, EDGE_PC), abs=1e-6, rel=1e-12)
        assert centre == pytest.approx(spherical.band_center_pc(shell, band, EDGE_PC), abs=1e-6, rel=1e-12)


def test_addresses_at_points_and_designations_match_python():
    rng = random.Random(5)
    points = []
    for _ in range(3000):
        radius = rng.choice([rng.uniform(0, 30), rng.uniform(0, 15000)])
        theta = math.acos(rng.uniform(-1, 1))
        phi = rng.uniform(-math.pi, math.pi)
        points.append((radius * math.sin(theta) * math.cos(phi), radius * math.sin(theta) * math.sin(phi),
                       radius * math.cos(theta)))
    got = _run(f"""
const pts = {json.dumps(points)};
console.log(JSON.stringify(pts.map(p => {{ const a = P.sectorAddressAt(...p, {EDGE_PC}); return [a.shell, a.band, a.wedge]; }})));
""")
    assert got == [list(spherical.sector_address_at(p, EDGE_PC)) for p in points]

    # The band field holds |band| < 2048, so the designation tests keep to the shells that fit it.
    sectors = [s for s in _sample_sectors(1500, 7) if abs(s[1]) < 2040 and s[0] < 4000]
    codes = [spherical.provisional_sector_designation(*s) for s in sectors]
    made = _run(f"""
const cells = {json.dumps(sectors)};
const codes = cells.map(([s, b, w]) => P.sectorDesignation(s, b, w));
const back = codes.map(c => {{ const a = P.parseSectorDesignation(c); return a ? [a.shell, a.band, a.wedge] : null; }});
console.log(JSON.stringify({{codes, back}}));
""")
    assert made["codes"] == codes
    assert made["back"] == [list(s) for s in sectors]
    for code in ("", "G1", "1" * 11, "0"):
        assert _run(f"console.log(JSON.stringify(P.parseSectorDesignation({json.dumps(code)})));") is None


def test_block_chains_and_bounds_match_python():
    sectors = _sample_sectors(600, 11)
    chains = [[(b.m, b.shell, b.band, b.wedge) for b in drill_chain_of(*s)] for s in sectors]
    got = _run(f"""
const secs = {json.dumps(sectors)};
console.log(JSON.stringify(secs.map(a => P.drillChainOf(...a).map(b => [b.m, b.shell, b.band, b.wedge]))));
""")
    assert got == [[list(b) for b in chain] for chain in chains]
    blocks = sorted({b for chain in chains for b in chain})
    got = _run(f"""
const blocks = {json.dumps(blocks)};
console.log(JSON.stringify(blocks.map(([m, s, b, w]) => {{
  const x = P.drillBlockBounds({{m, shell: s, band: b, wedge: w}}, {EDGE_PC});
  const parent = P.drillParent({{m, shell: s, band: b, wedge: w}});
  return [[x.r0, x.r1, x.c0, x.c1, x.t0, x.t1], parent ? [parent.m, parent.shell, parent.band, parent.wedge] : null,
    m > 1 ? P.drillBandCount(m, s) : null, m > 1 ? P.drillWedgeCount(m, s, b) : null];
}})));
""")
    for block, (bounds, parent, rows, wedges) in zip(blocks, got):
        (inner, outer), (north, south), (start, end) = drill_block_bounds(block, EDGE_PC)
        assert bounds == pytest.approx([inner, outer, north, south, start, end], abs=1e-9, rel=1e-12), block
        expected_parent = drill_parent(DrillBlock(*block))
        assert parent == (list(expected_parent) if expected_parent else None)
        if block[0] > 1:
            assert (rows, wedges) == (drill_band_count(block[0], block[1]), drill_wedge_count(block[0], block[1], block[2]))


@pytest.mark.parametrize("g", [3, 9, 11])
def test_windows_match_python(g):
    rng = random.Random(g)
    centres = []
    for m in DRILL_LEVELS:
        for _ in range(5):
            if m == 1:
                shell = rng.choice([rng.randint(2, 40), rng.randint(40, 3000)])
                low, high = spherical.band_range(shell)
                band = rng.randint(low, high)
                centres.append(DrillBlock(1, shell, band, rng.randrange(spherical.wedge_count(shell, band))))
            else:
                shell = rng.choice([0, 1, 2, 3, rng.randint(0, max(1, 3700 // m))])
                rows = drill_band_count(m, shell)
                band = rng.randint(-(rows // 2), rows // 2)
                centres.append(DrillBlock(m, shell, band, rng.randrange(drill_wedge_count(m, shell, band))))
    got = _run(f"""
const centres = {json.dumps([list(c) for c in centres])};
console.log(JSON.stringify(centres.map(([m, s, b, w]) =>
  P.windowCells({{m, shell: s, band: b, wedge: w}}, {EDGE_PC}, {g}).map(c => [c.m, c.shell, c.band, c.wedge]))));
""")
    for centre, cells in zip(centres, got):
        assert cells == [list(c) for c in window_cells(centre, EDGE_PC, g)], tuple(centre)


def test_the_outline_and_what_a_block_can_hold_match_python():
    """The browser works the galaxy's outline out from the density shape, and counts the sectors under it that a
    block can hold; both must be what the generator's skeleton and a walk of the block's sectors say."""
    threshold = 1.0 / expected_system_count_at_density_1(pc_to_ly(EDGE_PC))
    extents, outer, _confirmed = build_shell_extents(SHAPE, EDGE_PC, threshold)
    bounds = GalaxyBounds(extents, EDGE_PC)
    rng = random.Random(13)
    blocks = []
    for m, count in ((3, 12), (9, 5)):
        for _ in range(count):
            shell = rng.randint(0, max(1, outer // m))
            rows = drill_band_count(m, shell)
            band = rng.randint(-min(rows // 2, 12), min(rows // 2, 12))
            blocks.append(DrillBlock(m, shell, band, rng.randrange(drill_wedge_count(m, shell, band))))
    blocks.append(DrillBlock(27, 0, 0, 0))                  # the block holding the core
    blocks.append(DrillBlock(3, outer // 3, 0, 0))          # on the outer edge
    pole = drill_band_count(3, 667) // 2
    blocks.append(DrillBlock(3, 667, pole, 0))              # high above the disk
    expected = [sum(1 for a in block_addresses(b) if bounds.qualifies(a[0], a[1])) for b in blocks]
    got = _run(f"""
const outlineShape = {json.dumps(_shape_with_threshold())};
const outline = P.galaxyOutline({EDGE_PC}, outlineShape, {(outer + 1) * EDGE_PC});
const blocks = {json.dumps([list(b) for b in blocks])}.map(([m, s, b, w]) => ({{m, shell: s, band: b, wedge: w}}));
console.log(JSON.stringify({{
  high: Array.from(outline.high), outer: outline.outerShell,
  totals: blocks.map(b => P.drillBlockTotal(b, outline)), holds: blocks.map(b => P.drillBlockHolds(b, outline)),
}}));
""")
    assert got["high"] == [high for _low, high in extents]
    assert got["outer"] == outer
    assert got["totals"] == expected
    assert got["holds"] == [total > 0 for total in expected]
    assert any(total > 0 for total in expected) and any(total == 0 for total in expected), expected


def test_the_outline_decides_which_sectors_are_in_the_galaxy():
    threshold = 1.0 / expected_system_count_at_density_1(pc_to_ly(EDGE_PC))
    extents, outer, _confirmed = build_shell_extents(SHAPE, EDGE_PC, threshold)
    bounds = GalaxyBounds(extents, EDGE_PC)
    sectors = [(shell, band, 0) for shell in (0, 1, 2, 100, outer, outer + 1)
               for band in {max(-3, spherical.band_range(shell)[0]), 0, min(3, spherical.band_range(shell)[1]),
                            spherical.band_range(shell)[1]}]
    got = _run(f"""
const outlineShape = {json.dumps(_shape_with_threshold())};
const outline = P.galaxyOutline({EDGE_PC}, outlineShape, {(outer + 1) * EDGE_PC});
console.log(JSON.stringify({json.dumps(sectors)}.map(([s, b, w]) => P.sectorInOutline(s, b, w, outline))));
""")
    assert got == [bounds.qualifies(shell, band) for shell, band, _ in sectors]


@pytest.mark.parametrize("address", [(0, 0, 0), (1, 0, 3), (1, 0, 5), (30, 0, 10), (30, 5, 7), (500, -9, 100)])
def test_a_sector_prism_is_closed_wound_outward_and_has_the_cells_volume(address):
    """The faceted solid's volume by the divergence theorem equals the cell's exact volume to within the facets'
    bow (a closed surface wound counter-clockwise from outside gives a positive volume)."""
    got = _run(f"""
const b = P.sectorBounds(...{json.dumps(list(address))}, {EDGE_PC});
const g = P.buildPrismGeometry([b]);
let volume = 0;
for (let k = 0; k < g.indices.length; k += 3) {{
  const v = [0, 1, 2].map(n => Array.from(g.positions.slice(3 * g.indices[k + n], 3 * g.indices[k + n] + 3)));
  volume += (v[0][0] * (v[1][1] * v[2][2] - v[1][2] * v[2][1]) - v[0][1] * (v[1][0] * v[2][2] - v[1][2] * v[2][0])
    + v[0][2] * (v[1][0] * v[2][1] - v[1][1] * v[2][0])) / 6;
}}
console.log(JSON.stringify({{volume, count: g.owners.length, owners: Array.from(new Set(g.owners)),
  inRange: Array.from(g.indices).every(i => i < g.owners.length), normals: g.normals.length, uvs: g.uvs.length}}));
""")
    exact = spherical.cell_volume(*address, EDGE_PC)
    assert got["volume"] == pytest.approx(exact, rel=0.02)
    assert got["volume"] <= exact * (1 + 1e-9), "chords lie inside the arcs"
    assert got["owners"] == [0] and got["inRange"]
    assert got["normals"] == 3 * got["count"] and got["uvs"] == 2 * got["count"]


def test_many_cells_in_one_geometry_keep_their_owners():
    got = _run(f"""
const cells = [[30, 0, 1], [30, 1, 1], [31, 0, 1]].map(a => P.sectorBounds(...a, {EDGE_PC}));
const g = P.buildPrismGeometry(cells);
console.log(JSON.stringify({{owners: Array.from(new Set(g.owners)), count: g.owners.length,
  indexMax: Math.max(...g.indices)}}));
""")
    assert got["owners"] == [0, 1, 2]
    assert got["indexMax"] < got["count"]


# --- galaxyblocks.js -----------------------------------------------------------

PALETTE = {"dim": [0.01, 0.02, 0.05], "accent": [0.1, 0.2, 0.9], "hot": [0.8, 0.85, 1.0],
           "placedLow": [0.06, 0.05, 0.03], "placedHigh": [1.0, 0.92, 0.75]}

SCENE_CELLS = """
const palette = %s;
const centre = {m: 9, shell: 40, band: 2, wedge: 5};
const cells = P.windowCells(centre, %s, 5).map(c => Object.assign({m: c.m, shell: c.shell, band: c.band, wedge: c.wedge},
  P.drillBlockBounds(c, %s), {filled: 0, total: 729, stats: null}));
const plain = part => Object.fromEntries(Object.entries(part).map(([k, v]) => [k, ArrayBuffer.isView(v) ? Array.from(v) : v]));
""" % (json.dumps(PALETTE), EDGE_PC, EDGE_PC)


def _scene(extra="", shape_expr="shape"):
    """Packs a window of cells with galaxyblocks' scene; `extra` runs before it (to fill some cells in), and
    the result has the parts' arrays as lists, the cells each part holds and the eye used."""
    return _run(SCENE_CELLS + extra + f"""
const scene = B.createBlockScene({{shape: {shape_expr}, palette}});
const eye = [9000, 0, 500];
const built = scene.buildCells(cells, eye, null, "default");
console.log(JSON.stringify({{solid: plain(built.solid), glass: plain(built.glass), held: built.cells.map(
  list => list.map(c => [c.shell, c.band, c.wedge, c.filled, c.total, c.density === undefined ? null : c.density, c.opacity])),
  total: cells.length}}));
""")


def test_block_scene_packs_every_cell_once():
    out = _scene()
    held = [tuple(cell[:3]) for part in out["held"] for cell in part]
    assert len(held) == out["total"] > 50
    assert len(set(held)) == len(held)
    for name in ("solid", "glass"):
        part = out[name]
        count = part["vertexCount"]
        assert len(part["positions"]) == len(part["centers"]) == len(part["colors"]) == 3 * count
        assert len(part["uvs"]) == 2 * count
        assert len(part["alphas"]) == len(part["fills"]) == len(part["owners"]) == count
        assert all(0 <= index < count for index in part["indices"])
    # Nothing is filled, so everything is translucent glass, 10-30% opaque (MAP.37).
    assert out["solid"]["vertexCount"] == 0
    assert all(25 <= a <= 77 for a in out["glass"]["alphas"])
    assert set(out["glass"]["fills"]) == {0}
    owners = out["glass"]["owners"]
    assert all(0 <= owner < len(out["held"][1]) for owner in owners)


def test_block_scene_sorts_glass_far_to_near():
    out = _scene()
    glass = out["glass"]
    eye = [9000, 0, 500]
    seen = []
    for v, owner in enumerate(glass["owners"]):
        if not seen or seen[-1][0] != owner:
            seen.append((owner, math.dist(glass["centers"][3 * v:3 * v + 3], eye)))
    distances = [d for _, d in seen]
    assert distances == sorted(distances, reverse=True)


def test_block_scene_draws_a_filled_cell_at_the_cap_and_a_denser_one_warmer():
    out = _scene("""
const stats = (systems) => ({systems, expected_systems: systems, stars: systems, mean_age_gy: 4.5, luminosity_sol: 3000});
Object.assign(cells[3], {filled: 729, total: 729, stats: stats(1e9)});
Object.assign(cells[4], {filled: 5, total: 729, stats: stats(900)});
""")
    held = {tuple(cell[:3]): cell for part in out["held"] for cell in part}
    filled = [cell for cell in held.values() if cell[3] > 0]
    assert sorted(cell[3] for cell in filled) == [5, 729]
    # A full cell is at the 50% cap (never solid), a partly filled one between unfilled space and it.
    glass = out["glass"]
    full_alphas = {a for a, f in zip(glass["alphas"], glass["fills"]) if f == 255}
    assert full_alphas == {128}
    assert all(cell[6] <= 0.5 for cell in held.values())
    assert 0.1 <= min(cell[6] for cell in held.values() if cell[3] == 0)
    partial = [cell for cell in filled if cell[3] == 5][0]
    assert 0.1 < partial[6] < 0.5


def test_block_scene_without_a_shape_has_no_density():
    out = _scene("cells.length = 6;", shape_expr="null")
    held = [cell for part in out["held"] for cell in part]
    assert len(held) == 6
    assert all(cell[5] is None for cell in held)
    assert all(cell[6] == pytest.approx(0.1) for cell in held)


# --- MAP.128/129: blocks colored by their sectors' stats ----------------------

STATS_SETUP = """
const stats = (systems, stars, age, lum) => ({systems, expected_systems: systems, stars, mean_age_gy: age, luminosity_sol: lum});
const cell = (filled, total, s) => ({density: 1, filled, total, stats: s});
const young = cell(1, 1, stats(40, 60, 0.5, 5000));
const old = cell(1, 1, stats(40, 60, 11, 5000));
const sparse = cell(1, 1, stats(2, 3, 4.5, 3));
const dense = cell(1, 1, stats(900, 1200, 4.5, 3000));
const empty = cell(1, 1, stats(0, 0, null, 0));
const unknown = cell(1, 1, null);
const unfilled = {density: 1, filled: 0, total: 1, stats: null};
const all = [young, old, sparse, dense, empty, unknown, unfilled];
const ranges = B.statsRanges(all);
"""


def test_stellar_age_sets_the_hue_blue_young_slate_average_amber_old():
    out = _run(STATS_SETUP + """
console.log(JSON.stringify({young: B.ageColor(0), disk: B.ageColor(B.DISK_AGE_GY), old: B.ageColor(B.OLD_AGE_GY),
  older: B.ageColor(13), mid: B.ageColor(B.DISK_AGE_GY / 2), none: B.ageColor(null)}));
""")
    assert out["young"][2] > out["young"][0]  # blue
    assert out["old"][0] > out["old"][2]  # amber
    assert out["disk"][0] == pytest.approx(out["disk"][1], abs=0.08)  # slate: nearly grey
    assert out["older"] == out["old"]
    assert out["none"] == out["young"]
    assert out["mid"] == pytest.approx([(a + b) / 2 for a, b in zip(out["young"], out["disk"])])


def test_a_filled_cells_opacity_follows_what_it_holds_and_never_passes_the_cap():
    out = _run(STATS_SETUP + """
const one = cell(1, 1, stats(1, 1, 8, 1));
const huge = cell(1, 1, stats(1e9, 1e9, 4.5, 3000));
console.log(JSON.stringify({sparse: B.statsOpacity(sparse), dense: B.statsOpacity(dense), one: B.statsOpacity(one),
  huge: B.statsOpacity(huge), young: B.statsOpacity(young), empty: B.statsOpacity(empty), unknown: B.statsOpacity(unknown),
  unfilled: B.statsOpacity(unfilled), base: B.blockOpacity(unfilled), range: [B.FILLED_OPACITY_SPARSE,
  B.FILLED_OPACITY_DENSE], step: B.FILLED_EMPTY_STEP, cap: B.MAX_FILL_OPACITY}));
""")
    low, high = out["range"]
    assert out["cap"] == 0.5 and high == out["cap"]
    assert out["unfilled"] == pytest.approx(out["base"])
    assert out["one"] <= out["sparse"] < out["young"] < out["dense"] <= out["cap"]
    assert out["huge"] == pytest.approx(out["cap"])
    # One system, however alone in view, is faint: Boss's lone one-star sector.
    assert low <= out["one"] < 0.2
    assert out["empty"] == pytest.approx(out["base"] + out["step"])
    assert out["unknown"] == pytest.approx(out["empty"])


def test_no_fill_is_ever_more_opaque_than_half():
    """Boss (2026-10-08): "never under any frame be more opaque than 50%": not
    a block with every sector filled, nor a sector holding millions of stars."""
    out = _run(STATS_SETUP + """
const worst = [];
for (const total of [1, 10, 531441]) {
  for (const filled of [0, 1, total]) {
    for (const systems of [0, 1, 50, 5e5, 1e10]) {
      const c = cell(filled, total, stats(systems, systems, 5, 1e6));
      const plain = {density: 100, filled, total};
      worst.push(B.blockOpacity(plain), B.blockOpacity({density: 1e-6, filled, total}));
      worst.push(B.statsOpacity(c));
    }
  }
}
console.log(JSON.stringify({max: Math.max(...worst), cap: B.MAX_FILL_OPACITY}));
""")
    assert out["max"] <= 0.5 and out["cap"] == 0.5


def test_a_partly_filled_block_is_between_unfilled_space_and_its_stars_opacity():
    out = _run(STATS_SETUP + """
const one = cell(1, 531441, stats(900, 1200, 4.5, 3000));
const half = cell(5, 10, stats(4500, 6000, 4.5, 3000));
console.log(JSON.stringify({base: B.blockOpacity(unfilled), one: B.statsOpacity(one), half: B.statsOpacity(half),
  dense: B.statsOpacity(dense), min: B.FILLED_MIN_STEP}));
""")
    assert out["base"] < out["one"] < out["half"] <= out["dense"] <= 0.5
    assert out["one"] - out["base"] >= out["min"] * (out["dense"] - out["base"]) - 1e-9


def test_one_or_two_cells_put_their_mean_in_the_middle_of_each_scale():
    """Boss (2026-10-08): with one or two stars, the average goes in the middle of the
    scale and the ends are padded, instead of the lone value sitting at an end."""
    out = _run(STATS_SETUP + """
const unfilledColor = [0.1, 0.2, 0.4];
const lone = cell(1, 1, stats(40, 60, 4, 500));
const a = cell(1, 1, stats(40, 60, 4, 100));
const b = cell(1, 1, stats(40, 60, 4, 1000));
const r1 = B.statsRanges([lone]);
const r2 = B.statsRanges([a, b]);
const shade = (c, r) => B.statsColor(c, unfilledColor, r);
const bright = (c, r) => B.statsColor(c, unfilledColor, r, "luminosity");
console.log(JSON.stringify({
  ends1: B.scaleEnds(r1.luminosity, "luminosity"), mean1: r1.luminosity.mean, count1: r1.luminosity.count,
  ends2: B.scaleEnds(r2.luminosity, "luminosity"), mean2: r2.luminosity.mean,
  loneBrightness: shade(lone, r1), loneRamp: bright(lone, r1), aRamp: bright(a, r2), bRamp: bright(b, r2),
  hue: B.ageColor(4), midRamp: B.statRampColor(0.5), lowRamp: B.statRampColor(0), highRamp: B.statRampColor(1),
  floor: B.BRIGHTNESS_FLOOR, legend: B.legendFor("luminosity", r1), few: B.FEW_ITEMS,
  pad: B.SCALE_MIN_HALF_SPAN.luminosity,
}));
""")
    lone_log = math.log10(500)
    low, high = out["ends1"]
    assert out["count1"] == 1 and out["mean1"] == pytest.approx(lone_log)
    assert (low + high) / 2 == pytest.approx(lone_log) and high - low == pytest.approx(2 * out["pad"])
    # Two cells: still centred on their mean, the ends at least as far out as the cells themselves.
    low2, high2 = out["ends2"]
    assert (low2 + high2) / 2 == pytest.approx(out["mean2"])
    assert low2 <= math.log10(100) and high2 >= math.log10(1000)
    # A lone cell sits in the middle of the brightness range (and of a single-statistic ramp), not at its bright end.
    def linear(v):
        return v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4

    mid = out["floor"] + (1 - out["floor"]) * 0.5
    assert out["loneBrightness"] == pytest.approx([linear(v) * mid for v in out["hue"]], rel=1e-6)
    assert out["loneRamp"] == pytest.approx([linear(v) for v in out["midRamp"]], abs=1e-9)
    # The two cells sit either side of the middle (a decade apart, so about a fifth either way of
    # the padded two-decade scale or more), the brighter nearer the ramp's high end.
    a_share = [linear(v) for v in out["lowRamp"]]
    assert out["aRamp"] != out["bRamp"] and out["loneRamp"] not in (out["aRamp"], out["bRamp"])
    assert out["aRamp"] != a_share and out["bRamp"] != [linear(v) for v in out["highRamp"]]
    legend = out["legend"]
    assert legend["lowText"] != legend["highText"]


def test_luminosity_sets_the_brightness_of_the_age_hue():
    out = _run(STATS_SETUP + """
const dim = cell(1, 1, stats(40, 60, 11, 1));
const lit = cell(1, 1, stats(40, 60, 11, 100000));
const r = B.statsRanges([dim, lit]);
const unfilledColor = [0.1, 0.2, 0.4];
console.log(JSON.stringify({dim: B.statsColor(dim, unfilledColor, r), lit: B.statsColor(lit, unfilledColor, r),
  none: B.statsColor(unfilled, unfilledColor, r), empty: B.statsColor(empty, unfilledColor, ranges),
  floor: B.BRIGHTNESS_FLOOR}));
""")
    assert out["none"] == [0.1, 0.2, 0.4]
    assert all(l > d for l, d in zip(out["lit"], out["dim"]))
    assert out["dim"][0] / out["lit"][0] == pytest.approx(out["floor"])
    # Same hue, only the brightness differs.
    assert out["dim"][0] / out["dim"][2] == pytest.approx(out["lit"][0] / out["lit"][2])
    # No stars: the unfilled color a shade more saturated (same mean).
    grey = sum([0.1, 0.2, 0.4]) / 3
    assert sum(out["empty"]) / 3 == pytest.approx(grey)
    assert out["empty"][2] - grey > 0.4 - grey and out["empty"][0] < 0.1


def test_a_mixed_block_is_colored_by_its_star_weighted_mean_age():
    # The stage API weights each sector's mean age by its stars (MAP.129);
    # the page colors the block from that one mean, so it lands between a
    # young and an old sector's hue, nearer the one with more stars.
    out = _run(STATS_SETUP + """
const mixed = cell(2, 2, stats(100, 100, (80 * 0.5 + 20 * 11) / 100, 800));
const r = B.statsRanges([mixed]);
console.log(JSON.stringify({mixed: B.statsColor(mixed, [0, 0, 0], r), young: B.statsColor(young, [0, 0, 0], B.statsRanges([young])),
  old: B.statsColor(old, [0, 0, 0], B.statsRanges([old]))}));
""")
    assert out["young"][0] < out["mixed"][0] < out["old"][0]
    assert out["mixed"][0] - out["young"][0] < out["old"][0] - out["mixed"][0]


def test_density_matches_the_python_model_with_a_core_and_a_denser_arm():
    from planetgen.galaxy.density import arm_terms, relative_density
    amplitude, level = arm_terms(2.2, 0.3)
    shape = build_galaxy_shape(2800.0, 350.0, 200.0, 1.0, 2, math.radians(15.0), amplitude,
                               arm_level=level, core_amplitude=12.0)
    points = [(8000.0, 100.0, 50.0), (0.0, 0.0, 0.0), (15.0, -10.0, 3.0), (300.0, -200.0, 10.0), (-12000.0, 4000.0, -600.0)]
    got = _run(f"console.log(JSON.stringify({json.dumps(points)}.map(p => P.relativeDensity(...p, shape))));", shape)
    for point, value in zip(points, got):
        assert value == pytest.approx(relative_density(point, shape), rel=1e-9)
