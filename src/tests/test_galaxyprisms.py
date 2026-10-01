"""
Tests for `html/static/galaxyprisms.js`, the Galaxy Map's density prisms,
run under node (skipped where node isn't installed): its density must
match `stellarObjects.galaxyDensity.relative_density` exactly, its
mega-blocks must follow the pixel scale, hold whole sectors, list only the
solid's surface and stay within budget, and its geometry must be well
formed.
"""

import json
import math
import os
import shutil
import subprocess

import pytest

from stellarObjects.galaxyDensity import build_galaxy_shape, relative_density
from stellarObjects.galaxyGeometry import ring_sector_count

NODE = shutil.which("node")
MODULE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "html", "static", "galaxyprisms.js")

pytestmark = pytest.mark.skipif(NODE is None, reason="node is not installed")

SHAPE = build_galaxy_shape(2800.0, 350.0, 200.0, 1.0, 2, math.radians(15.0), 0.4)
EDGE_PC = 4.0
GALAXY_RADIUS_PC = 24000.0


def _run(script):
    """Runs `script` as an ES module with the prisms module imported as `P`
    and the test shape as `shape`; returns its JSON output."""
    source = (
        f"import * as P from {json.dumps('file://' + os.path.abspath(MODULE))};\n"
        f"const shape = {json.dumps(SHAPE._asdict())};\n"
        f"{script}\n"
    )
    result = subprocess.run([NODE, "--input-type=module", "-e", source], capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


def test_density_matches_the_python_model():
    points = [(8000.0, 100.0, 50.0), (0.0, 0.0, 0.0), (300.0, -200.0, 10.0), (-12000.0, 4000.0, -600.0),
              (5.0, 5.0, 900.0), (20000.0, -3.0, 0.0)]
    got = _run(f"console.log(JSON.stringify({json.dumps(points)}.map(p => P.relativeDensity(...p, shape))));")
    for point, value in zip(points, got):
        assert value == pytest.approx(relative_density(point, SHAPE), rel=1e-12)


# The page's camera: 50 degree field of view on a 600 px tall canvas, and
# the view ball 1.6 x the orbit radius (FETCH_RADIUS_FACTOR).
PC_PER_PX_PER_ORBIT_PC = 2 * math.tan(math.radians(25.0)) / 600
VIEW_RADIUS_PER_ORBIT = 1.6


def _threshold_shape():
    """The test shape with the galaxy's own sector threshold, as the page
    embeds it (lib/galaxymap3d._density_shape)."""
    from stellarObjects.galaxySkeleton import expected_system_count_at_density_1
    from stellarObjects.utils import pc_to_ly

    return {**SHAPE._asdict(), "sector_min_density": 1.0 / expected_system_count_at_density_1(pc_to_ly(EDGE_PC))}


@pytest.mark.parametrize("pc_per_px", [0.0, 0.5, 1.0, 2.9, 3.1, 50.0, 1000.0, 1e9])
def test_block_size_follows_the_pixel_scale(pc_per_px):
    """m is the smallest power of 3 at least BLOCK_MIN_PX pixels across,
    capped where the disk is a few block rings wide."""
    out = _run(f"""
console.log(JSON.stringify({{m: P.blockSizeForScale({pc_per_px}, {EDGE_PC}, 0, {GALAXY_RADIUS_PC}), minPx: P.BLOCK_MIN_PX}}));
""")
    m, min_px = out["m"], out["minPx"]
    assert min_px == 4
    assert 3 ** round(math.log(m, 3)) == m
    cap = GALAXY_RADIUS_PC / (3 * EDGE_PC)
    assert m <= cap
    if m * 3 <= cap:
        assert m * EDGE_PC >= min_px * pc_per_px
    if m > 1:
        assert (m / 3) * EDGE_PC < min_px * pc_per_px


@pytest.mark.parametrize("m", [1, 3, 9, 27])
def test_block_boundaries_fall_on_sector_boundaries(m):
    """A block's rings and layers are whole sector rings and layers, and
    block layer 0 straddles the plane like sector layer 0."""
    from stellarObjects.galaxyGeometry import layer_bounds_pc, ring_bounds_pc

    for ring, slab in ((0, 0), (2, -1), (5, 3)):
        ranges = _run(f"console.log(JSON.stringify(P.blockSectorRanges({ring}, {slab}, {m})));")
        size = m * EDGE_PC
        assert ring_bounds_pc(ranges["ringFirst"], EDGE_PC)[0] == pytest.approx(ring * size)
        assert ring_bounds_pc(ranges["ringLast"], EDGE_PC)[1] == pytest.approx((ring + 1) * size)
        assert layer_bounds_pc(ranges["layerFirst"], EDGE_PC)[0] == pytest.approx((slab - 0.5) * size)
        assert layer_bounds_pc(ranges["layerLast"], EDGE_PC)[1] == pytest.approx((slab + 0.5) * size)


@pytest.mark.parametrize("m", [1, 3, 9, 27, 81, 243])
def test_every_sector_is_in_exactly_one_block(m):
    """Around a block ring the wedges' slot ranges tile every member ring
    with no gap or overlap, so the counts sum to all its sectors."""
    rings = sorted({0, 1, 2, 5, 13, 3855 // m - 1})
    out = _run(f"""
const m = {m};
console.log(JSON.stringify({json.dumps(rings)}.map(ring => {{
  const w = P.blockWedgeCount(ring, m);
  const members = [];
  for (let i = ring * m; i < ring * m + m; i++) {{
    const ranges = [];
    for (let seg = 0; seg < w; seg++) ranges.push(P.blockSlotRange(ring, seg, m, i));
    members.push([P.ringSectorCount(i), ranges]);
  }}
  let total = 0;
  for (let seg = 0; seg < w; seg++) total += P.blockSectorCount(ring, seg, 0, m, {EDGE_PC}, null);
  return {{members, total}};
}})));
""")
    for ring, entry in zip(rings, out):
        expected = 0
        for n, ranges in entry["members"]:
            covered = [k for r in ranges for k in range(r["first"], r["last"] + 1)]
            assert covered == list(range(n))
            expected += n * m
        assert entry["total"] == expected


def test_block_counts_skip_sectors_the_skeleton_leaves_out():
    from stellarObjects.galaxySkeleton import bound_relative_density_at

    shape = _threshold_shape()
    threshold = shape["sector_min_density"]
    m = 9
    cases = [(ring, slab) for ring in (0, 3, 100, 380) for slab in (0, 1, 3)]
    out = _run(f"""
const s = {json.dumps(shape)};
console.log(JSON.stringify({json.dumps(cases)}.map(([ring, slab]) =>
  [P.blockSectorCount(ring, 0, slab, {m}, {EDGE_PC}, s), P.blockExists(ring, slab, {m}, {EDGE_PC}, s, {GALAXY_RADIUS_PC})])));
""")
    for (ring, slab), (count, exists) in zip(cases, out):
        slots = _run(f"""
const w = P.blockWedgeCount({ring}, {m});
const out = [];
for (let i = {ring * m}; i < {ring * m + m}; i++) {{ const r = P.blockSlotRange({ring}, 0, {m}, i); out.push(r.last - r.first + 1); }}
console.log(JSON.stringify(out));
""")
        expected = 0
        for offset, per_layer in enumerate(slots):
            i = ring * m + offset
            for layer in range(slab * m - (m - 1) // 2, slab * m + (m - 1) // 2 + 1):
                if bound_relative_density_at(SHAPE, (i + 0.5) * EDGE_PC, layer * EDGE_PC) >= threshold:
                    expected += per_layer
        assert count == expected
        assert exists == (expected > 0)


@pytest.mark.parametrize("m, low, high", [(3, 0.6, 1.4), (9, 0.75, 1.3), (27, 0.75, 1.3), (81, 0.75, 1.3), (243, 0.75, 1.3)])
def test_blocks_hold_about_m_cubed_sectors(m, low, high):
    """Out past the core (block ring 5), a block holds m^3 sectors give or
    take the wedge rounding, and its wedge count is within 20% of the
    sector rule's round(2 pi (I + 1/2))."""
    out = _run(f"""
const m = {m};
const counts = [];
const ratios = [];
for (let ring = 5; ring < Math.floor(3855 / m); ring++) {{
  const w = P.blockWedgeCount(ring, m);
  ratios.push(w / Math.round(2 * Math.PI * (ring + 0.5)));
  for (let seg = 0; seg < w; seg += Math.max(1, Math.floor(w / 20))) counts.push(P.blockSectorCount(ring, seg, 0, m, {EDGE_PC}, null) / m ** 3);
}}
console.log(JSON.stringify({{counts, ratios}}));
""")
    assert out["counts"]
    assert low <= min(out["counts"]) and max(out["counts"]) <= high
    assert all(1 / 1.2 - 1e-9 <= r <= 1.2 + 1e-9 for r in out["ratios"])


@pytest.mark.parametrize("m", [27, 81, 243])
def test_most_block_wedges_follow_master_lines(m):
    """A wedge count dividing the innermost member ring's master count puts
    every wedge side on a slot wall in every member ring; most rings get
    one once blocks are large."""
    from stellarObjects.galaxyGeometry import ring_master_count

    rings = list(range(3855 // m))
    counts = _run(f"console.log(JSON.stringify({json.dumps(rings)}.map(i => P.blockWedgeCount(i, {m}))))")
    aligned = 0
    for ring, w in zip(rings, counts):
        if ring_master_count(ring * m) % w == 0:
            aligned += 1
            assert all(ring_sector_count(i) % w == 0 for i in range(ring * m, ring * m + m))
    assert aligned >= 0.85 * len(rings)


@pytest.mark.parametrize("center, view_radius, m, slice_slab, viewer_z", [
    ((8000.0, 0.0, 0.0), 300.0, 3, 0, None),
    ((8000.0, 0.0, 600.0), 300.0, 3, None, None),
    ((200.0, 100.0, 40.0), 150.0, 1, 10, None),
    ((0.0, 0.0, 0.0), 6000.0, 81, 0, None),
    ((0.0, 0.0, 0.0), 6000.0, 81, 0, 5000.0),
    ((8000.0, 0.0, 600.0), 300.0, 3, None, -900.0),
    ((8000.0, 0.0, 600.0), 300.0, 3, None, 650.0),
])
def test_surface_listing_matches_a_brute_force_check(center, view_radius, m, slice_slab, viewer_z):
    """The surface listing is exactly the blocks in view with a missing
    ring or layer neighbour (the axis counts as filled; layers above the
    slice count as empty; a missing layer neighbour counts only when the
    face it exposes can face the viewer), and nothing above the slice."""
    out = _run(f"""
const s = {json.dumps(_threshold_shape())};
const options = {{slice: {json.dumps(slice_slab)}, viewerZ: {json.dumps(viewer_z)}}};
const all = P.blocksInView({json.dumps(center)}, {view_radius}, {m}, {EDGE_PC}, s, {GALAXY_RADIUS_PC}, null, {{...options, surfaceOnly: false}});
const surface = P.blocksInView({json.dumps(center)}, {view_radius}, {m}, {EDGE_PC}, s, {GALAXY_RADIUS_PC}, null, options);
const cut = options.slice === null ? Infinity : options.slice;
const exists = (ring, slab) => ring < 0 || (slab <= cut && P.blockExists(ring, slab, {m}, {EDGE_PC}, s, {GALAXY_RADIUS_PC}));
const v = options.viewerZ;
const expected = all.filter(b => !(exists(b.ring - 1, b.slab) && exists(b.ring + 1, b.slab)
  && (exists(b.ring, b.slab - 1) || (v !== null && v >= b.z0))
  && (exists(b.ring, b.slab + 1) || (v !== null && v <= b.z1))));
const key = b => b.ring + "/" + b.seg + "/" + b.slab;
console.log(JSON.stringify({{all: all.length, surface: surface.map(key).sort(), expected: expected.map(key).sort(),
  maxSlab: Math.max(...all.map(b => b.slab))}}));
""")
    assert out["surface"] == out["expected"]
    assert 0 < len(out["surface"]) < out["all"]
    if slice_slab is not None:
        assert out["maxSlab"] <= slice_slab


@pytest.mark.parametrize("focus", [(0.0, 0.0, 0.0), (8000.0, 0.0, 0.0), (8000.0, 0.0, 120.0)])
@pytest.mark.parametrize("sliced", [True, False])
def test_every_zoom_fits_the_budget_and_shows_the_focus(focus, sliced):
    """From the whole galaxy down to a few sectors: the blocks stay within
    budget, are never finer than the scale allows, and (sliced at the
    focus) the focus's own block is drawn when the skeleton allows it."""
    out = _run(f"""
const s = {json.dumps(_threshold_shape())};
const focus = {json.dumps(focus)};
const out = [];
for (let orbit = 32000; orbit >= 10; orbit /= 2) {{
  const pcPerPixel = orbit * {PC_PER_PX_PER_ORBIT_PC};
  const view = P.blocksForView(focus, orbit * {VIEW_RADIUS_PER_ORBIT}, {EDGE_PC}, {GALAXY_RADIUS_PC}, s, () => new Map(),
    {{pcPerPixel, sliceZ: {json.dumps(sliced)} ? focus[2] : null}});
  const size = view.m * {EDGE_PC};
  const ring = Math.floor(Math.hypot(focus[0], focus[1]) / size);
  const slab = Math.round(focus[2] / size);
  const own = view.blocks.filter(b => b.ring === ring && b.slab === slab
    && Math.atan2(focus[1], focus[0]) >= b.t0 - 1e-9 && Math.atan2(focus[1], focus[0]) < b.t1);
  out.push({{orbit, m: view.m, pcPerPixel, count: view.blocks.length, slice: view.slice, own: own.length,
    exists: P.blockExists(ring, slab, view.m, {EDGE_PC}, s, {GALAXY_RADIUS_PC}),
    minM: P.blockSizeForScale(pcPerPixel, {EDGE_PC}, 0, {GALAXY_RADIUS_PC}), budget: P.BLOCK_BUDGET}});
}}
console.log(JSON.stringify(out));
""")
    for view in out:
        assert view["count"] <= view["budget"]
        assert view["m"] >= view["minM"]
        if sliced:
            assert view["slice"] == round(focus[2] / (view["m"] * EDGE_PC))
            assert view["own"] == (1 if view["exists"] else 0)
        else:
            assert view["slice"] is None
    assert out[0]["m"] > out[-1]["m"]


def test_js_grid_matches_python():
    from stellarObjects.galaxyGeometry import (
        provisional_sector_designation, sector_address_at, sector_cell_vertices_pc,
    )

    points = [(0.0, 0.0, 0.0), (8000.0, 12.5, -40.0), (-3.2, 7.9, 1.7), (1234.5, -9876.5, 300.0), (0.1, -0.1, -1.76)]
    got = _run(f"""
const pts = {json.dumps(points)};
console.log(JSON.stringify(pts.map(p => {{
  const a = P.sectorAddressAt(...p, {EDGE_PC});
  return {{a, v: P.cellVertices(P.sectorCellBounds(a.ring, a.layer, a.slot, {EDGE_PC}))}};
}})));
""")
    for point, entry in zip(points, got):
        address = sector_address_at(point, EDGE_PC)
        assert (entry["a"]["ring"], entry["a"]["layer"], entry["a"]["slot"]) == address
        for js, py in zip(entry["v"], sector_cell_vertices_pc(*address, EDGE_PC)):
            assert js == pytest.approx(py, abs=1e-9)
        assert provisional_sector_designation(*address)


def test_blocks_skip_empty_space():
    out = _run(f"""
const blocks = P.blocksInView([8000, 0, 3000], 1000, 27, {EDGE_PC}, shape, {GALAXY_RADIUS_PC}, new Map(), {{surfaceOnly: false}});
console.log(JSON.stringify(blocks.length));
""")
    assert out == 0


def test_blocks_carry_their_azimuthal_mean_for_arm_shading():
    """Each block's `meanDensity` is its density with the arm factor held at
    1, so density / mean is the arm factor the page shades the spiral by:
    within 1 +/- arm_amplitude, and spanning most of that range around a
    ring. With no arms the two are equal."""
    out = _run(f"""
const armed = P.blocksInView([8000, 0, 0], 4000, 27, {EDGE_PC}, shape, {GALAXY_RADIUS_PC}, new Map(), {{surfaceOnly: false}});
const flat = P.blocksInView([8000, 0, 0], 4000, 27, {EDGE_PC}, {{...shape, arm_amplitude: 0}}, {GALAXY_RADIUS_PC}, new Map(), {{surfaceOnly: false}});
console.log(JSON.stringify({{
  armed: armed.map(p => [p.ring, p.density, p.meanDensity]),
  flat: flat.map(p => [p.density, p.meanDensity]),
}}));
""")
    amplitude = SHAPE.arm_amplitude
    ratios = [density / mean for ring, density, mean in out["armed"] if ring * 27 * EDGE_PC > 2000]
    assert ratios
    assert all(1 - amplitude - 1e-9 <= r <= 1 + amplitude + 1e-9 for r in ratios)
    assert max(ratios) - min(ratios) > amplitude
    assert out["flat"] and all(d == pytest.approx(m, rel=1e-12) for d, m in out["flat"])


def test_geometry_is_closed_and_faces_outward():
    out = _run("""
const prisms = [
  {r0: 0, r1: 10, t0: 0, t1: Math.PI / 3, z0: -5, z1: 5},
  {r0: 100, r1: 110, t0: 1, t1: 1.1, z0: 0, z1: 10},
];
const g = P.buildPrismGeometry(prisms);
console.log(JSON.stringify({positions: Array.from(g.positions), normals: Array.from(g.normals),
  owners: Array.from(g.owners), indices: Array.from(g.indices)}));
""")
    positions = out["positions"]
    normals = out["normals"]
    indices = out["indices"]
    assert len(indices) % 3 == 0
    assert max(indices) < len(out["owners"])
    assert set(out["owners"]) == {0, 1}
    for t in range(0, len(indices), 3):
        a, b, c = (positions[3 * i:3 * i + 3] for i in indices[t:t + 3])
        u = [b[k] - a[k] for k in range(3)]
        v = [c[k] - a[k] for k in range(3)]
        n = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
        if math.hypot(*n) < 1e-9:
            continue  # the center ring's wedges meet in a point at the axis
        normal = normals[3 * indices[t]:3 * indices[t] + 3]
        assert sum(n[k] * normal[k] for k in range(3)) > 0


def test_wedge_counts_follow_the_cylindrical_sector_rule():
    rings = list(range(0, 40)) + [100, 1234]
    got = _run(f"console.log(JSON.stringify({json.dumps(rings)}.map(i => P.azimuthSegments(i))));")
    assert got == [ring_sector_count(i) for i in rings]


def test_js_master_wedge_rule_matches_python_to_the_edge():
    """ringSectorCount and ringMasterCount mirror galaxyGeometry ring by
    ring out to ring 4,000, past the default galaxy's edge (3,855)."""
    from stellarObjects.galaxyGeometry import ring_master_count

    got = _run("""
const out = [];
for (let i = 0; i <= 4000; i++) out.push([P.ringSectorCount(i), P.ringMasterCount(i)]);
console.log(JSON.stringify(out));
""")
    assert got == [[ring_sector_count(i), ring_master_count(i)] for i in range(4001)]


@pytest.mark.parametrize("center, view_radius", [((0.0, 0.0, 1268.0), 30.0), ((15420.0, 0.0, 0.0), 30.0)])
def test_one_sector_blocks_outline_exactly_the_skeletons_layers(center, view_radius):
    # With the galaxy's sector threshold, a single-sector block exists
    # exactly when build_layer_extents' bound lets that sector exist.
    from stellarObjects.galaxySkeleton import bound_relative_density_at, expected_system_count_at_density_1
    from stellarObjects.utils import pc_to_ly

    threshold = 1.0 / expected_system_count_at_density_1(pc_to_ly(EDGE_PC))
    out = _run(f"""
const s = Object.assign({{}}, shape, {{sector_min_density: {threshold}}});
console.log(JSON.stringify(P.blocksInView({json.dumps(center)}, {view_radius}, 1, {EDGE_PC}, s, {GALAXY_RADIUS_PC}, new Map(), {{surfaceOnly: false}})));
""")
    drawn = {(p["ring"], p["slab"]) for p in out}
    everything = _run(f"""
const s = Object.assign({{}}, shape, {{sector_min_density: 1e-30}});
console.log(JSON.stringify(P.blocksInView({json.dumps(center)}, {view_radius}, 1, {EDGE_PC}, s, {GALAXY_RADIUS_PC}, new Map(), {{surfaceOnly: false}})));
""")
    candidates = {(p["ring"], p["slab"]) for p in everything}
    assert drawn and candidates - drawn  # the view straddles the galaxy's edge
    for ring, layer in candidates:
        qualifies = bound_relative_density_at(SHAPE, (ring + 0.5) * EDGE_PC, layer * EDGE_PC) >= threshold
        assert ((ring, layer) in drawn) == qualifies


@pytest.mark.parametrize("m", [1, 3, 27])
def test_block_address_and_block_at_match_the_listing(m):
    """blockAddressAt puts each listed block's own middle in that block,
    blockAt rebuilds it exactly, and at m = 1 the block holding a sector's
    center is that sector (the page sums filled sectors this way)."""
    out = _run(f"""
const blocks = P.blocksInView([8000, 0, 0], 60 * {m}, {m}, {EDGE_PC}, shape, {GALAXY_RADIUS_PC}, null, {{surfaceOnly: false}});
const checked = blocks.filter((b, n) => n % 7 === 0).map(b => {{
  const r = (b.r0 + b.r1) / 2, t = (b.t0 + b.t1) / 2, z = (b.z0 + b.z1) / 2;
  const a = P.blockAddressAt(r * Math.cos(t), r * Math.sin(t), z, {m}, {EDGE_PC});
  return {{b, a, again: P.blockAt(b.ring, b.seg, b.slab, {m}, {EDGE_PC}, shape, null)}};
}});
const sectors = [[2000, 0, 5], [3, -2, 1], [0, 0, 2]].map(([ring, layer, slot]) => {{
  const c = P.cellCoordinates(P.sectorCellBounds(ring, layer, slot, {EDGE_PC})).cartesian;
  return [[ring, layer, slot], P.blockAddressAt(...c, 1, {EDGE_PC})];
}});
console.log(JSON.stringify({{checked, sectors, bare: P.blockAt(4, 1, 0, {m}, {EDGE_PC}, null, null)}}));
""")
    assert out["checked"]
    for entry in out["checked"]:
        block, address, again = entry["b"], entry["a"], entry["again"]
        assert (address["ring"], address["seg"], address["slab"]) == (block["ring"], block["seg"], block["slab"])
        for field in ("r0", "r1", "t0", "t1", "z0", "z1", "density", "meanDensity"):
            assert again[field] == pytest.approx(block[field], rel=1e-12, abs=1e-12)
    if m == 1:
        for (ring, layer, slot), address in out["sectors"]:
            assert (address["ring"], address["slab"], address["seg"]) == (ring, layer, slot)
    assert out["bare"]["density"] is None
