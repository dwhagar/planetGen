"""
Tests for `html/static/galaxyprisms.js`, the Galaxy Map's density prisms,
run under node (skipped where node isn't installed): its density must
match `planetgen.galaxy.density.relative_density` exactly, its
mega-blocks must follow the pixel scale, hold whole sectors, list only the
solid's surface and stay within budget, and its geometry must be well
formed. Also `html/static/galaxyblocks.js`, which packs a view's blocks
for the GPU (in the page's Web Worker): every block once, filled sectors
counted into the right blocks, and arrays that agree with each other.
"""

import json
import math
import os
import shutil
import subprocess

import pytest

from planetgen.galaxy.density import build_galaxy_shape, relative_density, shape_with_terms
from planetgen.galaxy.geometry import ring_sector_count

NODE = shutil.which("node")
STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "html", "static")
MODULE = os.path.join(STATIC, "galaxyprisms.js")
BLOCKS_MODULE = os.path.join(STATIC, "galaxyblocks.js")

pytestmark = pytest.mark.skipif(NODE is None, reason="node is not installed")

SHAPE = build_galaxy_shape(2800.0, 350.0, 200.0, 1.0, 2, math.radians(15.0), 0.4)
EDGE_PC = 4.0
GALAXY_RADIUS_PC = 24000.0


def _run(script):
    """Runs `script` as an ES module with the prisms module imported as `P`,
    the block scene module as `B` and the test shape as `shape`; returns
    its JSON output."""
    source = (
        f"import * as P from {json.dumps('file://' + os.path.abspath(MODULE))};\n"
        f"import * as B from {json.dumps('file://' + os.path.abspath(BLOCKS_MODULE))};\n"
        f"const shape = {json.dumps(shape_with_terms(SHAPE))};\n"
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
    from planetgen.galaxy.skeleton import expected_system_count_at_density_1
    from planetgen.physics.units import pc_to_ly

    return {**shape_with_terms(SHAPE), "sector_min_density": 1.0 / expected_system_count_at_density_1(pc_to_ly(EDGE_PC))}


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
    from planetgen.galaxy.geometry import layer_bounds_pc, ring_bounds_pc

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
    from planetgen.galaxy.skeleton import bound_relative_density_at

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
    from planetgen.galaxy.geometry import ring_master_count

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
    ((8000.0, 0.0, 1640.0), 300.0, 3, None, None),
    ((200.0, 100.0, 40.0), 150.0, 1, 10, None),
    ((0.0, 0.0, 0.0), 6000.0, 81, 0, None),
    ((0.0, 0.0, 0.0), 6000.0, 81, 0, 5000.0),
    ((8000.0, 0.0, 1640.0), 300.0, 3, None, -900.0),
    ((8000.0, 0.0, 1640.0), 300.0, 3, None, 1690.0),
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
    from planetgen.galaxy.geometry import (
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


def test_js_grid_matches_python_one_ulp_inside_every_face():
    """GEN.31: a point one ulp under a layer's top face (layer 0's above
    all), or just inside a ring's face, gets the same cell in JavaScript
    as in Python. Slots rest on atan2, which node and libm may round an
    ulp apart, so right at a slot face each language is checked against
    its own angle bounds instead."""
    from planetgen.galaxy.geometry import (
        layer_bounds_pc, ring_bounds_pc, ring_radius_pc, sector_address_at, slot_angle_bounds,
    )

    points = []
    r = ring_radius_pc(5, EDGE_PC)
    for layer in (-3, -1, 0, 1, 2, 317):
        bottom, top = layer_bounds_pc(layer, EDGE_PC)
        points += [(r, 0.0, bottom), (r, 0.0, math.nextafter(top, -math.inf)), (r, 0.0, top)]
    for ring in (1, 2, 37, 3855):
        inner, outer = ring_bounds_pc(ring, EDGE_PC)
        points += [(inner, 0.0, 0.0), (math.nextafter(outer, 0.0), 0.0, 0.0), (0.0, inner, -2.0)]
    got = _run(f"""
const pts = {json.dumps(points)};
console.log(JSON.stringify(pts.map(p => P.sectorAddressAt(...p, {EDGE_PC}))));
""")
    for point, entry in zip(points, got):
        assert (entry["ring"], entry["layer"], entry["slot"]) == sector_address_at(point, EDGE_PC), point
    _bottom, top = layer_bounds_pc(0, EDGE_PC)
    assert sector_address_at((r, 0.0, math.nextafter(top, -math.inf)), EDGE_PC)[1] == 0

    faces = []
    for ring in (1, 2, 37, 3855):
        rr = ring_radius_pc(ring, EDGE_PC)
        for slot in (0, 1, 5):
            start, end = slot_angle_bounds(ring, slot)
            for theta in (start, math.nextafter(end, 0.0)):
                faces.append((rr * math.cos(theta), rr * math.sin(theta), 0.0))
    held = _run(f"""
const pts = {json.dumps(faces)};
console.log(JSON.stringify(pts.map(p => {{
  const a = P.sectorAddressAt(...p, {EDGE_PC});
  const b = P.sectorCellBounds(a.ring, a.layer, a.slot, {EDGE_PC});
  let t = Math.atan2(p[1], p[0]);
  if (t < 0) t += 2 * Math.PI;
  return b.t0 <= t && t < b.t1;
}})));
""")
    assert all(held)
    for point in faces:
        ring, _layer, slot = sector_address_at(point, EDGE_PC)
        low, high = slot_angle_bounds(ring, slot)
        assert low <= math.atan2(point[1], point[0]) % (2 * math.pi) < high


def test_blocks_skip_empty_space():
    out = _run(f"""
const blocks = P.blocksInView([8000, 0, 4500], 1000, 27, {EDGE_PC}, shape, {GALAXY_RADIUS_PC}, new Map(), {{surfaceOnly: false}});
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



def test_geometry_leaves_out_faces_shared_by_listed_neighbours():
    """Opaque blocks from one listing drop the faces they share: a stacked
    pair loses the top and bottom between them, wedge neighbours their
    common side, and rings with matching wedge counts the wall between
    them. Blocks without an address, or built without skipShared (the
    translucent ones), keep every face."""
    out = _run("""
const dt = 2 * Math.PI / 12;
function block(ring, seg, slab) {
  return {ring: ring, seg: seg, slab: slab, r0: 10 * ring, r1: 10 * ring + 10, t0: seg * dt, t1: (seg + 1) * dt,
          z0: 10 * slab - 5, z1: 10 * slab + 5};
}
function count(prisms) {
  return P.buildPrismGeometry(prisms, {skipShared: true}).owners.length;
}
function bare(prisms) {
  return prisms.map(p => ({r0: p.r0, r1: p.r1, t0: p.t0, t1: p.t1, z0: p.z0, z1: p.z1}));
}
const one = count([block(3, 0, 0)]);
const stacked = [block(3, 0, 0), block(3, 0, 1)];
const side = [block(3, 0, 0), block(3, 1, 0)];
const wrap = [block(3, 11, 0), block(3, 0, 0)];
const walls = [block(3, 0, 0), block(4, 0, 0)];
const glass = P.buildPrismGeometry(stacked).owners.length;
console.log(JSON.stringify({one: one, glass: glass, stacked: [count(stacked), count(bare(stacked))],
  side: [count(side), count(bare(side))], wrap: [count(wrap), count(bare(wrap))],
  walls: [count(walls), count(bare(walls))]}));
""")
    arcs_row = (out["one"] - 8) // 8  # 4 curved faces of 2 rows each, plus two 4-vertex sides
    assert out["stacked"] == [out["stacked"][1] - 2 * 2 * arcs_row, 2 * out["one"]]
    assert out["side"] == [2 * out["one"] - 8, 2 * out["one"]]
    assert out["wrap"] == [2 * out["one"] - 8, 2 * out["one"]]
    assert out["walls"] == [2 * out["one"] - 2 * 2 * arcs_row, 2 * out["one"]]
    assert out["glass"] == 2 * out["one"]

def test_wedge_counts_follow_the_cylindrical_sector_rule():
    rings = list(range(0, 40)) + [100, 1234]
    got = _run(f"console.log(JSON.stringify({json.dumps(rings)}.map(i => P.azimuthSegments(i))));")
    assert got == [ring_sector_count(i) for i in rings]


def test_js_master_wedge_rule_matches_python_to_the_edge():
    """ringSectorCount and ringMasterCount mirror galaxyGeometry ring by
    ring out to ring 4,000, past the default galaxy's edge (3,855)."""
    from planetgen.galaxy.geometry import ring_master_count

    got = _run("""
const out = [];
for (let i = 0; i <= 4000; i++) out.push([P.ringSectorCount(i), P.ringMasterCount(i)]);
console.log(JSON.stringify(out));
""")
    assert got == [[ring_sector_count(i), ring_master_count(i)] for i in range(4001)]


@pytest.mark.parametrize("center, view_radius", [((0.0, 0.0, 4760.0), 30.0), ((16212.0, 0.0, 0.0), 30.0)])
def test_one_sector_blocks_outline_exactly_the_skeletons_layers(center, view_radius):
    # With the galaxy's sector threshold, a single-sector block exists
    # exactly when build_layer_extents' bound lets that sector exist.
    from planetgen.galaxy.skeleton import bound_relative_density_at, expected_system_count_at_density_1
    from planetgen.physics.units import pc_to_ly

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


# --- galaxyblocks.js -----------------------------------------------------------

PALETTE = {"dim": [0.01, 0.02, 0.05], "accent": [0.1, 0.2, 0.9], "hot": [0.8, 0.85, 1.0],
           "placedLow": [0.06, 0.05, 0.03], "placedHigh": [1.0, 0.92, 0.75]}


def _build(view, points=(), shape_expr="shape"):
    """galaxyblocks' build for `view`, with filled `points` ([x, y, z,
    count, is_sector, relative]); returns the parts' arrays as lists plus
    what the same view's listing holds."""
    config = {"edgePc": EDGE_PC, "galaxyRadius": GALAXY_RADIUS_PC, "minPx": 4, "budget": 60000, "palette": PALETTE}
    return _run(f"""
const config = Object.assign({json.dumps(config)}, {{shape: {shape_expr}}});
const scene = B.createBlockScene(config);
scene.setFilled(Float64Array.from({json.dumps([v for p in points for v in p])}));
const view = {json.dumps(view)};
const built = scene.build(view);
const listed = config.shape ? P.blocksForView(view.center, view.viewRadius, {EDGE_PC}, {GALAXY_RADIUS_PC}, config.shape, () => null,
  {{pcPerPixel: view.pcPerPixel, sliceZ: view.sliceZ, viewerZ: view.viewerZ, minPx: 4, budget: 60000}}) : {{blocks: []}};
const plain = part => Object.fromEntries(Object.entries(part).map(([k, v]) => [k, ArrayBuffer.isView(v) ? Array.from(v) : v]));
console.log(JSON.stringify({{m: built.m, solid: plain(built.solid), glass: plain(built.glass), cellStride: B.CELL_STRIDE,
  pointStride: B.POINT_STRIDE, listed: listed.blocks.map(b => [b.ring, b.seg, b.slab]),
  transfers: B.transferablesOf(built).length}}));
""")


def _cells(part, stride):
    cells = part["cells"]
    return [cells[i:i + stride] for i in range(0, len(cells), stride)]


VIEW = {"center": [8000.0, 0.0, 0.0], "viewRadius": 120.0, "pcPerPixel": 2.0, "sliceZ": 0.0, "viewerZ": 300.0,
        "eye": [8200.0, 0.0, 300.0]}


def test_block_scene_packs_every_listed_block_once():
    out = _build(VIEW)
    stride = out["cellStride"]
    assert out["pointStride"] == 6
    keys = [tuple(int(v) for v in cell[:3]) for part in ("solid", "glass") for cell in _cells(out[part], stride)]
    assert sorted(keys) == sorted(tuple(b) for b in out["listed"])
    assert len(set(keys)) == len(keys)
    assert out["transfers"] == 18
    for name in ("solid", "glass"):
        part = out[name]
        count = part["vertexCount"]
        assert len(part["positions"]) == len(part["centers"]) == len(part["colors"]) == 3 * count
        assert len(part["uvs"]) == 2 * count
        assert len(part["alphas"]) == len(part["fills"]) == len(part["owners"]) == count
        assert all(0 <= owner < len(part["cells"]) // stride for owner in part["owners"])
        assert all(0 <= index < count for index in part["indices"])
    # Nothing is filled, so everything is translucent glass, 10-30% opaque (MAP.37).
    assert out["solid"]["vertexCount"] == 0
    assert all(25 <= a <= 77 for a in out["glass"]["alphas"])
    assert set(out["glass"]["fills"]) == {0}


def test_block_scene_sorts_glass_far_to_near():
    out = _build(VIEW)
    glass = out["glass"]
    eye = VIEW["eye"]
    seen = []
    for v, owner in enumerate(glass["owners"]):
        if not seen or seen[-1][0] != owner:
            center = glass["centers"][3 * v:3 * v + 3]
            seen.append((owner, math.dist(center, eye)))
    distances = [d for _, d in seen]
    assert distances == sorted(distances, reverse=True)


def test_block_scene_counts_filled_sectors_into_blocks():
    """At one sector per block: a generated sector's own block is drawn at
    the opacity cap (never solid), keeps its point index and takes the warm color; a coarse cell's count
    lands in the block holding its middle, unless the slice hides it."""
    view = dict(VIEW, pcPerPixel=0.1, viewRadius=40.0)
    ring, layer, slot = 2000, 0, 5
    center = _run(f"console.log(JSON.stringify(P.cellCoordinates(P.sectorCellBounds({ring}, {layer}, {slot}, {EDGE_PC})).cartesian));")
    points = [[*center, 1, 1, 1.0], [8010.0, 3.0, 0.0, 5, 0, 0], [8010.0, 3.0, 8.0, 2, 0, 0]]
    out = _build(view, points)
    stride = out["cellStride"]
    assert out["m"] == 1
    assert out["solid"]["vertexCount"] == 0
    solid = [cell for cell in _cells(out["glass"], stride) if cell[5] > 0]
    by_address = {tuple(int(v) for v in cell[:3]): cell for cell in solid}
    sector_cell = by_address.pop((ring, slot, layer))
    assert sector_cell[5] == sector_cell[6] == 1  # filled == total
    assert sector_cell[10] == 0  # the first point
    assert sector_cell[7:10] == pytest.approx(center)
    # Five in one sector's block counts as all of it.
    [(_, coarse)] = by_address.items()
    assert coarse[5] == coarse[6] == 5 and coarse[10] is None
    # Filled blocks sit at the cap, 50% (Boss, 2026-10-08): the fill alphas of the filled
    # cells are all at most half of 255.
    filled_owners = {i for i, cell in enumerate(_cells(out["glass"], stride)) if cell[5] > 0}
    filled_alphas = {a for a, owner, f in zip(out["glass"]["alphas"], out["glass"]["owners"], out["glass"]["fills"])
                     if owner in filled_owners and f == 255}
    assert filled_alphas == {128}
    # The sector's warm placed color: red above blue, unlike the bluish
    # density ramp the other block takes.
    cells = _cells(out["glass"], stride)
    owners = out["glass"]["owners"]
    colors = out["glass"]["colors"]
    sector_index = cells.index(sector_cell)
    coarse_index = cells.index(coarse)
    v_sector = owners.index(sector_index)
    v_coarse = owners.index(coarse_index)
    assert colors[3 * v_sector] > colors[3 * v_sector + 2]
    assert colors[3 * v_coarse] < colors[3 * v_coarse + 2]
    assert len(filled_owners) == 2


def test_block_scene_without_a_shape_draws_only_filled_blocks():
    out = _build(dict(VIEW, pcPerPixel=0.1), [[8000.0, 0.0, 0.0, 1, 0, 0]], shape_expr="null")
    stride = out["cellStride"]
    cells = _cells(out["solid"], stride) + _cells(out["glass"], stride)
    assert len(cells) == 1
    assert cells[0][3] is None  # no density without a shape (NaN, as JSON null)


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
