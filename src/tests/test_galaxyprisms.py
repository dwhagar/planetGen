"""
Tests for `html/static/galaxyprisms.js`, the Galaxy Map's density prisms,
run under node (skipped where node isn't installed): its density must
match `stellarObjects.galaxyDensity.relative_density` exactly, its grid
must stay within budget and cover the view, and its geometry must be
well formed.
"""

# TODO(galaxy-map #13/#14): add tests for
#   - blockSizeForScale against pcPerPixel;
#   - block sector counts within a few percent of m³, exact with aligned
#     wedges;
#   - surfaceBlocksInView against a brute-force exposed-block check;
#   - the slice hiding the right layers;
#   - the budget holding at every zoom.
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


@pytest.mark.parametrize("center, view_radius", [
    ((0.0, 0.0, 0.0), 32000.0),
    ((8000.0, 0.0, 0.0), 4000.0),
    ((8000.0, 0.0, 0.0), 800.0),
    ((8000.0, 0.0, 20.0), 51.0),
    ((0.0, 0.0, 0.0), 300.0),
])
def test_prisms_fit_the_budget_and_cover_the_view(center, view_radius):
    out = _run(f"""
const center = {json.dumps(center)};
const m = P.sectorsPerPrism(center, {view_radius}, {EDGE_PC}, {GALAXY_RADIUS_PC}, shape, 0);
const prisms = P.prismsInView(center, {view_radius}, m, {EDGE_PC}, shape, {GALAXY_RADIUS_PC}, new Map());
console.log(JSON.stringify({{m, budget: P.PRISM_BUDGET, prisms}}));
""")
    m = out["m"]
    size = m * EDGE_PC
    prisms = out["prisms"]
    # A power of three sectors a side, and finer detail as the view narrows.
    assert m >= 1 and 3 ** round(math.log(m, 3)) == m
    assert size <= view_radius
    assert 0 < len(prisms) <= 1.5 * out["budget"]
    for p in prisms:
        assert p["r0"] == pytest.approx(p["ring"] * size)
        assert p["r1"] - p["r0"] == pytest.approx(size)
        assert p["z1"] - p["z0"] == pytest.approx(size)
        assert p["density"] >= 0.02
    rc = math.hypot(center[0], center[1])
    theta = math.atan2(center[1], center[0]) % (2 * math.pi)
    assert any(
        p["r0"] <= rc <= p["r1"] and p["z0"] <= center[2] <= p["z1"]
        and (p["t0"] <= theta <= p["t1"] or p["r0"] == 0)
        for p in prisms
    )


def test_groups_stay_big_enough_on_screen():
    # At 50 pc per pixel a group must be at least MIN_PRISM_PX pixels wide.
    out = _run(f"""
const m = P.sectorsPerPrism([8000, 0, 0], 400, {EDGE_PC}, {GALAXY_RADIUS_PC}, shape, 50);
console.log(JSON.stringify({{m, minPx: P.MIN_PRISM_PX}}));
""")
    assert out["m"] * EDGE_PC >= 50 * out["minPx"]
    assert out["m"] * EDGE_PC / 3 < 50 * out["minPx"]


@pytest.mark.parametrize("m", [1, 3, 9, 27])
def test_group_boundaries_fall_on_sector_boundaries(m):
    """A group's rings and layers are whole sector rings and layers, and
    group layer 0 straddles the plane like sector layer 0."""
    from stellarObjects.galaxyGeometry import layer_bounds_pc, ring_bounds_pc

    for ring, slab in ((0, 0), (2, -1), (5, 3)):
        ranges = _run(f"console.log(JSON.stringify(P.groupSectorRanges({ring}, {slab}, {m})));")
        size = m * EDGE_PC
        assert ring_bounds_pc(ranges["ringFirst"], EDGE_PC)[0] == pytest.approx(ring * size)
        assert ring_bounds_pc(ranges["ringLast"], EDGE_PC)[1] == pytest.approx((ring + 1) * size)
        assert layer_bounds_pc(ranges["layerFirst"], EDGE_PC)[0] == pytest.approx((slab - 0.5) * size)
        assert layer_bounds_pc(ranges["layerLast"], EDGE_PC)[1] == pytest.approx((slab + 0.5) * size)


def test_group_sector_counts_add_up_to_every_sector():
    """Every sector's center falls in exactly one group, so a group ring's
    counts sum to all its member rings' slots times m layers."""
    from stellarObjects.galaxyGeometry import ring_sector_count

    m = 9
    for ring in (0, 1, 4):
        total = _run(f"""
const n = P.ringSectorCount({ring});
let total = 0;
for (let seg = 0; seg < n; seg++) {{
  const step = 2 * Math.PI / n;
  total += P.groupSectorCount({{ring: {ring}, slab: 0, t0: seg * step, t1: (seg + 1) * step}}, {m});
}}
console.log(JSON.stringify(total));
""")
        expected = sum(ring_sector_count(i) for i in range(ring * m, ring * m + m)) * m
        assert total == expected


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


def test_prisms_skip_empty_space():
    out = _run(f"""
const prisms = P.prismsInView([8000, 0, 3000], 1000, 27, {EDGE_PC}, shape, {GALAXY_RADIUS_PC}, new Map());
console.log(JSON.stringify(prisms.length));
""")
    assert out == 0


def test_prisms_carry_their_azimuthal_mean_for_arm_shading():
    """Each prism's `meanDensity` is its density with the arm factor held at
    1, so density / mean is the arm factor the page shades the spiral by:
    within 1 +/- arm_amplitude, and spanning most of that range around a
    ring. With no arms the two are equal."""
    out = _run(f"""
const armed = P.prismsInView([8000, 0, 0], 4000, 27, {EDGE_PC}, shape, {GALAXY_RADIUS_PC}, new Map());
const flat = P.prismsInView([8000, 0, 0], 4000, 27, {EDGE_PC}, {{...shape, arm_amplitude: 0}}, {GALAXY_RADIUS_PC}, new Map());
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
def test_one_sector_prisms_outline_exactly_the_skeletons_layers(center, view_radius):
    # With the galaxy's sector threshold, a single-sector prism is drawn
    # exactly when build_layer_extents' bound lets that sector exist.
    from stellarObjects.galaxySkeleton import bound_relative_density_at, expected_system_count_at_density_1
    from stellarObjects.utils import pc_to_ly

    threshold = 1.0 / expected_system_count_at_density_1(pc_to_ly(EDGE_PC))
    out = _run(f"""
const s = Object.assign({{}}, shape, {{sector_min_density: {threshold}}});
console.log(JSON.stringify(P.prismsInView({json.dumps(center)}, {view_radius}, 1, {EDGE_PC}, s, {GALAXY_RADIUS_PC}, new Map())));
""")
    drawn = {(p["ring"], p["slab"]) for p in out}
    everything = _run(f"""
const s = Object.assign({{}}, shape, {{sector_min_density: 1e-30}});
console.log(JSON.stringify(P.prismsInView({json.dumps(center)}, {view_radius}, 1, {EDGE_PC}, s, {GALAXY_RADIUS_PC}, new Map())));
""")
    candidates = {(p["ring"], p["slab"]) for p in everything}
    assert drawn and candidates - drawn  # the view straddles the galaxy's edge
    for ring, layer in candidates:
        qualifies = bound_relative_density_at(SHAPE, (ring + 0.5) * EDGE_PC, layer * EDGE_PC) >= threshold
        assert ((ring, layer) in drawn) == qualifies
