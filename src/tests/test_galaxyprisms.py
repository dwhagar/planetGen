"""
Tests for `html/static/galaxyprisms.js`, the Galaxy Map's density prisms,
run under node (skipped where node isn't installed): its density must
match `stellarObjects.galaxyDensity.relative_density` exactly, its grid
must stay within budget and cover the view, and its geometry must be
well formed.
"""

import json
import math
import os
import shutil
import subprocess

import pytest

from stellarObjects.galaxyDensity import build_galaxy_shape, relative_density

NODE = shutil.which("node")
MODULE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "html", "static", "galaxyprisms.js")

pytestmark = pytest.mark.skipif(NODE is None, reason="node is not installed")

SHAPE = build_galaxy_shape(2800.0, 350.0, 200.0, 1.0, 2, math.radians(15.0), 0.4)
EDGE_PC = 3.526
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
const shells = P.shellsPerPrism(center, {view_radius}, {EDGE_PC}, {GALAXY_RADIUS_PC}, shape);
const prisms = P.prismsInView(center, {view_radius}, shells * {EDGE_PC}, shape, {GALAXY_RADIUS_PC}, new Map());
console.log(JSON.stringify({{shells, budget: P.PRISM_BUDGET, prisms}}));
""")
    shells = out["shells"]
    size = shells * EDGE_PC
    prisms = out["prisms"]
    # A power of two shells, and the view gets finer detail as it narrows.
    assert shells >= 1 and shells & (shells - 1) == 0
    assert size <= view_radius
    assert 0 < len(prisms) <= 1.5 * out["budget"]
    for p in prisms:
        assert p["r0"] == pytest.approx(p["ring"] * size)
        assert p["r1"] - p["r0"] == pytest.approx(size)
        assert p["z1"] - p["z0"] == pytest.approx(size)
        assert p["density"] >= 0.02
    # The prism holding the view's center is there (every one of these
    # centers sits in the disk, so it's dense enough to draw).
    rc = math.hypot(center[0], center[1])
    theta = math.atan2(center[1], center[0]) % (2 * math.pi)
    assert any(
        p["r0"] <= rc <= p["r1"] and p["z0"] <= center[2] <= p["z1"]
        and (p["t0"] <= theta <= p["t1"] or p["r0"] == 0)
        for p in prisms
    )


def test_prisms_skip_empty_space():
    out = _run(f"""
const prisms = P.prismsInView([8000, 0, 3000], 1000, 32 * {EDGE_PC}, shape, {GALAXY_RADIUS_PC}, new Map());
console.log(JSON.stringify(prisms.length));
""")
    assert out == 0


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
    assert got == [4 * max(1, round(2 * math.pi * (i + 0.5) / 4)) for i in rings]
