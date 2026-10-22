# tests/test_gen_grid_seams.py

"""
`docs/TODO.md` TEST.30 (grid seams and the nucleus), on the real sector
grid (`tuning.DEFAULT_SECTOR_EDGE_PC`, 4 pc), now the spherical one:

- the outermost planned shell and band of a real plan (the Milky-Way
  defaults `planetgen plan` uses, and a small toy shape) against
  `GalaxyBounds`: the last shell and each shell's highest band are inside,
  one more is not, and the faces between them fall on the documented side;
- the stored plan: `planetgen plan` stores exactly the outline
  `build_shell_extents` computes;
- the core (shell 0, band 0, wedge 0): it holds the origin and every point
  on the axis out to the first cap, and its neighbors are the first shell's
  cells.

The angular and face seams of the grid itself (the `2*pi` wrap, one ulp
inside every face, the caps and the core) are pinned in
`test_galaxy_spherical.py`; `test_fuzz_galaxy_grid.py` fuzzes the same
functions over arbitrary edges.
"""

import math
import random

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import spherical
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.skeleton import (
    GalaxyBounds, bound_relative_density_at, build_shell_extents, expected_system_count_at_density_1,
)

from tests.bughunt_support import mysql_argv, run_cli

EDGE = float(tuning.DEFAULT_SECTOR_EDGE_PC)
THRESHOLD = 1.0 / expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY)

# planetgen plan's own defaults (a real Milky-Way-scale outline, ~3,800
# shells) and test_galaxy_gen.py's small toy shape.
_MILKY_WAY = dict(disk_scale_length_pc=2800.0, disk_scale_height_pc=350.0, bulge_scale_radius_pc=200.0,
                  bulge_amplitude=1.0, arm_count=2, pitch_angle_deg=15.0, arm_amplitude=0.4)
_TOY = dict(disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
            bulge_amplitude=2.0, arm_count=2, pitch_angle_deg=15.0, arm_amplitude=0.4)


def _shape(params):
    params = dict(params)
    return build_galaxy_shape(pitch_angle_rad=math.radians(params.pop("pitch_angle_deg")), **params)


_PLANS = {}


def _plan(name):
    if name not in _PLANS:
        shape = _shape(_MILKY_WAY if name == "milky_way" else _TOY)
        extents, outer, confirmed = build_shell_extents(shape, EDGE, THRESHOLD)
        assert confirmed and extents
        _PLANS[name] = shape, extents, outer
    return _PLANS[name]


# ---------------------------------------------------------------------------
# The outermost planned shell and band against GalaxyBounds
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("plan", ["milky_way", "toy"])
def test_each_shells_highest_band_is_inside_and_the_next_is_not(plan):
    shape, extents, outer = _plan(plan)
    bounds = GalaxyBounds(extents, EDGE)
    assert bounds.outer_shell_index == outer == extents[-1][0]
    for shell, high in extents[::max(1, len(extents) // 400)]:
        top = spherical.band_range(shell)[1]
        assert bounds.qualifies(shell, high) and bounds.qualifies(shell, -high) and bounds.qualifies(shell, 0)
        assert bounds.qualifies(shell, 0)
        plane, z = spherical.band_center_pc(shell, high, EDGE)
        assert bound_relative_density_at(shape, plane, z) >= THRESHOLD
        if high < top:
            assert not bounds.qualifies(shell, high + 1) and not bounds.qualifies(shell, -high - 1)
            plane, z = spherical.band_center_pc(shell, high + 1, EDGE)
            assert bound_relative_density_at(shape, plane, z) < THRESHOLD or high + 1 == top
            assert bounds.describe_miss(shell, high + 1) == \
                f"band {high + 1} is outside the galaxy (shell {shell} reaches bands -{high} to {high})"


@pytest.mark.parametrize("plan", ["milky_way", "toy"])
def test_the_shell_past_the_outermost_is_outside(plan):
    _shape_, extents, outer = _plan(plan)
    bounds = GalaxyBounds(extents, EDGE)
    assert not bounds.qualifies(outer + 1, 0)
    assert bounds.describe_miss(outer + 1, 0) == \
        f"shell {outer + 1} is outside the galaxy (shells run from 0 to {outer})"


@pytest.mark.parametrize("plan", ["milky_way", "toy"])
def test_outer_face_of_the_outermost_shell_belongs_to_the_shell_outside(plan):
    """Shells are half-open `[k - 1/2, k + 1/2) * edge`: a point just inside
    the outermost shell's outer face is planned, the face itself is the next
    shell and outside the galaxy."""
    _shape_, extents, outer = _plan(plan)
    bounds = GalaxyBounds(extents, EDGE)
    _inner, outer_face = spherical.shell_bounds_pc(outer, EDGE)

    inside = (math.nextafter(outer_face, 0.0), 0.0, 0.0)
    shell, band, wedge = spherical.sector_address_at(inside, EDGE)
    assert (shell, band) == (outer, 0)
    assert bounds.qualifies(shell, band)

    on_face = (outer_face, 0.0, 0.0)
    shell, band, _wedge = spherical.sector_address_at(on_face, EDGE)
    assert (shell, band) == (outer + 1, 0)
    assert not bounds.qualifies(shell, band)


@pytest.mark.parametrize("plan", ["milky_way", "toy"])
def test_random_addresses_reach_but_never_pass_the_outermost_shell(plan):
    _shape_, extents, outer = _plan(plan)
    bounds = GalaxyBounds(extents, EDGE)
    rng = random.Random(30)
    for max_shell in (outer, outer + 1, outer - 1, 0):
        for _ in range(200):
            shell, band, wedge = bounds.random_address(rng, max_shell=max_shell)
            assert bounds.qualifies(shell, band)
            assert shell <= min(outer, max_shell)
            assert 0 <= wedge < spherical.wedge_count(shell, band)


def test_stored_plan_bounds_match_the_computed_outline(mysql_config):
    """`planetgen plan` stores exactly the outline `build_shell_extents`
    computes, and `store.get_galaxy_bounds` hands it back on the same
    4 pc grid, outermost shell and band included."""
    _shape_, extents, outer = _plan("toy")
    argv = []
    for key, value in _TOY.items():
        if key != "arm_amplitude":
            flag = "--disk-scale-length" if key == "disk_scale_length_pc" else f"--{key.replace('_', '-')}"
            argv += [flag, f"{value}pc" if key == "disk_scale_length_pc" else str(value)]
    argv += ["--arm-density", "1.4", "--interarm-density", "0.6"]  # the arm amplitude 0.4 in the density terms
    run_cli("plan", argv + ["--no-bright-stars"] + mysql_argv(mysql_config))

    conn = store.get_connection(mysql_config)
    try:
        bounds = store.get_galaxy_bounds(conn)
        stored = store.get_galaxy_shells(conn)
    finally:
        conn.close()
    assert bounds.edge_pc == EDGE
    # The plan draws the globular clusters, whose bands extend the disk's outline (GEN.160), so the
    # expected outline is walked on the shape the plan stored, clusters and all.
    conn = store.get_connection(mysql_config)
    try:
        clustered = store.get_galaxy_shape(conn).shape
    finally:
        conn.close()
    assert clustered.components is not None
    extents, outer, confirmed = build_shell_extents(clustered, EDGE, THRESHOLD)
    assert confirmed
    assert bounds.high_band == dict(extents)
    assert bounds.outer_shell_index == outer
    assert [row[:2] for row in stored] == extents
    assert [row[2] for row in stored] == [spherical.cells_in_band_range(shell, high) for shell, high in extents]
    assert not bounds.qualifies(outer + 1, 0)


# ---------------------------------------------------------------------------
# The core: shell 0, band 0, wedge 0
# ---------------------------------------------------------------------------

CORE = (0, 0, 0)


@pytest.mark.parametrize("point", [(0.0, 0.0, 0.0), (-0.0, 0.0, -0.0), (1.9, 0.0, 0.0), (0.0, -1.9, 0.0),
                                   (0.0, 0.0, 1.9), (0.0, 0.0, -1.9), (1.0, 1.0, 1.0)])
def test_points_within_half_an_edge_of_the_origin_are_in_the_core(point):
    assert spherical.sector_address_at(point, EDGE) == CORE


def test_the_core_ends_at_half_an_edge():
    assert spherical.sector_address_at((EDGE / 2, 0.0, 0.0), EDGE)[0] == 1
    assert spherical.sector_address_at((math.nextafter(EDGE / 2, 0.0), 0.0, 0.0), EDGE) == CORE


@pytest.mark.parametrize("z", [2.0 + 1e-9, 4.0, 5.9, -4.0, -5.9])
def test_axis_points_outside_the_core_fall_in_a_pole_cap_of_the_first_shells(z):
    shell, band, wedge = spherical.sector_address_at((0.0, 0.0, z), EDGE)
    assert shell == 1 and wedge == 0
    assert band == (spherical.band_range(1)[1] if z > 0 else spherical.band_range(1)[0])


def test_the_core_is_the_neighbor_of_exactly_the_first_shells_cells_that_touch_it():
    neighbors = spherical.neighbor_addresses(*CORE)
    assert neighbors and all(shell == 1 for shell, _band, _wedge in neighbors)
    for neighbor in neighbors:
        assert CORE in spherical.neighbor_addresses(*neighbor)


def test_the_core_designation_round_trips_and_is_distinct():
    designation = spherical.provisional_sector_designation(*CORE)
    assert spherical.parse_sector_designation(designation) == CORE
    others = {spherical.provisional_sector_designation(*cell) for cell in spherical.neighbor_addresses(*CORE)}
    assert designation not in others


@pytest.mark.parametrize("plan", ["milky_way", "toy"])
def test_every_real_plan_holds_the_core(plan):
    _shape_, extents, _outer = _plan(plan)
    bounds = GalaxyBounds(extents, EDGE)
    assert extents[0] == (0, 0)
    assert bounds.qualifies(0, 0)
    assert bounds.cells_in_shell(0) == 1
