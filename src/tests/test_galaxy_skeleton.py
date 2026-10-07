# tests/test_galaxy_skeleton.py

"""
Tests for `planetgen.galaxy.skeleton` -- the per-layer outline the
galaxy-wide density skeleton is built from. See that module's docstring
for the reasoning: an exact upper bound over `theta` falls monotonically
with both `R` and `|z|`, so the galaxy is a stack of layers symmetric
about the plane, each running from ring 0 out to one outer ring.

The headline property is soundness: every sector whose exact
`galaxyDensity.relative_density` clears the threshold must fall inside its
layer's extent. An extent that is too wide only costs a few extra live
checks.
"""

import math

import pytest

from planetgen.galaxy.density import build_galaxy_shape, relative_density
from planetgen.galaxy.geometry import ring_radius_pc, ring_sector_count, sector_position_pc
from planetgen.galaxy.skeleton import (
    bound_relative_density_at,
    build_layer_extents,
    GalaxyBounds,
    candidate_sector_count,
    column_extents,
    expected_system_count_at_density_1,
)

SHAPE = build_galaxy_shape(
    disk_scale_length_pc=40.0,
    disk_scale_height_pc=12.0,
    bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0,
    arm_count=2,
    pitch_angle_rad=math.radians(15),
    arm_amplitude=0.4,
)
EDGE_PC = 4.0


def test_expected_system_count_at_density_1_matches_spacesector():
    from planetgen.galaxy.sector import SpaceSector
    expected = SpaceSector(name="x", edge_ly=11.5).expected_system_count()
    assert expected_system_count_at_density_1(11.5) == pytest.approx(expected)


def test_bound_is_symmetric_in_z_and_falls_with_height():
    r = ring_radius_pc(20, EDGE_PC)
    previous = bound_relative_density_at(SHAPE, r, 0.0)
    for z in (1.0, 5.0, 20.0, 60.0):
        assert bound_relative_density_at(SHAPE, r, z) == pytest.approx(bound_relative_density_at(SHAPE, r, -z))
        current = bound_relative_density_at(SHAPE, r, z)
        assert current < previous
        previous = current


def test_bound_falls_with_radius():
    previous = bound_relative_density_at(SHAPE, ring_radius_pc(0, EDGE_PC), 8.0)
    for ring in range(1, 60):
        current = bound_relative_density_at(SHAPE, ring_radius_pc(ring, EDGE_PC), 8.0)
        assert current < previous
        previous = current


def test_bound_is_never_below_the_exact_density():
    for ring in (0, 3, 15, 40):
        r = ring_radius_pc(ring, EDGE_PC)
        for z in (0.0, 4.0, -9.0):
            bound = bound_relative_density_at(SHAPE, r, z)
            for k in range(0, 360, 15):
                theta = math.radians(k)
                exact = relative_density((r * math.cos(theta), r * math.sin(theta), z), SHAPE)
                assert exact <= bound * (1 + 1e-12)


def _qualifies(ring, layer, threshold_rho):
    return bound_relative_density_at(SHAPE, ring_radius_pc(ring, EDGE_PC), layer * EDGE_PC) >= threshold_rho


@pytest.mark.parametrize("threshold_rho", [0.01, 0.05, 0.15, 0.5])
def test_layer_extents_never_exclude_a_true_qualifier(threshold_rho):
    extents = dict(build_layer_extents(SHAPE, EDGE_PC, threshold_rho)[0])
    for ring_index in range(0, 60, 3):
        for layer in range(-30, 31):
            for slot in range(0, ring_sector_count(ring_index), 2):
                position = sector_position_pc(ring_index, layer, slot, EDGE_PC)
                if relative_density(position, SHAPE) >= threshold_rho:
                    assert layer in extents and ring_index <= extents[layer], (
                        f"({ring_index}, {layer}, {slot}) qualifies but falls outside the layer's extent"
                    )


@pytest.mark.parametrize("threshold_rho", [0.01, 0.05, 0.15])
def test_layer_extents_match_a_brute_force_scan(threshold_rho):
    extents, outer, confirmed = build_layer_extents(SHAPE, EDGE_PC, threshold_rho)
    assert confirmed
    brute = {}
    for layer in range(-90, 91):
        rings = [ring for ring in range(0, 200) if _qualifies(ring, layer, threshold_rho)]
        if rings:
            assert rings == list(range(len(rings))), "a layer's rings must run unbroken from ring 0"
            brute[layer] = rings[-1]
    assert dict(extents) == brute
    assert outer == brute[0]


def test_layer_extents_run_top_to_bottom_symmetric_and_shrink_off_the_plane():
    extents, _outer, _confirmed = build_layer_extents(SHAPE, EDGE_PC, 0.05)
    layers = [layer for layer, _ring in extents]
    top = layers[0]
    assert layers == list(range(top, -top - 1, -1))
    rings = dict(extents)
    for layer in range(0, top):
        assert rings[layer] == rings[-layer]
        assert rings[layer + 1] <= rings[layer]


def test_layer_extents_are_tight_at_their_edges():
    extents, _outer, _confirmed = build_layer_extents(SHAPE, EDGE_PC, 0.05)
    top = extents[0][0]
    for layer, outer_ring in extents:
        assert _qualifies(outer_ring, layer, 0.05)
        assert not _qualifies(outer_ring + 1, layer, 0.05)
    assert not _qualifies(0, top + 1, 0.05)


def test_layer_extents_report_an_unconfirmed_edge_when_capped():
    extents, outer, confirmed = build_layer_extents(SHAPE, EDGE_PC, 1e-12, max_ring=10, max_layer=5)
    assert not confirmed
    assert outer == 10
    assert len(extents) == 11
    assert all(ring <= 10 for _layer, ring in extents)


def test_layer_extents_with_nothing_qualifying():
    assert build_layer_extents(SHAPE, EDGE_PC, 1e9) == ([], -1, True)


def test_candidate_sector_count_sums_every_ring_of_every_layer():
    extents = [(1, 2), (0, 3), (-1, 2)]
    by_ring = [ring_sector_count(r) for r in range(4)]
    assert candidate_sector_count(extents) == 2 * sum(by_ring[:3]) + sum(by_ring)
    assert candidate_sector_count([]) == 0


def test_column_extents_turn_the_layer_outline_on_its_side():
    extents, outer, _confirmed = build_layer_extents(SHAPE, EDGE_PC, 0.05)
    layers = dict(extents)
    columns = column_extents(extents)
    assert [ring for ring, _lo, _hi in columns] == list(range(outer + 1))
    for ring, lo, hi in columns:
        reached = [layer for layer, last in layers.items() if ring <= last]
        assert (lo, hi) == (min(reached), max(reached))
        assert lo == -hi


def test_galaxy_bounds_contains_exactly_the_outline():
    bounds = GalaxyBounds([(1, 1), (0, 3), (-1, 1)], EDGE_PC)
    assert (bounds.top_layer_index, bounds.outer_ring_index) == (1, 3)
    assert bounds.contains(3, 0) and not bounds.contains(4, 0)
    assert bounds.contains(1, -1) and not bounds.contains(2, -1)
    assert not bounds.contains(0, 2) and not bounds.contains(0, -2)
    assert "layer 0 ends at ring 3" in bounds.describe_miss(4, 0)
    assert "layers run from 1 down to -1" in bounds.describe_miss(0, 5)
    assert bounds.cell_count() == candidate_sector_count([(1, 1), (0, 3), (-1, 1)])
    assert not GalaxyBounds([], EDGE_PC)


def test_random_addresses_stay_inside_and_reach_every_ring_evenly_per_sector():
    import random
    from collections import Counter

    bounds = GalaxyBounds([(1, 1), (0, 3), (-1, 1)], EDGE_PC)
    rng = random.Random(7)
    draws = [bounds.random_address(rng) for _ in range(20000)]
    for ring, layer, slot in draws:
        assert bounds.contains(ring, layer)
        assert 0 <= slot < ring_sector_count(ring)
    # Every sector equally likely: each (ring, layer) gets draws in
    # proportion to its slot count.
    counts = Counter((ring, layer) for ring, layer, _slot in draws)
    total = bounds.cell_count()
    for (ring, layer), n in counts.items():
        assert n / len(draws) == pytest.approx(ring_sector_count(ring) / total, rel=0.15)
    assert len(counts) == 2 * 2 + 4

    capped = {bounds.random_address(rng, max_ring=0)[0] for _ in range(200)}
    assert capped == {0}
    assert GalaxyBounds([], EDGE_PC).random_address(rng) is None
