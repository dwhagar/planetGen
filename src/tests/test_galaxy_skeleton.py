# tests/test_galaxy_skeleton.py

"""
Tests for `stellarObjects.galaxySkeleton` -- the per-ring layer-band
finder the galaxy-wide density skeleton is built from. See that module's
docstring for the reasoning: a ring's qualifying region is found via an
exact upper bound over `theta`, which falls monotonically with `|z|`, so
it is one band of layers centered on the plane.

The headline property is soundness: every sector whose exact
`galaxyDensity.relative_density` clears the threshold must fall inside its
ring's band. A band that is too wide only costs a few extra live checks.
"""

import math

import pytest

from stellarObjects.galaxyDensity import build_galaxy_shape, relative_density
from stellarObjects.galaxyGeometry import ring_radius_pc, ring_sector_count, sector_position_pc
from stellarObjects.galaxySkeleton import (
    RingBand,
    bound_relative_density_at,
    build_ring_bands,
    expected_system_count_at_density_1,
    find_ring_band,
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
EDGE_PC = 3.526


def test_expected_system_count_at_density_1_matches_spacesector():
    from stellarObjects.spaceSector import SpaceSector
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


def test_bound_is_never_below_the_exact_density():
    for ring in (0, 3, 15, 40):
        r = ring_radius_pc(ring, EDGE_PC)
        for z in (0.0, 4.0, -9.0):
            bound = bound_relative_density_at(SHAPE, r, z)
            for k in range(0, 360, 15):
                theta = math.radians(k)
                exact = relative_density((r * math.cos(theta), r * math.sin(theta), z), SHAPE)
                assert exact <= bound * (1 + 1e-12)


@pytest.mark.parametrize("ring_index, threshold_rho", [
    (0, 0.05), (5, 0.05), (5, 0.2), (15, 0.05), (15, 0.15), (30, 0.02), (60, 0.01),
])
def test_find_ring_band_never_excludes_a_true_qualifier(ring_index, threshold_rho):
    band = find_ring_band(SHAPE, EDGE_PC, ring_index, threshold_rho)
    for layer in range(-40, 41):
        for slot in range(ring_sector_count(ring_index)):
            position = sector_position_pc(ring_index, layer, slot, EDGE_PC)
            if relative_density(position, SHAPE) >= threshold_rho:
                assert band is not None and band.layer_index_min <= layer <= band.layer_index_max, (
                    f"({ring_index}, {layer}, {slot}) qualifies but falls outside {band}"
                )


def test_find_ring_band_is_symmetric_and_tight_at_its_edge():
    band = find_ring_band(SHAPE, EDGE_PC, 10, 0.05)
    assert isinstance(band, RingBand)
    assert band.layer_index_min == -band.layer_index_max
    r = ring_radius_pc(10, EDGE_PC)
    assert bound_relative_density_at(SHAPE, r, band.layer_index_max * EDGE_PC) >= 0.05
    assert bound_relative_density_at(SHAPE, r, (band.layer_index_max + 1) * EDGE_PC) < 0.05


def test_find_ring_band_empty_far_beyond_the_edge():
    assert find_ring_band(SHAPE, EDGE_PC, 5000, threshold_rho=1e-6) is None


def test_find_ring_band_respects_max_layer():
    band = find_ring_band(SHAPE, EDGE_PC, 0, 1e-9, max_layer=3)
    assert band == RingBand(-3, 3)


def test_build_ring_bands_stops_after_an_empty_streak():
    bands, outer, confirmed = build_ring_bands(SHAPE, EDGE_PC, 0.05, empty_streak_to_stop=5)
    assert confirmed
    assert bands and bands[0][0] == 0
    assert outer == bands[-1][0]
    rings = [ring for ring, _lo, _hi in bands]
    assert rings == sorted(rings)
    for ring, lo, hi in bands:
        assert (lo, hi) == tuple(find_ring_band(SHAPE, EDGE_PC, ring, 0.05))
    for ring in range(outer + 1, outer + 6):
        assert find_ring_band(SHAPE, EDGE_PC, ring, 0.05) is None


def test_build_ring_bands_reports_an_unconfirmed_edge_when_capped():
    bands, outer, confirmed = build_ring_bands(SHAPE, EDGE_PC, 1e-12, max_ring=10)
    assert not confirmed
    assert outer == 10
    assert len(bands) == 11


def test_build_ring_bands_with_nothing_qualifying():
    bands, outer, confirmed = build_ring_bands(SHAPE, EDGE_PC, 1e9, empty_streak_to_stop=3)
    assert (bands, outer, confirmed) == ([], -1, True)
