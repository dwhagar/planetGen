"""
The phenomena scatter draws a layer's rings as arrays (PERF.63): the array
density and factors equal the scalar ones, the layer holds as many objects
of each kind as the ring-by-ring walk expects, they sit where it puts them,
and the rows are the same shape as before.
"""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import density, remnant_distribution
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.geometry import ring_sector_count, sector_address_at
from planetgen.generation import bright_stars
from planetgen.generation import phenomenon_scatter as scatter
from planetgen.util import draw
from planetgen.physics.units import ly_to_pc

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
SHAPE = build_galaxy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
E_VALUE = 2000.0
OUTER = 14


def _layer(layer_index=0, seed=5, outer=OUTER, **kwargs):
    kwargs.setdefault("min_mass_solar", 1.0)
    return list(scatter.scatter_layer(SHAPE, layer_index, outer, EDGE_PC, E_VALUE, seed, **kwargs))


def _points(count, seed=3):
    rng = draw.Stream(seed)
    return [(rng.uniform(-300, 300), rng.uniform(-300, 300), rng.uniform(-60, 60)) for _ in range(count)]


@pytest.mark.parametrize("bulge_amplitude, core_amplitude", [(2.0, 0.0), (2.0, 5.0), (-1.0, 0.0), (-1.0, -2.0)])
def test_the_array_density_is_the_sum_of_the_scalar_populations(bulge_amplitude, core_amplitude):
    shape = build_galaxy_shape(
        disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
        bulge_amplitude=bulge_amplitude, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
        core_amplitude=core_amplitude,
    )
    for z in (0.0, 7.0, -40.0, 5000.0):
        points = [(x, y, z) for x, y, _z in _points(40)]
        got = density.total_density_array(np.array([p[0] for p in points]), np.array([p[1] for p in points]), z, shape)
        want = [sum(bright_stars._densities(point, shape).values()) for point in points]
        assert got == pytest.approx(want, rel=1e-9, abs=1e-12)


def test_the_array_radial_factor_matches_the_scalar_one():
    radii = [0.0, 3.0, 40.0, 900.0, 25000.0]
    for kind in ("black-hole", "neutron-star", "planetary-nebula"):
        assert remnant_distribution.radial_factor_array(kind, radii).tolist() == pytest.approx(
            [remnant_distribution.radial_factor(kind, radius) for radius in radii])


def test_the_first_vector_ring_is_the_first_with_every_angle_bin():
    first = scatter._first_full_ring()
    assert ring_sector_count(first) >= bright_stars.ANGLE_BINS > ring_sector_count(first - 1)


def test_every_row_sits_in_its_own_cell_with_a_mass_above_the_cut():
    rows = _layer(0, min_mass_solar=10.0)
    assert len(rows) > 200
    assert max(row[0] for row in rows) > scatter._first_full_ring()
    for row in map(dict, (zip(scatter.PHENOMENON_SCATTER_COLUMNS, row) for row in rows)):
        point = (row["position_x_mpc"] / 1000, row["position_y_mpc"] / 1000, row["position_z_mpc"] / 1000)
        assert sector_address_at(point, EDGE_PC) == (row["ring_index"], row["layer_index"], row["ring_slot_index"])
        assert row["ring_index"] <= OUTER and 0 <= row["seed"] < 2 ** 63
        if row["kind"] == "black-hole":
            assert row["mass_solar"] >= 10.0
        elif row["kind"] == "neutron-star":
            assert row["mass_solar"] is None or row["mass_solar"] >= 10.0
        else:
            assert row["mass_solar"] is None


def test_skipped_cells_get_no_object():
    rows = _layer(0)
    skip = {tuple(row[:3]) for row in rows[::3]}
    again = _layer(0, skip_addresses=skip)
    assert skip.isdisjoint(tuple(row[:3]) for row in again)
    assert len(again) == len(rows) - sum(1 for row in rows if tuple(row[:3]) in skip)


def test_a_layer_is_the_same_from_the_same_seed_and_not_from_another():
    assert _layer(1, seed=9) == _layer(1, seed=9)
    assert _layer(1, seed=9) != _layer(1, seed=10)


def test_a_layer_inside_the_first_rings_uses_only_the_ring_walk():
    rows = _layer(0, outer=2)
    assert rows and max(row[0] for row in rows) <= 2


def test_the_counts_match_the_ring_by_ring_walk():
    """Per kind and per ring group, the layer holds what the walk's Poisson means add up to (Monte Carlo, 5 sigma)."""
    rates = scatter._scattered_rates(1.0)
    first = scatter._first_full_ring()
    seeds = 12
    rows = [row for seed in range(seeds) for row in _layer(0, seed=seed)]
    for kind, subtype, rate in rates:
        for low, high in ((first, 8), (9, OUTER)):
            expected = 0.0
            for ring in range(low, high + 1):
                slots, bins = bright_stars._ring_bins(ring, 0, SHAPE, E_VALUE, EDGE_PC)
                base = [sum(d.values()) for d in bins]
                weights = scatter._kind_weights(kind, base, ring, 0, EDGE_PC, SHAPE.disk_scale_height_pc)
                expected += rate * E_VALUE * slots / len(bins) * sum(weights)
            expected *= seeds
            got = sum(1 for row in rows if row[3] == kind and row[4] == subtype and low <= row[0] <= high)
            assert abs(got - expected) <= 5 * math.sqrt(max(expected, 1.0)) + 2, (kind, subtype, low, high, got, expected)
