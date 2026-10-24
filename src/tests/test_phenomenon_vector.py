"""
The phenomena scatter draws a shell's bands as arrays (PERF.63): the array
density and factors equal the scalar ones, the shell holds as many objects
of each kind as the band-by-band means expect, they sit where those put them,
and the rows are the same shape as before.
"""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import remnant_density, density, remnant_distribution
from planetgen.galaxy import spherical
from planetgen.galaxy.spherical import sector_address_at, wedge_count
from planetgen.generation import bright_stars
from planetgen.generation import phenomenon_scatter as scatter
from planetgen.util import draw
from planetgen.physics.units import ly_to_pc
from tests.toy_shape import toy_shape

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
SHAPE = toy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
E_VALUE = 2000.0
SHELL = 6
HIGH_BAND = 3


def _shell(shell=SHELL, high_band=HIGH_BAND, seed=5, **kwargs):
    kwargs.setdefault("min_mass_solar", 1.0)
    return list(scatter.scatter_shell(SHAPE, shell, high_band, EDGE_PC, E_VALUE, seed, **kwargs))


def _points(count, seed=3):
    rng = draw.Stream(seed)
    return [(rng.uniform(-300, 300), rng.uniform(-300, 300), rng.uniform(-60, 60)) for _ in range(count)]


@pytest.mark.parametrize("bulge_amplitude, nuclear_amplitude", [(2.0, 0.0), (2.0, 1.0), (-1.0, 0.0), (-1.0, 3.0)])
def test_the_array_density_is_the_sum_of_the_scalar_populations(bulge_amplitude, nuclear_amplitude):
    shape = toy_shape(
        disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
        bulge_amplitude=bulge_amplitude, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
        nuclear_amplitude=nuclear_amplitude,
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


def test_every_row_sits_in_its_own_cell_with_a_mass_above_the_cut():
    rows = _shell(min_mass_solar=10.0)
    assert len(rows) > 100
    assert len({row[1] for row in rows}) > 1, "the shell's rows should fall in more than one band"
    for row in map(dict, (zip(scatter.PHENOMENON_SCATTER_COLUMNS, row) for row in rows)):
        point = (row["position_x_mpc"] / 1000, row["position_y_mpc"] / 1000, row["position_z_mpc"] / 1000)
        assert sector_address_at(point, EDGE_PC) == (row["shell_index"], row["band_index"], row["wedge_index"])
        assert row["shell_index"] == SHELL and abs(row["band_index"]) <= HIGH_BAND and 0 <= row["seed"] < 2 ** 63
        if row["kind"] == "black-hole":
            assert row["mass_solar"] >= 10.0
        elif row["kind"] == "neutron-star":
            # A magnetar ignores the mass cut (GEN.278): it follows the young stars, whatever the cut leaves.
            assert row["subtype"] == scatter.MAGNETAR_SUBTYPE or row["mass_solar"] is None or row["mass_solar"] >= 10.0
        else:
            assert row["mass_solar"] is None


def test_skipped_cells_get_no_object():
    rows = _shell()
    skip = {tuple(row[:3]) for row in rows[::3]}
    again = _shell(skip_addresses=skip)
    assert skip.isdisjoint(tuple(row[:3]) for row in again)
    assert len(again) == len(rows) - sum(1 for row in rows if tuple(row[:3]) in skip)


def test_a_shell_is_the_same_from_the_same_seed_and_not_from_another():
    assert _shell(seed=9) == _shell(seed=9)
    assert _shell(seed=9) != _shell(seed=10)


def test_a_shell_with_no_bands_to_draw_has_no_rows():
    assert _shell(high_band=-1) == []


def _expected_in_shell(kind, rate, shell, high_band, subtype=None):
    """What the band-by-band means add up to for one kind in one shell (the mean the arrays draw its Poisson
    counts from, worked out here with the scalar density and factors)."""
    bins = bright_stars.ANGLE_BINS
    total = 0.0
    for band in range(-high_band, high_band + 1):
        plane, z = spherical.band_center_pc(shell, band, EDGE_PC)
        weight = 0.0
        for k in range(bins):
            angle = (k + 0.5) * 2.0 * math.pi / bins
            point = (plane * math.cos(angle), plane * math.sin(angle), z)
            stars = sum(bright_stars._densities(point, SHAPE).values())
            # A stellar-mass black hole weighs the stars by the failed share of their metals (GEN.563), a neutron star by
            # what is left, bound (GEN.554).
            weight += (stars * remnant_density.place_factor(point, SHAPE, kind)
                       if remnant_density.follows_metals(kind, subtype) else stars)
        if kind in tuning.REMNANT_SCALE_HEIGHT_RATIO:
            weight *= (remnant_distribution.vertical_factor(kind, z, SHAPE.disk_scale_height_pc)
                       * remnant_distribution.radial_factor(kind, plane))
        total += (rate * E_VALUE * wedge_count(shell, band) * bright_stars.volume_ratio(shell, band, EDGE_PC)
                  / bins * weight)
    return total


def test_the_counts_match_the_band_by_band_means():
    """Per kind and shell, the shell holds what the bands' Poisson means add up to (Monte Carlo, 5 sigma)."""
    rates = scatter._scattered_rates(1.0)
    seeds = 12
    for shell in (4, 8):
        rows = [row for seed in range(seeds) for row in _shell(shell, seed=seed)]
        for kind, subtype, rate in rates:
            expected = _expected_in_shell(kind, rate, shell, HIGH_BAND, subtype) * seeds
            got = sum(1 for row in rows if row[3] == kind and row[4] == subtype)
            assert abs(got - expected) <= 5 * math.sqrt(max(expected, 1.0)) + 2, (kind, subtype, shell, got, expected)
