# tests/test_king.py

"""
GEN.159: the [R2359] solver. Its concentrations match the textbook values, a cluster's expected number of systems is
its mass over the mean system mass, a radius drawn by inverse transform lands where the table says, and a cluster whose
King edge would pass its Jacobi radius is given a lower central potential.
"""

import math

import pytest

from planetgen import tuning
from planetgen.physics import king

TEXTBOOK_CONCENTRATIONS = {3: 0.67, 5: 1.03, 7: 1.53, 9: 2.12, 12: 2.74}


@pytest.mark.parametrize("w0, expected", TEXTBOOK_CONCENTRATIONS.items())
def test_the_solver_reproduces_the_textbook_concentration(w0, expected):
    assert king.king_table(w0).concentration == pytest.approx(expected, abs=0.01)


def test_every_table_in_the_range_is_monotonic_and_normalised():
    low, high = tuning.KING_W0_RANGE
    for w0 in range(low, high + 1):
        table = king.king_table(w0)
        assert len(table.x) == tuning.KING_TABLE_POINTS
        assert table.mass_fraction[0] == 0.0 and table.mass_fraction[-1] == 1.0
        assert all(b >= a for a, b in zip(table.mass_fraction, table.mass_fraction[1:]))
        assert table.density[0] == 1.0 and table.density[-1] < table.density[0]
        assert table.half_mass_x < table.tidal_x


def test_a_table_outside_the_range_is_refused():
    with pytest.raises(ValueError):
        king.king_table(tuning.KING_W0_RANGE[1] + 1)
    with pytest.raises(ValueError):
        king.king_table(5.5)


def test_the_expected_systems_are_the_mass_over_the_mean_system_mass():
    assert king.expected_systems(2.0e5, 0.4) == pytest.approx(5.0e5)
    assert king.expected_systems(8.0e5, 0.5) == pytest.approx(1.6e6)


@pytest.mark.parametrize("w0", (3, 7, 12))
def test_the_density_integrates_to_the_mass(w0):
    mass, half = 2.0e5, 5.0
    edge = king.tidal_radius_pc(w0, half)
    steps = 20000
    total = 0.0
    for i in range(steps):
        r = (i + 0.5) * edge / steps
        total += 4.0 * math.pi * r * r * king.density_per_pc3(w0, mass, half, r) * (edge / steps)
    assert total == pytest.approx(mass, rel=0.01)


@pytest.mark.parametrize("w0", (3, 7, 12))
def test_a_half_mass_radius_holds_half_the_mass(w0):
    assert king.mass_within_pc(w0, 6.0, 6.0) == pytest.approx(0.5, abs=0.002)
    assert king.radius_pc(w0, 6.0, 0.5) == pytest.approx(6.0, rel=0.002)


@pytest.mark.parametrize("fraction", (0.01, 0.1, 0.5, 0.9, 0.99))
def test_a_radius_drawn_by_inverse_transform_lands_on_the_tables_mass_share(fraction):
    r = king.radius_pc(7, 5.0, fraction)
    assert king.mass_within_pc(7, 5.0, r) == pytest.approx(fraction, abs=0.002)


def test_nothing_lies_beyond_the_tidal_radius():
    edge = king.tidal_radius_pc(7, 5.0)
    assert king.mass_within_pc(7, 5.0, edge * 2) == 1.0
    assert king.density_per_pc3(7, 1.0e5, 5.0, edge * 1.01) == 0.0
    assert king.radius_pc(7, 5.0, 1.0) == pytest.approx(edge, rel=1e-6)


def test_the_nearest_w0_inverts_a_concentration():
    for w0, concentration in TEXTBOOK_CONCENTRATIONS.items():
        assert king.nearest_w0(concentration) == w0
    assert king.nearest_w0(-1.0) == tuning.KING_W0_RANGE[0]
    assert king.nearest_w0(9.0) == tuning.KING_W0_RANGE[1]


def test_the_jacobi_radius_grows_with_mass_and_distance():
    # The galaxy's mass inside 8 and 16 kpc: 9e10 and 1.8e11 solar masses, an isothermal halo of 220 km/s (M = v^2 R / G).
    assert king.jacobi_radius_pc(2.0e5, 8000.0, 9.0e10) < king.jacobi_radius_pc(2.0e6, 8000.0, 9.0e10)
    assert king.jacobi_radius_pc(2.0e5, 8000.0, 9.0e10) < king.jacobi_radius_pc(2.0e5, 16000.0, 1.8e11)
    # 2e5 solar masses at 8 kpc keep about 80 pc, as globular clusters do (tens of parsecs); a heavier galaxy squeezes the cluster.
    assert 30.0 < king.jacobi_radius_pc(2.0e5, 8000.0, 9.0e10) < 120.0
    assert king.jacobi_radius_pc(2.0e5, 8000.0, 2.0e11) < king.jacobi_radius_pc(2.0e5, 8000.0, 9.0e10)


def test_limit_w0_lowers_the_potential_until_the_edge_fits():
    half = 5.0
    wide = king.tidal_radius_pc(12, half)
    fitted = king.limit_w0(12, half, wide / 2)
    assert fitted < 12
    assert king.tidal_radius_pc(fitted, half) <= wide / 2 or fitted == tuning.KING_W0_RANGE[0]
    assert king.limit_w0(7, half, 1.0e9) == 7
    assert king.limit_w0(12, half, 0.001) == tuning.KING_W0_RANGE[0]
