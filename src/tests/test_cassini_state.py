"""
GEN.1536: the first Cassini state of a locked or resonant planet from its real precession constant (`planetgen.physics.cassini_state`): Mercury's MESSENGER figure against the nodal regression of its
orbit gives its measured tilt, the hydrostatic figure is the one Table 2 of [R5243] prints, and the fossil figure of a cold despun body lies between the Moon's and Mercury's.
"""

import math
from types import SimpleNamespace

import pytest

from planetgen.physics import bulk_draw, cassini_state as cs, constants, obliquity_variation as ov

MERCURY_MASS_KG = 3.3011e23
MERCURY_RADIUS_M = 2439.7e3
MERCURY_N = 2.0 * math.pi / (87.969 * 86400.0)
MERCURY_E = 0.2056
MERCURY_MODE = ov.Mode(-2.298e-5, math.radians(8.6))           # the regression of Mercury's node and its inclination to the Laplace plane [R5242]


def test_mercurys_hydrostatic_figure_is_the_one_the_table_prints():
    j2, c22 = cs.hydrostatic_figure(MERCURY_MASS_KG, MERCURY_RADIUS_M, 0.346, 1.5 * MERCURY_N, MERCURY_N, MERCURY_E)
    assert j2 == pytest.approx(0.534e-6, rel=0.08)        # [R5243: Table 2], <J2> of Mercury, with its 3:2 spin
    assert c22 == pytest.approx(0.3 * j2)


def test_mercurys_real_figure_holds_its_axis_to_the_orbit_at_the_measured_tilt():
    precession = cs.precession_rad_yr(5.031e-5, 8.088e-6, 0.346, MERCURY_N, MERCURY_E, True)
    assert precession * 206264.8 == pytest.approx(1000.0, rel=0.02)        # arcsec per year, against the node's 4.7
    tilt = cs.tilt_deg([MERCURY_MODE], precession) * 60.0
    assert tilt == pytest.approx(2.04, abs=0.6)           # arcmin, [R5242] measure 2.04


def test_the_fossil_bracket_is_the_moons_and_mercurys_and_a_relaxed_body_keeps_the_hydrostatic_figure():
    j2_low, c22_low = cs.frozen_figure(1.0, False, 0.0)
    j2_high, _ = cs.frozen_figure(1.0, False, 1.0)
    assert j2_low == pytest.approx(22.7) and j2_high == pytest.approx(94.2) and c22_low == pytest.approx(cs.FOSSIL_C22_OVER_J2 * 22.7)
    assert cs.frozen_figure(1.0, True, 0.7) == (1.0, cs.HYDROSTATIC_C22_OVER_J2)


def test_the_tilt_falls_as_the_precession_rises_and_is_the_inclination_where_it_vanishes():
    modes = [MERCURY_MODE]
    assert cs.tilt_deg([], 1.0) == 0.0 and cs.tilt_deg([ov.Mode(0.0, 0.1)], 1.0) == 0.0
    assert cs.tilt_deg(modes, 0.0) == pytest.approx(8.6)
    assert cs.tilt_deg(modes, 1e-4) > cs.tilt_deg(modes, 1e-3) > cs.tilt_deg(modes, 1e-2)
    two = cs.tilt_deg([MERCURY_MODE, MERCURY_MODE], 1e-3)
    assert two == pytest.approx(math.sqrt(2.0) * cs.tilt_deg(modes, 1e-3))      # the modes add in quadrature


def test_a_one_to_one_body_has_no_c22_term_and_a_resonant_one_has():
    locked = cs.precession_rad_yr(5e-5, 8e-6, 0.346, MERCURY_N, MERCURY_E, False)
    resonant = cs.precession_rad_yr(5e-5, 8e-6, 0.346, MERCURY_N, MERCURY_E, True)
    assert resonant > locked > 0.0


def test_an_icy_body_and_a_hot_rocky_one_relax_and_a_cold_rocky_one_does_not():
    icy = SimpleNamespace(star=SimpleNamespace(snow_line_au=2.0, mass=constants.SOLAR_MASS_TO_KG), formation_distance_au=5.0, distance=5.0)
    rocky = SimpleNamespace(star=SimpleNamespace(snow_line_au=2.0, mass=constants.SOLAR_MASS_TO_KG), formation_distance_au=0.4, distance=0.4)
    assert bulk_draw.origin_of(icy) == bulk_draw.ICY
    assert cs.is_relaxed(icy, 0.0) and cs.is_relaxed(rocky, 2.5) and not cs.is_relaxed(rocky, 0.05)


def test_the_hydrostatic_figure_needs_a_body():
    with pytest.raises(ValueError):
        cs.hydrostatic_figure(0.0, MERCURY_RADIUS_M, 0.346, MERCURY_N, MERCURY_N, 0.0)


JUPITER = SimpleNamespace(mass=1.898e27, radius=71492.0, period=11.86, rotation_period_hours=9.925)
EARTH = SimpleNamespace(mass=5.972e24, radius=6371.0, period=1.0, rotation_period_hours=23.934)


def _moon(planet, mass, radius, a_km, inclination, formation_au):
    return SimpleNamespace(mass=mass, radius=radius, distance=a_km * 1e3 / constants.AU_TO_M, eccentricity=0.0, orbital_inclination_deg=inclination, body_type="t",
                           star=SimpleNamespace(snow_line_au=2.7, mass=constants.SOLAR_MASS_TO_KG), formation_distance_au=formation_au)


def test_the_fluid_figure_of_jupiter_is_an_upper_bound_of_its_j2_and_the_moons_node_regresses_in_18_6_years():
    j2, _ = cs.hydrostatic_figure(JUPITER.mass, JUPITER.radius * 1000.0, 0.254, 2.0 * math.pi / (JUPITER.rotation_period_hours * 3600.0), 0.0, 0.0)
    assert 0.0147 < j2 < 1.6 * 0.0147        # the Darwin-Radau figure of a giant is high (0.0228), the measured one is 0.0147
    n = 2.0 * math.pi / (27.32 * 86400.0)
    n_earth = 2.0 * math.pi / (365.25 * 86400.0)
    s = cs.moon_regression_rad_yr(n, 384400e3, 1.08e-3, 6371e3, n_earth)
    assert 2.0 * math.pi / s == pytest.approx(18.6, rel=0.06)        # the period of the node of the Moon's orbit, years: the lowest-order solar term gives 17.8


@pytest.mark.parametrize("name, planet, mass, radius, a_km, period_d, inclination, formation_au, flux, low, high", [
    ("Io", JUPITER, 8.93e22, 1821.6, 421800, 1.769, 0.04, 5.2, 2.5, 0.0, 0.04),
    ("Europa", JUPITER, 4.8e22, 1560.8, 671100, 3.551, 0.47, 5.2, 0.1, 0.05, 0.25),
    ("Ganymede", JUPITER, 1.48e23, 2634.1, 1070400, 7.155, 0.2, 5.2, 0.05, 0.0, 0.35),
    ("Callisto", JUPITER, 1.08e23, 2410.3, 1882700, 16.69, 0.19, 5.2, 0.01, 0.0, 0.4),
    ("Moon", EARTH, 7.346e22, 1737.4, 384400, 27.32, 5.145, 1.0, 0.02, 3.0, 6.7)])
def test_the_galilean_moons_and_the_moon_come_out_at_the_tilts_they_have(name, planet, mass, radius, a_km, period_d, inclination, formation_au, flux, low, high):
    """GEN.1556: the measured tilts are Io under 0.04 degrees, Europa 0.1, Ganymede and Callisto a few tenths at most, the Moon 6.7 (its Cassini state 2, which the first state only approaches)."""
    moon = _moon(planet, mass, radius, a_km, inclination, formation_au)
    for u in (0.0, 0.5, 1.0):
        assert low <= cs.locked_moon_tilt_deg(moon, planet, flux, period_d / 365.25, u) <= high


def test_a_moon_without_a_regression_or_an_inclination_has_no_tilt():
    with pytest.raises(ValueError):
        cs.moon_regression_rad_yr(0.0, 1e8, 1e-3, 6e6, 1e-7)
    assert cs.locked_moon_tilt_deg(_moon(EARTH, 7.346e22, 1737.4, 384400, 0.0, 1.0), EARTH, 0.02, 27.32 / 365.25, 0.5) == 0.0
