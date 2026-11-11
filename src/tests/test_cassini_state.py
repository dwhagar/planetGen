"""
GEN.1536: the first Cassini state of a locked or resonant planet from its real precession constant (`planetgen.physics.cassini_state`): Mercury's MESSENGER figure against the nodal regression of its
orbit gives its measured tilt, the hydrostatic figure is the one Table 2 of [R5228] prints, and the fossil figure of a cold despun body lies between the Moon's and Mercury's.
"""

import math
from types import SimpleNamespace

import pytest

from planetgen.physics import bulk_draw, cassini_state as cs, constants, obliquity_variation as ov

MERCURY_MASS_KG = 3.3011e23
MERCURY_RADIUS_M = 2439.7e3
MERCURY_N = 2.0 * math.pi / (87.969 * 86400.0)
MERCURY_E = 0.2056
MERCURY_MODE = ov.Mode(-2.298e-5, math.radians(8.6))           # the regression of Mercury's node and its inclination to the Laplace plane ([R5229])


def test_mercurys_hydrostatic_figure_is_the_one_the_table_prints():
    j2, c22 = cs.hydrostatic_figure(MERCURY_MASS_KG, MERCURY_RADIUS_M, 0.346, 1.5 * MERCURY_N, MERCURY_N, MERCURY_E)
    assert j2 == pytest.approx(0.534e-6, rel=0.08)        # [R5228] Table 2, <J2> of Mercury, with its 3:2 spin
    assert c22 == pytest.approx(0.3 * j2)


def test_mercurys_real_figure_holds_its_axis_to_the_orbit_at_the_measured_tilt():
    precession = cs.precession_rad_yr(5.031e-5, 8.088e-6, 0.346, MERCURY_N, MERCURY_E, True)
    assert precession * 206264.8 == pytest.approx(1000.0, rel=0.02)        # arcsec per year, against the node's 4.7
    tilt = cs.tilt_deg([MERCURY_MODE], precession) * 60.0
    assert tilt == pytest.approx(2.04, abs=0.6)           # arcmin, [R5229] measure 2.04


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
    icy = SimpleNamespace(star=SimpleNamespace(snow_line_au=2.0), formation_distance_au=5.0, distance=5.0)
    rocky = SimpleNamespace(star=SimpleNamespace(snow_line_au=2.0), formation_distance_au=0.4, distance=0.4)
    assert bulk_draw.origin_of(icy) == bulk_draw.ICY
    assert cs.is_relaxed(icy, 0.0) and cs.is_relaxed(rocky, 2.5) and not cs.is_relaxed(rocky, 0.05)


def test_the_hydrostatic_figure_needs_a_body():
    with pytest.raises(ValueError):
        cs.hydrostatic_figure(0.0, MERCURY_RADIUS_M, 0.346, MERCURY_N, MERCURY_N, 0.0)
