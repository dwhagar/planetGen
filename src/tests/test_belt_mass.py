# tests/test_belt_mass.py

"""GEN.477 E2: the mass left in a planetesimal belt after collisional grinding (`physics/belt_mass.py`; Wyatt et al. 2007)."""

import math

import pytest

from planetgen.physics import belt_mass as bm, constants, formation
from planetgen.util import draw

SUN_KG = constants.SOLAR_MASS_TO_KG
SUN_W = 3.828e26
KRASINSKY_2002_EARTH = 0.6e-3      # the inferred mass of the Solar System's belt, as Wyatt et al. 2007 quote it: a test point only


def disc_mass(lower, upper, star_kg=SUN_KG):
    return bm.primordial_mass_earth(lower, upper, star_kg, formation.snow_line_au(SUN_W))


def test_wyatts_equation_17_is_reproduced_at_its_stated_assumptions():
    """Eq. 17: M_max = 0.009 r^3.5 (dr/r) Dc M*^-0.5 t^-1 / G for rho = 2700 and e = I; with `STRENGTH_FLOOR_J_PER_KG` as the strength and Dc = 1 km."""
    r, dr, dc, age = 10.0, 5.0, 1.0, 100.0
    lower, upper = r - dr / 2, r + dr / 2
    a = bm.collisional_constant(lower, upper, 1.0, dc)
    e = (upper - lower) / (upper + lower)
    xc = 1e-3 * (r * bm.STRENGTH_FLOOR_J_PER_KG / e ** 2) ** (1 / 3)
    g = xc ** -0.5 + 0.67 * xc ** -1.5 + 0.2 * xc ** -2.5 - 1.87
    assert a / age == pytest.approx(0.009 * r ** 3.5 * (dr / r) * dc / age / g, rel=0.02)


def test_a_young_belt_keeps_what_formed_and_an_old_one_keeps_the_same_whatever_formed():
    light = bm.belt_mass_earth(1e-3, 2.2, 3.4, 1.0, 4500.0, 500.0)
    heavy = bm.belt_mass_earth(1e3, 2.2, 3.4, 1.0, 4500.0, 500.0)
    assert heavy == pytest.approx(bm.collisional_constant(2.2, 3.4, 1.0, 500.0) / 4500.0, rel=0.01)
    assert light < heavy < 1e-3
    assert bm.belt_mass_earth(1.0, 2.2, 3.4, 1.0, 0.0, 500.0) == pytest.approx(1.0)
    assert bm.belt_mass_earth(1.0, 2.2, 3.4, 1.0, 1.0, 500.0) > 0.4


def test_the_mass_falls_as_one_over_the_age_late_on_and_never_rises():
    ages = [1e2, 1e3, 4.5e3, 1e4]
    values = [bm.belt_mass_earth(7.0, 2.2, 3.4, 1.0, t, 500.0) for t in ages]
    assert values == sorted(values, reverse=True)
    assert values[-2] / values[-1] == pytest.approx(10.0 / 4.5, rel=0.05)


def test_the_solar_systems_belt_comes_out_within_a_few_times_the_inferred_mass():
    """A test and not an input: 2.2 to 3.4 AU, 4.5 Gyr, Ceres-size bodies at 500 to 1000 km, a primordial disc of a few Earth masses."""
    primordial = disc_mass(2.2, 3.4)
    assert 1.0 < primordial < 20.0
    for largest in (500.0, 1000.0):
        left = bm.belt_mass_earth(primordial, 2.2, 3.4, 1.0, 4500.0, largest)
        assert KRASINSKY_2002_EARTH / 5 < left < KRASINSKY_2002_EARTH * 5
    assert primordial / bm.belt_mass_earth(primordial, 2.2, 3.4, 1.0, 4500.0, 500.0) > 1000.0


def test_the_disruption_threshold_is_the_gravity_bound_one_for_big_bodies_and_the_strength_floor_for_small():
    assert bm.disruption_threshold_j_per_kg(1.0) == bm.STRENGTH_FLOOR_J_PER_KG
    assert bm.disruption_threshold_j_per_kg(1000.0) > 1e5
    assert bm.disruption_threshold_j_per_kg(1000.0) > bm.disruption_threshold_j_per_kg(100.0)


def test_a_heavier_star_leaves_less_and_a_wider_belt_more_per_the_law():
    base = bm.collisional_constant(2.2, 3.4, 1.0, 500.0)
    assert bm.collisional_constant(2.2, 3.4, 4.0, 500.0) == pytest.approx(base / 2.0)
    assert bm.collisional_constant(2.0, 3.6, 1.0, 500.0) != pytest.approx(base)


def test_the_density_word_follows_the_mass_per_area_against_its_reference():
    assert bm.density_word(10.0, 1.0) == "dense"
    assert bm.density_word(0.1, 1.0) == "sparse"
    assert bm.density_word(1.0, 1.0) == "typical"
    assert bm.density_word(bm.DENSE_FACTOR, 1.0) == "dense"
    assert math.isclose(bm.surface_density_earth_per_au2(math.pi * (9 - 1), 1.0, 3.0), 1.0)


def test_a_belt_state_gives_the_mass_in_kilograms_and_the_word_against_the_median_largest_body():
    """GEN.1000: the largest body comes from the size law and is settled with the mass (`belt_ledger.belt_state`)."""
    from planetgen.physics import belt_ledger as bl

    line = formation.snow_line_au(SUN_W)
    states = [bl.belt_state(2.2, 3.4, SUN_KG, line, 1.0, 4.5, bl.ROCKY_LAW, draw.Stream(f"belt:{i}")) for i in range(30)]
    masses = sorted(state[0] / constants.EARTH_MASS_TO_KG for state in states)
    assert 1e-4 < masses[len(masses) // 2] < 1e-2
    assert {state[1] for state in states} <= {"dense", "typical", "sparse"} and "typical" in {state[1] for state in states}
    mass_kg, _, settled = states[0]
    assert mass_kg == pytest.approx(settled.mass_earth * constants.EARTH_MASS_TO_KG) and settled.largest_km > 0.0
    assert bl.belt_state(2.2, 3.4, SUN_KG, line, 1.0, 4.5, bl.ROCKY_LAW, draw.Stream("belt:0"))[0] == mass_kg, "the belt's own stream: the same belt, the same mass"
    young = bl.belt_state(2.2, 3.4, SUN_KG, line, 1.0, 0.01, bl.ROCKY_LAW, draw.Stream("belt:0"))[0]
    assert young > mass_kg, "a young belt has not been ground down"
    assert bl.belt_state(2.2, 3.4, SUN_KG, line, 2.0, 4.5, bl.ROCKY_LAW, draw.Stream("belt:0"))[0] >= mass_kg * 0.999, "a metal-rich disc leaves no less"


def test_bad_inputs_raise():
    with pytest.raises(ValueError):
        bm.belt_mass_earth(0.0, 2.2, 3.4, 1.0, 100.0, 500.0)
    with pytest.raises(ValueError):
        bm.primordial_mass_earth(3.0, 2.0, SUN_KG, 2.7)
    with pytest.raises(ValueError):
        bm.disruption_threshold_j_per_kg(0.0)
    with pytest.raises(ValueError):
        bm.density_word(0.0, 1.0)


def test_a_belt_whose_largest_body_is_too_strong_to_disrupt_keeps_what_formed():
    """Past Xc = 1 the constant of Wyatt et al.'s law goes negative; the mass is then the primordial one, not a negative number (GEN.1061)."""
    assert bm.collisional_constant(10.0, 10.2, 1.0, 5000.0) < 0.0
    assert bm.belt_mass_earth(2.0, 10.0, 10.2, 1.0, 4500.0, 5000.0) == pytest.approx(2.0)
