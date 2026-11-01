# tests/test_body_status.py

"""
GEN.1063: the status of a body, a pure function of its mass, composition, orbit and the neighbours in the ledger. The Solar System bodies are test targets (the masses
below are read from the literature for the test only); the module takes none of them as input.
"""

import math

import pytest

from planetgen.physics import body_status as bs

EARTH = 5.972e24
SUN_NEIGHBOURS = {   # (name, a_au, mass_kg)
    "mercury": (0.387, 3.30e23), "venus": (0.723, 4.87e24), "earth": (1.0, EARTH), "mars": (1.524, 6.42e23),
    "jupiter": (5.203, 1.898e27), "saturn": (9.58, 5.68e26),
}


def _status(name, mass, ice, a, others):
    return bs.status(mass, ice, a, 1.0, others)


@pytest.mark.parametrize("name, mass, ice, a, others, word", [
    ("Earth", EARTH, 0.0, 1.0, 3.5e18, bs.PLANET),
    ("Mars", 6.42e23, 0.0, 1.524, 1.0e19, bs.PLANET),
    ("Jupiter", 1.898e27, 0.0, 5.203, 1.0e24, bs.PLANET),
    ("Ceres", 9.39e20, 0.1, 2.77, 2.8e21, bs.DWARF_PLANET),
    ("Pluto", 1.30e22, 0.3, 39.5, 1.2e23, bs.DWARF_PLANET),
    ("Eris", 1.66e22, 0.3, 67.7, 1.2e23, bs.DWARF_PLANET),
    ("Trojan", 1.0e19, 0.2, 5.203, 1.898e27, bs.MINOR_BODY),
    ("Juno", 2.7e19, 0.0, 2.67, 2.8e21, bs.MINOR_BODY),
])
def test_the_status_word_of_the_known_bodies(name, mass, ice, a, others, word):
    assert _status(name, mass, ice, a, others).word == word


def test_the_numbers_are_kept_beside_the_word_and_the_two_discriminants_agree():
    earth = _status("Earth", EARTH, 0.0, 1.0, 3.5e18)
    assert earth.mu == pytest.approx(1.7e6, rel=0.02)
    assert earth.pi == pytest.approx(816, rel=0.01)       # 1 / (1.9e-4 C^1.5) at 1 AU
    assert earth.lam > 1e5
    ceres = _status("Ceres", 9.39e20, 0.1, 2.77, 2.8e21)
    assert ceres.mu == pytest.approx(1 / 3, rel=0.05)
    assert ceres.pi < 1 and ceres.lam < 1
    assert all(_status(n, m, i, a, o).pi_agrees for n, m, i, a, o in [
        ("Earth", EARTH, 0.0, 1.0, 3.5e18), ("Ceres", 9.39e20, 0.1, 2.77, 2.8e21), ("Pluto", 1.3e22, 0.3, 39.5, 1.2e23)])


def test_margot_pi_of_jupiter_and_mars_is_the_printed_value():
    assert bs.margot_pi(1.898e27, 5.203, 1.0) == pytest.approx(4.0e4, rel=0.05)
    assert bs.margot_pi(6.42e23, 1.524, 1.0) == pytest.approx(54, rel=0.05)


def test_lambda_is_one_at_margots_clearing_mass_at_one_au_round_the_sun():
    mass = bs.margot_clearing_mass_earth(1.0, 1.0) * bs.constants.EARTH_MASS_TO_KG
    assert bs.stern_levison_lambda(mass, 1.0, 1.0) == pytest.approx(1.0)


def test_a_moon_carries_only_its_roundness():
    moon = bs.status(7.35e22, 0.0, 1.0, 1.0, 0.0, parent_is_planet=True)
    assert moon.word == bs.MOON and moon.roundness > 1.0
    assert moon.mu is None and moon.pi is None and moon.lam is None and moon.pi_agrees is None


def test_roundness_is_the_engines_ratio_and_rises_with_mass():
    assert bs.roundness(bs.round_mass_kg(0.0), 0.0) == pytest.approx(1.0, rel=1e-6)
    assert bs.roundness(bs.round_mass_kg(1.0), 1.0) == pytest.approx(1.0, rel=1e-6)
    assert 1.0e19 < bs.round_mass_kg(1.0) < 1.0e20 < bs.round_mass_kg(0.0) < 1.0e21          # ice turns round at a smaller mass than rock
    assert bs.roundness(1e19, 0.0) < 1.0 < bs.roundness(1e22, 0.0)
    masses = [1e18, 1e19, 3e19, 1e20, 3e20, 1e21, 1e22]
    for ice in (0.0, 0.5, 1.0):
        values = [bs.roundness(m, ice) for m in masses]
        assert values == sorted(values)
    assert bs.roundness(1e20, 1.0) > bs.roundness(1e20, 0.0)
    with pytest.raises(ValueError):
        bs.roundness(0.0, 0.5)


def test_the_zone_counts_the_neighbours_in_it_and_the_share_of_a_belt():
    half = bs.ZONE_HILL_RADII * bs.hill_radius_au(1.0, EARTH, 1.0)
    inside, outside = (1.0 + 0.5 * half, 1e20), (1.0 + 2 * half, 1e22)
    assert bs.zone_mass_kg(1.0, EARTH, 1.0, [(1.0, EARTH), inside, outside]) == 1e20          # the body itself is left out
    belt = (1.0 - 10 * half, 1.0 + 10 * half, 1e21)
    assert bs.zone_mass_kg(1.0, EARTH, 1.0, belts=[belt]) == pytest.approx(1e21 * 2 * half / (20 * half))
    assert bs.zone_mass_kg(1.0, EARTH, 1.0, belts=[(5.0, 6.0, 1e21)]) == 0.0
    assert bs.soter_mu(EARTH, 0.0) == math.inf


def test_the_solar_system_planets_clear_their_zones_from_the_planet_list_and_the_belt_remainder():
    belts = [(2.1, 3.3, 3.0e21), (30.0, 50.0, 2.4e23)]       # test-only: the main belt and the Kuiper belt
    for name, (a, mass) in SUN_NEIGHBOURS.items():
        others = bs.zone_mass_kg(a, mass, 1.0, list(SUN_NEIGHBOURS.values()), belts)
        assert bs.status(mass, 0.0, a, 1.0, others).word == bs.PLANET, name
