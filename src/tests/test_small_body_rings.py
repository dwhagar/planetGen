# tests/test_small_body_rings.py

"""GEN.1056: the resonance rings of small bodies (R4) and the gate against rings close to a star (R9), `physics/small_body_rings.py`.
Quaoar, Haumea and Chariklo are test targets only: nothing in the module is centred on one."""

import math

import pytest

from planetgen.physics import moon_composition, rings
from planetgen.physics import small_body_rings as sbr
from planetgen.util import draw

SUN_KG = 1.989e30

# name: (mass kg, mean radius km, spin period h, distance from the Sun AU, observed ring radius km, how far the ring may lie from the resonance, share)
QUAOAR = (1.2e21, 555.0, 17.68, 43.7, 4057.0, 0.05)
HAUMEA = (4.0e21, 780.0, 3.9155, 43.1, 2287.0, 0.05)
CHARIKLO = (8.0e18, 124.0, 7.004, 15.8, 398.0, 0.10)


class Always:
    """A stream whose chance draw always passes the rate (`random()` is 0) and whose other draws are the middle of their range."""

    def random(self):
        return 0.0

    def uniform(self, a, b):
        return 0.5 * (a + b)


@pytest.mark.parametrize("case", [QUAOAR, HAUMEA, CHARIKLO], ids=["Quaoar", "Haumea", "Chariklo"])
def test_the_ring_orbits_three_times_slower_than_the_body_turns_where_the_known_rings_are(case):
    mass, _radius, spin, _distance, observed, tolerance = case
    assert sbr.resonance_radius_km(mass, spin) == pytest.approx(observed, rel=tolerance)


@pytest.mark.parametrize("case", [QUAOAR, HAUMEA], ids=["Quaoar", "Haumea"])
def test_the_known_rings_lie_past_the_roche_limit_of_ice_and_inside_the_stable_zone(case):
    mass, radius, spin, distance, observed, _tolerance = case
    found = sbr.resonance_ring(mass, radius, spin, distance, SUN_KG, True, Always())
    assert found is not None
    assert found.radius_km > found.roche_km
    assert found.ring.ring_inner_km < found.radius_km < found.ring.ring_outer_km
    assert found.radius_km == pytest.approx(observed, rel=0.05)
    assert found.ring.ring_kind == "narrow" and found.ring.ring_origin in rings.ORIGINS


def test_the_ring_is_thin_light_and_dense_enough_to_see():
    mass, radius, spin, distance, _observed, _tolerance = QUAOAR
    found = sbr.resonance_ring(mass, radius, spin, distance, SUN_KG, True, Always())
    ring = found.ring
    assert ring.ring_mass_kg < rings.MAX_MASS_OVER_PLANET * mass
    assert (ring.ring_outer_km - ring.ring_inner_km) / found.radius_km < 0.06
    assert 1.0e-3 < ring.ring_optical_depth < 10.0


def test_a_fast_dense_spinner_whose_resonance_is_inside_the_roche_limit_has_no_such_ring():
    # a very fast rotation puts the resonance close in: the ordinary rings of rings.py are that case
    mass, radius, _spin, distance, _observed, _tolerance = QUAOAR
    assert sbr.resonance_radius_km(mass, 2.0) < 2.0 * radius
    assert sbr.resonance_ring(mass, radius, 2.0, distance, SUN_KG, True, Always()) is None


def test_a_resonance_outside_the_stable_zone_is_no_ring():
    mass, radius, _spin, _distance, _observed, _tolerance = QUAOAR
    # a slow spinner close to its star: the resonance lies far outside the half of the Hill radius that stays bound
    assert sbr.resonance_radius_km(mass, 300.0) > sbr.HILL_FRACTION_MAX * 0.2 * 1.495978707e8 * (mass / (3.0 * SUN_KG)) ** (1.0 / 3.0)
    assert sbr.resonance_ring(mass, radius, 300.0, 0.2, SUN_KG, True, Always()) is None


def test_a_planet_sized_body_takes_the_rule_of_rings_py_and_not_this_one():
    heavy = moon_composition.DWARF_PLANET_MAX_MASS_KG * 2.0
    assert sbr.resonance_ring(heavy, 2000.0, 20.0, 40.0, SUN_KG, True, Always()) is None


def test_the_chance_is_the_registered_rate():
    mass, radius, spin, distance, _observed, _tolerance = QUAOAR
    hits = 0
    with draw.bound(3):
        for number in range(4000):
            hits += sbr.resonance_ring(mass, radius, spin, distance, SUN_KG, True, draw.Stream(f"r4:{number}")) is not None
    assert hits / 4000 == pytest.approx(sbr.RESONANCE_RING_RATE, abs=0.02)


def test_the_stream_is_used_the_same_whether_or_not_a_ring_results():
    mass, radius, spin, distance, _observed, _tolerance = QUAOAR
    for u_stream in ("a", "b", "c"):
        passing, failing = draw.Stream(u_stream), draw.Stream(u_stream)
        sbr.resonance_ring(mass, radius, spin, distance, SUN_KG, True, passing)
        sbr.resonance_ring(mass, radius, spin, 0.05, SUN_KG, True, failing)     # R9 refuses, the draws are still taken
        assert passing.random() == failing.random()


def test_no_ring_survives_inside_a_tenth_of_an_au():
    assert not sbr.ring_survives(0.05) and not sbr.ring_survives(0.0999)
    assert sbr.ring_survives(0.1) and sbr.ring_survives(1.0)
    with pytest.raises(ValueError):
        sbr.ring_survives(0.0)


def test_a_body_inside_the_gate_has_no_ring_whatever_the_draw():
    mass, radius, spin, _distance, _observed, _tolerance = QUAOAR
    assert sbr.resonance_ring(mass, radius, spin, 0.05, SUN_KG, True, Always()) is None


def test_the_resonance_radius_follows_kepler_and_scales_with_the_order():
    first = sbr.resonance_radius_km(1.0e21, 10.0, order=1)
    assert sbr.resonance_radius_km(1.0e21, 10.0, order=8) == pytest.approx(first * 4.0)
    assert sbr.resonance_radius_km(8.0e21, 10.0, order=1) == pytest.approx(first * 2.0)
    for bad in ((0.0, 10.0), (1.0e21, 0.0)):
        with pytest.raises(ValueError):
            sbr.resonance_radius_km(*bad)


def test_the_rate_the_limits_and_the_ranges_are_registered_stand_ins():
    assert 0.0 < sbr.RESONANCE_RING_RATE < 1.0 and sbr.INNER_RING_LIMIT_AU == 0.1
    assert math.isclose(sbr.HILL_FRACTION_MAX, 0.5)
