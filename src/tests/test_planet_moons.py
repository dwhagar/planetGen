# tests/test_planet_moons.py

"""GEN.1157 slice 0: moons of planet size around the giants (physics/planet_moons.py), pure and unwired."""

import pytest

from planetgen import tuning
from planetgen.physics import moon_composition
from planetgen.physics import planet_moons as pm
from planetgen.physics.constants import EARTH_MASS_TO_KG, G
from planetgen.util import draw

JUPITER = moon_composition.JUPITER_MASS_KG


def test_a_disk_cannot_hold_more_than_canup_and_ward_total():
    assert pm.DISK_TOTAL_RATIO_MAX == pytest.approx(6.4e-4, rel=0.02)
    assert pm.disk_moon_ceiling_kg(JUPITER) == pytest.approx(pm.DISK_TOTAL_RATIO_MAX * JUPITER)
    # and no draw of the disk's own total passes it
    top = moon_composition.total_moon_mass_ratio(moon_composition.ALPHA_OVER_F[1])
    assert top == pytest.approx(pm.DISK_TOTAL_RATIO_MAX)


def test_a_disk_moon_of_planet_size_needs_a_super_jovian_share():
    # a Jupiter's whole disk holds about 0.2 Earth masses: planet size only if one moon took the lot
    assert pm.disk_moon_ceiling_kg(JUPITER) < 3 * pm.PLANET_MOON_MIN_MASS_KG
    assert pm.disk_moon_ceiling_kg(10 * JUPITER) > 1.0 * EARTH_MASS_TO_KG


JUPITER_RADIUS_M = 7.1492e7
AU_M = 1.495978707e11
SUN_KG = 1.98841e30

# The worked case of [R5229]: Kepler-1625 b, 10 Jupiter masses at 0.87 au round 1.079 solar masses.
CASE_PARENT = 10 * JUPITER
CASE_HILL_M = 0.87 * AU_M * (CASE_PARENT / (3 * 1.079 * SUN_KG)) ** (1.0 / 3.0)
CASE_SPEED_M_S = (G * 1.079 * SUN_KG / (0.87 * AU_M)) ** 0.5


def test_the_bound_of_the_binary_encounter_reproduces_the_partner_masses_of_the_worked_case():
    """[R5229] finds that a Neptune-mass (17 Earth-mass) moon of Kepler-1625 b needs an ejected partner of about 5 to 50 Earth masses for an encounter at 1 to 30 km/s; the bound has the
    partner mass that leaves 17 Earth masses in that range (about 4.6 at 1 km/s and 53 at the planet's own 30 km/s)."""
    assert CASE_HILL_M == pytest.approx(261 * JUPITER_RADIUS_M, rel=0.02)
    assert CASE_SPEED_M_S == pytest.approx(30e3, rel=0.12)      # [R5229] rounds the planet's Keplerian speed to 30 km/s
    needed = {}
    for speed_km_s in (1.0, 30.0):
        low, high = 1.0 * EARTH_MASS_TO_KG, 1.0e4 * EARTH_MASS_TO_KG
        for _ in range(100):
            middle = (low * high) ** 0.5
            if pm.encounter_limit_kg(CASE_PARENT, middle, CASE_HILL_M, speed_km_s * 1e3 / pm.EXCHANGE_V_INFINITY_FRACTION) < 17 * EARTH_MASS_TO_KG:
                low = middle
            else:
                high = middle
        needed[speed_km_s] = low / EARTH_MASS_TO_KG
    assert needed[1.0] == pytest.approx(4.6, rel=0.1)
    assert needed[30.0] == pytest.approx(53.0, rel=0.1)


def test_the_bound_rises_with_the_partner_falls_with_the_speed_and_never_passes_the_partner():
    bounds = [pm.exchange_bound_kg(CASE_PARENT, m * EARTH_MASS_TO_KG, CASE_HILL_M, CASE_SPEED_M_S) for m in (1, 3, 10, 30, 100, 300)]
    assert bounds == sorted(bounds)
    assert bounds[0] == 0.0 and bounds[-1] > 0.0
    assert all(bound <= m * EARTH_MASS_TO_KG for bound, m in zip(bounds, (1, 3, 10, 30, 100, 300)))
    slow = pm.exchange_bound_kg(CASE_PARENT, 30 * EARTH_MASS_TO_KG, CASE_HILL_M, 0.1 * CASE_SPEED_M_S)
    assert slow > pm.exchange_bound_kg(CASE_PARENT, 30 * EARTH_MASS_TO_KG, CASE_HILL_M, CASE_SPEED_M_S)
    assert pm.exchange_bound_kg(CASE_PARENT, 30 * EARTH_MASS_TO_KG, 0.0, CASE_SPEED_M_S) == 0.0


def test_the_exchange_chance_is_one_in_seventy_where_a_partner_can_leave_a_moon_and_zero_elsewhere():
    assert pm.exchange_chance(CASE_PARENT, CASE_HILL_M, CASE_SPEED_M_S) == pytest.approx(1.0 / 70.0)
    assert pm.exchange_chance(1.0 * EARTH_MASS_TO_KG, CASE_HILL_M, CASE_SPEED_M_S) == 0.0       # no partner lighter than the parent leaves a moon
    assert pm.exchange_chance(CASE_PARENT, 0.0, CASE_SPEED_M_S) == 0.0                           # no Hill sphere
    assert pm.exchange_chance(CASE_PARENT, CASE_HILL_M, 100.0 * CASE_SPEED_M_S) == 0.0            # an encounter too fast for any partner
    assert pm.minimum_partner_kg(CASE_PARENT, CASE_HILL_M, CASE_SPEED_M_S) > pm.PLANET_MOON_MIN_MASS_KG
    lightest = pm.minimum_partner_kg(CASE_PARENT, CASE_HILL_M, CASE_SPEED_M_S)
    assert pm.exchange_bound_kg(CASE_PARENT, lightest, CASE_HILL_M, CASE_SPEED_M_S) == pytest.approx(pm.PLANET_MOON_MIN_MASS_KG, rel=1e-6)


def test_the_mass_runs_from_the_floor_to_the_ceiling_by_the_generators_law():
    assert pm._law_mass_kg(1.0, 100.0, 0.0) == pytest.approx(1.0)
    assert pm._law_mass_kg(1.0, 100.0, 1.0) == pytest.approx(100.0)
    masses = [pm._law_mass_kg(1.0, 100.0, i / 100) for i in range(101)]
    assert masses == sorted(masses)
    assert tuning.MOON_MASS_FUNCTION_SLOPE > 0
    assert pm._law_mass_kg(5.0, 5.0, 0.3) == 5.0
    with pytest.raises(ValueError):
        pm._law_mass_kg(1.0, 100.0, 1.5)


def test_a_moon_never_passes_the_bound_of_its_encounter():
    seen = 0
    for seed in range(1200):
        stream = draw.Stream(f"pm:{seed}")
        moon = pm.draw_exchange_moon(CASE_PARENT, CASE_HILL_M, CASE_SPEED_M_S, CASE_PARENT / 10, stream)
        if moon is not None:
            seen += 1
            assert moon.origin == "exchange"
            assert pm.PLANET_MOON_MIN_MASS_KG <= moon.mass_kg <= CASE_PARENT / 10      # the generator's one-tenth limit holds the bound of the encounter (which no partner, the parent at most, passes)
    assert seen > 0


def test_the_draw_takes_three_numbers_whatever_the_parent():
    class Counting:
        def __init__(self):
            self.n = 0

        def random(self):
            self.n += 1
            return 0.001

    for parent, hill in ((EARTH_MASS_TO_KG, CASE_HILL_M), (JUPITER, CASE_HILL_M), (20 * JUPITER, CASE_HILL_M), (JUPITER, 0.0)):
        rng = Counting()
        pm.draw_exchange_moon(parent, hill, CASE_SPEED_M_S, parent / 10, rng)
        assert rng.n == 3


def test_about_one_giant_in_seventy_holds_one():
    stream = draw.Stream("pm:rate")
    hits = sum(pm.draw_exchange_moon(CASE_PARENT, CASE_HILL_M, CASE_SPEED_M_S, CASE_PARENT / 10, stream) is not None for _ in range(14000))
    assert hits == pytest.approx(200, rel=0.25)


def test_bad_inputs_are_refused():
    with pytest.raises(ValueError):
        pm.disk_moon_ceiling_kg(0.0)
    with pytest.raises(ValueError):
        pm.exchange_chance(0.0, CASE_HILL_M, CASE_SPEED_M_S)
