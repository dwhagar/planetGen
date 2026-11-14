# tests/test_planet_moons.py

"""GEN.1157 slice 0: moons of planet size around the giants (physics/planet_moons.py), pure and unwired."""

import pytest

from planetgen import tuning
from planetgen.physics import moon_composition
from planetgen.physics import planet_moons as pm
from planetgen.physics.constants import EARTH_MASS_TO_KG
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


def test_the_exchange_chance_is_one_in_seventy_for_a_giant_and_zero_for_a_pebble():
    assert pm.exchange_chance(JUPITER) == pytest.approx(1.0 / 70.0)
    assert pm.exchange_chance(1.0 * EARTH_MASS_TO_KG) == 0.0
    assert pm.exchange_chance(17 * EARTH_MASS_TO_KG) > 0.0     # a Neptune can hold a tenth of an Earth mass at 1e-2


def test_the_mass_runs_from_the_floor_to_the_ceiling_by_the_generators_law():
    ceiling = pm.exchange_ceiling_kg(JUPITER)
    assert pm.exchange_mass_kg(JUPITER, 0.0) == pytest.approx(pm.PLANET_MOON_MIN_MASS_KG)
    assert pm.exchange_mass_kg(JUPITER, 1.0) == pytest.approx(ceiling)
    masses = [pm.exchange_mass_kg(JUPITER, i / 100) for i in range(101)]
    assert masses == sorted(masses)
    assert tuning.MOON_MASS_FUNCTION_SLOPE > 0
    with pytest.raises(ValueError):
        pm.exchange_mass_kg(1.0 * EARTH_MASS_TO_KG, 0.5)
    with pytest.raises(ValueError):
        pm.exchange_mass_kg(JUPITER, 1.5)


def test_a_moon_never_passes_its_ceiling_or_the_generators_own_limit():
    for seed in range(300):
        moon = pm.draw_exchange_moon(10 * JUPITER, draw.Stream(f"pm:{seed}"))
        if moon is not None:
            assert moon.origin == "exchange"
            assert pm.PLANET_MOON_MIN_MASS_KG <= moon.mass_kg <= pm.exchange_ceiling_kg(10 * JUPITER)
            assert moon.mass_kg <= 10 * JUPITER / 10       # planets.moon_size_limits: a tenth of the parent


def test_the_draw_takes_two_numbers_whatever_the_parent():
    class Counting:
        def __init__(self):
            self.n = 0

        def random(self):
            self.n += 1
            return 0.9

    for parent in (EARTH_MASS_TO_KG, JUPITER, 20 * JUPITER):
        rng = Counting()
        pm.draw_exchange_moon(parent, rng)
        assert rng.n == 2


def test_about_one_giant_in_seventy_holds_one():
    stream = draw.Stream("pm:rate")
    hits = sum(pm.draw_exchange_moon(JUPITER, stream) is not None for _ in range(14000))
    assert hits == pytest.approx(200, rel=0.25)


def test_bad_inputs_are_refused():
    for call in (pm.disk_moon_ceiling_kg, pm.exchange_chance, pm.exchange_ceiling_kg):
        with pytest.raises(ValueError):
            call(0.0)
