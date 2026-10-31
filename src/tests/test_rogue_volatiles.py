# tests/test_rogue_volatiles.py

"""GEN.1112 RP3: the air and water a thrown-out planet carries (`planetgen.physics.rogue_volatiles`): the envelope it leaves with, its pressure, the giant class, the water and the frozen air."""

import random

import pytest

from planetgen import tuning
from planetgen.physics import hydrosphere, rogue_volatiles as rv


def test_a_planet_thrown_out_in_its_first_gigayear_keeps_the_envelope_formation_gave_it():
    kept = rv.envelope_at_ejection(5.0, 0.05, 100.0, 1.0, 1.0, 0.5)
    assert kept == pytest.approx(5.0 * 0.05 / 0.95)                  # the loss law starts at ONSET_GY; what is lost earlier is in the starting envelope


def test_a_planet_that_stayed_close_to_its_star_for_gigayears_loses_it_and_one_far_out_does_not():
    assert rv.envelope_at_ejection(5.0, 0.05, 100.0, 1.0, 1.0, 3.0) == 0.0
    far = rv.envelope_at_ejection(5.0, 0.05, 1.0, 1.0, 1.0, 3.0)
    assert far == pytest.approx(5.0 * 0.05 / 0.95, rel=0.02)
    assert rv.envelope_at_ejection(5.0, 0.0, 1.0, 1.0, 1.0, 3.0) == 0.0


def test_the_base_pressure_of_a_rock_with_a_thin_envelope_is_stevensons_hundreds_of_bar():
    pressure = rv.envelope_base_pressure_bar(1.0, 1.0, 0.001)
    assert 100.0 < pressure < 3000.0
    assert pressure == pytest.approx(1150.0, rel=0.01)
    assert rv.envelope_base_pressure_bar(1.0, 1.0, 0.0) == 0.0
    assert rv.envelope_base_pressure_bar(1.0, 1.0, 0.002) == pytest.approx(2.0 * pressure)


def test_the_giant_class_is_the_radius_the_generator_already_calls_a_giant():
    assert rv.envelope_radius_km(10.0, 0.1, 3.0) > tuning.GIANT_MIN_RADIUS_KM
    assert rv.is_gas_giant(10.0, 0.1, 3.0)
    assert not rv.is_gas_giant(5.0, 0.02, 3.0)
    assert not rv.is_gas_giant(1.0, 0.001, 3.0)


def test_the_mass_threshold_falls_as_the_planet_gets_heavier():
    fractions = [rv.smallest_giant_envelope_fraction(m, 3.0) for m in (5.0, 10.0, 17.0)]
    assert fractions == sorted(fractions, reverse=True)
    assert fractions[0] == pytest.approx(0.05, abs=0.01) and fractions[2] == pytest.approx(0.015, abs=0.005)
    assert rv.smallest_giant_envelope_fraction(1.0, 3.0) == pytest.approx(0.06, abs=0.01)
    assert rv.smallest_giant_envelope_fraction(0.01, 3.0) == rv.smallest_giant_envelope_fraction(1.0, 3.0)     # held at the radius law's lower edge
    assert rv.envelope_radius_km(100.0, 0.5, 3.0) == rv.envelope_radius_km(30.0, 0.5, 3.0)
    for m, f in zip((5.0, 10.0, 17.0), fractions):
        assert rv.is_gas_giant(m, f + 1e-3, 3.0) and not rv.is_gas_giant(m, f - 1e-3, 3.0)


def test_water_follows_the_snow_line_the_planet_formed_at():
    inside = [rv.water_mass_fraction(1.0, 1.0, 2.7, random.Random(i)) for i in range(60)]
    outside = [rv.water_mass_fraction(1.0, 6.0, 2.7, random.Random(i)) for i in range(60)]
    assert sum(inside) / 60 < sum(outside) / 60
    assert sum(w >= hydrosphere.ICY_MIN_WATER_FRACTION for w in outside) > sum(w >= hydrosphere.ICY_MIN_WATER_FRACTION for w in inside)


def test_the_water_draw_is_the_in_system_law_and_takes_two_numbers():
    class Counting(random.Random):
        draws = 0

        def random(self):
            Counting.draws += 1
            return super().random()
    value = rv.water_mass_fraction(2.0, 4.0, 2.7, Counting(3))
    assert value > 0.0 and Counting.draws == 2
    assert rv.is_water_rich(0.05) and not rv.is_water_rich(1e-4)


def test_an_earth_held_an_air_and_the_moon_did_not():
    assert rv.had_air(1.0, 1.0, 1.0)
    assert not rv.had_air(0.0123, 0.273, 1.0)
    assert rv.had_air(0.107, 0.53, 0.43)                              # Mars, thinly
    assert rv.had_air(0.1, 0.5, 0.01)                                 # far from a star even a small rock keeps it


def test_bad_inputs_are_refused():
    with pytest.raises(rv.RogueVolatilesError):
        rv.envelope_at_ejection(0.0, 0.05, 1.0, 1.0, 1.0, 1.0)
    with pytest.raises(rv.RogueVolatilesError):
        rv.envelope_at_ejection(5.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    with pytest.raises(rv.RogueVolatilesError):
        rv.water_mass_fraction(1.0, 0.0, 2.7)
    with pytest.raises(rv.RogueVolatilesError):
        rv.had_air(1.0, 0.0, 1.0)
