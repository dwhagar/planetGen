# tests/test_sulfuric_vapor.py

"""GEN.1424: the vapor pressure of sulfuric acid and water over aqueous sulfuric acid and the window the liquid condenses in
(`planetgen.physics.sulfuric_vapor`): Zeleznik's Table 7, the pure-acid and water curves, Tabazadeh's Table 1, and a Venus-like and a cold air."""

import json
import math

import pytest

from planetgen.physics import saturation
from planetgen.physics import sulfuric_vapor as sv


# [R4877: Table 7], read from the paper: (T, x acid, -G(r)/RT, -mu_acid(r)/RT, -mu_water(r)/RT).
ZELEZNIK_TABLE_7 = (
    (200.0, 0.2, 7.0325, 22.2093, 3.2383),
    (200.0, 0.5, 8.8609, 6.4547, 11.2671),
    (200.0, 0.8, 5.5743, 2.5271, 17.7633),
    (250.0, 0.2, 4.9699, 16.5640, 2.0713),
    (250.0, 0.5, 6.2077, 3.8505, 8.5649),
    (250.0, 0.8, 3.3074, 0.7680, 13.4652),
    (300.0, 0.4, 5.0548, 5.3346, 4.8682),
    (300.0, 0.8, 2.5047, 0.4121, 10.8755),
)

# Tabazadeh and others 1997 (Geophys. Res. Lett. 24, 1931) Table 1: ln p(H2O, mbar) = a + b/T + c/T^2 over a solution of the weight percent.
TABAZADEH_TABLE_1 = {
    20: (19.761, -4414.7, -142940.0),
    40: (20.379, -4828.5, -112550.0),
    60: (20.405, -5252.1, -100840.0),
    75: (21.169, -6310.6, -48396.0),
    80: (21.808, -6985.9, -12170.0),
}


def water_pressure_over_solution_pa(weight, temperature_k):
    return math.exp(sv.ln_water_activity(sv.mole_fraction(weight), temperature_k)) * sv.water_pressure_pa(temperature_k)


def test_the_rule_file_is_complete():
    document = json.load(open(sv.DATA_FILE, encoding="utf-8"))
    assert document["version"] == sv.VERSION
    assert set(document["zeleznik_1991"]["parameters"]) == {"mu_111", "mu_121", "mu_221", "mu_122", "eps_111", "eps_121", "eps_211", "eps_221", "eps_122", "eps_212"}
    for name, entry in document["numbers"].items():
        assert entry["tier"].startswith("T") and entry["tweak"].startswith("L") and entry["basis"] and "value" in entry, name


def test_a_rule_file_entry_without_a_tier_is_refused():
    with pytest.raises(sv.SulfuricVaporError):
        sv._check("x", {"tweak": "L0", "basis": "b"})


@pytest.mark.parametrize("temperature_k,x,g,mu_acid,mu_water", ZELEZNIK_TABLE_7)
def test_the_model_reproduces_zeleznik_table_7(temperature_k, x, g, mu_acid, mu_water):
    value, _ = sv.relative_gibbs(x, temperature_k)
    acid, water = sv.relative_potentials(x, temperature_k)
    assert -value == pytest.approx(g, abs=2.0e-4)
    assert -acid == pytest.approx(mu_acid, abs=2.0e-4)
    assert -water == pytest.approx(mu_water, abs=2.0e-4)


def test_the_activities_are_one_for_the_pure_liquids_and_fall_with_the_other_component():
    for temperature_k in (200.0, 260.0, 340.0):
        assert sv.ln_water_activity(0.0, temperature_k) == pytest.approx(0.0, abs=1.0e-9)
        assert sv.ln_acid_activity(1.0, temperature_k) == pytest.approx(0.0, abs=1.0e-9)
        grid = [i / 20.0 for i in range(1, 20)]
        water = [sv.ln_water_activity(x, temperature_k) for x in grid]
        acid = [sv.ln_acid_activity(x, temperature_k) for x in grid]
        assert all(b < a for a, b in zip(water, water[1:])), "more acid, less water activity"
        assert all(b > a for a, b in zip(acid, acid[1:])), "more acid, more acid activity"


def test_the_two_activities_obey_the_gibbs_duhem_relation():
    h = 1.0e-5
    for temperature_k in (210.0, 250.0, 330.0):
        for x in (0.1, 0.3, 0.6, 0.85):
            slope_acid = (sv.ln_acid_activity(x + h, temperature_k) - sv.ln_acid_activity(x - h, temperature_k)) / (2.0 * h)
            slope_water = (sv.ln_water_activity(x + h, temperature_k) - sv.ln_water_activity(x - h, temperature_k)) / (2.0 * h)
            assert x * slope_acid + (1.0 - x) * slope_water == pytest.approx(0.0, abs=1.0e-6)


def test_the_weight_percent_and_the_mole_fraction_are_each_others_inverse():
    for weight in (10.0, 50.0, 75.0, 98.0):
        assert sv.weight_percent(sv.mole_fraction(weight)) == pytest.approx(weight, rel=1.0e-12)
    assert sv.weight_percent(0.0) == 0.0 and sv.weight_percent(1.0) == pytest.approx(100.0)


def test_the_pure_acid_curve_is_the_published_one_and_the_refit_is_a_constant_factor_above():
    # ln p(atm) at 360 K is 16.259 - 10156/360 by construction; the other values are the form's (the three open papers print it)
    assert sv.pure_acid_pressure_pa(360.0) / sv.ATMOSPHERE_PA == pytest.approx(math.exp(16.259 - 10156.0 / 360.0), rel=1.0e-9)
    for temperature_k, atmospheres in ((190.0, 1.15e-17), (250.0, 1.54e-11), (300.0, 2.02e-8), (400.0, 1.04e-4)):
        assert sv.pure_acid_pressure_pa(temperature_k) / sv.ATMOSPHERE_PA == pytest.approx(atmospheres, rel=5.0e-3)
    ratios = [sv.pure_acid_pressure_pa(t, True) / sv.pure_acid_pressure_pa(t) for t in (190.0, 250.0, 360.0, 445.0)]
    assert all(1.27 < r < 1.29 for r in ratios) and max(ratios) - min(ratios) < 0.01


def test_the_liquid_water_curve_joins_the_liquid_curve_of_saturation_at_the_triple_point():
    below = sv.water_pressure_pa(saturation.WATER_TRIPLE_POINT_K - 1.0e-6)
    above = sv.water_pressure_pa(saturation.WATER_TRIPLE_POINT_K)
    assert below == pytest.approx(above, rel=2.0e-3)
    assert sv.water_pressure_pa(250.0) > saturation.pressure_pa("H2O", 250.0), "supercooled liquid is above the ice curve"
    assert sv.water_pressure_pa(300.0) == saturation.pressure_pa("H2O", 300.0)


def test_a_cold_dilute_solution_agrees_with_tabazadeh_and_a_strong_one_is_within_the_models_spread():
    """Tabazadeh's fits are to the Clegg and Brimblecombe model, Zeleznik's to its own data: they agree to 2 percent up to 40 percent acid at 190 to 260 K and
    part by up to 40 percent in the water pressure at 80 percent acid (the registered spread, SV3)."""
    for temperature_k in (190.0, 200.0, 230.0, 260.0):
        for weight, (a, b, c) in TABAZADEH_TABLE_1.items():
            theirs = 100.0 * math.exp(a + b / temperature_k + c / temperature_k ** 2)
            ratio = water_pressure_over_solution_pa(weight, temperature_k) / theirs
            if weight <= 40:
                assert ratio == pytest.approx(1.0, abs=0.02), (temperature_k, weight)
            else:
                assert 0.55 < ratio < 1.02, (temperature_k, weight)


def test_an_air_with_more_water_makes_a_weaker_solution_and_the_limits_are_pure_water_and_pure_acid():
    temperature_k = 280.0
    p0 = sv.water_pressure_pa(temperature_k)
    weights = [sv.composition(temperature_k, share * p0).weight_percent for share in (1.0e-4, 1.0e-3, 1.0e-2, 0.1, 0.5)]
    assert all(b < a for a, b in zip(weights, weights[1:]))
    saturated = sv.composition(temperature_k, p0)
    assert saturated.status == sv.WATER_SATURATED and saturated.weight_percent == 0.0
    assert sv.acid_pressure_over_solution_pa(temperature_k, 2.0 * p0) == 0.0
    dry = sv.composition(temperature_k, 0.0)
    assert dry.weight_percent == 100.0 and dry.acid_fraction == 1.0


def test_the_composition_solves_the_water_activity_it_was_given():
    for temperature_k, share in ((220.0, 1.0e-3), (290.0, 3.0e-2), (340.0, 2.0e-4)):
        p = share * sv.water_pressure_pa(temperature_k)
        solution = sv.composition(temperature_k, p)
        assert math.exp(sv.ln_water_activity(solution.acid_fraction, temperature_k)) == pytest.approx(solution.water_activity, rel=1.0e-9)


def test_the_status_says_when_the_solution_model_is_carried_past_its_window():
    assert sv.status_at(250.0) == sv.OK and sv.status_at(350.0) == sv.OK
    assert sv.status_at(361.0) == sv.EXTRAPOLATED and sv.status_at(190.0) == sv.EXTRAPOLATED
    assert sv.composition(361.0, 30.0e-6 * 1.3e5).status == sv.EXTRAPOLATED


def test_a_venus_like_air_condenses_the_acid_near_the_published_base_and_the_strength_falls_going_up():
    """The published lower cloud base is 345 to 367 K at 0.98 to 1.39 bar with the acid 95 percent and over, the top 75 to 96 percent ([R0157: Table 1]
    through `sulfuric_cloud.json`, the strengths of the Mode 2 table in the Venus acid-cloud report): 30 ppm of water and 3 to 10 ppm of acid vapor
    at 1.3 bar put the base inside that, and the strength at the levels above falls from 98 to the 80s."""
    pressure = 1.3e5
    for acid_ppm in (3.0, 10.0):
        base = sv.condensation_temperature_k(acid_ppm * 1.0e-6 * pressure, 30.0e-6 * pressure)
        assert 345.0 <= base <= 372.0, acid_ppm
        assert sv.composition(base, 30.0e-6 * pressure).weight_percent >= 95.0
    levels = ((232.0, 0.037e5, 2.0e-6), (243.0, 0.098e5, 2.0e-6), (263.0, 0.235e5, 3.0e-6), (300.0, 0.53e5, 5.0e-6), (348.0, 1.07e5, 20.0e-6), (361.0, 1.3e5, 30.0e-6))
    strengths = [sv.composition(t, w * p).weight_percent for t, p, w in levels]
    assert all(b > a for a, b in zip(strengths, strengths[1:]))
    assert 75.0 <= strengths[0] <= 96.0 and strengths[-1] >= 95.0
    assert sv.condensation_temperature_k(10.0e-6 * pressure, 30.0e-6 * pressure, alternative=True) < sv.condensation_temperature_k(10.0e-6 * pressure, 30.0e-6 * pressure)


def test_the_window_over_a_column_is_the_levels_cooler_than_the_base():
    levels = [(1.3e5, 361.0), (1.0e5, 348.0), (0.5e5, 300.0), (0.1e5, 245.0)]
    flags = sv.window(levels, 2.0e-6, 30.0e-6)
    assert not flags[0] and flags[-1] and flags == sorted(flags), "warm below, cloud above, no gap between"
    assert sv.window(levels, 0.0, 30.0e-6) == [False] * 4


def test_a_cold_dry_air_makes_a_strong_acid_aerosol_and_next_to_no_acid_condenses_at_the_dew_point():
    stratosphere = sv.composition(200.0, 5.0e-6 * 5000.0)
    assert 55.0 < stratosphere.weight_percent < 70.0 and stratosphere.status == sv.OK
    assert sv.condensation_temperature_k(0.0, 30.0e-6 * 1.3e5) is None
    water = 30.0e-6 * 1.3e5
    dew = sv.condensation_temperature_k(1.0e-40, water)
    assert sv.water_pressure_pa(dew) == pytest.approx(water, rel=1.0e-3), "next to no acid: the liquid forms where the water itself does (the dew point)"
    assert not sv.condenses(361.0, 1.0e-8, 30.0e-6 * 1.3e5)


def test_the_pressures_must_have_a_positive_temperature():
    with pytest.raises(ValueError):
        sv.water_pressure_pa(0.0)
    with pytest.raises(ValueError):
        sv.pure_acid_pressure_pa(-1.0)
