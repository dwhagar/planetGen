# tests/test_venus_spin.py

"""
GEN.1038: the Venus-type spin (`planetgen.physics.venus_spin`): Venus's 243 day retrograde rotation and its fast prograde twin, the Earth that has not reached an
equilibrium, Mercury without a thermal tide, the torque balance at every state, the retrograde threshold, the end state by birth tilt, and the inputs refused.
"""

import math

import pytest

from planetgen.physics import constants
from planetgen.physics import venus_spin as vs

VENUS = {"mass_kg": 4.8675e24, "radius_m": 6.0518e6, "surface_pressure_pa": 92.0e5, "distance_au": 0.723332, "eccentricity": 0.0068,
             "star_mass_kg": constants.SOLAR_MASS_TO_KG, "star_luminosity_w": constants.SOLAR_LUMINOSITY}
EARTH = {"mass_kg": 5.972e24, "radius_m": 6.371e6, "surface_pressure_pa": 1.013e5, "distance_au": 1.0, "eccentricity": 0.0167,
             "star_mass_kg": constants.SOLAR_MASS_TO_KG, "star_luminosity_w": constants.SOLAR_LUMINOSITY}
MERCURY = {"mass_kg": 3.3011e23, "radius_m": 2.4397e6, "surface_pressure_pa": 0.0, "distance_au": 0.387098, "eccentricity": 0.2056,
               "star_mass_kg": constants.SOLAR_MASS_TO_KG, "star_luminosity_w": constants.SOLAR_LUMINOSITY}


def _balance(ratio, strength, e):
    """The paper's equation 54 residual, with the signs of the three ranges, independent of the closed forms in the module."""
    sign = lambda x: math.copysign(1.0, x)
    return ratio - (1.0 + 6.0 * e * e) - strength * ((1.0 - 10.5 * e * e) * sign(ratio - 1.0) - e * e * sign(2.0 * ratio - 1.0)
                                                      + 9.0 * e * e * sign(2.0 * ratio - 3.0))


def test_venus_comes_out_at_its_own_strength_and_damping_time():
    assert vs.strength_over_orbital_rate(**{k: VENUS[k] for k in ("mass_kg", "radius_m", "distance_au", "star_mass_kg", "star_luminosity_w",
                                                                  "surface_pressure_pa")}) == pytest.approx(1.9255, rel=1e-9)
    assert vs.damping_time_gy(VENUS["mass_kg"], VENUS["radius_m"], VENUS["distance_au"], VENUS["star_mass_kg"]) == pytest.approx(2.3, rel=1e-9)


def test_venus_rotates_243_days_backwards_or_77_days_forwards():
    result = vs.venus_type_spin(**VENUS, age_gy=4.5, birth_obliquity_deg=100.0)
    by_label = {s.label: s for s in result.equilibria}
    assert set(by_label) == {"2+", "2-"}
    n_days = 224.7
    assert n_days / by_label["2-"].ratio == pytest.approx(-242.9, abs=0.3)
    assert n_days / by_label["2+"].ratio == pytest.approx(76.8, abs=0.3)
    assert result.reached and result.retrograde and result.final_obliquity_deg == 180.0
    assert result.period_days == pytest.approx(242.9, abs=0.5)


def test_venus_born_with_a_low_tilt_would_have_gone_to_the_fast_prograde_state():
    result = vs.venus_type_spin(**VENUS, age_gy=4.5, birth_obliquity_deg=20.0)
    assert not result.retrograde and result.end_state.label == "2+" and result.final_obliquity_deg == 0.0
    assert result.period_days == pytest.approx(76.8, abs=0.5)


def test_the_earth_has_not_reached_an_equilibrium_in_its_age_and_has_no_retrograde_state():
    result = vs.venus_type_spin(**EARTH, age_gy=4.5, birth_obliquity_deg=120.0)
    assert result.damping_time_gy > 4.5 and not result.reached and result.period_days is None
    assert result.strength < 1.0 and not result.retrograde
    assert all(not s.retrograde for s in result.equilibria)
    assert result.final_obliquity_deg == 120.0           # the tilt it was born with, nothing has moved it


def test_mercury_has_no_thermal_tide_and_its_one_state_is_the_eccentric_gravitational_one():
    result = vs.venus_type_spin(**MERCURY, age_gy=4.5, birth_obliquity_deg=100.0)
    assert result.thin_air and result.strength == 0.0 and result.reached
    assert [s.label for s in result.equilibria] == ["1+"]
    assert result.end_state.ratio == pytest.approx(1.0 + 6.0 * 0.2056 ** 2, rel=1e-12)
    assert not result.retrograde
    # the 3:2 spin-orbit resonance that Mercury is actually in is tidal_spin's, not this module's
    assert result.period_days < 88.0


def test_a_round_orbit_without_a_thermal_tide_is_synchronous():
    states = vs.equilibria(0.0, 0.0)
    assert [(s.label, s.ratio) for s in states] == [("1+", 1.0)]


def test_every_state_balances_the_torques_and_sits_in_its_range():
    for strength in [0.0, 0.05, 0.3, 0.5, 0.55, 0.9, 1.0, 1.1, 1.9255, 2.5, 3.0]:
        for e in [0.0, 0.05, 0.1, 0.2, 0.3]:
            states = vs.equilibria(strength, e)
            assert 1 <= len(states) <= 3
            assert [s.ratio for s in states] == sorted((s.ratio for s in states), reverse=True)
            for s in states:
                assert _balance(s.ratio, strength, e) == pytest.approx(0.0, abs=1e-9), (strength, e, s)
            labels = {s.label for s in states}
            assert not ({"1-", "1+"} <= labels and "2-" in labels)


def test_no_state_is_synchronous_once_there_is_a_thermal_tide():
    for strength in [0.05, 0.5, 1.9255, 3.0]:
        for e in [0.0, 0.1, 0.2]:
            assert all(abs(s.ratio - 1.0) > 1e-6 for s in vs.equilibria(strength, e))


def test_a_retrograde_state_appears_exactly_at_the_published_threshold():
    for e in [0.0, 0.1, 0.2, 0.3]:
        threshold = (1.0 + 6.0 * e * e) / (1.0 - 2.5 * e * e)
        assert vs.retrograde_possible(threshold * 1.001, e) and not vs.retrograde_possible(threshold * 0.999, e)
        above, below = vs.equilibria(threshold * 1.01, e), vs.equilibria(threshold * 0.99, e)
        assert above[-1].retrograde and above[-1].label == "2-"
        assert not any(s.retrograde for s in below)
    # the paper's first-order form, 1 + 17 e^2 / 2 (1.09 at e 0.1), is this threshold's expansion
    assert vs.retrograde_possible(1.08, 0.1) is False and vs.retrograde_possible(1.09, 0.1) is True


def test_the_retrograde_period_is_always_long_in_the_design_notes_sense():
    # weird-worlds.md section 5: W2 spin period over 20 days and sign negative
    for strength in [1.1, 1.5, 1.9255, 2.4]:
        slow = vs.equilibria(strength, 0.0)[-1]
        assert slow.retrograde
        result_days = 224.7 / abs(slow.ratio)
        assert result_days > 20.0


def test_the_birth_tilt_picks_the_state():
    two = vs.equilibria(1.9255, 0.0)
    assert vs.end_state(two, 10.0).label == "2+" and vs.end_state(two, 59.9).label == "2+"
    assert vs.end_state(two, 60.0).label == "2-" and vs.end_state(two, 179.0).label == "2-"
    three = vs.equilibria(0.55, 0.2)
    assert [s.label for s in three] == ["2+", "1+", "1-"]
    assert [vs.end_state(three, d).label for d in (10.0, 45.0, 99.0, 100.0, 170.0)] == ["2+", "1+", "1+", "1-", "1-"]
    one = vs.equilibria(0.0, 0.1)
    assert vs.end_state(one, 150.0) is one[0]


def test_the_strength_scales_with_light_and_distance_as_the_paper_has_it():
    base = vs.strength_over_orbital_rate(4.8675e24, 6.0518e6, 0.723332, constants.SOLAR_MASS_TO_KG, constants.SOLAR_LUMINOSITY, 92.0e5)
    assert vs.strength_over_orbital_rate(4.8675e24, 6.0518e6, 0.723332, constants.SOLAR_MASS_TO_KG, 2 * constants.SOLAR_LUMINOSITY, 92.0e5) == pytest.approx(2 * base)
    assert vs.strength_over_orbital_rate(4.8675e24, 6.0518e6, 2 * 0.723332, constants.SOLAR_MASS_TO_KG, constants.SOLAR_LUMINOSITY, 92.0e5) == pytest.approx(
        base * 2 ** 2.5)
    assert vs.strength_over_orbital_rate(2 * 4.8675e24, 6.0518e6, 0.723332, constants.SOLAR_MASS_TO_KG, constants.SOLAR_LUMINOSITY, 92.0e5) == pytest.approx(base / 2)
    # a thin air has a smaller tide, in proportion to its pressure below 10 bar, and none at all without air
    assert vs.strength_over_orbital_rate(4.8675e24, 6.0518e6, 0.723332, constants.SOLAR_MASS_TO_KG, constants.SOLAR_LUMINOSITY, 1.0e5) == pytest.approx(base / 10)
    assert vs.strength_over_orbital_rate(4.8675e24, 6.0518e6, 0.723332, constants.SOLAR_MASS_TO_KG, constants.SOLAR_LUMINOSITY, 0.0) == 0.0
    assert vs.strength_over_orbital_rate(4.8675e24, 6.0518e6, 0.723332, constants.SOLAR_MASS_TO_KG, constants.SOLAR_LUMINOSITY, 20.0e5) == pytest.approx(base)


def test_the_damping_time_grows_steeply_with_distance_and_falls_with_the_star_mass():
    venus = vs.damping_time_gy(4.8675e24, 6.0518e6, 0.723332, constants.SOLAR_MASS_TO_KG)
    assert vs.damping_time_gy(4.8675e24, 6.0518e6, 2 * 0.723332, constants.SOLAR_MASS_TO_KG) == pytest.approx(64 * venus)
    assert vs.damping_time_gy(4.8675e24, 6.0518e6, 0.723332, 2 * constants.SOLAR_MASS_TO_KG) == pytest.approx(venus / 4)


def test_an_eccentricity_past_the_papers_limit_is_computed_at_the_limit_and_flagged():
    flagged = vs.venus_type_spin(**{**VENUS, "eccentricity": 0.5}, age_gy=4.5, birth_obliquity_deg=100.0)
    limit = vs.venus_type_spin(**{**VENUS, "eccentricity": 0.3}, age_gy=4.5, birth_obliquity_deg=100.0)
    assert flagged.eccentricity_clamped and not limit.eccentricity_clamped
    assert flagged.equilibria == limit.equilibria


def test_inputs_the_model_cannot_use_are_refused():
    with pytest.raises(vs.VenusSpinError):
        vs.strength_over_orbital_rate(0.0, 6.0e6, 1.0, constants.SOLAR_MASS_TO_KG, constants.SOLAR_LUMINOSITY, 1.0e5)
    with pytest.raises(vs.VenusSpinError):
        vs.equilibria(1.0, 0.4)
    with pytest.raises(vs.VenusSpinError):
        vs.equilibria(-1.0, 0.1)
    with pytest.raises(vs.VenusSpinError):
        vs.end_state(vs.equilibria(1.9255, 0.0), 200.0)
