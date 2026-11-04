# tests/test_agb_carbon.py

"""GEN.1026: carbon stars and S stars from the thermal pulses of the asymptotic giant branch (`planetgen/physics/agb_carbon.py`)."""

import pytest

from planetgen.physics import agb_carbon, post_main_sequence

AGB = "asymptotic-giant"


def _kind_at(mass, fe_h, share):
    line = post_main_sequence.timeline(mass, fe_h)
    window = line.window(post_main_sequence.ASYMPTOTIC_GIANT)
    age = window.start_gy + share * (window.end_gy - window.start_gy)
    return agb_carbon.kind_of(AGB, mass, line.state(age).mass_sol, fe_h, age)


def test_the_pulse_and_dredge_tables_are_the_published_rows():
    assert agb_carbon.pulses_in_branch(2.75) == 25.0 and agb_carbon.pulses_in_branch(3.25) == 22.0
    assert agb_carbon.pulses_in_branch(1.0) == 17.0    # held below the table
    assert agb_carbon.pulses_in_branch(2.25) == pytest.approx(22.0)    # linear between 1.75 (19) and 2.75 (25)


def test_the_dredged_mass_grows_with_the_pulses_and_is_zero_before_the_first():
    grown = [agb_carbon.dredged_core_mass_sol(2.5, n) for n in range(0, 25)]
    assert grown[0] == 0.0 and all(b > a for a, b in zip(grown, grown[1:]))


def test_a_star_that_is_too_light_or_has_run_no_pulse_keeps_the_envelope_c_to_o():
    ratio = agb_carbon.surface_c_to_o(1.2, 0.0, 1.1, 0.9, 12)
    assert ratio == pytest.approx(10.0 ** (8.46 - 8.69), rel=1.0e-9)
    assert agb_carbon.surface_c_to_o(2.5, 0.0, 2.3, 1.5, 0) == ratio


def test_c_to_o_rises_with_the_pulses_run():
    ratios = [agb_carbon.surface_c_to_o(2.5, 0.0, 2.3, 1.2, n) for n in range(0, 24, 3)]
    assert all(b > a for a, b in zip(ratios, ratios[1:]))


def test_carbon_stars_are_the_dredging_stars_below_the_hot_bottom_burning_mass_at_solar_metals():
    assert _kind_at(1.2, 0.0, 0.99) is None    # too light to dredge
    for mass in (1.6, 2.0, 2.5, 3.0, 3.5):
        assert _kind_at(mass, 0.0, 0.99) == agb_carbon.CARBON_STAR, mass
    assert _kind_at(4.5, 0.0, 0.99) is None    # hot bottom burning


def test_a_star_is_oxygen_rich_early_in_the_branch_and_carbon_rich_late():
    assert _kind_at(2.5, 0.0, 0.02) is None
    assert _kind_at(2.5, 0.0, 0.99) == agb_carbon.CARBON_STAR


def test_an_s_star_comes_between_the_oxygen_rich_phase_and_the_carbon_star_phase():
    seen = [_kind_at(2.5, 0.0, share / 100.0) for share in range(5, 100)]
    first_s, first_c = seen.index(agb_carbon.S_STAR), seen.index(agb_carbon.CARBON_STAR)
    assert seen[0] is None and first_s < first_c
    assert all(kind in (None, agb_carbon.S_STAR, agb_carbon.CARBON_STAR) for kind in seen)


def test_a_metal_poor_star_dredges_at_lower_mass_and_turns_carbon_rich_sooner():
    assert agb_carbon.tdu_minimum_mass_sol(-1.0) < agb_carbon.tdu_minimum_mass_sol(0.0)
    assert agb_carbon.tdu_minimum_mass_sol(-5.0) == agb_carbon.TDU_MIN_MASS_FLOOR_SOL
    assert agb_carbon.tdu_minimum_mass_sol(0.5) == agb_carbon.TDU_MIN_MASS_AT_SOLAR_SOL    # no help from metal-rich
    assert _kind_at(1.3, -1.0, 0.99) == agb_carbon.CARBON_STAR and _kind_at(1.3, 0.0, 0.99) is None


def test_a_star_off_the_branch_or_without_its_values_is_neither():
    assert agb_carbon.kind_of("red-giant", 2.5, 1.5, 0.0, 0.5) is None
    assert agb_carbon.kind_of(AGB, 2.5, 1.5, None, 0.5) is None
    assert agb_carbon.kind_of(AGB, None, 1.5, 0.0, 0.5) is None


def test_the_initial_mass_is_found_from_the_lifespan_a_star_stores():
    for mass in (1.5, 2.5, 3.5):
        for fe_h in (0.0, -1.0):
            lifespan = post_main_sequence.timeline(mass, fe_h).end_gy
            assert agb_carbon.initial_mass_from_lifespan_sol(lifespan, fe_h) == pytest.approx(mass, rel=1.0e-6)
    assert agb_carbon.initial_mass_from_lifespan_sol(None, 0.0) is None
    assert agb_carbon.initial_mass_from_lifespan_sol(1.0e-3, 0.0) is None    # no star of 8 solar masses or less lives so briefly
