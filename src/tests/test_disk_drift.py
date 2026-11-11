"""Pebble drift and evaporation at the ice lines of a disc (GEN.1473, physics/disk_drift.py)."""

import math

import pytest

from planetgen import tuning
from planetgen.physics import disk_condensation as dc, disk_drift as dd, mixture_giant

SUN_SNOW_LINE_AU = 2.3


def test_the_front_and_its_inverse_agree_and_the_sun_ends_its_pebbles_at_the_edge():
    assert dd.front_au(1.0) == pytest.approx(dd.FRONT_AU_AT_1_MYR)
    for distance in (0.5, 5.0, 35.0, 120.0):
        assert dd.front_au(dd.arrival_myr(distance)) == pytest.approx(distance)
    assert dd.front_au(1.0, 8.0) == pytest.approx(2.0 * dd.FRONT_AU_AT_1_MYR)           # the star's mass to the 1/3
    edge = tuning.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER * SUN_SNOW_LINE_AU
    assert dd.pebble_end_myr(SUN_SNOW_LINE_AU) == pytest.approx(dd.arrival_myr(edge))
    assert 0.1 < dd.pebble_end_myr(SUN_SNOW_LINE_AU) < 0.4


def test_the_enrichment_ratio_is_the_design_notes_two_to_three_for_the_sun():
    for age in (0.05, 0.1):
        assert 1.8 < dd.instant_ratio(age) < 2.8
    assert dd.instant_ratio(0.05) == pytest.approx(2.27, rel=0.03)


def test_the_disc_mass_cancels_in_the_ratio_so_the_ratio_does_not_depend_on_the_star_mass():
    assert dd.instant_ratio(0.3, mass_solar=0.3) == pytest.approx(dd.instant_ratio(0.3, mass_solar=1.0))
    assert dd.instant_ratio(0.3, mass_solar=1.5) == pytest.approx(dd.instant_ratio(0.3, mass_solar=1.0))


def test_a_richer_disc_has_more_pebbles_per_gas_to_the_two_thirds():
    assert dd.instant_ratio(0.1, metal_ratio=2.0) == pytest.approx(dd.instant_ratio(0.1) * 2.0 ** (2.0 / 3.0), rel=1e-6)


def test_a_species_has_no_pebbles_before_the_front_arrives_or_after_they_end():
    arrival = dd.arrival_myr(dd.line_au_of("CO2", SUN_SNOW_LINE_AU))
    end = dd.pebble_end_myr(SUN_SNOW_LINE_AU)
    assert 0.0 < arrival < end
    assert dd.enrichment_ratio("CO2", 0.5 * arrival, SUN_SNOW_LINE_AU) == 0.0
    assert dd.enrichment_ratio("CO2", 0.5 * (arrival + end), SUN_SNOW_LINE_AU) > 1.0
    assert dd.enrichment_ratio("CO2", end, SUN_SNOW_LINE_AU) == 0.0
    assert dd.enrichment_ratio("CO2", 5.0, SUN_SNOW_LINE_AU) == 0.0
    assert dd.enrichment_ratio("H2O", 0.01, SUN_SNOW_LINE_AU) > 1.0                  # the water line is crossed first


def test_a_line_outside_the_disc_edge_gets_no_pebbles():
    assert dd.line_au_of("CO", SUN_SNOW_LINE_AU) > dd.edge_au(SUN_SNOW_LINE_AU)
    assert dd.enrichment_ratio("CO", 0.1, SUN_SNOW_LINE_AU) == 0.0
    with pytest.raises(ValueError):
        dd.enrichment_ratio("Unobtainium", 0.1, SUN_SNOW_LINE_AU)


def test_the_lag_is_a_boxcar_of_the_ratio_that_never_looks_before_the_start_of_the_disc():
    line = dd.line_au_of("CO2", SUN_SNOW_LINE_AU)
    place = 5.0
    lag = dd.lag_myr(line, place)
    assert lag == pytest.approx(dd.LAG_MYR_PER_AU * (line - place))
    assert dd.lag_myr(5.0, 35.0) == 0.0                                              # a place outside the line has no lag
    first, end = dd.arrival_myr(line), dd.pebble_end_myr(SUN_SNOW_LINE_AU)
    # before the lagged signal could come (window start at the start of the disc, pebbles not yet at the line): zero
    assert dd.enrichment_ratio("CO2", 0.5 * first, SUN_SNOW_LINE_AU, to_au=place) == 0.0
    # after the pebbles end, the window still holds them: the lagged value is nonzero for as long as the window overlaps the epoch
    late = end + 0.5 * lag
    assert dd.enrichment_ratio("CO2", late, SUN_SNOW_LINE_AU, to_au=place) > 0.0
    assert dd.enrichment_ratio("CO2", end + 2.0 * lag, SUN_SNOW_LINE_AU, to_au=place) == 0.0
    # the mean never exceeds the largest instant ratio in the epoch
    top = max(dd.instant_ratio(t) for t in (first, 0.5 * (first + end), end))
    assert 0.0 < dd.enrichment_ratio("CO2", end, SUN_SNOW_LINE_AU, to_au=place) < top


def test_the_lagged_mean_of_a_constant_window_is_the_instant_value():
    # a place so near the line that the lag is tiny: the mean over the window equals the instant ratio
    line = dd.line_au_of("CO2", SUN_SNOW_LINE_AU)
    age = 0.5 * (dd.arrival_myr(line) + dd.pebble_end_myr(SUN_SNOW_LINE_AU))
    assert dd.enrichment_ratio("CO2", age, SUN_SNOW_LINE_AU, to_au=line - 1e-4) == pytest.approx(dd.enrichment_ratio("CO2", age, SUN_SNOW_LINE_AU), rel=1e-3)


@pytest.mark.parametrize("temperature", [30.0, 45.0, 60.0, 100.0, 150.0, 182.0, 250.0, 450.0])
def test_weights_of_one_give_the_static_gas_and_the_static_c_over_o_exactly(temperature):
    ones = {species: 1.0 for species in dc.LINES_K}
    assert dc._gas_and_solid_atoms(temperature, 1.0, 1.0, 1.0, 1.0, 1.0, ones) == dc._gas_and_solid_atoms(temperature, 1.0, 1.0, 1.0, 1.0, 1.0)
    for enrichment in (0.0, 1.0, 5.0, 40.0):
        assert dc.envelope_c_over_o(temperature, enrichment, ice_weights=ones) == dc.envelope_c_over_o(temperature, enrichment)


def test_no_weights_is_the_static_disc():
    assert dc.envelope_c_over_o(120.0, 3.0, ice_weights=None) == dc.envelope_c_over_o(120.0, 3.0, ice_weights={})


def test_weights_of_zero_leave_the_supply_gas_and_more_than_one_adds_vapor_inside_the_line():
    temperature = 300.0                                                                # warmer than the water and ammonia lines
    static = dc._gas_and_solid_atoms(temperature, 1.0, 1.0, 1.0, 1.0, 1.0)
    none_carried = dc._gas_and_solid_atoms(temperature, 1.0, 1.0, 1.0, 1.0, 1.0, {"H2O": 0.0})
    more = dc._gas_and_solid_atoms(temperature, 1.0, 1.0, 1.0, 1.0, 1.0, {"H2O": 2.0})
    water = dc.inventory().moles["H2O"]
    assert none_carried[1] == pytest.approx(static[1] - water)                         # no water pebbles: the oxygen of the water is not in the gas
    assert more[1] == pytest.approx(static[1] + water)
    assert none_carried[0] == static[0] == more[0]                                     # water carries no carbon
    assert none_carried[2:] == static[2:] == more[2:]                                  # the solids stay the static ones


def test_the_carbon_and_oxygen_of_a_species_in_the_gas_follow_its_molecule():
    temperature = 100.0                                                                # warmer than the CO2 line, colder than the water line
    static = dc._gas_and_solid_atoms(temperature, 1.0, 1.0, 1.0, 1.0, 1.0)
    weighted = dc._gas_and_solid_atoms(temperature, 1.0, 1.0, 1.0, 1.0, 1.0, {"CO2": 3.0})
    co2 = dc.inventory().moles["CO2"]
    assert weighted[0] - static[0] == pytest.approx(2.0 * co2)
    assert weighted[1] - static[1] == pytest.approx(4.0 * co2)


def test_a_frozen_ice_is_not_put_in_the_gas_by_a_weight():
    temperature = 100.0                                                                # water is frozen here
    assert dc._gas_and_solid_atoms(temperature, 1.0, 1.0, 1.0, 1.0, 1.0, {"H2O": 5.0}) == dc._gas_and_solid_atoms(temperature, 1.0, 1.0, 1.0, 1.0, 1.0)


def test_the_carbon_and_oxygen_in_the_gas_and_the_solids_stay_the_budget_for_the_static_disc_and_the_envelope_holds_the_enrichment():
    held = dc.inventory()
    gas_c, gas_o, solid_c, solid_o = dc._gas_and_solid_atoms(100.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    assert gas_c + solid_c == pytest.approx(held.carbon) and gas_o + solid_o == pytest.approx(held.oxygen)
    # drift adds vapor from pebbles that came from outside: the envelope still holds exactly the enrichment of the budget when it is above the gas's share
    weights = {"CO2": 2.0}
    enriched = dc._gas_and_solid_atoms(100.0, 1.0, 1.0, 1.0, 1.0, 1.0, weights)
    sigma = dc._planetesimal_share(10.0, *enriched, metals_star=held.carbon + held.oxygen)
    assert enriched[0] + enriched[1] + sigma * (enriched[2] + enriched[3]) == pytest.approx(10.0 * (held.carbon + held.oxygen))


def test_inside_the_water_line_the_pebble_epoch_gas_is_oxygen_rich_and_the_late_gas_is_the_supply():
    early = dd.ice_weights(300.0, 0.05, SUN_SNOW_LINE_AU)
    assert early["H2O"] > 1.0
    late = dd.ice_weights(300.0, 5.0, SUN_SNOW_LINE_AU)
    assert late["H2O"] == 0.0
    c_static = dc.envelope_c_over_o(300.0, 0.0)
    c_early = dc.envelope_c_over_o(300.0, 0.0, ice_weights=early)
    c_late = dc.envelope_c_over_o(300.0, 0.0, ice_weights=late)
    assert c_early < c_static < c_late                                                 # the notes' direction: low early, high late
    assert c_early < 0.35


def test_ice_weights_hold_only_the_lines_between_the_edge_and_the_place():
    weights = dd.ice_weights(120.0, 0.1, SUN_SNOW_LINE_AU)
    assert "H2O" not in weights and "NH3" not in weights                               # still frozen at 120 K
    assert "CO2" in weights and "H2S" in weights
    assert "CH4" not in weights and "CO" not in weights                                # colder than the edge: already in the supply gas
    assert dd.ice_weights(dd.edge_temperature_k() - 1.0, 0.1, SUN_SNOW_LINE_AU) == {}


def test_a_giant_gathers_its_gas_when_the_pebbles_have_built_the_critical_core_and_not_later_than_the_pebbles_end():
    age = dd.gathering_age_myr(150.0, SUN_SNOW_LINE_AU)
    place = dd.place_au(150.0, SUN_SNOW_LINE_AU)
    assert dd.arrival_myr(place) < age < dd.pebble_end_myr(SUN_SNOW_LINE_AU)
    # the pebble mass that passed the place by then is the core
    coefficient = (dd.PEBBLE_FLUX_EARTH_MASSES_PER_YR * (dd.solids_over_gas() / dd.REFERENCE_SOLIDS_OVER_GAS) ** dd.PEBBLE_FLUX_Z_EXPONENT * dd.MYR_YR * 1.5)
    built = coefficient * (age ** (2.0 / 3.0) - dd.arrival_myr(place) ** (2.0 / 3.0))
    assert built == pytest.approx(dd.CRITICAL_CORE_EARTH, rel=1e-6)


def test_a_core_the_pebbles_cannot_build_gathers_gas_as_the_pebbles_end_and_sees_no_pebbles():
    # a far place where the front arrives only at the end of the epoch
    temperature = dd.edge_temperature_k() + 0.5
    place = dd.place_au(temperature, SUN_SNOW_LINE_AU)
    assert dd.arrival_myr(place) < dd.pebble_end_myr(SUN_SNOW_LINE_AU)
    age = dd.gathering_age_myr(temperature, SUN_SNOW_LINE_AU)
    assert age == pytest.approx(dd.pebble_end_myr(SUN_SNOW_LINE_AU))
    assert dd.enrichment_ratio("H2O", age, SUN_SNOW_LINE_AU) == 0.0


def test_the_gathering_age_is_capped_at_the_disc_lifetime_and_is_a_function_of_the_inputs_only():
    from planetgen.physics import scatter
    for mass in (0.1, 0.3, 1.0, 1.5, 2.5):
        age = dd.gathering_age_myr(200.0, SUN_SNOW_LINE_AU * mass, mass)
        assert 0.0 < age <= scatter.DISK_LIFETIME_YR / dd.MYR_YR
        assert age == dd.gathering_age_myr(200.0, SUN_SNOW_LINE_AU * mass, mass)


def test_the_epochs_of_a_dwarf_and_of_a_heavier_star_follow_the_front():
    dwarf = dd.pebble_end_myr(0.3, 0.3)
    heavy = dd.pebble_end_myr(4.0, 1.5)
    assert dwarf < 0.05 < dd.pebble_end_myr(SUN_SNOW_LINE_AU) < heavy


def test_the_cold_finger_enrichment_is_one_over_one_minus_the_leak_cut_by_the_ice_share():
    assert dd.cold_finger_enrichment(1.0) == pytest.approx(1.0 / (1.0 - dd.COLD_FINGER_LEAK))
    assert dd.cold_finger_enrichment(0.0) == 1.0
    assert 1.0 < dd.cold_finger_enrichment(0.4) < dd.cold_finger_enrichment(1.0)
    with pytest.raises(ValueError):
        dd.cold_finger_enrichment(1.5)


def test_a_giant_subject_with_a_snow_line_takes_the_drift_weights_and_one_without_keeps_the_static_disc():
    subject = {"formation_temperature_k": 300.0, "disk_snow_line_au": SUN_SNOW_LINE_AU, "host_mass_solar": 1.0}
    assert mixture_giant._disk_ice_weights(subject, 300.0)["H2O"] > 1.0
    assert mixture_giant._disk_ice_weights({"formation_temperature_k": 300.0}, 300.0) is None
    drift = mixture_giant._disk_c_over_o(subject, 0.3)
    static = mixture_giant._disk_c_over_o({"formation_temperature_k": 300.0}, 0.3)
    assert drift < static
    assert math.isfinite(drift)


def test_a_giant_subject_with_a_snow_line_gives_its_nitrogen_sulfur_and_noble_gases_the_drift_weights_too():
    # GEN.1516: `_disk_enrichments` passes the same weights `_disk_c_over_o` does; a subject without a snow line keeps the static factors.
    subject = {"formation_temperature_k": 250.0, "disk_snow_line_au": SUN_SNOW_LINE_AU, "host_mass_solar": 1.0}
    drift = mixture_giant._disk_enrichments(subject, 0.6)
    static = mixture_giant._disk_enrichments({"formation_temperature_k": 250.0}, 0.6)
    assert static == dc.envelope_enrichments(250.0, 10.0 ** 0.6, 1.0, 1.0, 1.0, 1.0, 1.0)
    assert drift != static
    weights = mixture_giant._disk_ice_weights(subject, 250.0)
    assert drift["Ar"] == pytest.approx(weights["Ar"]) and drift["Xe"] == pytest.approx(weights["Xe"])
    assert all(math.isfinite(value) and value >= 0.0 for value in drift.values())


def test_a_giant_whose_gas_holds_no_argon_still_makes_a_mixture():
    from planetgen.physics import gas_equilibrium
    budget = mixture_giant.element_budget(0.5, disk_factors={"N": 1.5, "S": 1.0, "Ar": 0.0, "Kr": 0.0, "Xe": 0.0})
    assert budget["Ar"] == 0.0 and budget["H"] > 0.0
    result = gas_equilibrium.solve({symbol: value for symbol, value in budget.items() if value > 0.0}, 800.0, 1.0)
    assert sum(result.mole_fractions.values()) == pytest.approx(1.0)
