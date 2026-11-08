"""The species-by-species disc condensation (GEN.1363 part A, physics/disc_condensation.py)."""

import pytest

from planetgen.physics import condensate_shares, disc_condensation as dc

SUN_L = 3.828e26


def test_the_suns_budget_per_ten_thousand_hydrogen_is_the_asplund_one():
    assert dc.solar_per_1e4_h("C") == pytest.approx(2.884, rel=1e-3)
    assert dc.solar_per_1e4_h("O") == pytest.approx(4.898, rel=1e-3)
    assert dc.solar_per_1e4_h("N") == pytest.approx(0.676, rel=1e-3)
    assert dc.solar_per_1e4_h("Kr") == pytest.approx(1.32e-5, rel=1e-2) and dc.solar_per_1e4_h("Xe") == pytest.approx(1.66e-6, rel=1e-2)
    with pytest.raises(ValueError):
        dc.solar_per_1e4_h("Unobtainium")


def test_the_carbon_is_split_between_grains_and_the_gases_and_none_is_lost():
    held = dc.inventory()
    moles = held.moles
    assert held.refractory_carbon == pytest.approx(0.5 * held.carbon)
    assert held.refractory_carbon + moles["CO"] + moles["CO2"] + moles["CH4"] == pytest.approx(held.carbon)
    assert moles["CO"] / moles["CO2"] == pytest.approx(5.0)                         # Oberg's CO 1.5 against CO2 0.3
    assert 0.03 < moles["CH4"] / (held.carbon - held.refractory_carbon) < 0.06


def test_the_oxygen_of_the_rock_the_co_and_the_co2_is_out_of_the_water_and_none_is_lost():
    held = dc.inventory()
    moles = held.moles
    assert held.rock_oxygen + moles["CO"] + 2.0 * moles["CO2"] + moles["H2O"] == pytest.approx(held.oxygen)
    short = dc.inventory(oxygen=0.3)                                                # barely more than the rock's share: the water is gone, not negative
    assert short.moles["H2O"] == 0.0 and dc.inventory(oxygen=0.0).moles["H2O"] == 0.0


def test_the_nitrogen_and_sulfur_are_split_and_the_noble_gases_follow_their_ratio():
    held = dc.inventory(nitrogen=2.0, sulfur=0.5, noble=3.0)
    n = dc.solar_per_1e4_h("N") * 2.0
    assert 2.0 * held.moles["N2"] + held.moles["NH3"] + dc.ammonium_salts(nitrogen=2.0) == pytest.approx(n)
    assert held.moles["H2S"] == pytest.approx(0.3 * dc.solar_per_1e4_h("S") * 0.5)
    assert held.moles["Ar"] == pytest.approx(3.0 * dc.solar_per_1e4_h("Ar"))


def test_each_ice_freezes_below_its_line_and_the_lines_run_warm_to_cold():
    lines = list(dc.LINES_K.values())
    assert lines == sorted(lines, reverse=True) and dc.ORDER[0] == "H2O" and dc.ORDER[-1] == "CO"
    assert dc.LINES_K["H2O"] == condensate_shares.ICE_CONDENSATION_K["H2O"]
    assert set(dc.condensed(600.0)) == set()                                         # whole rock, nothing frozen above the grains' line
    assert set(dc.condensed(300.0)) == {"C"}
    assert set(dc.condensed(150.0)) == {"C", "H2O"}
    assert {"NH3", "NH4"} <= set(dc.condensed(100.0)) and "H2S" not in dc.condensed(100.0)
    assert "H2S" in dc.condensed(70.0) and "Xe" not in dc.condensed(70.0)
    for temperature, frozen in ((60.0, "Xe"), (50.0, "Kr"), (47.5, "Ar"), (45.0, "CO2"), (35.0, "CH4"), (21.0, "N2"), (19.0, "CO")):
        assert frozen in dc.condensed(temperature)
    assert "CO" not in dc.condensed(21.0) and "N2" not in dc.condensed(23.0)
    assert set(dc.condensed(10.0)) == set(dc.ORDER) | {"NH4", "C"}


def test_a_species_a_star_lacks_does_not_freeze():
    assert "H2O" not in dc.condensed(100.0, oxygen=0.2)
    assert "CO" not in dc.condensed(10.0, carbon=0.0)


def test_masses_are_percent_of_the_nebula_on_the_water_of_table_11():
    held = dc.inventory(carbon=0.0, oxygen=1.0, nitrogen=0.0, sulfur=0.0, noble=0.0)
    masses = dc.masses_percent(150.0, carbon=0.0, nitrogen=0.0, sulfur=0.0, noble=0.0)
    assert masses["rock"] == condensate_shares.TABLE_11_MASS_PERCENT["rock"]
    assert masses["H2O"] == pytest.approx(condensate_shares.TABLE_11_MASS_PERCENT["H2O"])        # all the oxygen the rock leaves, no CO to take any
    assert held.moles["H2O"] > 3.7
    full = dc.masses_percent(150.0)
    assert full["H2O"] < masses["H2O"]                                                        # the CO and CO2 took oxygen from the water
    assert full["C"] == pytest.approx(dc.inventory().refractory_carbon * 12.011 / 18.015 * dc.PERCENT_PER_H2O_MOLE)


def test_the_ammonia_of_a_solar_disc_is_near_table_11_less_the_nitrogen_that_is_n2():
    ammonia = dc.masses_percent(100.0)["NH3"]
    assert ammonia == pytest.approx(condensate_shares.TABLE_11_MASS_PERCENT["NH3"] * (1.0 - dc.NITROGEN_AS_N2 - dc.NITROGEN_AS_SALTS), rel=0.05)


def test_the_gas_c_over_o_rises_outward_past_the_snow_lines_and_the_solids_stay_oxygen_rich():
    sun_total = dc.solar_per_1e4_h("C") / dc.solar_per_1e4_h("O")
    hot_gas, hot_solid = dc.gas_and_solid_c_over_o(600.0)
    assert hot_gas > sun_total and hot_solid == 0.0                                          # the rock took oxygen and nothing took carbon
    warm_gas, warm_solid = dc.gas_and_solid_c_over_o(150.0)                                  # water frozen: the gas is richer in carbon again
    assert warm_gas > hot_gas > sun_total
    co2_gas, _solid = dc.gas_and_solid_c_over_o(40.0)
    assert co2_gas > warm_gas                                                                # CO2 frozen out too
    for temperature in (600.0, 150.0, 40.0, 25.0):
        gas, solid = dc.gas_and_solid_c_over_o(temperature)
        assert solid is None or solid < 1.0                                                   # solids are oxygen-rich
    assert dc.gas_and_solid_c_over_o(10.0)[0] is None                                        # everything is frozen: no gas left


def test_the_c_over_o_budget_closes_at_every_temperature():
    held = dc.inventory()
    for temperature in (600.0, 300.0, 150.0, 100.0, 60.0, 45.0, 35.0, 25.0, 15.0):
        frozen = dc.condensed(temperature)
        carbon_solid = frozen.get("C", 0.0) + frozen.get("CO", 0.0) + frozen.get("CO2", 0.0) + frozen.get("CH4", 0.0)
        gas, solid = dc.gas_and_solid_c_over_o(temperature)
        oxygen_solid = held.rock_oxygen + frozen.get("H2O", 0.0) + frozen.get("CO", 0.0) + 2.0 * frozen.get("CO2", 0.0)
        if solid is not None:
            assert solid == pytest.approx(carbon_solid / oxygen_solid)
        if gas is not None:
            assert gas == pytest.approx((held.carbon - carbon_solid) / (held.oxygen - oxygen_solid))


def test_the_lines_in_au_move_out_with_the_star_and_are_ordered():
    water, carbon_monoxide = dc.line_au("H2O", SUN_L), dc.line_au("CO", SUN_L)
    assert water == pytest.approx(condensate_shares.ice_line_au(SUN_L, "H2O"), rel=1e-9)
    assert carbon_monoxide > dc.line_au("CO2", SUN_L) > water
    assert dc.line_au("H2O", 4 * SUN_L) == pytest.approx(2.0 * water, rel=1e-9)
    with pytest.raises(ValueError):
        dc.line_au("Ne", SUN_L)


def test_a_negative_ratio_or_temperature_is_refused():
    with pytest.raises(ValueError):
        dc.inventory(carbon=-1.0)
    with pytest.raises(ValueError):
        dc.condensed(0.0)


def test_the_envelope_of_an_unenriched_giant_is_the_star_and_an_enriched_one_takes_the_oxygen_rich_solids():
    star = dc.envelope_c_over_o(100.0, 1.0)
    gas, solid = dc.gas_and_solid_c_over_o(100.0)
    assert star == pytest.approx(dc.inventory().carbon / dc.inventory().oxygen)       # gas and solids in the disc's own proportion: the star's C/O
    rich = dc.envelope_c_over_o(100.0, 60.0)
    assert rich == pytest.approx(solid, rel=0.1) and rich < star
    assert dc.envelope_c_over_o(100.0, 0.05) == pytest.approx(gas)                   # less than the gas's own share of the metals: gas alone


def test_the_envelope_c_over_o_follows_the_star_and_refuses_a_negative_enrichment():
    assert dc.envelope_c_over_o(100.0, 1.0, carbon=1.2) > dc.envelope_c_over_o(100.0, 1.0) > dc.envelope_c_over_o(100.0, 1.0, carbon=0.8)
    with pytest.raises(ValueError):
        dc.envelope_c_over_o(100.0, -1.0)
    with pytest.raises(ValueError):
        dc.envelope_c_over_o(0.0, 1.0)
