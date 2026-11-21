# tests/test_sealed_fuel.py

"""GEN.1825: the chemical life of an ocean under an ice lid is held by the moles of its reductants, not by its heat flow's share of Earth's chemosynthesis."""

import pytest

from planetgen import tuning
from planetgen.physics import life_classes, life_supply, ocean_oxygen

ICE = {"hydrosphere": "ice-covered ocean", "ocean_depth_km": None, "surface_dose_msv_yr": 0.0, "internal_heat_flux_w_m2": 0.0614, "tidal_flux_w_m2": 0.0,
       "h2_supply": 1.0, "donor_supply": 1.0}
"""A sealed ocean with Earth's heat flow and Earth's hydrothermal reductants (both shares 1), no water depth counted for the radiolysis."""

CHEMICAL = ("3a", "3b", "4a", "4b", "4c", "8")


def production(class_id, values):
    return life_classes.net_primary_production(life_classes.describe(class_id), {}, None, values)


def test_the_fuel_of_earths_vents_fixes_3e_minus_5_of_earths_production():
    expected = 0.3 * ocean_oxygen.H2_O2_EQUIVALENT_MOL_M2_YR / 11.5
    assert life_classes.sealed_fuel_npp(ICE) == pytest.approx(expected)
    assert expected == pytest.approx(3.61e-5, rel=0.01)
    for class_id in CHEMICAL:
        assert production(class_id, ICE) == pytest.approx(expected)               # the heat share's 1e-3 is the larger, so the fuel is the cap


def test_the_cap_follows_the_reductants_and_the_radiolysis_adds_to_it():
    half = {**ICE, "h2_supply": 0.5, "donor_supply": 0.5}
    assert life_classes.sealed_fuel_npp(half) == pytest.approx(life_classes.sealed_fuel_npp(ICE) / 2.0)
    none = {**ICE, "h2_supply": 0.0, "donor_supply": 0.0}
    assert life_classes.sealed_fuel_npp(none) == 0.0
    deep = {**none, "ocean_depth_km": 100.0}
    assert life_classes.sealed_fuel_npp(deep) == pytest.approx(0.3 * ocean_oxygen.radiolysis(deep) / 11.5)
    assert life_classes.sealed_fuel_npp(deep) > 0.0                             # the water's own radiolysis fuels a sealed ocean with no vents


def test_a_sealed_ocean_with_more_heat_than_earths_is_not_given_more_fuel_by_it():
    hot = {**ICE, "internal_heat_flux_w_m2": 0.6}
    assert production("3a", hot) == pytest.approx(production("3a", ICE))       # the share of reductants is capped at Earth's (`life_supply.h2_supply`), and so is the cap
    assert tuning.LIFE_NPP["chemosynthesis_share"] * life_supply.heat_share(hot) > production("3a", hot)


def test_an_open_ocean_and_a_dry_body_keep_the_heat_flows_share():
    for state in ("surface ocean", "dry", "brine seep", "hycean"):
        values = {**ICE, "hydrosphere": state}
        assert production("3a", values) == pytest.approx(tuning.LIFE_NPP["chemosynthesis_share"] * life_supply.heat_share(values))


def test_a_class_that_is_not_chemical_does_not_read_the_fuel():
    assert life_classes.describe("1a").energy["source"] != "chemical"
    assert "sealed_fuel_efficiency" in life_supply.NUMBERS and "ocean_production_mol_c_m2_yr" in life_supply.NUMBERS
