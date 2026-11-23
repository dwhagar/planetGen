"""GEN.1821: the carbon a biosphere buries and what it does to the CO2 of its world (`physics.carbon_burial`, `surface_pressure.burial_rate`, the levelled weathering of `continental_weathering` and `waterworld_carbon`)."""

from types import SimpleNamespace

import pytest

from planetgen.physics import biogenic_oxygen, carbon_burial, continental_weathering as cw, life_air, life_supply, surface_pressure as sp, waterworld_carbon as wc

EARTH_LAW = wc.SeafloorLaw(activation_j_mol=75.0e3, ph_exponent=0.27, spreading_exponent=0.47)
Q = sp.HEAT_FLUX_EARTH_W_M2


def _balance(burial=0.0, supply_bar=10.0, **kwargs):
    return cw.pco2_bar(supply_bar, Q, sp.LAND_FRACTION_EARTH, sp.TS_MODERN_K, cw.EARTH_ABSORBED_FLUX_W_M2, (1.0, 1.5), EARTH_LAW, burial)


def test_the_earths_net_burial_share_is_its_reductant_sink_over_its_outgassing():
    assert sp.EARTH_NET_BURIAL_SHARE == pytest.approx(3.0 / 6.6)
    assert sp.EARTH_NET_BURIAL_SHARE == pytest.approx(biogenic_oxygen.EARTH_REDUCTANT_MOL_YR / carbon_burial.OUTGASSING_EARTH_MOL_YR)
    assert sp.LIFELESS_WEATHERING_SHARE == pytest.approx(1.0 - sp.EARTH_NET_BURIAL_SHARE)


def test_the_earth_with_its_biosphere_keeps_its_modern_pressure_and_the_earth_without_life_holds_more_carbon_dioxide():
    assert _balance(burial=sp.EARTH_NET_BURIAL_SHARE) == pytest.approx(sp.CO2_MODERN_BAR, rel=1e-6)
    assert _balance(burial=0.0) > 1.3 * sp.CO2_MODERN_BAR                      # weathering has to take 1.83 times its modern rate
    assert _balance(burial=0.0) > _balance(burial=0.2) > _balance(burial=sp.EARTH_NET_BURIAL_SHARE)    # the more a biosphere buries, the less carbon dioxide is left


def test_a_biosphere_that_buries_more_than_earths_lowers_the_carbon_dioxide_under_the_earths_and_the_burial_never_drains_it_under_the_floor():
    assert _balance(burial=0.9) < sp.CO2_MODERN_BAR
    huge = _balance(burial=5.0)                                                # more than the volcanoes give: the floor of the producers holds
    assert huge >= sp.BURIAL_FLOOR_BAR * (1.0 - 1e-6)
    assert huge < 2.0 * sp.BURIAL_FLOOR_BAR * sp.BURIAL_FULL_FACTOR


def test_the_burial_rises_in_a_line_from_nothing_at_half_a_pascal_to_the_whole_at_five_pascal():
    floor, full = sp.BURIAL_FLOOR_BAR, sp.BURIAL_FLOOR_BAR * sp.BURIAL_FULL_FACTOR
    assert floor == pytest.approx(0.5e-5)                                      # 0.5 Pa, the CO2 gate of class 1a
    assert sp.burial_rate(0.4, 0.5 * floor) == 0.0 and sp.burial_rate(0.4, floor) == 0.0
    assert sp.burial_rate(0.4, 0.5 * (floor + full)) == pytest.approx(0.2)
    assert sp.burial_rate(0.4, full) == pytest.approx(0.4) and sp.burial_rate(0.4, 10.0 * full) == pytest.approx(0.4)
    assert sp.burial_rate(0.0, full) == 0.0 and sp.burial_rate(-1.0, full) == 0.0


def test_the_seafloor_sink_of_a_world_with_no_land_takes_the_burial_too():
    def sea(burial):
        return wc.seafloor_pco2_bar(0.01, Q, lambda p_bar: sp.TS_MODERN_K, EARTH_LAW, 1.5, burial)
    assert sea(0.0) == pytest.approx(0.01)                                       # the seafloor alone cannot keep up with the volcanoes of an Earth, whatever the temperature: it keeps all of them
    assert sea(1.5) < sea(0.0)                                                 # a biosphere that buries more than the volcanoes give takes the rest


def test_the_ledger_gives_classes_2_to_4_a_burial_row_and_no_other_class_has_one():
    assert {number for number, row in life_air.CLASSES.items() if row.get("burial")} == {"2", "3", "4"}
    assert life_air.number("anoxic_burial_share") == pytest.approx(0.02)


def _row(number):
    return life_air.CLASSES[number]


def test_an_earth_sized_world_with_earths_donors_fixes_the_20_tmol_of_carbon_the_photoferrotrophs_do_and_buries_a_fiftieth():
    body = SimpleNamespace(radius=6371.0)
    supplies = {"donor_supply": 1.0, "h2_supply": 1.0}
    area_ratio = biogenic_oxygen.area_m2(body) / biogenic_oxygen.EARTH_AREA_M2
    assert carbon_burial.fixed_carbon_mol_yr(body, _row("2"), supplies) == pytest.approx(20.0e12 * area_ratio)         # 80 Tmol of Fe2+ over 4 electrons a carbon [R3393]
    assert carbon_burial.fixed_carbon_mol_yr(body, _row("4"), supplies) == pytest.approx(20.0e12 * area_ratio * life_supply.number("sealed_fuel_efficiency"))
    hydrogen = carbon_burial.fixed_carbon_mol_yr(body, _row("3"), supplies)
    assert hydrogen == pytest.approx(life_air.number("earth_h2_tmol_yr") * 1e12 * 2.0 / 4.0 * area_ratio * life_supply.number("sealed_fuel_efficiency"))
    assert 0.0 < hydrogen < 1.0e12


def test_a_class_with_no_row_or_a_supply_of_nothing_fixes_no_carbon():
    body = SimpleNamespace(radius=6371.0)
    assert carbon_burial.fixed_carbon_mol_yr(body, _row("5"), {"donor_supply": 1.0}) == 0.0
    assert carbon_burial.fixed_carbon_mol_yr(body, _row("2"), {"donor_supply": 0.0}) == 0.0


def test_classes_on_one_fuel_are_counted_once_and_a_class_that_has_not_begun_buries_nothing(monkeypatch):
    body = SimpleNamespace(radius=6371.0, star=SimpleNamespace(age=4.0), life_age_photosynthesis_gyr=None)
    monkeypatch.setattr(life_supply, "supplies", lambda values: {"donor_supply": 1.0, "h2_supply": 0.0})
    monkeypatch.setattr(carbon_burial.habitability_world, "values_from_body", lambda body: {})

    def living(number, started):
        return life_air.Living(SimpleNamespace(number=int(number)), 0.5, started)
    area_ratio = biogenic_oxygen.area_m2(body) / biogenic_oxygen.EARTH_AREA_M2
    both = carbon_burial.anoxic_net_mol_yr(body, [living("2", 1.0), living("4", 1.0)])
    assert both == pytest.approx(0.02 * 20.0e12 * area_ratio)                  # the larger of the two (class 2's), not their sum
    assert carbon_burial.anoxic_net_mol_yr(body, [living("2", 5.0)]) == 0.0    # its producers begin after the star's age


def test_a_giant_a_world_with_no_life_and_a_body_with_no_radius_take_no_burial():
    assert carbon_burial.rate_of(SimpleNamespace(body_type="g", radius=60000.0), [object()]) == 0.0
    assert carbon_burial.rate_of(SimpleNamespace(radius=6371.0), []) == 0.0
    assert carbon_burial.rate_of(SimpleNamespace(radius=None), [object()]) == 0.0
