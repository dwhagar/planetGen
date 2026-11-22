# tests/test_ocean_oxygen.py

"""GEN.1824: the oxygen of a sealed ocean, from the lid's oxidant, the water's radiolysis and the faults' radicals over the reductants of the floor, and the animals of every class put to the oxygen floors with it."""

import pytest

from planetgen import tuning
from planetgen.physics import biogenic_oxygen, life_classes, life_complexity, life_supply, ocean_oxygen as oo

ICE = {"hydrosphere": "ice-covered ocean", "ocean_depth_km": 100.0, "surface_dose_msv_yr": 0.0, "internal_heat_flux_w_m2": 0.0, "tidal_flux_w_m2": 0.0}
"""A sealed ocean with no heat and no dose: only its own water's radiolysis makes oxygen."""


def sealed(**fields):
    return {**ICE, **fields}


def test_only_an_ocean_under_a_lid_has_an_oxygen_of_its_own():
    for state in ("surface ocean", "brine seep", "ice", "dry", "vapor", "hycean", None):
        assert oo.ocean_o2_kpa({**ICE, "hydrosphere": state}) is None
        assert oo.balance({**ICE, "hydrosphere": state}) is None
    assert oo.ocean_o2_kpa(ICE) is not None


def test_the_lid_delivers_the_share_of_europas_surface_oxidant_that_the_dose_gives_to_the_published_top():
    at_reference = oo.lid_delivery(sealed(surface_dose_msv_yr=oo.REFERENCE_DOSE_MSV_YR))
    assert at_reference == pytest.approx(0.85 * 1.2e10 / 3.1e13, rel=0.05)           # 12 kg a second of oxygen over Europa's area, the share that gets through
    assert oo.lid_delivery(sealed(surface_dose_msv_yr=oo.REFERENCE_DOSE_MSV_YR / 4.0)) == pytest.approx(at_reference / 4.0)
    assert oo.lid_delivery(sealed()) == 0.0 and oo.lid_delivery(sealed(surface_dose_msv_yr=None)) == 0.0
    top = oo.lid_delivery(sealed(surface_dose_msv_yr=1.0e12))
    assert top == oo.LID_DELIVERY_MAX_MOL_M2_YR == pytest.approx(3.0e11 / 3.1e13, rel=0.05)      # the Greenberg end of the range


def test_the_waters_radiolysis_goes_with_its_depth_to_a_hundred_kilometres():
    assert oo.radiolysis(sealed(ocean_depth_km=3.7)) == pytest.approx(3.2e-6, rel=0.1)        # the report's Earth column of 3.7 km
    assert oo.radiolysis(sealed(ocean_depth_km=100.0)) == pytest.approx(8.5e-5, rel=0.1)
    assert oo.radiolysis(sealed(ocean_depth_km=400.0)) == oo.radiolysis(sealed(ocean_depth_km=100.0))
    assert oo.radiolysis(sealed(ocean_depth_km=None)) == 0.0


def test_the_faults_go_with_the_heat_to_earths_and_no_further():
    earth = life_supply.surface_pressure.HEAT_FLUX_EARTH_W_M2
    assert oo.fracture(sealed(internal_heat_flux_w_m2=earth)) == pytest.approx(oo.FRACTURE_MAX_MOL_M2_YR)
    assert oo.fracture(sealed(internal_heat_flux_w_m2=earth / 10.0)) == pytest.approx(oo.FRACTURE_MAX_MOL_M2_YR / 10.0)
    assert oo.fracture(sealed(internal_heat_flux_w_m2=earth * 8.0)) == pytest.approx(oo.FRACTURE_MAX_MOL_M2_YR)
    assert oo.fracture(sealed()) == 0.0


def test_the_sink_is_the_mean_of_the_two_supplies_in_oxygen_terms_and_the_water_keeps_its_hydrogen():
    shares = {"h2_supply": 0.5, "donor_supply": 2.0}
    assert oo.sink(ICE, supplies={"h2_supply": 0.0, "donor_supply": 0.0}) == pytest.approx(oo.radiolysis(ICE))      # only the water's own hydrogen
    expected = oo.H2_O2_EQUIVALENT_MOL_M2_YR * 1.25 + oo.radiolysis(ICE)
    assert oo.sink(ICE, supplies=shares) == pytest.approx(expected)
    assert oo.H2_O2_EQUIVALENT_MOL_M2_YR == pytest.approx(1.4e-3, rel=0.02)
    assert oo.sink({**ICE, "h2_supply": 0.5, "donor_supply": 2.0}) == pytest.approx(expected)          # a stored share replaces the one worked out


def test_the_waters_radiolysis_is_net_zero_until_its_hydrogen_leaves(monkeypatch):
    assert oo.balance(ICE).source == pytest.approx(oo.radiolysis(ICE)) and oo.balance(ICE).sink >= oo.radiolysis(ICE)
    assert oo.ocean_o2_kpa(ICE) == 0.0                                      # the water's oxygen is taken back by its own hydrogen
    monkeypatch.setattr(oo, "RADIOLYSIS_H2_ESCAPE", 1.0)
    assert oo.ocean_o2_kpa(ICE) > tuning.LIFE_COMPLEXITY["oxygen_floor_kpa"]["A1"]      # counted gross, a quiet 100 km sea would clear the first floor


def test_a_world_with_heat_holds_no_oxygen_for_the_reductants_take_it_all():
    warm = sealed(internal_heat_flux_w_m2=0.06, surface_dose_msv_yr=oo.REFERENCE_DOSE_MSV_YR)
    held = oo.balance(warm)
    assert held.sink > held.source and held.kpa == 0.0 and held.um == 0.0


def test_the_bank_is_the_surplus_over_the_pulse_interval_in_the_water_the_ocean_holds():
    # the research report's arithmetic: Europa's 1.2e10 mol a year over 3.06e13 m2, an ocean of 100 km, the sink of 1e9 mol of H2 a year, 50 Myr
    source, reductants = 1.18e10 / 3.06e13, 0.5 * 1.0e9 / 3.06e13
    values = sealed(ocean_depth_km=100.0)
    assert oo.held_um(source, reductants, values) == pytest.approx(185.0, rel=0.03)
    assert oo.held_um(source, 0.5 * 1.6e10 / 3.06e13, values) == pytest.approx(62.0, rel=0.05)
    assert oo.held_um(source, source * 2.0, values) == 0.0
    assert oo.held_um(source, reductants, sealed(ocean_depth_km=0.0)) == 0.0
    assert oo.held_um(source, reductants, sealed(ocean_depth_km=200.0)) == pytest.approx(oo.held_um(source, reductants, values) / 2.0)   # the same oxygen in twice the water
    assert oo.held_kpa(source, reductants, values) == pytest.approx(185.0 / oo.HENRY_UM_PER_KPA, rel=0.03)


def test_the_bank_is_held_at_the_pressure_the_air_law_allows_and_no_higher():
    from planetgen.physics import mixture_oxygen
    rich = oo.held_kpa(1.0, 0.0, sealed(ocean_depth_km=1.0))
    assert rich == pytest.approx(mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA / 1000.0)


def test_the_comet_oxygen_is_a_stock_the_reductants_drain_and_it_is_off_by_default(monkeypatch):
    assert oo.INHERITED_O2_RETAINED == 0.0 and oo.held_um(0.0, 0.0, sealed()) == 0.0
    monkeypatch.setattr(oo, "INHERITED_O2_RETAINED", 0.01)
    rich = oo.held_um(0.0, 0.0, sealed(ocean_depth_km=100.0))
    assert rich == pytest.approx(0.01 * 0.038 * 55.5 * 1.0e6)                        # about 21 millimolar, untouched with no reductants
    drained = oo.held_um(0.0, 1.0e-3, sealed(ocean_depth_km=100.0))                  # 4.5 Gyr of a strong sink in 1e8 litres: 4.5e4 micromolar, more than it had
    assert drained == 0.0
    mild = oo.held_um(0.0, 1.0e-6, sealed(ocean_depth_km=100.0))
    assert 0.0 < mild < rich
    assert oo.held_um(1.0e-3, 0.0, sealed(ocean_depth_km=100.0)) > oo.held_um(0.0, 0.0, sealed(ocean_depth_km=100.0))     # the bank adds to it


def test_a_quiet_world_with_a_heavy_belt_dose_reaches_the_floors_and_a_busy_one_does_not():
    floor = tuning.LIFE_COMPLEXITY["oxygen_floor_kpa"]
    heavy = sealed(surface_dose_msv_yr=1.0e7, ocean_depth_km=60.0, internal_heat_flux_w_m2=0.002)
    assert oo.ocean_o2_kpa(heavy) >= floor["A1"]
    assert oo.ocean_o2_kpa({**heavy, "surface_dose_msv_yr": 100.0}) < floor["A1"]            # the cosmic rays alone make too little against the sink
    assert oo.ocean_o2_kpa({**heavy, "internal_heat_flux_w_m2": 0.06, "tidal_flux_w_m2": 0.5}) < floor["A1"]
    assert oo.ocean_o2_kpa(sealed(ocean_depth_km=50.0)) < floor["A1"]


def _row(class_id):
    return life_classes.BY_ID[class_id]


INDEX_A1 = tuning.LIFE_COMPLEXITY["thresholds"]["A1"]


def test_a_nitrate_class_in_a_sealed_ocean_needs_the_oceans_oxygen_for_its_animals():
    row = _row("4b")
    with_air_only = life_complexity.assess(row, 1.0e-2, 0.0, 273.15)
    assert row.oxidant != "o2" and with_air_only.grade >= 1                    # the rule of every other place: a class that breathes nitrate needs no oxygen floor
    dry = life_complexity.assess(row, 1.0e-2, 0.0, 273.15, ocean_o2_kpa=0.0)
    assert dry.grade == 0 and dry.index == pytest.approx(with_air_only.index)     # the index is the class's own; only the grade is read against the floor
    first = life_complexity.assess(row, 1.0e-2, 0.0, 273.15, ocean_o2_kpa=tuning.LIFE_COMPLEXITY["oxygen_floor_kpa"]["A1"])
    assert first.grade == 1
    second = life_complexity.assess(row, 1.0e-2, 0.0, 273.15, ocean_o2_kpa=tuning.LIFE_COMPLEXITY["oxygen_floor_kpa"]["A2"])
    assert second.grade == 2


def test_the_animals_that_lost_aerobic_respiration_are_the_way_round_the_floor():
    row = _row("3b")
    index_npp = 1.0e-3                                                      # an index past the multicellular threshold, under the A1 one at the anaerobic yield
    plain = life_complexity.assess(row, index_npp, 0.0, 273.15, ocean_o2_kpa=0.0)
    anoxic = life_complexity.assess(row, index_npp, 0.0, 273.15, anoxic_animals=True, ocean_o2_kpa=0.0)
    assert plain.grade == 0 and anoxic.grade == 1


def test_an_oxygen_breather_reads_the_larger_of_its_air_and_the_ocean():
    row = _row("4a")
    air = life_complexity.assess(row, 1.0e-2, 5.0, 273.15)
    assert life_complexity.assess(row, 1.0e-2, 5.0, 273.15, ocean_o2_kpa=0.0).index == pytest.approx(air.index)
    assert life_complexity.assess(row, 1.0e-2, 0.0, 273.15, ocean_o2_kpa=5.0).index == pytest.approx(air.index)
    assert life_complexity.assess(row, 1.0e-2, 5.0, 273.15, ocean_o2_kpa=1.0).oxidant_yield == pytest.approx(air.oxidant_yield)


def test_a_body_with_no_lid_is_assessed_as_before():
    for class_id in ("3a", "3b", "4a", "4b", "4c"):
        row = _row(class_id)
        for o2 in (0.0, 0.5, 5.0):
            assert life_complexity.assess(row, 3.0e-3, o2, 280.0) == life_complexity.assess(row, 3.0e-3, o2, 280.0, ocean_o2_kpa=None)
