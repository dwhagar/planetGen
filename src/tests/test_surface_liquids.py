# tests/test_surface_liquids.py

"""The liquids a surface can hold (GEN.1035): the table, the windows against the draft's computed table, what a Titan-like, a warm and an outside-every-window
body have, a sweep that shows each liquid in its window, and the refusal of a bad row."""

import copy
import math

import pytest

from planetgen.physics import saturation, surface_liquids as sl


def test_the_table_loads_with_every_row_sourced_and_graded():
    assert len(sl.LIQUIDS) >= 50
    assert len({row.id for row in sl.LIQUIDS}) == len(sl.LIQUIDS)
    for row in sl.LIQUIDS:
        assert row.kind in sl.KINDS and row.tier in sl.TIERS and row.tweak in sl.TWEAKS and row.source
        assert row.density_kg_m3 > 0 and row.viscosity_mpa_s and (row.saturation or row.boiling_k > row.melting_k)
    assert {row.kind for row in sl.LIQUIDS} == set(sl.KINDS)
    for formula in ("H2", "Ne", "N2", "CO", "Ar", "O2", "CH4", "C2H6", "C2H4", "Kr", "Xe", "H2S", "CO2", "NH3", "SO2", "HCN", "CH3OH", "C6H6", "H2O"):
        assert any(row.saturation == formula for row in sl.LIQUIDS), formula      # every species of saturation.py that has a liquid curve


@pytest.mark.parametrize("liquid_id,pressure_pa,low,high", [
    ("nitrogen", 1e5, 63.15, 77.2), ("methane", 1e5, 90.7, 111.5), ("ethane", 1e5, 90.4, 184.3), ("ammonia", 1e5, 195.5, 239.6),
    ("water", 1e5, 273.16, 372.8), ("carbon_dioxide", 10e5, 216.6, 233.0), ("methane", 10e5, 90.7, 149.1), ("hydrogen", 1e5, 13.96, 20.3),
])
def test_the_window_of_a_species_with_a_saturation_curve_is_the_drafts_computed_one(liquid_id, pressure_pa, low, high):
    found_low, found_high = sl.window_k(liquid_id, pressure_pa)
    assert found_low == pytest.approx(low, abs=0.3) and found_high == pytest.approx(high, abs=1.0)


def test_carbon_dioxide_is_no_liquid_under_its_triple_pressure():
    assert sl.window_k("carbon_dioxide", 1e5) is None and sl.window_k("carbon_dioxide", 4e5) is None
    assert sl.window_k("carbon_dioxide", 6e5) is not None


def test_a_row_without_a_curve_follows_clausius_clapeyron_from_its_boiling_point():
    assert sl.window_k("propane", 1e5) == pytest.approx((85.5, 231.1), abs=0.5)
    assert sl.boiling_k("propane", 10e5) > sl.boiling_k("propane", 1e5) > sl.boiling_k("propane", 0.1e5)
    assert sl.window_k("sulfur", 1e5) == pytest.approx((388.4, 717.8), abs=0.5)
    assert sl.window_k("mercury", 1e5) == pytest.approx((234.3, 629.9), abs=0.5)
    assert sl.boiling_k("mercury", 1e12) == 1750.0           # held at the critical temperature


def test_every_liquid_is_stable_inside_its_own_window():
    for row in sl.LIQUIDS:
        pressure = 10e5 if row.saturation == "CO2" else 1e5
        low, high = sl.window_k(row.id, pressure)
        middle = 0.5 * (low + high)
        assert sl.is_stable(row.id, middle, pressure) and row.id in sl.stable_liquids(middle, pressure), row.id
        assert not sl.is_stable(row.id, low - 1.0, pressure) and not sl.is_stable(row.id, high + 1.0, pressure), row.id


def test_a_body_outside_every_window_holds_no_liquid():
    assert sl.stable_liquids(4000.0, 1e5) == ()
    assert sl.stable_liquids(2.0, 1e5) == ("helium4",)
    assert sl.dominant(4000.0, 1e5, {"CH4": 1e3}, {"water": 1.0}) is None


def test_a_titan_like_world_draws_a_hydrocarbon_sea():
    stable = sl.stable_liquids(94.0, 1.5e5)
    assert {"methane", "ethane", "propane"} <= set(stable) and "nitrogen" not in stable and "water" not in stable
    pick = sl.dominant(94.0, 1.5e5, {"CH4": 7.5e3, "N2": 1.4e5}, {"methane": 0.05, "water": 0.4})
    assert pick.liquid_id == "methane" and pick.basis == "reservoir"
    sea = sl.ideal_sea(94.0, {"CH4": 7.5e3, "N2": 1.4e5})
    assert sum(sea.values()) == pytest.approx(1.0)
    assert set(sea) == {"CH4", "N2", "C2H6"} and sea["CH4"] > sea["N2"] > 0.05 and 0.01 < sea["C2H6"] < 0.9
    density = sl.sea_density_kg_m3(sea, 94.0)
    assert 440.0 < density < 640.0                          # methane 0.45 g/cm3, ethane 0.65 and nitrogen 0.85 at those temperatures


def test_a_warm_world_still_draws_water():
    assert "water" in sl.stable_liquids(288.0, 1e5, kinds=("aqueous",)) and "nitrogen" not in sl.stable_liquids(288.0, 1e5)
    pick = sl.dominant(288.0, 1e5, {"H2O": 1.7e3}, {"water": 0.3})
    assert pick.liquid_id == "water"
    assert sl.dominant(288.0, 1e5, {"CH4": 1e3, "N2": 7e4}, {"methane": 0.05}) is None        # methane is gas at 288 K


def test_a_saturated_volatile_condenses_and_ranks_above_a_reservoir():
    cold = sl.candidates(95.0, 1.5e5, {"CH4": 2.5e4, "N2": 1.2e5}, {"water": 1.0, "ethane": 0.2})
    assert cold[0].basis == "condensing" and cold[0].liquid_id == "methane" and cold[0].score >= 1.0
    assert [c.liquid_id for c in cold if c.basis == "reservoir"] == ["ethane"]      # water is ice at 95 K, not stable
    assert sl.condensing(300.0, {"CH4": 1e4}) == {}                                   # over its critical temperature


def test_density_and_viscosity_follow_the_row():
    assert sl.density_kg_m3("water", 293.0) == pytest.approx(998.2)
    assert sl.density_kg_m3("methane", 111.7) == pytest.approx(422.6)
    assert sl.density_kg_m3("methane", 91.0) > sl.density_kg_m3("methane", 111.7)
    assert sl.density_kg_m3("methane", 1000.0) == pytest.approx(0.3 * 422.6)
    assert sl.viscosity_pa_s("water", 293.0) == pytest.approx(1.0e-3)
    assert sl.viscosity_pa_s("water", 373.0) < sl.viscosity_pa_s("water", 330.0) < sl.viscosity_pa_s("water", 293.0)
    assert sl.viscosity_pa_s("sulfur", 460.0) == pytest.approx(93.0)         # the lambda peak
    assert sl.viscosity_pa_s("sulfur", 430.0) < sl.viscosity_pa_s("sulfur", 392.0)
    assert sl.viscosity_pa_s("basalt_melt", 1473.0) > 1e3 * sl.viscosity_pa_s("komatiite_melt", 1900.0) * 0.01
    assert sl.viscosity_pa_s("rhyolite_melt", 1250.0) > 1e4 * sl.viscosity_pa_s("basalt_melt", 1473.0)


def test_a_refused_row_names_what_is_wrong():
    entry = {key: copy.deepcopy(value) for key, value in sl.liquid("water").__dict__.items()}
    good = {"id": "x", "name": "x", "kind": "metal", "melting_k": 300.0, "boiling_k": 600.0, "vaporization_kj_mol": 50.0,
            "density_kg_m3": {"at_k": 300.0, "value": 1000.0}, "expansion_per_k": 0.0001, "viscosity_mpa_s": [[300.0, 1.0]],
            "surface_tension_mn_m": 100.0, "tier": "T2", "tweak": "L0", "source": "test"}
    assert sl._row(good).id == "x"
    for change in ({"tier": "T10"}, {"tweak": "L9"}, {"kind": "plasma"}, {"source": ""}, {"boiling_k": 200.0}, {"saturation": "Unobtainium"}):
        with pytest.raises(sl.SurfaceLiquidError):
            sl._row({**good, **change})
    bad = dict(good)
    del bad["tier"]
    with pytest.raises(sl.SurfaceLiquidError, match="tier"):
        sl._row(bad)
    with pytest.raises(sl.SurfaceLiquidError):
        sl.liquid("no such liquid")
    del entry
