# tests/test_atmosphere_mixture_evolved.py

"""The evolved, cold-body and comet mixtures (GEN.479): Venus-like, Mars-like, abiotic oxygen, Titan-like, Pluto-like, coma."""

import pytest

from planetgen.generation.possibility_table import Subject
from planetgen.physics import gases, mixture_draw as md
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.util import draw

BASE = dict(mass_earth=1.0, ts_k=288.0, pressure_bar=1.0, age_gy=4.5, long_lived_geology=True, ocean_present=True,
            snow_ratio=0.5, zone="e", host_type="G", is_moon=False, insolation=1.0, luminosity_solar=1.0, water_lost=False,
            runaway_ratio=0.5, teq_k=255.0, xuv_index=1.0, v_esc_kms=11.2, radius_km=6371.0, host_mass_solar=1.0)

MARS = dict(mass_earth=0.1, ts_k=215.0, pressure_bar=0.006, v_esc_kms=5.0, radius_km=3400.0, ocean_present=False, long_lived_geology=False)
TITAN = dict(mass_earth=0.0225, ts_k=94.0, pressure_bar=1.5, is_moon=True, radius_km=2575.0, insolation=0.001, zone="c", v_esc_kms=2.6)
# the shoreline gives the retention rows d4 and d5 that j3 and j4 are gated on only to a body with an escape speed this high
PLUTO = dict(mass_earth=0.002, ts_k=40.0, pressure_bar=1e-5, radius_km=1188.0, zone="c", v_esc_kms=3.0, insolation=0.0006)
ABIOTIC = dict(host_type="M", age_gy=5.0, ts_k=290.0)


def draws(body_class, changes, family, count=60):
    found = []
    for seed in range(count):
        with draw.bound(seed):
            result = md.draw_mixture(Subject(body_class, {**BASE, **changes}))
        if result.family == family:
            found.append(result)
    return found


def test_a_mars_like_body_has_a_thin_carbon_dioxide_nitrogen_air_without_water():
    found = draws("terrestrial", MARS, "mars")
    assert found
    for result in found:
        x = result.mole_fractions
        assert x["CO2"] > 0.9 and 0.0 < x["N2"] < 0.1 and x.get("H2O", 0.0) < 1e-3
        assert x.get("CO", 0.0) < 0.1 * x["CO2"]


def test_a_venus_like_body_is_carbon_dioxide_with_a_few_percent_nitrogen():
    changes = dict(water_lost=True, runaway_ratio=1.2, ts_k=735.0, pressure_bar=90.0, mass_earth=0.8, ocean_present=False)
    found = draws("terrestrial", changes, "venus", count=120)
    assert found
    for result in found:
        assert result.mole_fractions["CO2"] > 0.6 and result.mole_fractions.get("H2O", 0.0) < 0.01


def test_an_abiotic_oxygen_world_keeps_oxygen_in_its_air():
    found = draws("terrestrial", ABIOTIC, "abiotic-oxygen")
    assert found
    assert max(result.mole_fractions.get("O2", 0.0) for result in found) > 0.1


def test_a_titan_like_air_is_nitrogen_with_methane_at_the_saturation_and_a_hydrogen_trace():
    found = draws("terrestrial", TITAN, "titan")
    assert found
    for result in found:
        x = result.mole_fractions
        assert x["N2"] > 0.85 and 0.02 < x["CH4"] < 0.12 and 1e-4 < x["H2"] < 5e-3
        assert x["CH4"] <= result.constants["methane_saturation_bar"] / TITAN["pressure_bar"] + 1e-12


def test_a_pluto_like_air_follows_raoults_law_over_the_ice():
    found = draws("terrestrial", PLUTO, "frost")
    assert found
    for result in found:
        x = result.mole_fractions
        assert x["N2"] > 0.9
        assert result.pressure_bar == pytest.approx(sum(
            result.inventory[formula] * __import__("planetgen.physics.saturation", fromlist=["x"]).saturation(formula, 40.0).pressure_pa / 1e5
            for formula in ("N2", "CH4", "CO")), rel=1e-6)


def test_a_coma_is_water_with_the_drawn_ratios_and_the_class_ordering_of_carbon_dioxide_and_monoxide():
    subject = Subject("comet", {name: 1.0 for name in md.required_inputs("comet")})
    seen = set()
    for seed in range(40):
        with draw.bound(seed):
            result = md.draw_mixture(subject)
        x = result.mole_fractions
        assert result.family == "coma" and sum(x.values()) == pytest.approx(1.0, abs=1e-12)
        assert all(gases.has(formula) for formula in x)
        row = result.drawn["k"]
        seen.add(row)
        if row == "k1":
            assert x["CO2"] >= x["CO"]
        elif row == "k2":
            assert x["CO"] >= x["CO2"]
        else:
            assert x["CO"] == pytest.approx(x["CO2"])
        assert x["H2O"] > 0.25
    assert seen == {"k1", "k2", "k3"}


def test_a_white_dwarf_planet_without_a_vapour_exosphere_has_no_gas():
    subject = Subject("terrestrial", {**BASE, "host_type": "WD", "ts_k": 400.0, "tp_k": 400.0, "zone": "h"})
    kinds = set()
    for seed in range(30):
        with draw.bound(seed):
            result = md.draw_mixture(subject)
        kinds.add(result.family)
        if result.family == "airless":
            assert result.mole_fractions == {}
    assert "airless" in kinds


def test_every_evolved_draw_is_reproducible():
    for body_class, changes in (("terrestrial", MARS), ("terrestrial", TITAN), ("terrestrial", PLUTO), ("terrestrial", ABIOTIC)):
        for seed in (3, 4):
            with draw.bound(seed):
                a = md.draw_mixture(Subject(body_class, {**BASE, **changes}))
            with draw.bound(seed):
                b = md.draw_mixture(Subject(body_class, {**BASE, **changes}))
            assert a.mole_fractions == b.mole_fractions


def test_the_trace_gases_follow_the_rows_drawn_after_the_life_stage():
    from planetgen.physics import mixture_trace as mt
    life = dict(o2_pal=1.0, life_chemical=True, sep_reducing=False, lightning=True, co2_rich=False, h2_poor=True,
                has_solid_ground=True, h2o_present=True, carbon_bearing=True, gas_has_co_h2s=False, t_layer_k=280.0)
    needed = mt.required_trace_inputs("terrestrial")
    assert {"o2_pal", "life_chemical", "lightning", "co2_rich", "h2_poor", "has_solid_ground"} <= set(needed)
    seen = set()
    for seed in range(40):
        subject = Subject("terrestrial", {**BASE, **life})
        with draw.bound(seed):
            first = md.draw_mixture(subject)
            result = mt.draw_trace_gases(subject, first)
        assert sum(result.mole_fractions.values()) == pytest.approx(1.0, abs=1e-12)
        assert all(gases.has(formula) for formula in result.mole_fractions)
        for gas in ("N2O", "PH3", "HCl", "HF", "NO2", "CH3Cl"):
            if gas in result.mole_fractions and gas not in first.mole_fractions:
                seen.add(gas)
        assert set(first.mole_fractions) <= set(result.mole_fractions)
        assert "m1" in result.recipe.binaries                     # the ozone rows are in the result for the ozone code
    assert {"N2O", "PH3"} & seen


TITAN_MASS_KG = 0.0225 * 5.972e24
ICY = dict(water_mass_fraction=0.3, rock_mass_fraction=0.55)


def _titan_n2_kg(result, rock_fraction=ICY["rock_mass_fraction"]):
    """The N2 of a titan-family draw in kg for a body of Titan's mass: its moles per 100 g of mantle times the mantle."""
    return result.moles * rock_fraction * TITAN_MASS_KG / 0.1 * 0.028


def test_a_titan_like_air_takes_the_nitrogen_of_the_ammonia_ice_and_comes_out_of_the_order_of_titans():
    """GEN.821 ([R1282], section 4.2): the nitrogen is the share epsilon of the ammonia that is outgassed times the share xi that
    ends as N2; with the ammonia 0.2 to 1.5 percent of the water of a Titan-sized ice moon the air holds 1e17 to 1e20 kg of N2 (Titan's 1.5 bar is
    about 9e18 kg)."""
    found = draws("terrestrial", {**TITAN, **ICY}, "titan", count=80)
    assert len(found) > 20
    masses = sorted(_titan_n2_kg(result) for result in found)
    assert masses[0] > 3e16 and masses[-1] < 3e20
    assert 1e17 < masses[len(masses) // 2] < 3e19


def test_the_ammonia_nitrogen_is_the_drawn_shares_of_the_ammonia_added_to_the_nitrogen_of_the_rock():
    seen = 0
    for seed in range(40):
        with draw.bound(seed):
            with_ice = md.draw_mixture(Subject("terrestrial", {**BASE, **TITAN, **ICY}))
        with draw.bound(seed):
            without = md.draw_mixture(Subject("terrestrial", {**BASE, **TITAN}))
        if with_ice.family != "titan" or without.family != "titan":
            continue
        seen += 1
        constants = with_ice.constants
        ammonia = constants["ammonia_per_water"] * ICY["water_mass_fraction"] * 17.03052 / 18.01528      # grams per gram of body
        nitrogen_atoms = 100.0 * constants["ammonia_outgassed"] * constants["ammonia_to_n2"] * ammonia / (17.03052 * ICY["rock_mass_fraction"])
        n2_with = with_ice.moles * with_ice.mole_fractions["N2"]
        n2_without = without.moles * without.mole_fractions["N2"]
        assert n2_with - n2_without == pytest.approx(0.5 * nitrogen_atoms, rel=1e-6)
    assert seen > 10


def test_a_stored_ammonia_ice_is_used_as_it_is_and_the_ratio_to_the_water_is_not_drawn():
    found = draws("terrestrial", {**TITAN, **ICY, "nh3_mass_fraction": 0.004}, "titan", count=30)
    assert found
    for result in found:
        assert "ammonia_per_water" not in result.constants and "ammonia_outgassed" in result.constants
