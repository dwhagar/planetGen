# tests/test_atmosphere_mixture_lava.py

"""The vapour over a molten surface (GEN.499, GEN.456, GEN.479) against Kite et al. 2016's vapour pressure of a BSE melt and the
Miguel et al. 2011 sequence of dominant gases."""

import math

import pytest

from planetgen.generation.possibility_table import Subject
from planetgen.physics import gases, mixture_draw as md, mixture_lava as ml, species_thermo
from planetgen.util import draw

HOT = dict(mass_earth=2.0, ts_k=2200.0, tp_k=2200.0, pressure_bar=1e-3, age_gy=4.5, long_lived_geology=False, ocean_present=False,
           snow_ratio=0.5, zone="h", host_type="G", is_moon=False, insolation=1000.0, luminosity_solar=1.0, water_lost=False,
           runaway_ratio=2.0, teq_k=1800.0, xuv_index=1.0, v_esc_kms=15.0, radius_km=8000.0, host_mass_solar=1.0)


def vapour(temperature, redox, seed=1, gamma=1e-8):
    with draw.bound(seed):
        constants = {}
        fractions = ml.melt_composition(constants)
    activities = {oxide: value * (gamma if oxide in ("Na2O", "K2O") else 1.0) for oxide, value in fractions.items()}
    return ml.vapour_pressures(ml.element_potentials(activities, temperature, redox), temperature)


def mean_molar_mass(pressures):
    return sum(species_thermo.SPECIES[name].molar_mass * 1e3 * value for name, value in pressures.items()) / sum(pressures.values())


def test_the_vapour_pressure_over_the_melt_is_kites_within_a_factor_of_four_from_2000_to_2500_k():
    for temperature in (2000.0, 2500.0):
        kite = 10.0 ** (-34.7 * math.exp(-0.00112 * temperature))
        assert sum(vapour(temperature, 0.5).values()) == pytest.approx(kite, rel=3.0)


def test_sodium_and_potassium_lead_the_cool_vapour_and_silicon_monoxide_the_hot_one():
    cool, hot = vapour(1800.0, 0.0), vapour(3200.0, 0.0)
    assert max(cool, key=cool.get) in ("Na", "K")
    assert max(hot, key=hot.get) == "SiO"
    assert hot["SiO"] > hot["O2"]                           # type IV of Miguel et al.: SiO above O2
    assert hot["Mg"] > hot["Fe"] > 0.0


def test_a_more_oxidised_melt_has_more_oxygen_in_its_vapour_and_less_iron_and_magnesium():
    reduced, oxidised = vapour(2800.0, -2.0), vapour(2800.0, 2.0)
    assert oxidised["O2"] / sum(oxidised.values()) > reduced["O2"] / sum(reduced.values())
    assert oxidised["Fe"] < reduced["Fe"] and oxidised["Mg"] < reduced["Mg"]


def test_the_vapour_obeys_the_equilibria_of_its_own_species():
    # SiO2(L) = SiO + 1/2 O2 at the held oxygen fugacity: p(SiO) fO2^(1/2) = a(SiO2) K.
    temperature, redox = 2500.0, 0.0
    with draw.bound(4):
        fractions = ml.melt_composition({})
    pressures = vapour(temperature, redox, seed=4)
    log_k = species_thermo.log10_k({"SiO": 1.0, "O2": 0.5, "SiO2(L)": -1.0}, temperature)
    from planetgen.physics import oxygen_fugacity as fo2
    left = math.log10(pressures["SiO"]) + 0.5 * fo2.log_fo2(redox, temperature)
    assert left == pytest.approx(log_k + math.log10(fractions["SiO2"]), abs=1e-6)


def test_the_mole_fractions_name_gas_table_species_and_sum_to_one_for_a_lava_world():
    subject = Subject("A", HOT)
    seen = set()
    for seed in range(30):
        with draw.bound(seed):
            result = md.draw_mixture(subject)
        seen.add(result.family)
        if result.family == "lava":
            assert sum(result.mole_fractions.values()) == pytest.approx(1.0, abs=1e-12)
            assert all(gases.has(formula) for formula in result.mole_fractions)
            assert result.constants["vapour_pressure_bar"] > 0.0
    assert "lava" in seen


def test_the_volatile_bearing_branch_adds_carbon_dioxide_and_monoxide_at_the_surface_pressure():
    subject = Subject("A", {**HOT, "pressure_bar": 10.0})
    bearing = []
    for seed in range(60):
        with draw.bound(seed):
            result = md.draw_mixture(subject)
        if result.family == "lava" and result.constants.get("volatile_bearing"):
            bearing.append(result)
    assert bearing
    for result in bearing:
        assert result.mole_fractions["CO2"] + result.mole_fractions["CO"] > 0.5


def test_a_white_dwarf_planet_with_a_vapour_exosphere_is_a_lava_mixture_over_the_debris_redox():
    subject = Subject("A", {**HOT, "host_type": "WD", "ts_k": 2000.0, "tp_k": 2000.0})
    families = {}
    for seed in range(40):
        with draw.bound(seed):
            result = md.draw_mixture(subject)
        families[result.family] = families.get(result.family, 0) + 1
        if result.drawn.get("l") == "l1":
            assert result.family == "lava" and result.redox_diw == pytest.approx(-1.0, abs=2.5)        # a4: N(-1, 0.5)
    assert "lava" in families


def test_a_surface_that_is_not_molten_is_refused():
    with draw.bound(1):
        with pytest.raises(ml.MixtureError):
            ml.vapour_mixture(Subject("A", {**HOT, "tp_k": 900.0}), md.RECIPES.draw_all(Subject("A", HOT), categories=md.CATEGORIES))


def test_a_stored_melt_composition_replaces_the_drawn_one():
    constants = {}
    fractions = ml.melt_composition(constants, {"SiO2": 50.0, "MgO": 30.0, "FeO": 20.0})
    assert sum(fractions.values()) == pytest.approx(1.0)
    assert fractions["Na2O"] < 1e-9 and fractions["FeO"] > fractions["Al2O3"]
    # no draw was made for it: the same call twice agrees without a seed
    assert ml.melt_composition({}, {"SiO2": 50.0, "MgO": 30.0, "FeO": 20.0}) == fractions
