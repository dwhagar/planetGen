# tests/test_atmosphere_mixture_lava.py

"""The vapour over a molten surface (GEN.499, GEN.456, GEN.479) against [R2380]'s vapour pressure of a BSE melt and the
[R3059] sequence of dominant gases."""

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
    subject = Subject("terrestrial", HOT)
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
    subject = Subject("terrestrial", {**HOT, "pressure_bar": 10.0})
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
    subject = Subject("terrestrial", {**HOT, "host_type": "WD", "ts_k": 2000.0, "tp_k": 2000.0})
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
            ml.vapour_mixture(Subject("terrestrial", {**HOT, "tp_k": 900.0}), md.RECIPES.draw_all(Subject("terrestrial", HOT), categories=md.CATEGORIES))


def test_a_stored_melt_composition_replaces_the_drawn_one():
    constants = {}
    fractions = ml.melt_composition(constants, {"SiO2": 50.0, "MgO": 30.0, "FeO": 20.0})
    assert sum(fractions.values()) == pytest.approx(1.0)
    assert fractions["Na2O"] < 1e-9 and fractions["FeO"] > fractions["Al2O3"]
    # no draw was made for it: the same call twice agrees without a seed
    assert ml.melt_composition({}, {"SiO2": 50.0, "MgO": 30.0, "FeO": 20.0}) == fractions


@pytest.mark.parametrize("temperature", [2000.0, 2500.0])
def test_the_alkali_oxides_and_monoxide_of_the_vapour_are_in_the_gas_table(temperature):
    """GEN.751: K2, KO, K2O and Na2O were species the gas table lacked, so a draw dropped them as `unlisted`; they are listed now."""
    pressures = vapour(temperature, 0.0)
    for name in ("K2", "KO", "K2O", "Na2O"):
        assert pressures[name] > 0.0, name
    subject = Subject("terrestrial", {**HOT, "ts_k": temperature, "tp_k": temperature})
    with draw.bound(draw.Stream("gen751")):
        recipe = md.RECIPES.draw_all(subject, categories=md.CATEGORIES)
        result = ml.vapour_mixture(subject, recipe, temperature_name="tp_k")
    assert result.unlisted == {}
    for formula in ("K2", "KO", "K2O", "Na2O"):
        assert gases.has(formula)
    assert result.mole_fractions.get("KO", 0.0) > 0.0 and result.mole_fractions.get("Na2O", 0.0) > 0.0


def test_aluminium_is_in_the_vapour_over_the_melt_in_the_proportion_of_its_activity():
    """GEN.818: Al2O3 has an activity (its mole fraction), so the vapour lists Al, AlO, Al2O, AlO2 and Al2O2."""
    assert "Al" in ml.ELEMENTS and ml.ALUMINA_ACTIVITY_COEFFICIENT == 1.0
    pressures = vapour(2500.0, 0.0)
    aluminium = {name: pressures[name] for name in ("AL", "ALO", "AL2O", "ALO2", "AL2O2")}
    assert all(value > 0.0 for value in aluminium.values())
    share = sum(aluminium.values()) / sum(pressures.values())
    assert 1e-9 < share < 1e-4                                  # a refractory element: far below Mg and SiO
    assert aluminium["ALO"] > aluminium["AL"] > aluminium["AL2O"]
    assert sum(aluminium.values()) < pressures["Mg"] * 1e-3


def test_the_aluminium_vapour_obeys_the_equilibrium_of_its_own_species():
    # 1/2 Al2O3(L) = AlO + 1/4 O2 at the held oxygen fugacity: p(AlO) fO2^(1/4) = a(Al2O3)^(1/2) K.
    temperature, redox = 2500.0, 0.0
    with draw.bound(4):
        fractions = ml.melt_composition({})
    pressures = vapour(temperature, redox, seed=4)
    log_k = species_thermo.log10_k({"ALO": 1.0, "O2": 0.25, "AL2O3(L)": -0.5}, temperature)
    from planetgen.physics import oxygen_fugacity as fo2
    left = math.log10(pressures["ALO"]) + 0.25 * fo2.log_fo2(redox, temperature)
    assert left == pytest.approx(log_k + 0.5 * math.log10(fractions["Al2O3"]), abs=1e-6)


def test_a_more_aluminous_melt_has_more_aluminium_in_its_vapour():
    """The Al2O3 share of a stored melt composition moves the aluminium vapour with the square root of the activity."""
    low = ml.vapour_pressures(ml.element_potentials({"SiO2": 0.5, "MgO": 0.4, "FeO": 0.08, "Al2O3": 0.01, "CaO": 0.01, "Na2O": 1e-9, "TiO2": 0.001, "K2O": 1e-9}, 2500.0, 0.0), 2500.0)
    high = ml.vapour_pressures(ml.element_potentials({"SiO2": 0.5, "MgO": 0.4, "FeO": 0.08, "Al2O3": 0.04, "CaO": 0.01, "Na2O": 1e-9, "TiO2": 0.001, "K2O": 1e-9}, 2500.0, 0.0), 2500.0)
    assert high["ALO"] / low["ALO"] == pytest.approx(2.0, rel=1e-6)
