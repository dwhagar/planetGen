# tests/test_atmosphere_mixture_giant.py

"""The air of giants, brown dwarfs and hot Jupiters (GEN.454): element budget, clouds, quench and the regimes of part 05."""

import dataclasses
import math

import pytest

from planetgen.generation.possibility_table import Subject
from planetgen.physics import gas_equilibrium as ge, mixture_giant as mg, saturation, species_thermo as st
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.physics.mixture_base import MixtureError
from planetgen.util import draw

JUPITER = dict(mass_earth=318.0, mass_jupiter=1.0, radius_km=69911.0, teq_k=120.0, t_layer_k=120.0, tday_k=120.0, teff_k=120.0,
               zone="c", host_type="G", age_gy=4.5, is_moon=False, insolation=0.04, luminosity_solar=1.0, host_mass_solar=1.0,
               core_radius_earth=3.0, snow_ratio=3.0, xuv_index=0.04, v_esc_kms=59.0, grain_um=1.0, photosphere_bar=0.1)


def giant(body_class="giant", **changes):
    return Subject(body_class, {**JUPITER, **changes})


def mixture(seed, body_class="giant", **changes):
    subject = giant(body_class, **changes)
    with draw.bound(seed):
        recipe = RECIPES.draw_all(subject, categories="abdcefghijkl")
        return mg.giant_mixture(subject, recipe)


def test_the_element_budget_is_the_reference_scaled_in_the_metals_with_the_helium_fixed():
    reference = mg.element_budget(0.0)
    enriched = mg.element_budget(1.0)
    assert reference["H"] == 1.0 and enriched["He"] == reference["He"]
    assert reference["He"] == pytest.approx(0.0821, rel=2e-3)
    assert reference["C"] / reference["O"] == pytest.approx(0.589, abs=0.005)
    for symbol in mg.METALS:
        assert enriched[symbol] == pytest.approx(10.0 * reference[symbol])
    split = mg.element_budget(0.0, c_over_o=1.2)
    assert split["C"] / split["O"] == pytest.approx(1.2)
    assert split["C"] + split["O"] == pytest.approx(reference["C"] + reference["O"])


def test_the_reference_metal_mass_fraction_is_the_photospheric_z_and_inverts():
    assert mg.metal_mass_fraction(0.0) == pytest.approx(0.0139, rel=0.1)
    for z in (0.01, 0.1, 0.3, 0.6):
        assert mg.metal_mass_fraction(mg.enrichment_of_mass_fraction(z)) == pytest.approx(z, rel=1e-6)


def test_the_cool_solar_gas_is_the_hydride_table_of_the_research():
    # Part 05 section 2.2: all O in H2O, C in CH4, N in NH3, S in H2S, H as H2 (the table's x at 300 K and 0.3 bar).
    budget = {symbol: mg.element_budget(0.0)[symbol] for symbol in ("H", "He", "C", "N", "O", "Ne", "S", "Ar")}
    fractions = ge.solve(budget, 300.0, 0.3).mole_fractions
    for name, expected in {"H2": 0.857, "He": 0.141, "H2O": 8.4e-4, "CH4": 4.96e-4, "NH3": 1.16e-4, "Ne": 1.97e-4, "H2S": 2.3e-5,
                           "Ar": 4.1e-6}.items():
        assert fractions[name] == pytest.approx(expected, rel=0.02), name


def test_the_condensation_temperatures_are_the_fits_of_visscher_lodders_and_fegley():
    cond = mg.condensation_temperatures(1.0, 0.0)
    assert cond["iron"] == pytest.approx(1e4 / 5.44)
    assert cond["forsterite"] == pytest.approx(1e4 / 5.89)
    assert cond["enstatite"] == pytest.approx(1e4 / 6.26)
    richer = mg.condensation_temperatures(1.0, 1.0)
    assert all(richer[name] > cond[name] for name in cond)          # more metal condenses higher
    assert mg.condensation_temperatures(10.0, 0.0)["iron"] > cond["iron"]


def test_clouds_take_iron_silicates_and_about_a_fifth_of_the_oxygen():
    budget = mg.element_budget(0.0)
    cold, clouds = mg.remove_clouds(budget, 1000.0, 0.3, 0.0)
    assert {"iron", "forsterite", "enstatite", "oxides"} <= set(clouds)
    for symbol in ("Fe", "Mg", "Si", "Ca", "Ti", "V", "Cr"):
        assert cold[symbol] == 0.0
    assert 0.15 < 1.0 - cold["O"] / budget["O"] < 0.3               # [R4561]: about 20 percent
    hot, none = mg.remove_clouds(budget, 2500.0, 0.3, 0.0)
    assert hot == budget and not none
    trapped, _ = mg.remove_clouds(budget, 2500.0, 0.3, 0.0, keep_titanium_vanadium=False)
    assert trapped["Ti"] == 0.0 and trapped["V"] == 0.0 and trapped["Fe"] == budget["Fe"]
    no_alkali, _ = mg.remove_clouds(budget, 1500.0, 0.3, 0.0, keep_sodium=False, keep_potassium=False)
    assert no_alkali["Na"] == 0.0 and no_alkali["K"] == 0.0 and no_alkali["S"] < budget["S"]


def test_ammonium_hydrosulfide_forms_where_the_fit_says_and_leaves_the_gas_at_its_boundary():
    budget = mg.element_budget(0.0)
    warm, taken = mg.remove_ammonium_hydrosulfide(budget, 250.0, 0.3, 0.0)
    assert taken == 0.0 and warm == budget
    cold, taken = mg.remove_ammonium_hydrosulfide(budget, 150.0, 0.3, 0.0)
    assert taken > 0.0 and cold["N"] < budget["N"] and cold["S"] < budget["S"]
    boundary = (48.91 - 1e4 / 150.0) / 4.15
    nitrogen_over, sulfur_over = cold["N"] / budget["N"], cold["S"] / budget["S"]
    assert math.log10(0.3) + 0.5 * math.log10(nitrogen_over) + 0.5 * math.log10(sulfur_over) == pytest.approx(boundary, abs=1e-6)


def carbon_monoxide_share(temperature_k, pressure_bar, metallicity, gravity, kzz):
    quench = mg.quench_level(temperature_k, pressure_bar, metallicity, 2.3e-3, gravity, kzz)
    if quench is None:
        return None
    moles = ge.solve(mg.element_budget(math.log10(metallicity)), max(quench[0], mg.CHEMISTRY_FLOOR_K), quench[1]).moles
    return moles["CO"] / (moles["CO"] + moles["CH4"])


def test_the_methane_to_carbon_monoxide_crossover_of_a_planet_is_near_500_k():
    # [R4863]: about 500 K, 400 K with strong mixing, up to 600 K with sluggish mixing (Jupiter's gravity, 0.3 bar).
    def crossover(metallicity, kzz):
        for temperature in range(300, 1100, 10):
            share = carbon_monoxide_share(float(temperature), 0.3, metallicity, 25.0, kzz)
            if share is not None and share > 0.5:
                return temperature
        return None
    assert 400 <= crossover(1.0, 1e8) <= 650
    assert crossover(1.0, 1e10) < crossover(1.0, 1e8) < crossover(1.0, 1e6)      # stronger mixing quenches hotter, so more CO
    assert crossover(10.0, 1e8) < crossover(1.0, 1e8)                            # a richer gas converts faster


def test_the_methane_to_carbon_monoxide_crossover_of_a_brown_dwarf_is_near_1100_k():
    # Gravity of a 30 Jupiter-mass dwarf, its 1 bar layer: 1100 +/- 100 K at the lowest mixing, lower for stronger mixing.
    def crossover(kzz):
        for temperature in range(500, 1600, 10):
            share = carbon_monoxide_share(float(temperature), mg.LAYER_PRESSURE_BAR["brown-dwarf"], 1.0, 776.0, kzz)
            if share is not None and share > 0.5:
                return temperature
        return None
    assert 1000 <= crossover(1e2) <= 1250
    assert crossover(1e5) < crossover(1e2)


def test_a_layer_that_equilibrates_has_no_quench_level():
    assert mg.quench_level(2000.0, 0.3, 1.0, 2.3e-3, 25.0, 1e8) is None


def test_the_ices_of_a_cold_giant_are_capped_at_their_saturation_pressure():
    result = mixture(0, t_layer_k=70.0, teq_k=70.0)
    for formula, fraction in result.mole_fractions.items():
        if saturation.has(formula):
            curve = saturation.saturation(formula, 70.0)
            if curve.status != saturation.ABOVE_CRITICAL:
                assert fraction * result.pressure_bar * 1e5 <= curve.pressure_pa * 1.001, formula
    assert result.mole_fractions.get("H2O", 0.0) < 1e-10 and result.mole_fractions.get("NH3", 0.0) < 1e-6
    assert result.mole_fractions["CH4"] > 1e-4                         # the methane that is not cold enough to freeze
    assert {"H2O", "NH3"} & set(result.condensed) or "ammonium_hydrosulfide" in result.condensed


def test_a_warm_giant_has_methane_water_and_ammonia_and_a_little_quenched_carbon_monoxide():
    result = mixture(1, t_layer_k=350.0, teq_k=350.0, zone="e")
    fractions = result.mole_fractions
    assert fractions["H2"] > 0.8 and 0.1 < fractions["He"] < 0.17
    assert fractions["H2O"] > 1e-3 and fractions["CH4"] > 1e-3 and fractions["NH3"] > 1e-4
    assert 1e-10 < fractions["CO"] < 0.5 * fractions["CH4"]


def test_a_hot_giant_holds_its_carbon_as_carbon_monoxide_and_its_nitrogen_as_n2():
    result = mixture(2, t_layer_k=1300.0, teq_k=1300.0, tday_k=1300.0, zone="e")
    fractions = result.mole_fractions
    assert fractions["CO"] > 10 * fractions.get("CH4", 0.0)
    assert fractions["N2"] > 10 * fractions.get("NH3", 0.0)
    assert result.constants["carbon_monoxide_frozen"] > 0.9 or "quench_k" not in result.constants


def test_the_giant_mixture_follows_the_drawn_metallicity_and_carbon_to_oxygen():
    # same seed, same regime: only the budget differs, so a richer atmosphere has more of every heavy carrier
    poor = mixture(4, t_layer_k=350.0, teq_k=350.0, mass_jupiter=0.3, mass_earth=95.0)
    rich = mixture(4, t_layer_k=350.0, teq_k=350.0, mass_jupiter=3.0, mass_earth=950.0)
    assert poor.constants["log_enrichment"] > rich.constants["log_enrichment"]   # Thorngren: lighter giants are richer in metals
    assert poor.mole_fractions["H2O"] / poor.mole_fractions["H2"] > rich.mole_fractions["H2O"] / rich.mole_fractions["H2"] * 0.5
    assert set(poor.inventory) >= {"H", "He", "C", "O", "N", "S"}
    assert poor.inventory["C"] / poor.inventory["O"] == pytest.approx(poor.constants["c_over_o"])


def test_a_carbon_rich_giant_makes_hydrogen_cyanide_and_acetylene_where_it_is_hot():
    subject = giant(t_layer_k=1000.0, teq_k=1000.0, zone="e")
    with draw.bound(5):
        recipe = RECIPES.draw_all(subject, categories="abdcefghijkl")
    recipe = dataclasses.replace(recipe, values={**recipe.values, "f7": 1.4})
    with draw.bound(6):
        rich = mg.giant_mixture(subject, recipe)
    recipe = dataclasses.replace(recipe, values={**recipe.values, "f7": 0.4})
    with draw.bound(6):
        poor = mg.giant_mixture(subject, recipe)
    assert rich.mole_fractions.get("HCN", 0.0) > 100 * poor.mole_fractions.get("HCN", 1e-30)
    assert rich.mole_fractions.get("H2O", 0.0) < poor.mole_fractions["H2O"]


def test_a_leaky_cold_trap_keeps_titanium_and_vanadium_oxides_in_the_gas_of_a_hot_jupiter():
    subject = giant(zone="h", t_layer_k=2000.0, teq_k=2000.0, tday_k=2000.0, insolation=1500.0)
    with draw.bound(7):
        recipe = RECIPES.draw_all(subject, categories="abdcefghijkl")
    kept = dataclasses.replace(recipe, binaries={**recipe.binaries, "h4": (True, None)})
    lost = dataclasses.replace(recipe, binaries={**recipe.binaries, "h4": (False, None)})
    with draw.bound(8):
        with_oxides = mg.giant_mixture(subject, kept)
    with draw.bound(8):
        without = mg.giant_mixture(subject, lost)
    assert with_oxides.mole_fractions.get("TiO", 0.0) > 1e-9
    assert without.mole_fractions.get("TiO", 0.0) < 1e-12 and without.mole_fractions.get("VO", 0.0) < 1e-12


def test_a_brown_dwarf_mixture_is_the_stars_metals_not_a_planets_enrichment():
    subject = Subject("brown-dwarf", {**JUPITER, "mass_jupiter": 30.0, "mass_earth": 9500.0, "radius_km": 70000.0, "teff_k": 900.0,
                                      "t_layer_k": 900.0, "zone": "r"})
    with draw.bound(9):
        recipe = RECIPES.draw_all(subject, categories="abdcefghijkl")
        result = mg.giant_mixture(subject, recipe)
    assert result.family == "brown-dwarf"
    assert result.constants["log_enrichment"] == pytest.approx(0.0)
    assert result.pressure_bar == mg.LAYER_PRESSURE_BAR["brown-dwarf"]
    assert result.mole_fractions["H2"] > 0.8


def test_the_host_metallicity_enriches_every_metal():
    low = mixture(3, t_layer_k=350.0, teq_k=350.0, fe_h=-0.5)
    high = mixture(3, t_layer_k=350.0, teq_k=350.0, fe_h=0.3)
    assert high.constants["log_enrichment"] - low.constants["log_enrichment"] == pytest.approx(0.8)
    assert high.mole_fractions["CH4"] > low.mole_fractions["CH4"]


def test_a_giant_without_a_radius_is_refused_by_name():
    subject = Subject("giant", {k: v for k, v in JUPITER.items() if k != "radius_km"})
    with draw.bound(0):
        recipe = RECIPES.draw_all(subject, categories="abdcefghijkl")
        with pytest.raises(MixtureError, match="radius_km"):
            mg.giant_mixture(subject, recipe)


def test_the_giant_mixture_is_deterministic_and_sums_to_one():
    first, second = mixture(11), mixture(11)
    assert first.mole_fractions == second.mole_fractions
    assert sum(first.mole_fractions.values()) == pytest.approx(1.0, abs=1e-12)
    assert all(st.SPECIES and value > 0 for value in first.mole_fractions.values())


def made_median_keeps_the_disk(subject, recipe, disk):
    from planetgen.physics import disk_condensation
    median = dataclasses.replace(recipe, values={**recipe.values, "f7": disk_condensation.RECIPE_C_OVER_O_MEDIAN})
    with draw.bound(6):
        return mg.giant_mixture(subject, median).constants["c_over_o"] == pytest.approx(disk, rel=1e-6)


def test_a_giant_with_a_disk_takes_its_carbon_to_oxygen_from_the_disk_spread_by_the_recipe():
    from planetgen.physics import disk_condensation
    subject = giant(t_layer_k=350.0, teq_k=350.0, formation_temperature_k=100.0, disk_carbon=1.0, disk_oxygen=1.0)
    with draw.bound(5):
        recipe = RECIPES.draw_all(subject, categories="abdcefghijkl")
    disk = disk_condensation.envelope_c_over_o(100.0, disk_condensation.ENVELOPE_METAL_SHARE * recipe.values["f5"])      # GEN.1471: the envelope carries a share of the bulk metals
    for f7 in (disk_condensation.RECIPE_C_OVER_O_MEDIAN, 0.4, 1.2):
        shifted = dataclasses.replace(recipe, values={**recipe.values, "f7": f7})
        with draw.bound(6):
            made = mg.giant_mixture(subject, shifted)
        assert made.constants["c_over_o"] == pytest.approx(disk * (f7 / disk_condensation.RECIPE_C_OVER_O_MEDIAN) ** disk_condensation.RECIPE_C_OVER_O_SPREAD_EXPONENT, rel=1e-6)
    assert made_median_keeps_the_disk(subject, recipe, disk)
    bare = dataclasses.replace(recipe, values={**recipe.values, "f7": 0.9})
    with draw.bound(6):
        no_disk = mg.giant_mixture(giant(t_layer_k=350.0, teq_k=350.0), bare)
    assert no_disk.constants["c_over_o"] == 0.9          # a giant with no disk keeps the recipe's value


def test_a_giant_with_a_disk_enriches_its_nitrogen_sulfur_and_noble_gases_element_by_element():
    from planetgen.physics import disk_condensation
    cold = giant(t_layer_k=350.0, teq_k=350.0, formation_temperature_k=40.0, disk_carbon=1.0, disk_oxygen=1.0)
    with draw.bound(5):
        recipe = RECIPES.draw_all(cold, categories="abdcefghijkl")
    with draw.bound(6):
        made = mg.giant_mixture(cold, recipe)
    factors = disk_condensation.envelope_enrichments(40.0, recipe.values["f5"])
    for symbol, factor in factors.items():
        assert made.constants[f"enrichment_{symbol.lower()}"] == pytest.approx(factor)
    inventory = made.inventory
    assert inventory["Ar"] / mg.SOLAR_PER_HYDROGEN["Ar"] == pytest.approx(factors["Ar"])
    assert inventory["N"] / mg.SOLAR_PER_HYDROGEN["N"] == pytest.approx(factors["N"])
    assert inventory["Fe"] / mg.SOLAR_PER_HYDROGEN["Fe"] == pytest.approx(recipe.values["f5"])     # the rock-forming metals keep the one scale
    with draw.bound(6):
        bare = mg.giant_mixture(giant(t_layer_k=350.0, teq_k=350.0), recipe)
    assert "enrichment_n" not in bare.constants                                   # a subject with no disk keeps the one scale for every metal
    assert bare.inventory["N"] / mg.SOLAR_PER_HYDROGEN["N"] == pytest.approx(recipe.values["f5"])


def test_the_host_metallicity_scales_the_stars_own_amount_of_the_enriched_elements():
    budget = mg.element_budget(0.5, disk_factors={"N": 2.0}, star_log_scale=0.3)
    assert budget["N"] == pytest.approx(mg.SOLAR_PER_HYDROGEN["N"] * 10.0 ** 0.3 * 2.0)
    assert budget["Mg"] == pytest.approx(mg.SOLAR_PER_HYDROGEN["Mg"] * 10.0 ** 0.5)


# --- GEN.1374: the lithium test, the alkali chlorides and the oxide to hydride transition of a brown dwarf's air ---------------------------------

def one_bar_air(temperature_k, **changes):
    """The equilibrium gas of a solar brown dwarf layer at 1 bar and `temperature_k` (the clouds of `remove_clouds`, the alkali clouds where the rows put them)."""
    budget = mg.element_budget(0.0, lithium_burned=changes.pop("lithium_burned", False))
    budget, _ = mg.remove_clouds(budget, temperature_k, 1.0, 0.0, keep_sodium=temperature_k >= mg.ALKALI_ONSET_K["Na"],
                                 keep_potassium=temperature_k >= mg.ALKALI_ONSET_K["K"], **changes)
    return ge.solve(budget, temperature_k, 1.0).mole_fractions


def crossing(temperature_a, temperature_b, first, second):
    """The temperature (K) between two where the mole fractions of `first` and `second` at 1 bar are equal (bisection on the log ratio)."""
    def ratio(temperature):
        fractions = one_bar_air(temperature)
        return math.log(fractions[first] / fractions[second])
    low, high = temperature_a, temperature_b
    assert ratio(low) * ratio(high) < 0.0
    for _ in range(40):
        middle = 0.5 * (low + high)
        low, high = (middle, high) if ratio(middle) * ratio(low) > 0.0 else (low, middle)
    return 0.5 * (low + high)


def test_the_titanate_and_chromium_clouds_stand_where_the_review_puts_them_at_one_bar():
    cond = mg.condensation_temperatures(1.0, 0.0)
    assert cond["titanate"] == pytest.approx(1963.0)                 # [R2725]: Ca4Ti3O10 first
    assert cond["chromium"] == pytest.approx(1518.0)                  # [R2731]: chromium metal
    assert cond["iron"] == pytest.approx(1838.0, abs=10.0)            # the review's 1843 K
    assert cond["chromium"] < cond["forsterite"] < cond["iron"] < cond["titanate"]
    richer = mg.condensation_temperatures(1.0, 1.0)
    assert richer["titanate"] > cond["titanate"] and richer["chromium"] > cond["chromium"]


def test_the_oxides_leave_before_the_iron_and_the_iron_before_the_chromium_so_the_hydrides_outlast_the_oxides():
    budget = mg.element_budget(0.0)
    window, clouds = mg.remove_clouds(budget, 1900.0, 1.0, 0.0)               # between the titanate (1963 K) and the iron (1838 K)
    assert window["Ti"] == 0.0 and window["V"] == 0.0 and window["Ca"] == 0.0
    assert window["Fe"] == budget["Fe"] and window["Cr"] == budget["Cr"] and set(clouds) == {"oxides"}
    after_iron, clouds = mg.remove_clouds(budget, 1750.0, 1.0, 0.0)             # iron gone, forsterite (about 1700 K) and chromium not yet
    assert after_iron["Fe"] == 0.0 and after_iron["Mg"] == budget["Mg"] and after_iron["Cr"] == budget["Cr"] and "chromium" not in clouds
    below, clouds = mg.remove_clouds(budget, 1450.0, 1.0, 0.0)
    assert below["Cr"] == 0.0 and "chromium" in clouds
    metal, _ = mg.remove_clouds(budget, 1450.0, 1.0, 0.0)
    assert metal["O"] > mg.remove_clouds(budget, 1800.0, 1.0, 0.0)[0]["O"] * 0.7   # chromium metal takes no oxygen


def test_titanium_and_vanadium_oxides_fade_where_the_hydrides_of_iron_and_chromium_take_over():
    hot, window, cooler, cold = (one_bar_air(t) for t in (2200.0, 1900.0, 1700.0, 1400.0))
    assert hot["TiO"] > 1e-8 and hot["VO"] > 1e-10 and hot["FeH"] > 1e-9 and hot["CrH"] > 1e-10
    assert window.get("TiO", 0.0) == 0.0 and window.get("VO", 0.0) == 0.0         # the M to L transition: the oxides are gone first
    assert window["FeH"] > 1e-9 and window["CrH"] > 1e-10                          # while the hydrides are in the gas
    assert cooler.get("FeH", 0.0) == 0.0 and cooler["CrH"] > 1e-10                 # the iron cloud takes FeH, then the chromium cloud takes CrH
    assert cold.get("FeH", 0.0) == 0.0 and cold.get("CrH", 0.0) == 0.0


def test_the_lithium_gas_changes_from_the_atom_to_the_chloride_where_the_review_puts_it():
    # [R2731: section 1.5.5]: atomic Li is the main lithium gas down to 1520 K at 1 bar, then LiCl
    assert crossing(1400.0, 1700.0, "Li", "LiCL") == pytest.approx(1520.0, abs=40.0)
    hot, cold = one_bar_air(1800.0), one_bar_air(1100.0)
    assert hot["Li"] > hot["LiCL"] and cold["LiCL"] > cold["Li"]


def test_the_potassium_gas_changes_from_the_atom_to_the_chloride_near_980_k():
    # [R2731: section 1.5.4]: K and KCl are equal at about 980 K at 1 bar
    assert crossing(900.0, 1100.0, "K", "KCL") == pytest.approx(980.0, abs=40.0)


def test_a_brown_dwarf_above_the_lithium_burning_mass_has_no_lithium_and_one_below_keeps_it():
    def subject_of(mass_jupiter):
        return Subject("brown-dwarf", {**JUPITER, "mass_jupiter": mass_jupiter, "mass_earth": mass_jupiter * 318.0, "radius_km": 70000.0,
                                       "teff_k": 1600.0, "t_layer_k": 1600.0, "zone": "r"})
    results = {}
    for mass in (30.0, 75.0):
        subject = subject_of(mass)
        with draw.bound(9):
            results[mass] = mg.giant_mixture(subject, RECIPES.draw_all(subject, categories="abdcefghijkl"))
    assert not results[30.0].constants.get("lithium_burned") and results[30.0].inventory["Li"] > 1e-10
    assert results[75.0].constants["lithium_burned"] is True and results[75.0].inventory["Li"] == 0.0
    assert any(results[30.0].mole_fractions.get(name, 0.0) > 0.0 for name in ("Li", "LiCl", "LiOH", "LiH"))
    assert not any(results[75.0].mole_fractions.get(name, 0.0) > 0.0 for name in ("Li", "LiCl", "LiOH", "LiH"))
    assert mg.burned_lithium(75.0) and not mg.burned_lithium(30.0) and not mg.burned_lithium(None)


def test_a_giant_planet_never_burns_its_lithium_and_the_potassium_cloud_takes_chlorine_with_it():
    budget = mg.element_budget(0.0)
    assert budget["Li"] > 0.0 and budget["Cl"] == pytest.approx(10.0 ** (5.31 - 12.0))
    gone, clouds = mg.remove_clouds(budget, 600.0, 0.3, 0.0, keep_potassium=False)
    assert gone["K"] == 0.0 and "potassium_chloride" in clouds
    assert gone["Cl"] == pytest.approx(budget["Cl"] - budget["K"])


def test_lithium_chloride_condenses_under_about_680_k_at_one_bar_and_takes_its_chlorine_and_lithium():
    for temperature, condensed in ((900.0, False), (760.0, False), (620.0, True), (300.0, True)):
        budget, _ = mg.remove_clouds(mg.element_budget(0.0), temperature, 1.0, 0.0, keep_sodium=False, keep_potassium=False)
        result, _, ices = mg.cap_ices(budget, temperature, 1.0, lambda trial, t=temperature: ge.solve(trial, t, 1.0))
        assert ("LiCl" in ices) is condensed, temperature
    assert mg.nasa_saturation_pa("LiCl", 700.0) == pytest.approx(5.6e-9 * 1e5, rel=0.1)
    assert mg.nasa_saturation_pa("H2O", 300.0) is None
    assert mg.nasa_saturation_pa("LiCl", 100.0) == mg.nasa_saturation_pa("LiCl", 300.0)      # under the fit read at its lower end


def test_the_ice_caps_settle_when_the_oxygen_left_sits_in_a_species_that_does_not_condense():
    """
    A cold giant layer whose water, methane and ammonia have all condensed leaves 1e-37 of oxygen per hydrogen, and 1.5e-7 of lithium holds most of it as lithium hydroxide at the
    chemistry's 300 K. Water is over its cap at 26 K by 40 orders, so each round took a sliver of the water and the equilibrium made it again from the hydroxide: 1000 rounds, then
    `MixtureError`. Seed 533 of a G2V system with PLANETS and MAX_PLANETS (Lane 21's giant-density read, 2026-11-09) is this layer; the budget, 26.49 K, 0.3 bar and 300 K are its own.
    What stays over the cap is under `MIN_CONDENSED` per hydrogen atom, the amount a condensate is dropped from the record at.
    """
    budget = {"H": 1.0, "He": 0.08203515443298177, "C": 0.007923251589436131, "N": 0.0022446996159850695, "O": 0.04823973997996018, "Ne": 0.009541950492118954,
              "S": 2.168404344971009e-19, "Ar": 2.3988329190194897e-06, "P": 2.1361766269879255e-05, "Kr": 1.31825673855641e-09, "Xe": 2.0859830264783805e-08,
              "Cl": 7.204043020470654e-06, "Li": 1.4778731499033949e-07}
    without_co = tuple(name for name in st.GASES if name != "CO")      # the species of a layer whose carbon monoxide quenched (`mixture_giant.draw`)
    result, left, ices = mg.cap_ices(budget, 26.487872842214287, 0.3, lambda trial: ge.solve(trial, 300.0, 0.3, species=without_co))
    assert set(ices) >= {"H2O", "CH4", "NH3"} and ices["H2O"] == pytest.approx(0.04823973997996018, rel=1e-6)
    assert left["O"] < mg.MIN_CONDENSED and result.mole_fractions.get("H2O", 0.0) * 0.3 < 1e-30


def test_the_ice_caps_still_refuse_a_layer_that_cycles_with_oxygen_to_spare():
    """The settling above is only for a species the budget cannot make more of than `MIN_CONDENSED`; a cycle with more to condense still raises."""
    class Cycling:
        mole_fractions = {"H2O": 1e-3}
        moles = {"H2O": 1e-3}
        total_moles = 1.0
    with pytest.raises(MixtureError, match="did not settle"):
        mg.cap_ices({"H": 1.0, "O": 0.5}, 100.0, 1.0, lambda trial: Cycling(), iterations=5)
