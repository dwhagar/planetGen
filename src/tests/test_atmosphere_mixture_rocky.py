# tests/test_atmosphere_mixture_rocky.py

"""The rocky outgassed mixture (GEN.451, GEN.479): chemistry from the drawn regimes, not a template per class."""

import math

import pytest

from planetgen.generation.possibility_table import Subject
from planetgen.physics import gases, mixture_rocky as md
from planetgen.physics.atmosphere_recipes import RECIPES, RecipeDraw
from planetgen.physics.mixture_base import MixtureError
from planetgen.util import draw

EARTH = dict(mass_earth=1.0, ts_k=288.0, pressure_bar=1.0, age_gy=4.5, long_lived_geology=True, ocean_present=True,
             snow_ratio=0.5, zone="e", host_type="G", is_moon=False, insolation=1.0, luminosity_solar=1.0, water_lost=False,
             runaway_ratio=0.5, teq_k=255.0, xuv_index=1.0, v_esc_kms=11.2, radius_km=6371.0, host_mass_solar=1.0)


def subject(body_class="terrestrial", **changes):
    return Subject(body_class, {**EARTH, **changes})


def mixture(seed, **changes):
    with draw.bound(seed):
        return md.outgassed_mixture(subject(**changes))


def test_the_mole_fractions_sum_to_one_and_name_only_gases_of_the_table():
    for seed in range(12):
        result = mixture(seed)
        assert sum(result.mole_fractions.values()) == pytest.approx(1.0, abs=1e-12)
        assert all(gases.has(formula) for formula in result.mole_fractions)
        assert all(value > 0 for value in result.mole_fractions.values())


def test_a_draw_is_reproducible_and_the_seed_changes_it():
    assert mixture(5).mole_fractions == mixture(5).mole_fractions
    assert mixture(5).mole_fractions != mixture(6).mole_fractions


def test_the_redox_of_a_young_or_inert_body_is_the_magma_ocean_stage_and_of_an_old_tectonic_one_the_later_mantle():
    for seed in range(8):
        with draw.bound(seed):
            recipe = RECIPES.draw_all(subject(long_lived_geology=False), categories=md.ROCKY_CATEGORIES)
            inert = md.draw_redox(subject(long_lived_geology=False), recipe, {})
        assert inert == recipe.values["a1"]
    differences = []
    for seed in range(8):
        with draw.bound(seed):
            recipe = RECIPES.draw_all(subject(), categories=md.ROCKY_CATEGORIES)
            constants = {}
            old = md.draw_redox(subject(age_gy=40.0), recipe, constants)
            young = md.draw_redox(subject(age_gy=0.0), recipe, {})
        assert young == recipe.values["a1"]
        assert old == pytest.approx(recipe.values["a2"], abs=0.01)
        differences.append(old - young)
    assert max(differences) > 1.0


def test_a_white_dwarf_debris_body_takes_its_own_redox_row_and_a_small_body_the_building_block_row():
    with draw.bound(1):
        recipe = RECIPES.draw_all(subject(host_type="WD"), categories=md.ROCKY_CATEGORIES)
    if "a4" in recipe.values:
        assert md.draw_redox(subject(host_type="WD"), recipe, {}) == recipe.values["a4"]
    assert md.draw_redox(subject(), RecipeDraw(values={"a3": -0.7}), {}) == -0.7
    assert md.draw_redox(subject(), RecipeDraw(values={"a2": 3.0}), {}) == 3.0
    assert md.draw_redox(subject(), RecipeDraw(), {}) is None


def test_the_later_mantle_is_never_more_reduced_than_the_magma_ocean_surface_of_the_same_body():
    """GEN.482b: [R1075, R1951] find that the later processes only oxidise the mantle."""
    old = subject(long_lived_geology=True, age_gy=9.0)
    recipe = RecipeDraw(values={"a1": 1.2, "a2": -0.5})
    assert md.draw_redox(old, recipe, {}) == pytest.approx(1.2)
    recipe = RecipeDraw(values={"a1": 0.0, "a2": 3.0})
    assert md.draw_redox(old, recipe, {}) == pytest.approx(3.0, abs=0.05)


def test_a_reduced_mantle_gives_carbon_monoxide_and_an_oxidised_one_carbon_dioxide_and_sulfur_dioxide():
    reduced, oxidised = [], []
    for seed in range(40):
        result = mixture(seed, long_lived_geology=False, ts_k=1000.0, pressure_bar=1.0)
        fractions = result.mole_fractions
        (reduced if result.redox_diw < 0.0 else oxidised if result.redox_diw > 3.0 else []).append(fractions)
    assert reduced
    assert all(f.get("CO", 0.0) > f.get("CO2", 0.0) * 0.3 for f in reduced)
    for seed in range(40):
        with draw.bound(seed):
            result = md.outgassed_mixture(subject(mass_earth=3.0, ts_k=1000.0))
        if result.redox_diw > 3.0:
            assert result.mole_fractions.get("CO", 0.0) < result.mole_fractions.get("CO2", 0.0)
            return
    pytest.fail("no oxidised draw in 40 seeds")


def test_the_inventory_follows_the_delivered_water_class():
    # TEST.401: GEN.802 made the rows follow the Bern population in the formation distance, so at the Earth's snow ratio of 0.5 the draw
    # holds b1, b3 and b4 and b2 (the thin water class) only inside the snow line; the sweep is at a ratio of 0.3 where all of b1 to b3 appear.
    seen = {}
    for seed in range(300):     # b2 is 6 percent at this formation distance (GEN.802), so 60 seeds could miss it
        near = subject(ocean_present=False, snow_ratio=0.3)
        with draw.bound(seed):
            recipe = RECIPES.draw_all(near, categories=md.ROCKY_CATEGORIES)    # a world with an ocean also weathers (GEN.674)
            inventory = md.rocky_inventory(near, recipe, {})
        seen.setdefault(recipe.groups["b"].name, []).append(inventory)
    assert {"b1", "b2"} <= set(seen)
    dry_h = max(i["H"] for i in seen["b1"])
    wet_h = min(i["H"] for i in seen["b2"])
    assert dry_h < wet_h
    # CO2 per H2O in the wet classes is Ortenzi's 0.2 (by mass share 0.2 of the water)
    for inventory in seen["b2"]:
        water_moles = inventory["H"] / 2
        assert inventory["C"] * 44.0095 / (water_moles * 18.01528) == pytest.approx(0.2, rel=1e-6)


def test_argon_grows_with_the_age_of_the_body():
    with draw.bound(3):
        recipe = RECIPES.draw_all(subject(), categories=md.ROCKY_CATEGORIES)
    young, old = [], []
    for age, bucket in ((0.5, young), (4.5, old)):
        with draw.bound(9):
            bucket.append(md.rocky_inventory(subject(age_gy=age), recipe, {})["Ar"])
    assert old[0] > young[0] > 0


def test_the_equilibrium_temperature_is_the_surface_one_above_the_quench_and_the_quench_below_it():
    hot = mixture(2, ts_k=1100.0, long_lived_geology=False)
    assert hot.temperature_k == pytest.approx(1100.0)
    cold = mixture(2, ts_k=150.0, long_lived_geology=False)
    low, high, _ = md.QUENCH_TEMPERATURE_K
    assert low <= cold.temperature_k <= high
    assert cold.constants["quench_temperature_k"] == cold.temperature_k


def test_an_earth_like_oxidised_volcano_gas_is_water_carbon_dioxide_and_sulfur_dioxide_at_the_vent():
    # Part 05 3.4: dIW+4, 1200 K, 1 bar: H2O .78, CO2 .16, SO2 .021 for a wet inventory; the vent temperature is the first stage.
    for seed in range(80):
        with draw.bound(seed):
            result = md.outgassed_mixture(subject(ts_k=1200.0, pressure_bar=1.0))
        if result.drawn["b"] in ("b2", "b3") and result.redox_diw > 3.0 and result.constants["magmatic_temperature_k"] < 1200.0 + 1.0:
            fractions = result.mole_fractions
            assert next(iter(fractions)) == "H2O"
            assert fractions["CO2"] > max(fractions.get("CO", 0.0), fractions.get("CH4", 0.0)), "CO2 is the carbon gas of an oxidised vent"
            return
    pytest.skip("no wet, oxidised draw below 1200 K in 80 seeds")


def test_a_subject_without_a_quantity_the_tables_read_is_refused_by_name():
    with pytest.raises(MixtureError, match="pressure_bar"):
        with draw.bound(0):
            md.outgassed_mixture(Subject("terrestrial", {k: v for k, v in EARTH.items() if k != "pressure_bar"}))
    with pytest.raises(MixtureError):
        with draw.bound(0):
            md.outgassed_mixture(Subject("terrestrial", {k: v for k, v in EARTH.items() if k != "ocean_present"}))


def test_every_quantity_of_the_mixture_draw_is_in_the_result_for_the_record():
    result = mixture(7)
    assert {"magmatic_temperature_k", "quench_temperature_k", "nitrogen_per_carbon", "sulfur_per_carbon"} <= set(result.constants)
    assert set(result.inventory) == {"H", "C", "N", "S", "Ar", "Cl", "F"}       # the halogens of GEN.524 part 2
    assert {"chlorine_wt_percent", "fluorine_wt_percent"} <= set(result.constants)
    assert result.drawn["a"] in ("a5", "a6", "a7")
    assert math.isfinite(result.redox_diw)
    assert result.family == "rocky"
