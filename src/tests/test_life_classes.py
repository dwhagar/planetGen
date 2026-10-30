"""
The life classes as a function of a finished body (GEN.886, `planetgen.physics.life_classes`, docs/design/habitability/life-classes.md).
"""

import copy
import json
import math

import pytest

from planetgen import tuning
from planetgen.physics import atmosphere, habitability, habitability_world, life_classes, life_scores

EARTH = {
    "body_type": "t", "planet_class": "M", "atmospheric_pressure_pa": 101325.0, "surface_temperature_k": 288.0, "gravity_g": 1.0,
    "radius_km": 6371.0, "mass_kg": 5.972e24, "hydrosphere": "surface ocean", "ocean_class": "ocean", "ocean_ph": 8.1,
    "water_activity": 0.98, "phosphorus": "limited", "water_mass_fraction": 2.3e-4, "surface_dose_msv_yr": 2.4,
    "dose_ground_msv_yr": 2.4, "uv_surface_index": 5.0, "energy_flux_w_m2": 80.0, "internal_heat_flux_w_m2": 0.09,
    "tidal_flux_w_m2": None, "star_teff_k": 5772.0, "star_luminosity_w": 3.828e26, "star_distance_au": 1.0,
    atmosphere.GAS_FIELD: {"O2": 21278.0, "N2": 78084.0, "CO2": 40.0, "Ar": 934.0, "H2O": 1000.0},
}


def world(**over):
    values = copy.deepcopy(EARTH)
    values.update(over)
    return values


def test_the_rule_file_has_every_class_of_the_design():
    assert life_classes.IDS == ("1a", "1b", "1c", "2", "3a", "3b", "4a", "4b", "4c", "5a", "5b", "6", "7", "8", "9")
    assert sorted({c.number for c in life_classes.CLASSES}) == list(range(1, 10))


def test_every_number_of_the_rule_file_carries_a_tier_a_tweak_and_a_basis():
    seen = 0

    def walk(node):
        nonlocal seen
        if isinstance(node, dict):
            if any(isinstance(v, (int, float)) and not isinstance(v, bool) for v in node.values()) and "quantity" in node or "medium" in node \
                    or "source" in node:
                assert node["tier"] in life_classes.TIERS and node["tweak"] in life_classes.TWEAKS and node["basis"], node
                seen += 1
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for value in node:
                walk(value)

    with open(life_classes.DATA_FILE, encoding="utf-8") as handle:
        walk(json.load(handle)["classes"])
    assert seen >= 60


@pytest.mark.parametrize("damage", [
    lambda row: row.pop("tier"),
    lambda row: row["solvent"].update(medium="mercury"),
    lambda row: row["solvent"].pop("tweak"),
    lambda row: row["energy"].update(source="gravity"),
    lambda row: row["gates"].append({"quantity": "pressure_pa", "tier": "T9", "tweak": "L3", "basis": "no operator"}),
    lambda row: row["windows"].append({"quantity": "temperature_k", "low": 300.0, "high": 200.0, "sigma_low": 1.0, "sigma_high": 1.0,
                                       "tier": "T9", "tweak": "L3", "basis": "low over high"}),
    lambda row: row.update(supply={"quantity": "x", "placeholder": 2.0, "tier": "T9", "tweak": "L4", "basis": "too big"}),
])
def test_a_rule_file_with_a_missing_or_broken_number_is_refused(damage):
    with open(life_classes.DATA_FILE, encoding="utf-8") as handle:
        row = copy.deepcopy(json.load(handle)["classes"][0])
    damage(row)
    with pytest.raises(life_classes.LifeClassError):
        life_classes._check(row)


def test_class_1a_on_an_earthlike_world_is_todays_phi_bio():
    values = world()
    index = habitability_world.compute(copy.deepcopy(values))
    scores = life_classes.score_classes(values)
    assert scores["1a"].likelihood == pytest.approx(index["phi_bio"], abs=1e-3)
    assert scores["1a"].score == life_scores._from_likelihood(index["phi_bio"])
    assert scores["1a"].failed == [] and scores["1a"].missing == []


def test_the_oxygen_a_world_has_makes_the_anoxygenic_class_unsupported_and_none_makes_it_supported():
    assert life_classes.score_classes(world())["2"].score == 0
    anoxic = world(**{atmosphere.GAS_FIELD: {"N2": 78000.0, "CO2": 20000.0, "CH4": 1000.0}})
    scores = life_classes.score_classes(anoxic)
    assert scores["2"].score > 0 and scores["2"].placeholders == []  # the donors are worked out from the heat (GEN.887)
    assert scores["1a"].score > 0  # Class 1 makes the oxygen; it does not need it (GEN.855)


def test_a_world_under_ice_has_no_light_class_but_has_the_chemical_ones():
    scores = life_classes.score_classes(world(hydrosphere="ice-covered ocean"))
    assert all(scores[c].score == 0 for c in ("1a", "1b", "1c", "2"))
    assert scores["3a"].score > 0 and scores["4b"].score > 0


def test_a_dry_airless_rock_supports_nothing():
    scores = life_classes.score_classes(world(hydrosphere="dry", atmospheric_pressure_pa=0.0, surface_temperature_k=250.0,
                                              ocean_ph=None, water_activity=None, water_mass_fraction=0.0, **{atmosphere.GAS_FIELD: {}}))
    assert all(result.score == 0 for result in scores.values())
    assert life_classes.best(scores) is None


def test_a_methane_sea_world_supports_class_7_and_not_class_1a():
    temperature = 94.0
    methane_kpa = atmosphere.vapour_pressure_kpa("ch4", temperature)
    titan = world(hydrosphere="dry", surface_temperature_k=temperature, atmospheric_pressure_pa=146700.0, ocean_ph=None, water_activity=None,
                  water_mass_fraction=0.3, internal_heat_flux_w_m2=0.01, star_luminosity_w=3.828e26, star_distance_au=9.5,
                  energy_flux_w_m2=1e-3, **{atmosphere.GAS_FIELD: {"N2": 140000.0, "CH4": methane_kpa * 1000.0}})
    scores = life_classes.score_classes(titan)
    assert scores["7"].score > 0
    assert scores["1a"].score == 0 and scores["7"].likelihood > 0.0 and scores["7"].placeholders == []


def test_the_phase_of_carbon_dioxide_follows_its_triple_and_critical_points():
    assert life_classes.co2_phase(280.0, 6.0e6) == "liquid"
    assert life_classes.co2_phase(280.0, 1.0e6) == "gas"
    assert life_classes.co2_phase(320.0, 9.0e6) == "supercritical"
    assert life_classes.co2_phase(320.0, 1.0e6) == "gas"
    assert life_classes.co2_phase(200.0, 1.0e6) == "solid"
    assert life_classes.co2_saturation_pa(life_classes.CO2_TRIPLE_K) == pytest.approx(life_classes.CO2_TRIPLE_PA)
    assert life_classes.co2_saturation_pa(life_classes.CO2_CRITICAL_K) == pytest.approx(life_classes.CO2_CRITICAL_PA)


def test_a_dense_carbon_dioxide_world_supports_class_6():
    dense = world(hydrosphere="dry", surface_temperature_k=280.0, atmospheric_pressure_pa=6.0e6, ocean_ph=None, water_activity=None,
                  water_mass_fraction=0.0, h2_supply=0.5, **{atmosphere.GAS_FIELD: {"CO2": 5.7e6, "N2": 3.0e5}})
    result = life_classes.score_classes(dense)["6"]
    assert result.score > 0 and result.likelihood <= 0.2  # the 0.1 solvent factor, lifted by the blend of the other terms
    thin = world(hydrosphere="dry", surface_temperature_k=280.0, atmospheric_pressure_pa=1.0e5, ocean_ph=None, water_activity=None,
                 water_mass_fraction=0.0, h2_supply=0.5, **{atmosphere.GAS_FIELD: {"CO2": 9.5e4, "N2": 5.0e3}})
    assert life_classes.score_classes(thin)["6"].score == 0


def test_the_acid_and_ammonia_classes_wait_for_their_quantities_and_score_once_they_are_given():
    dry = dict(hydrosphere="dry", surface_temperature_k=330.0, ocean_ph=None, water_activity=None, water_mass_fraction=0.0,
               internal_heat_flux_w_m2=0.1, **{atmosphere.GAS_FIELD: {"CO2": 80000.0, "N2": 10000.0, "SO2": 10000.0}})
    # the body's own air holds no water, so no acid forms: the quantity is worked out (GEN.888) and is 0
    waiting = life_classes.score_classes(world(**dry))
    assert waiting["5b"].score == 0 and "acid_medium" in waiting["5b"].failed
    given = life_classes.score_classes(world(acid_medium=2.0, **dry))
    assert given["5b"].score > 0 and given["5b"].missing == []
    # a cloud layer alone is not a pool: Class 5b needs level 2, Class 5a level 1
    cloud = life_classes.score_classes(world(acid_medium=1.0, **dry))
    assert cloud["5b"].score == 0 and cloud["5a"].score > 0
    cold = dict(dry, surface_temperature_k=220.0)
    # a body that stores no ice share has nothing to work the ammonia out from
    assert life_classes.score_classes(world(**cold))["8"].missing == ["ammonia_water_mass_fraction"]
    # a dry rock has no water in touch with its rock to make hydrogen (GEN.887), so the hydrogen of this one is given
    assert life_classes.score_classes(world(ammonia_water_mass_fraction=0.2, h2_supply=0.5, **cold))["8"].score > 0


def test_the_ammonia_and_the_silicon_of_a_body_are_worked_out_from_its_bulk_composition():
    cold = dict(hydrosphere="dry", surface_temperature_k=220.0, ocean_ph=None, water_activity=None, internal_heat_flux_w_m2=0.1,
                h2_supply=0.5, **{atmosphere.GAS_FIELD: {"N2": 10000.0}})
    icy = life_classes.score_classes(world(water_mass_fraction=0.2, ice_mass_fraction=0.4, **cold))["8"]
    assert icy.missing == [] and icy.score > 0           # 0.1 of the other ice over the water: ammonia at 0.2 of the water
    drier = life_classes.score_classes(world(water_mass_fraction=0.2, ice_mass_fraction=0.1, **cold))["8"]
    assert drier.score == 0 and "ammonia_water_mass_fraction" in drier.failed
    assert life_classes.derived_quantities({"rock_mg_si_ratio": 0.8})["silicon_rich_crust"] == 1.0
    assert life_classes.derived_quantities({"rock_mg_si_ratio": 1.3})["silicon_rich_crust"] == 0.0
    assert "silicon_rich_crust" not in life_classes.derived_quantities({})
    assert "silicon_rich_crust" not in life_classes.derived_quantities({"rock_mg_si_ratio": 0.8, "silicon_rich_crust": 0.0})


def test_a_stored_supply_replaces_the_one_worked_out_from_the_body():
    from planetgen.physics import life_supply
    under_ice = world(hydrosphere="ice-covered ocean")
    derived = life_supply.supplies(under_ice)["h2_supply"]
    ref = life_classes.describe("3a").supplies[0]["ref"]
    worked_out = life_classes.score_classes(under_ice)["3a"]
    stored = life_classes.score_classes(world(hydrosphere="ice-covered ocean", h2_supply=0.25))["3a"]
    assert worked_out.placeholders == [] and stored.placeholders == []
    assert stored.likelihood == pytest.approx(worked_out.likelihood * min(1.0, 0.25 / ref) / min(1.0, derived / ref))


def test_a_cool_red_star_favours_the_far_red_pigment_relative_to_the_visible_one_more_than_the_sun_does():
    def far_red_against_visible(teff):
        scores = life_classes.score_classes(world(star_teff_k=teff))
        return scores["1c"].band_ratio / scores["1a"].band_ratio

    assert far_red_against_visible(3000.0) > far_red_against_visible(5772.0)
    hot = life_classes.score_classes(world(star_teff_k=9000.0))
    assert hot["1b"].band_ratio > hot["1c"].band_ratio


def test_the_share_of_photons_in_a_band_follows_the_blackbody():
    sun = life_classes.band_fraction(5772.0, (400.0, 700.0))
    assert 0.15 < sun < 0.35
    # a cool star puts more of its visible-and-red photons in the far red, a hot one more in the blue
    assert (life_classes.band_fraction(3000.0, (700.0, 800.0)) / life_classes.band_fraction(3000.0, (400.0, 700.0))
            > life_classes.band_fraction(5772.0, (700.0, 800.0)) / sun)
    assert life_classes.band_fraction(9000.0, (300.0, 500.0)) > life_classes.band_fraction(4000.0, (300.0, 500.0))
    assert life_classes.band_fraction(None, (400.0, 700.0)) == 0.0
    assert life_classes.band_fraction(5772.0, (100.0, 20000.0)) == pytest.approx(1.0, abs=0.01)


def test_a_gas_giant_or_a_body_without_a_temperature_is_scored_nowhere():
    assert life_classes.score_classes(world(body_type="g")) == {}
    assert life_classes.score_classes(world(surface_temperature_k=None)) == {}


def test_the_gate_operators():
    found = {"a": 5.0, "kind": "x"}
    for gate, ok in ((dict(quantity="a", min=5.0), True), (dict(quantity="a", lt=5.0), False), (dict(quantity="kind", equals="x"), True),
                     (dict(quantity="kind", not_equals="x"), False), (dict(quantity="kind", **{"in": ["x", "y"]}), True),
                     (dict(quantity="kind", not_in=["x"]), False), (dict(quantity="gone", min=0.0), False)):
        assert life_classes.holds([gate], found)[0] is ok
    assert life_classes.holds([dict(quantity="gone", min=0.0)], found)[2] == ["gone"]


def test_the_score_steps_like_the_life_scores():
    for likelihood in (0.0, 0.004, 0.005, 0.11, 0.5, 0.97, 1.0):
        assert life_classes.score_of(likelihood) == life_scores._from_likelihood(likelihood)
    assert life_classes.score_of(None) == 0 and life_classes.score_of(1.0) == tuning.LIFE_SCORE_MAX


def test_the_best_class_is_the_largest_likelihood_with_the_lower_number_winning_a_tie():
    scores = life_classes.score_classes(world())
    top = life_classes.best(scores)
    assert top.id == "1a"  # 1a, 1b and 1c tie on a Sun-like star: the lower class wins
    assert all(top.likelihood >= s.likelihood for s in scores.values())


def test_every_class_is_a_likelihood_between_zero_and_one_for_odd_worlds():
    for extra in (dict(surface_temperature_k=5.0), dict(surface_temperature_k=900.0), dict(atmospheric_pressure_pa=1e9),
                  dict(surface_dose_msv_yr=1e9, dose_ground_msv_yr=1e9), dict(star_teff_k=None, star_luminosity_w=None)):
        for result in life_classes.score_classes(world(**extra)).values():
            assert 0.0 <= result.likelihood <= 1.0 and math.isfinite(result.likelihood)
