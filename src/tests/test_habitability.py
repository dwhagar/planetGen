# tests/test_habitability.py

"""
The habitability index's reference maths (GEN.84,
`planetgen.physics.habitability`): the worked examples in
docs/design/habitability-index.md section 2 reproduce, and the bands and
derived constants hold.
"""

import math

import pytest

from planetgen.physics import habitability as hab
from planetgen.physics.habitability import World, scores


EARTH = World(101.3, 15, {"N2": 0.78, "O2": 0.21, "CO2": 0.0004}, ph=8.1, relative_humidity=0.7)

EXAMPLES = {
    "Earth": (EARTH, ("Blue",) * 6, 1.0, 0.969, 0.968, 1.0, "ideal"),
    "Dense CO2": (World(300, 40, {"CO2": 0.95, "N2": 0.05}, ph=6.5, relative_humidity=0.8, energy_flux_w_m2=40),
                  ("Green", "Yellow", "Red", "Blue", "Blue", "Red"), 0.6905, 0.969, 0.0, 0.479, "sealed suit"),
    "Hycean": (World(2000, 60, {"H2": 0.9, "He": 0.1}, ph=7, relative_humidity=0.99, energy_flux_w_m2=20,
                     life_stage="microbial"),
               ("Yellow", "Yellow", "Green", "Blue", "Blue", "Yellow"), 0.7265, 0.969, 0.0, 0.022, "sealed suit"),
    "Photochemical CO": (World(100, 10, {"N2": 0.9, "CO2": 0.05, "CO": 0.05}, ph=7.5),
                         ("Blue", "Blue", "Red", "Blue", "Blue", "Red"), 0.7368, 0.969, 0.0, 1.0, "sealed suit"),
    "Ice-sealed ocean": (World(1e-6, -170, {}, gravity_ms2=1.3, ph=8, energy_flux_w_m2=1e-3, surface_dose_sv_yr=0.05,
                               habitat_dose_sv_yr=0.002, water_source="ice", relative_humidity=0.05,
                               life_stage="microbial", life_on_surface=False),
                         ("Red", "Red", "Green", "Blue", "Green", "Red"), 0.0, 0.533, 0.0, 0.468, "sealed suit"),
    "Mars-analog brine": (World(0.6, -60, {"CO2": 0.95, "N2": 0.03}, gravity_ms2=3.71, water_activity=0.5,
                                chaotropicity_kj_kg=90, ph=8, inventories={"C": 0.05, "N": 0.01, "H": 0.05},
                                phosphate_umol_l=0.5, energy_flux_w_m2=40, water_source="ice",
                                relative_humidity=0.05),
                          ("Red", "Red", "Green", "Red", "Yellow", "Red"), 0.3999, 0.006, 0.0, 0.458, "sealed suit"),
    "Dune world": (World(90, 45, {"N2": 0.78, "O2": 0.21, "CO2": 0.001}, relative_humidity=0.05, water_activity=0.7,
                         ph=8, water_source="vapour", inventories={"C": 0.5, "N": 0.5, "H": 0.05},
                         life_stage="macroscopic"),
                   ("Blue", "Green", "Blue", "Yellow", "Blue", "Red"), 0.7499, 0.638, 0.637, 0.795, "ideal"),
    "Europan ocean": (World(1e-9, -160, {}, gravity_ms2=1.31, ph=9, energy_flux_w_m2=1e-4, surface_dose_sv_yr=2000,
                            habitat_dose_sv_yr=0.002, water_source="ice", relative_humidity=0.05,
                            life_stage="microbial", life_on_surface=False),
                      ("Red", "Red", "Green", "Green", "Red", "Red"), 0.0, 0.395, 0.0, 0.276,
                      "full life support with radiation hardening"),
}
"""dict: `(world, the six tiers, PHI-4, PHI_bio, PHI_cpx, Phi_tech, kit)`. The Earth world has a
biosphere (its `life_stage` is set below); the others carry the stage the example needs."""
EARTH.life_stage = "biosphere"


@pytest.mark.parametrize("name", list(EXAMPLES))
def test_the_worked_examples_reproduce(name):
    world, tiers, phi4, bio, cpx, tech, kit = EXAMPLES[name]
    got = scores(world)
    assert tuple(got[d][1] for d in hab.DOMAINS) == tiers
    assert got["PHI-4"] == pytest.approx(phi4, abs=1e-3)
    assert (got["PHI_bio"], got["PHI_cpx"], got["Phi_tech"]) == pytest.approx((bio, cpx, tech), abs=2e-3)
    assert got["equipment"] == kit


def test_the_mask_limit_is_derived_from_the_armstrong_limit():
    assert hab.MASK_MIN_KPA == pytest.approx(14.3)
    pure_o2 = World(hab.MASK_MIN_KPA, 20, {"O2": 1.0})
    assert hab.inspired_o2_kpa(pure_o2) == pytest.approx(hab.INSPIRED_O2_MIN_KPA)


def test_band_edges_score_their_tier_boundaries():
    band = hab.TEMPERATURE_BAND
    assert hab.band_score(10, band) == 1.0
    assert hab.band_score(45, band) == pytest.approx(0.75)
    assert hab.band_score(122, band) == pytest.approx(0.4)
    assert hab.band_score(200, band) == 0.0
    assert hab.band_score(-50, band) == pytest.approx(0.4)
    assert [hab.tier(s) for s in (1.0, 0.8, 0.5, 0.2)] == ["Blue", "Green", "Yellow", "Red"]


def test_surface_dose_fits_the_moon_mars_and_earth():
    assert hab.gcr_surface_dose_sv_yr(0) == pytest.approx(0.502)
    mars = hab.gcr_surface_dose_sv_yr(hab.column_mass_g_cm2(0.6, 3.71))
    assert 0.2 < mars < 0.3
    earth = hab.gcr_surface_dose_sv_yr(hab.column_mass_g_cm2(101.3, 9.81))
    assert earth == pytest.approx(0.0024, abs=2e-4)


def test_an_uncompensable_wet_bulb_is_never_better_than_yellow():
    humid = World(40, 38, {"N2": 0.79, "O2": 0.21}, relative_humidity=0.9)
    assert hab.wet_bulb_c(38, 0.9) > hab.WET_BULB_CRIT_C
    assert hab.tier(hab.temperature_score(humid)) == "Yellow"


def test_no_solvent_means_no_microbial_life():
    dry = World(101.3, 15, {"N2": 0.79, "O2": 0.21}, solvent=None)
    assert hab.phi_bio(dry) == 0.0
    assert hab.phi_cpx(dry) == 0.0


def test_dose_tiers_are_blue_to_20_msv_green_to_100_msv_yellow_to_1_sv():
    """GEN.148: Boss's defaults (2026-10-10), with no gap between 100 mSv/yr and the Red edge."""
    def colour(sv):
        return hab.tier(hab.radiation_score(World(101.3, 15, {"O2": 0.21}, surface_dose_sv_yr=sv)))
    assert [colour(sv) for sv in (0.002, 0.02, 0.03, 0.1, 0.236, 0.52, 0.99, 1.1, 50.0)] == [
        "Blue", "Blue", "Green", "Green", "Yellow", "Yellow", "Yellow", "Red", "Red"]
    assert hab.band_score(1000.0, hab.DOSE_BAND, log=True) == 0.0
    assert hab.SINGLE_YEAR_DOSE_SV == 0.05


def test_l_rad_reads_absorbed_dose_in_gy_not_effective_dose_in_sv():
    """Mars's surface is 0.236 Sv/yr in 0.077 Gy/yr, 130 times under the 10 Gy/yr line (B8)."""
    mars = World(0.6, -60, {"CO2": 0.95}, gravity_ms2=3.71, surface_dose_sv_yr=0.236)
    assert mars.habitat_dose_gy_yr == pytest.approx(0.077, abs=1e-3)
    assert hab.l_rad(mars) == 1.0
    # 31 Sv/yr is 10 Gy/yr: the threshold is in Gy, and a dose in Sv reaches it 3 times later.
    edge = World(0.6, -60, {}, surface_dose_sv_yr=10.0 * hab.SV_PER_GY)
    assert hab.l_rad(edge) == 1.0
    assert hab.l_rad(World(0.6, -60, {}, surface_dose_sv_yr=(10.0 + 100.0) * hab.SV_PER_GY)) == pytest.approx(
        math.exp(-0.5), abs=1e-6)


def test_phi_bio_is_the_liebig_blend_of_the_four_likelihoods():
    world = World(101.3, 15, {"O2": 0.21}, water_activity=0.64, ph=8.1)
    terms = [hab.l_solv(world), hab.l_chem(world), hab.l_ener(world), hab.l_rad(world)]
    geometric = math.prod(terms) ** 0.25
    assert min(terms) < 0.9 * geometric  # one weak term
    expected = min(terms) ** 0.65 * math.prod(terms) ** (0.35 / 4)
    assert hab.phi_bio(world) == pytest.approx(expected)
    assert math.prod(terms) < hab.phi_bio(world) < geometric  # kinder than the product, stricter than the mean
    # A term at 0 is a gate; with one term short and the rest at 1 the exponents add to 0.65 + 0.35 / 4.
    assert hab.phi_bio(World(101.3, 15, {}, solvent=None)) == 0.0
    one_short = World(101.3, 15, {"O2": 0.21}, energy_flux_w_m2=hab.ENERGY_FLUX_REF_W_M2)
    assert hab.phi_bio(one_short) == pytest.approx(hab.l_chem(one_short) ** (0.65 + 0.35 / 4))


# --- GEN.207: six domains, the colour gates and Suitability ---------------------------------

BLUE_ALL = {name: 1.0 for name in hab.DOMAINS}
SCORE_OF = {"Blue": 1.0, "Green": 0.8, "Yellow": 0.5, "Red": 0.2}


def _overall(**colours):
    """The overall score and colour with the named domains at the given colours and the rest Blue."""
    domains = {name: SCORE_OF[colours.get(name, "Blue")] for name in hab.DOMAINS}
    score = hab.overall_score(domains)
    return score, hab.tier(score)


def test_radiation_and_temperature_cap_the_overall_colour_at_their_own():
    assert _overall(Radiation="Yellow")[1] == "Yellow"
    assert _overall(Radiation="Red")[1] == "Red"
    assert _overall(Temperature="Yellow")[1] == "Yellow"
    assert _overall(Temperature="Red")[1] == "Red"
    assert _overall(Radiation="Green")[1] == "Green"  # a Green one only stops Blue


def test_pressure_red_is_red_and_yellow_is_yellow_unless_a_full_suit_is_needed_anyway():
    assert _overall(Pressure="Red")[1] == "Red"
    assert _overall(Pressure="Yellow")[1] == "Yellow"
    # Air in Red already needs a sealed suit, so the pressure adds nothing: Air Red alone gives Yellow.
    assert _overall(Air="Red")[1] == "Yellow"
    assert _overall(Pressure="Red", Air="Red")[1] == "Yellow"
    assert _overall(Pressure="Yellow", Air="Red")[1] == "Yellow"
    # Temperature or Radiation Red is Red whatever the pressure does.
    assert _overall(Pressure="Red", Temperature="Red")[1] == "Red"


def test_air_and_water_in_yellow_or_red_cap_the_overall_at_yellow_never_lower():
    for colour in ("Yellow", "Red"):
        assert _overall(Air=colour)[1] == "Yellow"
        assert _overall(Water=colour)[1] == "Yellow"
    assert _overall(Air="Green", Water="Green")[1] == "Green"


def test_a_red_suitability_drops_the_overall_to_yellow_but_never_to_red():
    assert _overall(Suitability="Red")[1] == "Yellow"
    assert _overall(Suitability="Yellow")[1] == "Green"


def test_blue_overall_needs_all_six_blue():
    assert hab.overall_score(BLUE_ALL) == 1.0
    for name in hab.DOMAINS:
        domains = {**BLUE_ALL, name: 0.99}
        assert hab.overall_score(domains) < 1.0 and hab.tier(hab.overall_score(domains)) == "Green"


def test_the_gates_come_first_then_the_geometric_mean_inside_them():
    # A Yellow temperature among Blues has a mean of 0.5 ** (1/6), 0.89, which the Yellow gate pulls under 0.75.
    domains = {**BLUE_ALL, "Temperature": 0.5}
    assert math.prod(domains.values()) ** (1 / 6) > 0.85
    assert hab.overall_score(domains) == pytest.approx(0.7499)
    # A mean already under its gate stays the mean.
    low = {name: 0.6 for name in hab.DOMAINS}
    assert hab.overall_score(low) == pytest.approx(0.6)
    # A zero in Radiation, Temperature or Pressure makes it 0; one in Air, Water or Suitability cannot.
    assert hab.overall_score({**BLUE_ALL, "Radiation": 0.0}) == 0.0
    for name in ("Air", "Water", "Suitability"):
        assert hab.tier(hab.overall_score({**BLUE_ALL, name: 0.0})) == "Yellow"


def _core(**colours):
    return {name: SCORE_OF[colours.get(name, "Blue")] for name in hab.CORE_DOMAINS}


def test_suitability_is_red_if_radiation_temperature_pressure_or_water_is_red():
    for name in ("Radiation", "Temperature", "Pressure", "Water"):
        assert hab.suitability_ceiling(_core(**{name: "Red"})) == "Red", name


def test_suitability_is_at_most_yellow_if_any_of_the_first_four_has_a_yellow_and_for_red_air():
    for name in ("Radiation", "Temperature", "Pressure", "Air", "Water"):
        assert hab.suitability_ceiling(_core(**{name: "Yellow"})) == "Yellow", name
    assert hab.suitability_ceiling(_core(Air="Red")) == "Yellow"


def test_suitability_is_green_only_if_the_first_five_are_green_or_blue_and_blue_only_if_all_blue():
    assert hab.suitability_ceiling(_core()) == "Blue"
    assert hab.suitability_ceiling(_core(Pressure="Green", Air="Green")) == "Green"


def test_suitability_is_the_lowest_of_water_food_and_shelter():
    world = World(101.3, 15, {"O2": 0.21}, ph=8.1, life_stage="biosphere")
    assert hab.suitability_score(world) == 1.0
    assert hab.suitability_score(World(101.3, 15, {"O2": 0.21}, ph=8.1, life_stage="macroscopic")) == 0.85
    assert hab.suitability_score(World(101.3, 15, {"O2": 0.21}, ph=8.1, life_stage="microbial")) == 0.5
    assert hab.suitability_score(World(101.3, 15, {"O2": 0.21}, ph=8.1, life_stage="none")) == 0.2
    # water within reach: ice is Yellow, a surface lake is Green when the water is not Blue
    ice = World(101.3, 15, {"O2": 0.21}, ph=8.1, life_stage="biosphere", water_source="ice")
    assert hab.suitability_score(ice) == pytest.approx(0.6)
    brine = World(101.3, 15, {"O2": 0.21}, ph=5.5, life_stage="biosphere")
    assert 0.75 <= hab.suitability_score(brine) < 1.0
    # shelter follows the radiation dose
    irradiated = World(101.3, 15, {"O2": 0.21}, ph=8.1, life_stage="biosphere", surface_dose_sv_yr=0.5)
    assert hab.suitability_score(irradiated) == pytest.approx(hab.radiation_score(irradiated))


def test_life_in_an_ocean_under_ice_is_no_food_on_the_surface():
    buried = World(101.3, 15, {"O2": 0.21}, ph=8.1, life_stage="biosphere", life_on_surface=False)
    assert hab.food_score(buried) == hab.FOOD_SCORE["none"]


def test_a_world_with_no_liquid_water_has_red_water_and_red_suitability_so_yellow_overall():
    dry = World(101.3, 15, {"O2": 0.21}, solvent=None, life_stage="biosphere")
    got = scores(dry)
    assert got["Water"][1] == "Red"
    assert got["Suitability"][1] == "Red"
    assert got["PHI-4"] <= 0.7499 and hab.tier(got["PHI-4"]) == "Yellow"
    assert got["equipment"] == "ideal"  # drinking and eating are not the kit


def test_the_equipment_does_not_depend_on_water_chemistry_or_suitability():
    acid = World(101.3, 15, {"O2": 0.21}, ph=1.5, life_stage="none")
    assert hab.tier(hab.water_score(acid)) in ("Yellow", "Red")
    assert hab.equipment(acid) == "ideal"
