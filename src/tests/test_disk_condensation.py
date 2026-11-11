"""The species-by-species disk condensation (GEN.1363 part A, physics/disk_condensation.py)."""

import math

import pytest

from planetgen.physics import condensate_shares, disk_condensation as dc

SUN_L = 3.828e26


def test_the_suns_budget_per_ten_thousand_hydrogen_is_the_asplund_one():
    assert dc.solar_per_1e4_h("C") == pytest.approx(2.884, rel=1e-3)
    assert dc.solar_per_1e4_h("O") == pytest.approx(4.898, rel=1e-3)
    assert dc.solar_per_1e4_h("N") == pytest.approx(0.676, rel=1e-3)
    assert dc.solar_per_1e4_h("Kr") == pytest.approx(1.32e-5, rel=1e-2) and dc.solar_per_1e4_h("Xe") == pytest.approx(1.66e-6, rel=1e-2)
    with pytest.raises(ValueError):
        dc.solar_per_1e4_h("Unobtainium")


def test_the_carbon_is_split_between_grains_and_the_gases_and_none_is_lost():
    held = dc.inventory()
    moles = held.moles
    assert held.refractory_carbon == pytest.approx(0.5 * held.carbon)
    assert held.refractory_carbon + moles["CO"] + moles["CO2"] + moles["CH4"] == pytest.approx(held.carbon)
    assert moles["CO"] / moles["CO2"] == pytest.approx(5.0)                         # Oberg's CO 1.5 against CO2 0.3
    assert 0.03 < moles["CH4"] / (held.carbon - held.refractory_carbon) < 0.06


def test_the_oxygen_of_the_rock_the_co_and_the_co2_is_out_of_the_water_and_none_is_lost():
    held = dc.inventory()
    moles = held.moles
    assert held.rock_oxygen + moles["CO"] + 2.0 * moles["CO2"] + moles["H2O"] == pytest.approx(held.oxygen)
    short = dc.inventory(oxygen=0.3)                                                # barely more than the rock's share: the water is gone, not negative
    assert short.moles["H2O"] == 0.0 and dc.inventory(oxygen=0.0).moles["H2O"] == 0.0


def test_the_nitrogen_and_sulfur_are_split_and_the_noble_gases_follow_their_ratio():
    held = dc.inventory(nitrogen=2.0, sulfur=0.5, noble=3.0)
    n = dc.solar_per_1e4_h("N") * 2.0
    assert 2.0 * held.moles["N2"] + held.moles["NH3"] + dc.ammonium_salts(nitrogen=2.0) == pytest.approx(n)
    assert held.moles["H2S"] == pytest.approx(0.3 * dc.solar_per_1e4_h("S") * 0.5)
    assert held.moles["Ar"] == pytest.approx(3.0 * dc.solar_per_1e4_h("Ar"))


def test_each_ice_freezes_below_its_line_and_the_lines_run_warm_to_cold():
    lines = list(dc.LINES_K.values())
    assert lines == sorted(lines, reverse=True) and dc.ORDER[0] == "H2O" and dc.ORDER[-1] == "CO"
    assert dc.LINES_K["H2O"] == condensate_shares.ICE_CONDENSATION_K["H2O"]
    assert set(dc.condensed(600.0)) == set()                                         # whole rock, nothing frozen above the grains' line
    assert set(dc.condensed(300.0)) == {"C"}
    assert set(dc.condensed(150.0)) == {"C", "H2O"}
    assert {"NH3", "NH4"} <= set(dc.condensed(100.0)) and "H2S" not in dc.condensed(100.0)
    assert "H2S" in dc.condensed(70.0) and "Xe" not in dc.condensed(70.0)
    for temperature, frozen in ((60.0, "Xe"), (50.0, "Kr"), (47.5, "Ar"), (45.0, "CO2"), (35.0, "CH4"), (21.0, "N2"), (19.0, "CO")):
        assert frozen in dc.condensed(temperature)
    assert "CO" not in dc.condensed(21.0) and "N2" not in dc.condensed(23.0)
    assert set(dc.condensed(10.0)) == set(dc.ORDER) | {"NH4", "C"}


def test_a_species_a_star_lacks_does_not_freeze():
    assert "H2O" not in dc.condensed(100.0, oxygen=0.2)
    assert "CO" not in dc.condensed(10.0, carbon=0.0)


def test_masses_are_percent_of_the_nebula_on_the_water_of_table_11():
    held = dc.inventory(carbon=0.0, oxygen=1.0, nitrogen=0.0, sulfur=0.0, noble=0.0)
    masses = dc.masses_percent(150.0, carbon=0.0, nitrogen=0.0, sulfur=0.0, noble=0.0)
    assert masses["rock"] == condensate_shares.TABLE_11_MASS_PERCENT["rock"]
    assert masses["H2O"] == pytest.approx(condensate_shares.TABLE_11_MASS_PERCENT["H2O"])        # all the oxygen the rock leaves, no CO to take any
    assert held.moles["H2O"] > 3.7
    full = dc.masses_percent(150.0)
    assert full["H2O"] < masses["H2O"]                                                        # the CO and CO2 took oxygen from the water
    assert full["C"] == pytest.approx(dc.inventory().refractory_carbon * 12.011 / 18.015 * dc.PERCENT_PER_H2O_MOLE)


def test_the_ammonia_of_a_solar_disk_is_near_table_11_less_the_nitrogen_that_is_n2():
    ammonia = dc.masses_percent(100.0)["NH3"]
    assert ammonia == pytest.approx(condensate_shares.TABLE_11_MASS_PERCENT["NH3"] * (1.0 - dc.NITROGEN_AS_N2 - dc.NITROGEN_AS_SALTS), rel=0.05)


def test_the_gas_c_over_o_rises_outward_past_the_snow_lines_and_the_solids_stay_oxygen_rich():
    sun_total = dc.solar_per_1e4_h("C") / dc.solar_per_1e4_h("O")
    hot_gas, hot_solid = dc.gas_and_solid_c_over_o(600.0)
    assert hot_gas > sun_total and hot_solid == 0.0                                          # the rock took oxygen and nothing took carbon
    warm_gas, warm_solid = dc.gas_and_solid_c_over_o(150.0)                                  # water frozen: the gas is richer in carbon again
    assert warm_gas > hot_gas > sun_total
    co2_gas, _solid = dc.gas_and_solid_c_over_o(40.0)
    assert co2_gas > warm_gas                                                                # CO2 frozen out too
    for temperature in (600.0, 150.0, 40.0, 25.0):
        gas, solid = dc.gas_and_solid_c_over_o(temperature)
        assert solid is None or solid < 1.0                                                   # solids are oxygen-rich
    assert dc.gas_and_solid_c_over_o(10.0)[0] is None                                        # everything is frozen: no gas left


def test_the_c_over_o_budget_closes_at_every_temperature():
    held = dc.inventory()
    for temperature in (600.0, 300.0, 150.0, 100.0, 60.0, 45.0, 35.0, 25.0, 15.0):
        frozen = dc.condensed(temperature)
        carbon_solid = frozen.get("C", 0.0) + frozen.get("CO", 0.0) + frozen.get("CO2", 0.0) + frozen.get("CH4", 0.0)
        gas, solid = dc.gas_and_solid_c_over_o(temperature)
        oxygen_solid = held.rock_oxygen + frozen.get("H2O", 0.0) + frozen.get("CO", 0.0) + 2.0 * frozen.get("CO2", 0.0)
        if solid is not None:
            assert solid == pytest.approx(carbon_solid / oxygen_solid)
        if gas is not None:
            assert gas == pytest.approx((held.carbon - carbon_solid) / (held.oxygen - oxygen_solid))


def test_the_lines_in_au_move_out_with_the_star_and_are_ordered():
    water, carbon_monoxide = dc.line_au("H2O", SUN_L), dc.line_au("CO", SUN_L)
    assert water == pytest.approx(condensate_shares.ice_line_au(SUN_L, "H2O"), rel=1e-9)
    assert carbon_monoxide > dc.line_au("CO2", SUN_L) > water
    assert dc.line_au("H2O", 4 * SUN_L) == pytest.approx(2.0 * water, rel=1e-9)
    with pytest.raises(ValueError):
        dc.line_au("Ne", SUN_L)


def test_a_negative_ratio_or_temperature_is_refused():
    with pytest.raises(ValueError):
        dc.inventory(carbon=-1.0)
    with pytest.raises(ValueError):
        dc.condensed(0.0)


def test_the_envelope_of_an_unenriched_giant_is_the_star_and_an_enriched_one_takes_the_oxygen_rich_solids():
    star = dc.envelope_c_over_o(100.0, 1.0)
    gas, solid = dc.gas_and_solid_c_over_o(100.0)
    assert star == pytest.approx(dc.inventory().carbon / dc.inventory().oxygen)       # gas and solids in the disk's own proportion: the star's C/O
    rich = dc.envelope_c_over_o(100.0, 60.0)
    assert rich == pytest.approx(solid, rel=0.1) and rich < star
    assert dc.envelope_c_over_o(100.0, 0.05) == pytest.approx(gas)                   # less than the gas's own share of the metals: gas alone


def test_the_envelope_c_over_o_follows_the_star_and_refuses_a_negative_enrichment():
    assert dc.envelope_c_over_o(100.0, 1.0, carbon=1.2) > dc.envelope_c_over_o(100.0, 1.0) > dc.envelope_c_over_o(100.0, 1.0, carbon=0.8)
    with pytest.raises(ValueError):
        dc.envelope_c_over_o(100.0, -1.0)
    with pytest.raises(ValueError):
        dc.envelope_c_over_o(0.0, 1.0)


def test_the_water_share_of_the_solids_is_zero_inside_the_line_and_never_rises_outward():
    assert dc.water_share_of_solids(300.0) == 0.0
    shares = [dc.water_share_of_solids(t) for t in (181.9, 150.0, 100.0, 60.0, 45.0, 30.0, 18.0)]
    assert shares[0] == pytest.approx(0.341, abs=0.002) and all(a >= b for a, b in zip(shares, shares[1:]))
    assert shares[-1] == pytest.approx(0.229, abs=0.003)


def test_the_oxygen_held_by_carbon_monoxide_and_dioxide_leaves_less_water_than_the_equilibrium_table():
    equilibrium = 0.571 / (0.489 + 0.571)                               # [R2726: Table 11], every spare oxygen atom in water
    assert dc.water_share_of_solids(150.0) < equilibrium - 0.1
    assert dc.water_share_of_solids(150.0, oxygen=1.5) > dc.water_share_of_solids(150.0) > dc.water_share_of_solids(150.0, oxygen=0.6)
    assert dc.water_share_of_solids(150.0, carbon=2.0) < dc.water_share_of_solids(150.0)       # a carbon-rich star: more of the oxygen sits in CO and CO2


# --- GEN.1372: the nitrogen, sulfur and noble gases of a giant's envelope, element by element ---------------------------------------------------

def test_a_jupiter_like_envelope_takes_its_nitrogen_and_sulfur_near_three_times_the_star_at_four_times_its_heavy_elements():
    # The Galileo probe: NH3 2.96 +/- 1.13 and H2S 2.75 +/- 0.66 times solar ([R5163: Table 2]), at a formation temperature between the ammonia and hydrogen sulfide lines.
    factors = dc.envelope_enrichments(110.0, 4.0)
    assert factors["N"] == pytest.approx(2.96, abs=0.5)
    assert factors["S"] == pytest.approx(2.75, abs=0.66 + 0.4)
    assert factors["Ar"] == factors["Kr"] == factors["Xe"] == 1.0              # warmer than their lines: gas only, the model's known miss of the probe's 3 times solar


def test_every_factor_is_one_in_a_pure_gas_and_follows_the_ice_lines_as_the_giant_forms_colder():
    assert dc.envelope_enrichments(300.0, 1.0) == pytest.approx({"N": 1.0, "S": 1.0, "Ar": 1.0, "Kr": 1.0, "Xe": 1.0})      # enrichment 1: gas and solids in the disk's proportion
    warm, cold = dc.envelope_enrichments(300.0, 4.0), dc.envelope_enrichments(60.0, 4.0)
    assert warm["N"] == 1.0                                                      # all the nitrogen still gas above the ammonia line
    assert cold["N"] > 2.0 and cold["S"] > warm["S"]
    assert dc.envelope_enrichments(40.0, 4.0)["Ar"] > 1.0 == dc.envelope_enrichments(60.0, 4.0)["Ar"]      # argon freezes below 48 K
    assert dc.envelope_enrichments(40.0, 4.0)["Xe"] == pytest.approx(dc.envelope_enrichments(40.0, 4.0)["Ar"])            # one noble ratio and one sigma


def test_a_barely_enriched_envelope_holds_only_the_gas_and_so_misses_the_frozen_nitrogen():
    low = dc.envelope_enrichments(100.0, 0.2)
    assert low["N"] == pytest.approx(dc.NITROGEN_AS_N2)                                                    # only the N2 (the salts and ammonia are ice, none taken)
    assert low["S"] == pytest.approx(dc.SULFUR_AS_H2S + (1.0 - dc.SULFUR_AS_H2S) * 0.2)                    # the H2S gas plus troilite at the rock's scale


def test_the_stars_nitrogen_sulfur_and_noble_ratios_scale_the_split_not_the_factor_of_a_pure_gas_or_a_pure_ice():
    assert dc.envelope_enrichments(300.0, 4.0, nitrogen=2.0, noble=3.0)["N"] == 1.0
    assert dc.envelope_enrichments(30.0, 4.0, noble=3.0)["Ar"] == pytest.approx(dc.envelope_enrichments(30.0, 4.0)["Ar"])
    with pytest.raises(ValueError):
        dc.envelope_enrichments(0.0, 1.0)
    with pytest.raises(ValueError):
        dc.envelope_enrichments(100.0, -1.0)


def test_the_jump_of_the_solids_at_the_snow_line_is_the_solids_over_the_rock_alone():
    """GEN.1396: rock + grains + water over the rock, read between the water and ammonia lines: 1.97 for the Sun, under the
    equilibrium table's 2.17 (CO and CO2 hold some of the oxygen) and under the classical nebula's 30 / 7."""
    masses = dc.masses_percent(dc.SNOW_LINE_SOLID_JUMP_K)
    assert set(masses) == {"rock", "C", "H2O"}
    assert dc.SNOW_LINE_SOLID_JUMP == pytest.approx((masses["rock"] + masses["C"] + masses["H2O"]) / masses["rock"])
    assert dc.SNOW_LINE_SOLID_JUMP == pytest.approx(1.967, abs=0.001)
    assert dc.SNOW_LINE_SOLID_JUMP < 1.0 + 0.571 / 0.489 < 30.0 / 7.0
    assert dc.SNOW_LINE_SOLID_JUMP > 1.0 + masses["H2O"] / (masses["rock"] + masses["C"])    # the divisor is not the rock plus the grains


def test_a_star_with_more_oxygen_has_a_bigger_jump_and_one_with_no_spare_oxygen_only_the_grains():
    assert dc.snow_line_solid_jump(oxygen=1.5) > dc.SNOW_LINE_SOLID_JUMP > dc.snow_line_solid_jump(oxygen=0.6)
    dry = dc.snow_line_solid_jump(oxygen=condensate_shares.ROCK_OXYGEN_SHARE)
    masses = dc.masses_percent(dc.SNOW_LINE_SOLID_JUMP_K, oxygen=condensate_shares.ROCK_OXYGEN_SHARE)
    assert "H2O" not in masses and dry == pytest.approx((masses["rock"] + masses["C"]) / masses["rock"])


def test_the_mass_law_reads_the_jump_beyond_the_snow_line_only():
    from planetgen.physics import formation
    inside = formation.mmsn_surface_density_gcm2(1.0, 2.0)
    assert formation.mmsn_surface_density_gcm2(1.0, 0.5) == pytest.approx(inside * dc.SNOW_LINE_SOLID_JUMP)


def test_the_recipe_median_is_the_median_of_the_drawn_f7_and_not_the_center_of_its_core():
    """GEN.1471: the core of `f7` is cut at 0.8, so the median of the draw (0.51) is below the core's center (0.57); dividing the draw by the center, as GEN.1363 part C did, lowered every giant."""
    from planetgen.generation.possibility_table import Subject
    from planetgen.physics.atmosphere_recipes import RECIPES
    s = Subject("giant", dict(mass_jupiter=0.3, teq_k=950.0, zone="h"))
    below = RECIPES.probability_in("f7", s, None, dc.RECIPE_C_OVER_O_MEDIAN)
    assert below == pytest.approx(0.5, abs=1e-3), below
    assert dc.RECIPE_C_OVER_O_MEDIAN < dc.RECIPE_C_OVER_O_CENTER - 0.04
    d = RECIPES._dist_of(RECIPES.entry("f7"), s)
    assert d["core"]["sigma"] == pytest.approx(dc.RECIPE_C_OVER_O_CORE_SIGMA_DEX) and d["core"]["range"] == pytest.approx([math.log10(x) for x in dc.RECIPE_C_OVER_O_CORE_RANGE])
    assert (d["tail"]["lo"], d["tail"]["hi"]) == pytest.approx(dc.RECIPE_C_OVER_O_TAIL_RANGE)


def test_a_giant_with_the_medians_f7_keeps_its_disks_c_over_o_and_a_lightly_mixed_envelope_is_the_gas_leaning_one():
    from planetgen.physics import disk_drift, mixture_giant
    subject = {"formation_temperature_k": 100.0, "disk_carbon": 1.0, "disk_oxygen": 1.0}
    bulk = mixture_giant._disk_c_over_o(subject, math.log10(18.0))             # 18 times the star's heavy elements in the planet
    assert bulk == pytest.approx(dc.envelope_c_over_o(100.0, dc.ENVELOPE_METAL_SHARE * 18.0, solids_edge_k=None))      # GEN.1520: a subject without a snow line and disk edge places no solids' edge
    assert bulk > dc.envelope_c_over_o(100.0, 18.0) + 0.04                     # the envelope carries a share of the bulk, so less of the oxygen-rich solids


# --- GEN.1516: the ice weights of the pebble drift reach the nitrogen, sulfur and noble gases -------------------------------------------------

ALL_LINES = tuple(dc.LINES_K)


def test_weights_of_one_give_the_static_envelope_factors_at_every_temperature():
    for temperature in (30.0, 60.0, 110.0, 150.0, 300.0, 700.0):
        weights = {species: 1.0 for species in ALL_LINES if species in dc.LINES_K and dc.LINES_K[species] <= temperature}
        for enrichment in (0.2, 1.0, 4.0):
            assert dc.envelope_enrichments(temperature, enrichment, ice_weights=weights) == pytest.approx(dc.envelope_enrichments(temperature, enrichment), rel=1e-12)
    assert dc.envelope_enrichments(110.0, 4.0, ice_weights={}) == dc.envelope_enrichments(110.0, 4.0)


def test_a_noble_gas_ice_that_the_pebbles_carried_in_is_the_gas_times_its_weight_and_one_that_was_not_is_absent():
    got = dc.envelope_enrichments(300.0, 4.0, ice_weights={"Ar": 3.0, "Kr": 0.0, "Xe": 1.5})
    assert got["Ar"] == pytest.approx(3.0) and got["Kr"] == 0.0 and got["Xe"] == pytest.approx(1.5)


def test_the_nitrogen_and_the_sulfur_take_the_weights_of_their_ices():
    n2 = 2.0 * dc.inventory().moles["N2"]
    nh3 = dc.inventory().moles["NH3"]
    total = n2 + nh3 + dc.ammonium_salts()
    got = dc.envelope_enrichments(300.0, 4.0, ice_weights={"NH3": 3.0, "N2": 2.0, "H2S": 2.0})
    assert got["N"] == pytest.approx((total + 2.0 * nh3 + n2) / total)
    sulfur = dc.solar_per_1e4_h("S")
    h2s = sulfur * dc.SULFUR_AS_H2S
    troilite = sulfur * (1.0 - dc.SULFUR_AS_H2S)
    assert got["S"] == pytest.approx((sulfur - troilite + h2s + 4.0 * troilite) / sulfur)


def test_a_weight_on_an_ice_the_giant_formed_colder_than_the_line_of_is_ignored_because_that_ice_is_frozen_there():
    assert dc.envelope_enrichments(40.0, 4.0, ice_weights={"H2O": 2.0, "NH3": 0.0, "Ar": 0.0}) == dc.envelope_enrichments(40.0, 4.0, ice_weights={"H2O": 1.0, "NH3": 1.0, "Ar": 1.0})
    assert dc.envelope_enrichments(40.0, 4.0, ice_weights={"H2O": 2.0})["Ar"] > 1.0                    # argon still freezes at 40 K: the solids' share, whatever the water's weight


# --- GEN.1519: the planetesimals a giant accretes formed between its place and the disk's edge ----------------------------------------------

EDGE_K = dc.LINES_K["H2O"] / 18.0 ** 0.5


def test_without_a_provenance_or_at_the_edge_the_solids_are_the_ones_that_froze_at_the_place():
    assert dc.provenance_condensed(30.0, EDGE_K) == dc.condensed(30.0)                 # at or beyond the edge: nothing farther out
    for temperature in (30.0, 110.0):
        assert dc.envelope_enrichments(temperature, 4.0, solids_edge_k=None) == dc.envelope_enrichments(temperature, 4.0)
    assert dc.envelope_enrichments(30.0, 4.0, solids_edge_k=EDGE_K) == pytest.approx(dc.envelope_enrichments(30.0, 4.0), rel=1e-12)
    assert dc.envelope_c_over_o(30.0, 4.0, solids_edge_k=EDGE_K) == pytest.approx(dc.envelope_c_over_o(30.0, 4.0), rel=1e-12)


def test_the_share_of_the_annuli_that_froze_an_ice_follows_the_distance_of_its_line():
    temperature = 300.0
    mixed = dc.provenance_condensed(temperature, EDGE_K)
    held = dc.inventory()
    edge_ratio = (temperature / EDGE_K) ** 2                                            # r_edge / r_place
    for species in ("Ar", "Kr", "Xe"):
        line_ratio = (temperature / dc.LINES_K[species]) ** 2                           # r_line / r_place
        share = (edge_ratio - line_ratio) / (edge_ratio - 1.0)                          # the annuli beyond the line, equal per AU
        assert mixed[species] == pytest.approx(share * held.moles[species], abs=2.0 / dc.PROVENANCE_SAMPLES * held.moles[species])
    assert mixed["Xe"] / held.moles["Xe"] > mixed["Kr"] / held.moles["Kr"] > mixed["Ar"] / held.moles["Ar"] > 0.0      # the warmer the line, the larger the share of the disk beyond it
    assert "CO" not in mixed and "N2" not in mixed                                     # colder than the edge: not solid anywhere in the disk


def test_a_giant_that_formed_warm_takes_its_noble_gases_from_the_planetesimals_that_formed_beyond_their_lines():
    warm = dc.envelope_enrichments(110.0, 4.0, solids_edge_k=EDGE_K)
    local = dc.envelope_enrichments(110.0, 4.0)
    assert local["Ar"] == local["Kr"] == local["Xe"] == 1.0                              # no solid noble gas at its own place
    assert 1.0 < warm["Ar"] < warm["Kr"] < warm["Xe"]
    assert warm["N"] == pytest.approx(local["N"], rel=0.1)                              # the nitrogen the probe measured is unchanged within a tenth
    assert all(math.isfinite(value) and value > 0.0 for value in warm.values())


def test_the_provenance_needs_positive_temperatures():
    with pytest.raises(ValueError):
        dc.provenance_condensed(0.0, EDGE_K)
    with pytest.raises(ValueError):
        dc.provenance_condensed(100.0, 0.0)
