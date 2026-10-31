# tests/test_roundness.py

"""
GEN.1158: the roundness engine and the binding energy of every kind of mass object. The bodies below are test targets (radii and densities from Lineweaver and Norman 2010 and
Probsthain 2018, Proteus from the Neptunian-moon literature); the module takes none of them as input.
"""

import copy
import json
import math

import pytest

from planetgen.physics import binding_energy as be
from planetgen.physics import constants, roundness as rn

G = constants.G
EARTH = constants.EARTH_MASS_TO_KG
SUN = constants.SOLAR_MASS_TO_KG
C = constants.SPEED_OF_LIGHT_M_S


def test_the_round_radius_is_exactly_the_calibration_at_the_reference_densities():
    materials = _document()["materials"]
    assert 150e3 < rn.round_radius_m(1000.0, 1.0) < 260e3          # ice rounds near Lineweaver and Norman's 200 km
    assert rn.round_radius_m(materials["rock"]["reference_density_kg_m3"], 0.0) == pytest.approx(materials["rock"]["round_radius_km"] * 1e3)


def test_the_strengths_are_a_few_megapascals_and_rock_is_stronger_than_ice():
    assert 0.4e6 < rn.ICE_STRENGTH_PA < 3.0e6             # Lineweaver and Norman's range for ice
    assert 1.5e7 < rn.ROCK_STRENGTH_PA < 1.8e7
    assert 5.0e7 <= rn.METAL_STRENGTH_PA <= 2.5e8


def test_lineweaver_and_normans_ten_megapascal_example_rounds_near_their_240_km():
    sigma, rho = 10.0e6, 2500.0
    radius = math.sqrt(sigma / ((2 * math.pi / 3) * (rn.SPHERICITY_FACTOR ** 2 - 1) * G * rho ** 2))
    assert radius / 1e3 == pytest.approx(233.0, abs=3.0)


def test_the_round_radius_goes_as_one_over_the_density_at_one_strength():
    assert rn.round_radius_m(2000.0, 0.0) == pytest.approx(2.0 * rn.round_radius_m(4000.0, 0.0))


def _ice_fraction_with_solid_volume_share(share):
    """The ice mass fraction of a rock and ice body whose solids are this share of the volume."""
    ice, rock = rn.GRAIN_DENSITIES_KG_M3["ice"], rn.GRAIN_DENSITIES_KG_M3["rock"]
    solid_mass = share * rock
    return (1.0 - share) * ice / ((1.0 - share) * ice + solid_mass)


def test_strength_runs_from_ice_to_rock_by_the_volume_share_of_the_solids():
    rule = rn.TABLE.mixing
    assert rn.yield_strength_pa(0.0) == pytest.approx(rn.ROCK_STRENGTH_PA)
    assert rn.yield_strength_pa(1.0) == pytest.approx(rn.ICE_STRENGTH_PA)
    for share in (0.1, 0.3, rule["tested_to"]):          # the measured law: ice strength x exp(2 f)
        assert rn.yield_strength_pa(_ice_fraction_with_solid_volume_share(share)) == pytest.approx(rn.ICE_STRENGTH_PA * math.exp(rule["exponent"] * share))
    assert rn.yield_strength_pa(_ice_fraction_with_solid_volume_share(rule["close_packing"])) == pytest.approx(rn.ICE_STRENGTH_PA * math.exp(rule["exponent"] * rule["close_packing"]))
    values = [rn.yield_strength_pa(x / 20.0) for x in range(21)]
    assert values == sorted(values, reverse=True)


def test_half_the_volume_in_rock_is_about_2_7_times_ice_not_the_8_5_of_a_straight_blend():
    assert rn.yield_strength_pa(_ice_fraction_with_solid_volume_share(0.5)) / rn.ICE_STRENGTH_PA == pytest.approx(2.72, rel=0.01)


def test_the_mixing_rule_uses_volume_not_mass():
    # 50 percent rock by mass is about 25 percent by volume: a factor exp(0.5) of 1.6 over ice
    ratio = rn.yield_strength_pa(0.5) / rn.ICE_STRENGTH_PA
    assert 1.5 < ratio < 1.8


def test_the_creep_law_of_ice_follows_the_research_table():
    # strength research section 11a, rate 1e-15 per second (MPa)
    for kelvin, mpa in ((100.0, 34.0), (120.0, 9.2), (140.0, 3.6), (160.0, 0.90), (200.0, 0.091), (250.0, 0.011)):
        assert rn.yield_strength_pa(1.0, temperature_k=kelvin) / 1e6 == pytest.approx(mpa, rel=0.05)


def test_the_default_ice_is_the_creep_law_near_158_kelvin():
    assert rn.yield_strength_pa(1.0, temperature_k=158.0) == pytest.approx(rn.ICE_STRENGTH_PA, rel=0.1)
    assert rn.ICE_STRENGTH_PA == pytest.approx(1.0e6)


def test_colder_ice_is_stronger_and_warmer_ice_weaker():
    values = [rn.yield_strength_pa(1.0, temperature_k=kelvin) for kelvin in (80.0, 100.0, 130.0, 160.0, 200.0, 240.0)]
    assert values == sorted(values, reverse=True)
    assert rn.round_radius_m(1000.0, 1.0, temperature_k=100.0) > rn.round_radius_m(1000.0, 1.0) > rn.round_radius_m(1000.0, 1.0, temperature_k=220.0)
    assert rn.yield_strength_pa(0.0, temperature_k=100.0) == rn.yield_strength_pa(0.0)        # rock's strength does not follow the ice law


def test_ice_above_the_melting_point_is_held_at_it():
    assert rn.yield_strength_pa(1.0, temperature_k=400.0) == rn.yield_strength_pa(1.0, temperature_k=constants.CELSIUS_ZERO_K)


def test_cold_ice_in_a_mix_is_stronger_than_cold_ice_alone_until_the_solids_jam():
    assert rn.yield_strength_pa(0.5, temperature_k=140.0) > rn.yield_strength_pa(1.0, temperature_k=140.0)


def test_a_temperature_at_or_below_zero_kelvin_is_refused():
    with pytest.raises(ValueError):
        rn.yield_strength_pa(1.0, temperature_k=0.0)


# name, mean radius km, density kg/m3, ice fraction, is round
BODIES = [
    ("Ceres", 469.7, 2161.0, 0.18, True),
    ("Vesta", 258.6, 3456.0, 0.0, True),
    ("Pallas", 272.5, 2761.0, 0.0, True),
    ("Juno", 115.0, 2700.0, 0.0, False),
    ("Enceladus", 252.1, 1609.0, 0.39, True),
    ("Miranda", 235.0, 1200.0, 0.67, True),
    ("Mimas", 198.2, 1149.0, 0.72, True),
    ("Proteus", 210.0, 700.0, 0.7, False),          # porous: the low density is why it is not round at Mimas' radius
    ("Hyperion", 140.0, 540.0, 0.5, False),
]


@pytest.mark.parametrize("name, radius_km, density, ice, round_", BODIES)
def test_the_known_bodies_fall_on_the_right_side(name, radius_km, density, ice, round_):
    ratio = rn.roundness_ratio(radius_km * 1e3, density, ice)
    assert rn.is_round(ratio) == round_, (name, ratio)


def test_the_ratio_says_how_far_from_round_a_body_is():
    mimas = rn.roundness_ratio(198.2e3, 1149.0, 0.72)
    hyperion = rn.roundness_ratio(140e3, 540.0, 0.5)
    ceres = rn.roundness_ratio(469.7e3, 2161.0, 0.18)
    assert hyperion < 0.5 < 1.0 <= mimas < 1.5 < ceres


def test_the_ratio_from_mass_is_the_ratio_of_the_sphere_they_make():
    mass = 4.0 / 3.0 * math.pi * (198.2e3) ** 3 * 1149.0
    assert rn.roundness_ratio_from_mass(mass, 1149.0, 0.72) == pytest.approx(rn.roundness_ratio(198.2e3, 1149.0, 0.72))
    assert rn.bulk_density_kg_m3(mass, 198.2e3) == pytest.approx(1149.0)


def test_the_round_mass_is_a_sphere_of_the_round_radius():
    for density, ice in ((1000.0, 1.0), (2500.0, 0.0), (1600.0, 0.4)):
        assert rn.roundness_ratio_from_mass(rn.round_mass_kg(density, ice), density, ice) == pytest.approx(1.0)


def test_ice_turns_round_at_a_smaller_mass_than_rock_and_the_bands_agree_with_the_report():
    ice = rn.round_mass_kg(1000.0, 1.0)
    rock = rn.round_mass_kg(2500.0, 0.0)
    assert 2e19 < ice < 4.5e19         # about Mimas (3.7e19 kg), which is round
    assert 2e20 < rock < 5e20


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -1.0, 0.0])
def test_a_bad_density_or_radius_is_refused(bad):
    with pytest.raises(ValueError):
        rn.round_radius_m(bad, 0.5)
    with pytest.raises(ValueError):
        rn.roundness_ratio(bad, 1000.0, 0.5)


def test_an_ice_fraction_outside_zero_to_one_is_refused():
    for bad in (-0.1, 1.1, float("nan")):
        with pytest.raises(ValueError):
            rn.yield_strength_pa(bad)


# --- binding energy ---------------------------------------------------------------------------------------------------------------------

def test_a_uniform_sphere_is_three_fifths_of_gm2_over_r():
    assert be.uniform_j(EARTH, 6.371e6) == pytest.approx(0.6 * G * EARTH ** 2 / 6.371e6)
    assert be.binding_energy_j("rocky", EARTH, 6.371e6) == be.uniform_j(EARTH, 6.371e6)


def test_earth_with_a_core_binds_more_than_a_uniform_one():
    uniform = be.uniform_j(EARTH, 6.371e6)
    differentiated = be.binding_energy_j("rocky", EARTH, 6.371e6, core_mass_kg=0.32 * EARTH, core_radius_m=3.48e6)
    assert 1.05 < differentiated / uniform < 1.25               # the real Earth is about 1.11 times the uniform value
    assert be.two_shell_j(EARTH, 6.371e6, 0.0, 0.0) == pytest.approx(uniform)


def test_polytropes_follow_three_over_five_minus_n():
    assert be.polytrope_j(SUN, 6.957e8, 0.0) == pytest.approx(0.6 * G * SUN ** 2 / 6.957e8)
    assert be.polytrope_j(SUN, 6.957e8, 3.0) == pytest.approx(1.5 * G * SUN ** 2 / 6.957e8)
    assert be.polytrope_j(SUN, 6.957e8, 1.5) == pytest.approx(6.0 / 7.0 * G * SUN ** 2 / 6.957e8)
    with pytest.raises(ValueError):
        be.polytrope_j(SUN, 6.957e8, 5.0)


def test_the_sun_is_a_few_1e41_joules():
    sun = be.binding_energy_j("star", SUN, 6.957e8)
    assert 4e41 < sun < 8e41


def test_the_star_index_switches_at_the_fully_convective_mass():
    assert be.star_index(0.2) == be.POLYTROPE_INDEX_CONVECTIVE
    assert be.star_index(0.5) == be.POLYTROPE_INDEX_RADIATIVE
    low = be.binding_energy_j("star", 0.2 * SUN, 0.23 * 6.957e8)
    assert low == be.polytrope_j(0.2 * SUN, 0.23 * 6.957e8, 1.5)


def test_jupiter_as_an_n_one_polytrope_is_three_quarters_of_gm2_over_r():
    jupiter = be.binding_energy_j("gas-giant", 1.898e27, 6.99e7)
    assert jupiter == pytest.approx(0.75 * G * 1.898e27 ** 2 / 6.99e7)


def test_brown_dwarfs_and_white_dwarfs_use_the_three_halves_polytrope():
    assert be.binding_energy_j("brown-dwarf", 0.05 * SUN, 7.0e7) == be.polytrope_j(0.05 * SUN, 7.0e7, 1.5)
    assert be.binding_energy_j("white-dwarf", 0.6 * SUN, 8.7e6) == be.polytrope_j(0.6 * SUN, 8.7e6, 1.5)


def test_a_neutron_star_is_lattimer_and_prakashs_fit():
    mass, radius = 1.4 * SUN, 12.0e3
    beta = G * mass / (radius * C ** 2)
    expected = mass * C ** 2 * 0.6 * beta / (1 - 0.5 * beta)
    assert be.binding_energy_j("neutron-star", mass, radius) == pytest.approx(expected)
    assert 0.1 < be.binding_energy_j("neutron-star", mass, radius) / (mass * C ** 2) < 0.2          # about 0.12 to 0.17 of the mass-energy
    with pytest.raises(ValueError):
        be.neutron_star_j(mass, 1.0)


def test_a_black_hole_stores_minus_one_meaning_infinite():
    assert be.binding_energy_j("black-hole", 10 * SUN) == be.INFINITE == -1.0
    assert be.is_infinite(be.INFINITE) and not be.is_infinite(1.0e30)
    assert be.binding_energy_j("black-hole", 10 * SUN, 3.0e4) == -1.0


def test_every_real_binding_energy_is_positive_so_minus_one_is_never_one():
    values = [be.uniform_j(1e19, 1e5), be.polytrope_j(SUN, 7e8, 3.0), be.neutron_star_j(1.4 * SUN, 12e3), be.two_shell_j(EARTH, 6.4e6, 0.3 * EARTH, 3.5e6)]
    assert all(v > 0.0 for v in values)


def test_an_unknown_kind_or_a_missing_radius_is_refused():
    with pytest.raises(ValueError):
        be.binding_energy_j("galaxy", SUN, 1e9)
    with pytest.raises(ValueError):
        be.binding_energy_j("star", SUN)
    with pytest.raises(ValueError):
        be.uniform_j(0.0, 1.0)


# --- the strength table: where real data drops in ---------------------------------------------------------------------------------------

def _document():
    with open(rn.DATA_FILE, encoding="utf-8") as handle:
        return json.load(handle)


def test_the_file_loads_to_the_table_in_use():
    table = rn.load_table(_document())
    assert table.strengths_pa == rn.TABLE.strengths_pa
    assert table.strengths_pa["metal"] == pytest.approx(1.0e8)                     # the research's weak-source stand-in for metal
    assert table.strengths_pa["ice"] < table.strengths_pa["rock"] < table.strengths_pa["metal"]


def test_a_material_can_take_the_strength_of_another_until_it_has_its_own():
    document = copy.deepcopy(_document())
    document["materials"]["metal"] = {"same_as": "rock"}
    assert rn.load_table(document).strengths_pa["metal"] == rn.load_table(document).strengths_pa["rock"]


def test_a_measured_strength_in_the_file_replaces_a_stand_in_without_changing_the_interface(monkeypatch):
    document = copy.deepcopy(_document())
    document["materials"]["rock"] = {"strength_pa": 4.0e7}
    document["materials"]["metal"] = {"strength_pa": 2.0e8}
    monkeypatch.setattr(rn, "TABLE", rn.load_table(document))
    assert rn.yield_strength_pa(0.0) == pytest.approx(4.0e7)
    assert rn.yield_strength_pa(0.0, metal_fraction=1.0) == pytest.approx(2.0e8)
    grains = rn.GRAIN_DENSITIES_KG_M3
    metal_share = (0.5 / grains["metal"]) / (0.5 / grains["metal"] + 0.5 / grains["rock"])
    assert rn.yield_strength_pa(0.0, metal_fraction=0.5) == pytest.approx(4.0e7 ** (1.0 - metal_share) * 2.0e8 ** metal_share)
    stronger = rn.round_radius_m(2500.0, 0.0)
    monkeypatch.undo()
    assert stronger > rn.round_radius_m(2500.0, 0.0)          # a stronger rock needs a larger body to round


def test_metal_with_its_own_strength_and_density_rounds_where_it_should(monkeypatch):
    document = copy.deepcopy(_document())
    document["materials"]["metal"] = {"strength_pa": 1.0e8}
    monkeypatch.setattr(rn, "TABLE", rn.load_table(document))
    iron = rn.round_radius_m(7800.0, 0.0, metal_fraction=1.0)
    assert iron == pytest.approx(math.sqrt(1.0e8 / (rn.stress_coefficient(rn.TABLE.sphericity_factor) * G * 7800.0 ** 2)))


def test_the_exponential_porosity_rule_weakens_a_porous_body(monkeypatch):
    document = copy.deepcopy(_document())
    document["porosity"] = {"rule": "exponential", "coefficient": 5.0}
    monkeypatch.setattr(rn, "TABLE", rn.load_table(document))
    assert rn.yield_strength_pa(0.0, porosity=0.0) == pytest.approx(rn.ROCK_STRENGTH_PA)
    assert rn.yield_strength_pa(0.0, porosity=0.5) == pytest.approx(rn.ROCK_STRENGTH_PA * math.exp(-2.5))


def test_the_stand_in_porosity_rule_leaves_the_strength_alone():
    assert rn.yield_strength_pa(0.3, porosity=0.4) == rn.yield_strength_pa(0.3)


@pytest.mark.parametrize("change", [
    lambda d: d.pop("sphericity_factor"),
    lambda d: d.__setitem__("sphericity_factor", 0.9),
    lambda d: d["materials"].pop("ice"),
    lambda d: d["materials"]["rock"].pop("round_radius_km"),
    lambda d: d["materials"]["metal"].pop("strength_pa"),
    lambda d: d["materials"].__setitem__("metal", {"same_as": "plasma"}),
    lambda d: d["mixing"].__setitem__("rule", "voigt"),
    lambda d: d["porosity"].__setitem__("rule", "sponge"),
    lambda d: d["temperature"]["ice"].__setitem__("rule", "balmy"),
    lambda d: d["temperature"].__setitem__("plasma", {"rule": "none"}),
    lambda d: d["temperature"]["rock"]["columns"].pop(),
    lambda d: d["temperature"]["rock"].pop("rate_per_s"),
    lambda d: d["temperature"]["rock"]["columns"][0].pop("e_kj_mol"),
    lambda d: d["mixing"].__setitem__("close_packing", 0.3),
    lambda d: d["temperature"]["ice"].__setitem__("regimes", []),
    lambda d: d["temperature"]["ice"]["regimes"][0].pop("n"),
])
def test_a_table_that_is_wrong_is_refused(change):
    document = copy.deepcopy(_document())
    change(document)
    with pytest.raises(rn.RoundnessError):
        rn.load_table(document)


def test_the_metal_and_ice_fractions_cannot_exceed_the_whole_body():
    with pytest.raises(ValueError):
        rn.yield_strength_pa(0.7, metal_fraction=0.6)


# --- the rubble-pile branch (GEN.1176) --------------------------------------------------------------------------------------------------

def _with_rubble(monkeypatch, **changes):
    document = copy.deepcopy(_document())
    document["rubble"].update(changes)
    monkeypatch.setattr(rn, "TABLE", rn.load_table(document))


def test_the_rubble_numbers_are_the_researchs():
    rule = rn.TABLE.rubble
    assert rule["friction_angle_deg"] == 35.0
    assert 44.0 <= rule["cohesion_pa"] <= 210.0
    assert 1.0e6 <= rule["hot_strength_pa"] <= 4.0e6


def test_the_failure_stress_is_mohr_coulomb_with_cohesion_and_the_centre_pressure():
    phi, cohesion = math.radians(35.0), rn.TABLE.rubble["cohesion_pa"]
    radius, density = 50.0e3, 2000.0
    pressure = (2.0 * math.pi / 3.0) * G * density ** 2 * radius ** 2
    expected = (2.0 * cohesion * math.cos(phi) + 2.0 * pressure * math.sin(phi)) / (1.0 - math.sin(phi))
    assert rn.rubble_strength_pa(radius, density) == pytest.approx(expected)


def test_a_cold_rubble_pile_is_never_round_at_any_size():
    ratios = [rn.roundness_ratio(radius, 2500.0, 0.0, structure="rubble") for radius in (1.0e2, 1.0e4, 3.0e5, 1.0e7)]
    assert max(ratios) < 0.3
    assert all(not rn.is_round(ratio) for ratio in ratios)
    assert rn.roundness_ratio(1.0e7, 2500.0, 0.0, structure="rubble") == pytest.approx(rn.roundness_ratio(1.0e6, 2500.0, 0.0, structure="rubble"), rel=1e-3)   # size-free at large size


def test_a_cold_rubble_pile_of_ceres_size_is_not_round_but_the_same_monolith_is():
    assert not rn.is_round(rn.roundness_ratio(469.7e3, 2161.0, 0.18, structure="rubble"))
    assert rn.is_round(rn.roundness_ratio(469.7e3, 2161.0, 0.18))


def test_a_hot_rubble_pile_rounds_at_the_hot_strength():
    radius = rn.rubble_round_radius_m(2500.0, heated=True)
    assert radius == pytest.approx(math.sqrt(rn.TABLE.rubble["hot_strength_pa"] / (rn.stress_coefficient(rn.TABLE.sphericity_factor) * G * 2500.0 ** 2)))
    assert rn.is_round(rn.roundness_ratio(radius * 1.01, 2500.0, 0.0, structure="rubble", heated=True))
    assert not rn.is_round(rn.roundness_ratio(radius * 0.99, 2500.0, 0.0, structure="rubble", heated=True))
    assert radius < rn.round_radius_m(2500.0, 0.0)                         # weaker than the cold monolith, so it rounds at a smaller size


def test_an_icy_rubble_pile_that_took_its_shape_warm_rounds_like_the_ice_creeps():
    warm = rn.rubble_round_radius_m(1000.0, 1.0, temperature_k=200.0)
    assert warm == pytest.approx(rn.round_radius_m(1000.0, 1.0, temperature_k=200.0))
    assert rn.rubble_round_radius_m(1000.0, 1.0, temperature_k=100.0) > warm      # colder ice is stronger, so a larger body is needed
    with pytest.raises(ValueError):
        rn.rubble_round_radius_m(1000.0, 1.0)                                   # no shape-set temperature: friction only


def test_a_cold_rubble_pile_has_no_round_radius():
    with pytest.raises(ValueError):
        rn.rubble_round_radius_m(2500.0)


def test_rock_with_a_shape_set_temperature_and_no_ice_is_not_made_weak_by_it():
    assert rn.rubble_strength_pa(100.0e3, 2500.0, 0.0, temperature_k=150.0) == rn.rubble_strength_pa(100.0e3, 2500.0)


def test_a_shallow_enough_friction_angle_does_let_a_pile_round_by_size(monkeypatch):
    _with_rubble(monkeypatch, friction_angle_deg=3.0)
    assert rn.is_round(rn.roundness_ratio(1.0e6, 2500.0, 0.0, structure="rubble"))
    assert math.isfinite(rn.rubble_round_radius_m(2500.0))


def test_the_ratio_from_mass_takes_the_structure():
    mass = 4.0 / 3.0 * math.pi * (100.0e3) ** 3 * 2500.0
    assert rn.roundness_ratio_from_mass(mass, 2500.0, 0.0, structure="rubble", heated=True) == pytest.approx(
        rn.roundness_ratio(100.0e3, 2500.0, 0.0, structure="rubble", heated=True))


def test_an_unknown_structure_or_a_bad_rubble_input_is_refused():
    with pytest.raises(ValueError):
        rn.roundness_ratio(1.0e5, 2500.0, 0.0, structure="sponge")
    for bad in (0.0, -1.0, float("nan")):
        with pytest.raises(ValueError):
            rn.rubble_strength_pa(bad, 2500.0)
    with pytest.raises(ValueError):
        rn.rubble_strength_pa(1.0e5, 2500.0, ice_fraction=1.5)


@pytest.mark.parametrize("change", [
    lambda d: d.pop("rubble"),
    lambda d: d["rubble"].__setitem__("rule", "sand"),
    lambda d: d["rubble"].__setitem__("friction_angle_deg", 95.0),
    lambda d: d["rubble"].__setitem__("friction_angle_deg", 0.0),
    lambda d: d["rubble"].pop("cohesion_pa"),
    lambda d: d["rubble"].pop("hot_strength_pa"),
])
def test_a_wrong_rubble_table_is_refused(change):
    document = copy.deepcopy(_document())
    change(document)
    with pytest.raises(rn.RoundnessError):
        rn.load_table(document)


# --- the strength draw (GEN.1178) -------------------------------------------------------------------------------------------------------

def test_the_median_draw_is_the_stand_in_and_none_is_the_same():
    assert rn.drawn_strengths_pa(0.5) == pytest.approx(rn.TABLE.strengths_pa)
    for call in (lambda draw: rn.yield_strength_pa(0.4, draw=draw), lambda draw: rn.round_radius_m(2000.0, 0.4, draw=draw),
                 lambda draw: rn.roundness_ratio(2.0e5, 2000.0, 0.4, draw=draw), lambda draw: rn.rubble_strength_pa(2.0e5, 2000.0, heated=True, draw=draw)):
        assert call(0.5) == pytest.approx(call(None))


def test_the_ends_of_the_draw_are_the_ends_of_each_range():
    low, high = rn.drawn_strengths_pa(0.0), rn.drawn_strengths_pa(1.0)
    for name, (lo, hi) in rn.RANGES_PA.items():
        assert low[name] == pytest.approx(lo) and high[name] == pytest.approx(hi)
    assert (low["ice"], low["rock"], low["metal"]) == pytest.approx((0.4e6, 4.0e6, 50.0e6))
    assert (high["ice"], high["rock"], high["metal"]) == pytest.approx((3.0e6, 30.0e6, 250.0e6))


def test_a_higher_draw_is_a_stronger_body_and_a_bigger_round_radius():
    draws = [x / 10.0 for x in range(11)]
    for values in ([rn.yield_strength_pa(0.0, draw=d) for d in draws], [rn.yield_strength_pa(1.0, draw=d) for d in draws],
                   [rn.round_radius_m(2500.0, 0.0, draw=d) for d in draws], [rn.rubble_round_radius_m(2500.0, heated=True, draw=d) for d in draws]):
        assert values == sorted(values) and values[0] < values[-1]


def test_the_round_radius_moves_about_one_and_a_half_either_way():
    low, mid, high = (rn.round_radius_m(2500.0, 0.0, draw=d) for d in (0.0, 0.5, 1.0))
    assert 1.3 < mid / low < 2.2 and 1.0 < high / mid < 1.8


def test_a_drawn_body_is_round_on_one_side_of_its_own_boundary():
    boundary = rn.round_radius_m(2500.0, 0.0, draw=0.2)
    assert rn.is_round(rn.roundness_ratio(boundary * 1.01, 2500.0, 0.0, draw=0.2))
    assert not rn.is_round(rn.roundness_ratio(boundary * 0.99, 2500.0, 0.0, draw=0.2))
    assert not rn.is_round(rn.roundness_ratio(boundary * 1.01, 2500.0, 0.0, draw=0.9))      # a tougher body of that size is still a potato


def test_the_draw_scales_the_ice_creep_law_and_the_rubble_numbers():
    assert rn.yield_strength_pa(1.0, temperature_k=140.0, draw=1.0) / rn.yield_strength_pa(1.0, temperature_k=140.0) == pytest.approx(3.0)
    assert rn.yield_strength_pa(1.0, temperature_k=140.0, draw=0.0) / rn.yield_strength_pa(1.0, temperature_k=140.0) == pytest.approx(0.4)
    assert rn.rubble_strength_pa(1.0e4, 2500.0, draw=1.0) > rn.rubble_strength_pa(1.0e4, 2500.0, draw=0.0)
    hot = rn.TABLE.rubble
    assert rn.rubble_strength_pa(1.0e6, 2500.0, heated=True, draw=1.0) == pytest.approx(hot["hot_range_pa"][1])


def test_a_cold_rubble_pile_stays_unround_at_every_draw():
    for draw in (0.0, 0.5, 1.0):
        assert rn.roundness_ratio(1.0e6, 2500.0, 0.0, structure="rubble", draw=draw) < 0.5


def test_the_draw_does_not_change_the_table_in_use():
    before = dict(rn.TABLE.strengths_pa), dict(rn.TABLE.rubble)
    rn.drawn_strengths_pa(0.9)
    rn.rubble_strength_pa(1.0e5, 2500.0, draw=0.9)
    assert (dict(rn.TABLE.strengths_pa), dict(rn.TABLE.rubble)) == before


@pytest.mark.parametrize("bad", [-0.1, 1.1, float("nan"), float("inf"), "high", True])
def test_a_draw_outside_zero_to_one_is_refused(bad):
    with pytest.raises(ValueError):
        rn.drawn_strengths_pa(bad)
    with pytest.raises(ValueError):
        rn.yield_strength_pa(0.5, draw=bad)


def test_a_material_with_no_range_does_not_move(monkeypatch):
    document = copy.deepcopy(_document())
    del document["materials"]["metal"]["range_pa"]
    monkeypatch.setattr(rn, "TABLE", rn.load_table(document))
    monkeypatch.setattr(rn, "RANGES_PA", rn.TABLE.ranges_pa)       # the draw reads the module's ranges (GEN.1178)
    assert rn.drawn_strengths_pa(0.0)["metal"] == rn.drawn_strengths_pa(1.0)["metal"] == rn.TABLE.strengths_pa["metal"]


@pytest.mark.parametrize("change", [
    lambda d: d["materials"]["rock"].__setitem__("range_pa", [20.0e6, 30.0e6]),         # the strength is outside its range
    lambda d: d["materials"]["rock"].__setitem__("range_pa", [4.0e6]),
    lambda d: d["rubble"].__setitem__("hot_range_pa", [4.0e6, 1.0e6]),
])
def test_a_wrong_range_is_refused(change):
    document = copy.deepcopy(_document())
    change(document)
    with pytest.raises(rn.RoundnessError):
        rn.load_table(document)


# --- the rock creep law (GEN.1184) -------------------------------------------------------------------------------------------------------

def _rock(kelvin, draw=None):
    return rn.yield_strength_pa(0.0, temperature_k=kelvin, draw=draw) / 1e6


def test_cold_rock_keeps_its_whole_body_strength_up_to_the_creep_line():
    for kelvin in (100.0, 300.0, 900.0, 1000.0, 1100.0):
        assert _rock(kelvin) == pytest.approx(rn.ROCK_STRENGTH_PA / 1e6)
    assert _rock(None) == pytest.approx(rn.ROCK_STRENGTH_PA / 1e6)


def test_hot_rock_weakens_by_the_olivine_creep_law():
    # strength research section 15 (MPa at 1e-15 per second): Hirth and Kohlstedt / Korenaga and Karato columns
    rule = rn.TABLE.temperature["rock"]
    for kelvin, hk, kk in ((1200.0, 7.3, 12.8), (1300.0, 2.3, 4.9), (1400.0, 0.84, 2.2)):
        low, high = (rn._creep_stress_mpa(column, rule["rate_per_s"], kelvin) for column in rule["columns"])
        assert low == pytest.approx(hk, rel=0.05) and high == pytest.approx(kk, rel=0.05)
        assert _rock(kelvin) == pytest.approx(math.sqrt(hk * kk), rel=0.07)        # the median is between the columns


def test_a_drawn_body_is_capped_by_its_own_cold_strength_and_otherwise_follows_its_column():
    assert _rock(1200.0, 0.0) == pytest.approx(4.0)                  # weak rock (4 MPa cold) is already under the 7.3 MPa creep stress at 1200 K
    assert _rock(1400.0, 0.0) == pytest.approx(0.84, rel=0.05)
    assert _rock(1300.0, 1.0) == pytest.approx(4.9, rel=0.05)
    assert [_rock(1300.0, d) for d in (0.0, 0.5, 1.0)] == sorted(_rock(1300.0, d) for d in (0.0, 0.5, 1.0))


def test_the_rock_creep_columns_match_at_1000_kelvin():
    rule = rn.TABLE.temperature["rock"]
    low, high = (rn._creep_stress_mpa(column, rule["rate_per_s"], 1000.0) for column in rule["columns"])
    assert low == pytest.approx(153.0, rel=0.03) and high == pytest.approx(152.0, rel=0.03)


def test_rock_strength_never_rises_as_it_warms_and_a_hotter_body_rounds_smaller():
    values = [_rock(kelvin) for kelvin in range(300, 1600, 50)]
    assert values == sorted(values, reverse=True)
    assert rn.round_radius_m(2500.0, 0.0, temperature_k=1300.0) < rn.round_radius_m(2500.0, 0.0) / 2.0


def test_the_hot_rubble_strength_is_olivine_creep_near_1330_kelvin():
    hot = rn.TABLE.rubble["hot_strength_pa"] / 1e6
    assert _rock(1330.0) == pytest.approx(hot, rel=0.45)
    assert _rock(1250.0) > 4.0 * 0.7 and _rock(1400.0) < 2.5          # the 1 to 4 MPa range is about 1250 to 1400 K


def test_a_hot_rubble_pile_with_a_shape_set_temperature_rounds_at_the_creep_strength():
    radius = rn.rubble_round_radius_m(2500.0, temperature_k=1300.0)
    assert radius == pytest.approx(rn.round_radius_m(2500.0, 0.0, temperature_k=1300.0))
    assert rn.is_round(rn.roundness_ratio(radius * 1.02, 2500.0, 0.0, temperature_k=1300.0, structure="rubble"))
    with pytest.raises(ValueError):
        rn.rubble_round_radius_m(2500.0, temperature_k=300.0)        # a cold rock pile has nothing to creep: friction only


def test_rock_in_an_icy_mix_is_unaffected_at_ice_temperatures():
    assert rn.yield_strength_pa(0.5, temperature_k=150.0) == pytest.approx(rn.yield_strength_pa(0.5, temperature_k=150.0, draw=0.5))
    assert rn.yield_strength_pa(0.2, temperature_k=100.0) > rn.yield_strength_pa(0.2)       # the cold-ice side of the mix is still stronger


def test_the_creep_rule_follows_the_table_rate(monkeypatch):
    document = copy.deepcopy(_document())
    document["temperature"]["rock"]["rate_per_s"] = 1.0e-14
    monkeypatch.setattr(rn, "TABLE", rn.load_table(document))
    assert _rock(1300.0) > 1.5 * 3.35            # a faster rate takes more stress (rate^(1/n)): 3.35 MPa at the stand-in rate


def test_a_table_may_leave_rock_without_a_temperature_rule(monkeypatch):
    document = copy.deepcopy(_document())
    del document["temperature"]["rock"]
    monkeypatch.setattr(rn, "TABLE", rn.load_table(document))
    assert _rock(1300.0) == pytest.approx(rn.ROCK_STRENGTH_PA / 1e6)


def test_a_creep_law_at_a_few_kelvin_does_not_overflow():
    rule = rn.TABLE.temperature["rock"]
    assert all(rn._creep_stress_mpa(column, rule["rate_per_s"], 18.0) == float("inf") for column in rule["columns"])
    assert rn.yield_strength_pa(0.0, temperature_k=18.0) == rn.yield_strength_pa(0.0)       # cold rock keeps its whole-body strength
