"""
The best-fit class table and scorer (GEN.1205, GEN.1206): `physics/data/planet_class_fit.json` gives every class its own positive
conditions and `physics/class_fit.py` scores a body against them. `planet_classes.classify` decides the stored letter from these rows (GEN.1207).
"""

import json
import math
import os

import pytest

from planetgen.physics import class_fit as cf, planet_classes as pc
from tests.test_planet_classes import CASES, EARTH, airless, body


def test_every_class_letter_has_its_own_complete_conditions():
    assert set(cf.LETTERS) == set(pc.LETTERS)
    by_kind = {(c.letter, c.kind) for c in cf.CLASSES}
    for cls in pc.CLASSES:
        assert (cls.letter, cls.kind) in by_kind, f"{cls.letter} ({cls.kind}) has no row in planet_class_fit.json"
    for cls in cf.CLASSES:
        assert cls.when and cls.name


def test_no_class_relies_on_an_earlier_one_saying_no():
    # a terrestrial class with air has at least two conditions of its own except the single-line ones the tree also had (A, V, E);
    # the point checked: no class has an empty condition list, and the data file has no ordering field.
    document = json.load(open(cf.DATA_FILE, encoding="utf-8"))
    assert all(row["when"] for row in document["classes"])
    assert "order" not in document and all("order" not in row for row in document["classes"])


def test_every_edge_a_condition_reads_has_a_width_in_the_fudge_table_and_the_hand_shoulders_are_gone():
    for cls in cf.CLASSES:
        for condition in cls.when:
            if "equals" not in condition:
                assert condition.get("min", condition.get("lt")) in cf.WIDTHS[condition["quantity"]]["edges"]
    document = json.load(open(cf.DATA_FILE, encoding="utf-8"))
    assert "shoulders" not in document and "shoulder_note" not in document      # GEN.1314: replaced by the table, not kept beside it
    assert hasattr(cf, "WIDTHS") and not hasattr(cf, "SHOULDERS")


def test_the_fudge_rows_are_registered():
    text = open(os.path.join(os.path.dirname(__file__), "..", "..", "docs", "design", "research", "gap-filling-plan.md"), encoding="utf-8").read()
    for row in ("FT1", "FT2", "FT3", "FT4", "CP1", "CP2"):
        assert f"| {row} |" in text


def test_membership_is_a_straight_line_across_the_shoulder():
    width = cf.WIDTHS["temperature_k"]["edges"][273.15]
    edge_min = {"quantity": "temperature_k", "min": 273.15}
    assert cf.membership(edge_min, {"temperature_k": 273.15}) == 1.0
    assert cf.membership(edge_min, {"temperature_k": 273.15 - width / 2}) == pytest.approx(0.5)
    assert cf.membership(edge_min, {"temperature_k": 273.15 - width}) == 0.0
    edge_lt = {"quantity": "temperature_k", "lt": 320.0}
    width = cf.WIDTHS["temperature_k"]["edges"][320.0]
    assert cf.membership(edge_lt, {"temperature_k": 319.0}) == 1.0
    assert cf.membership(edge_lt, {"temperature_k": 320.0 + width / 2}) == pytest.approx(0.5)
    assert cf.membership(edge_lt, {"temperature_k": 320.0 + 2 * width}) == 0.0
    assert cf.membership(edge_lt, {"temperature_k": None}) == 0.0


def test_a_pressure_edge_is_measured_in_the_log_of_a_ratio():
    edge = {"quantity": "pressure_pa", "min": 50000.0}
    width = cf.WIDTHS["pressure_pa"]["edges"][50000.0]
    assert cf.overshoot(edge, {"pressure_pa": 50000.0 * math.exp(-width / 2)}) == pytest.approx(0.5)
    assert cf.overshoot(edge, {"pressure_pa": 0.0}) == math.inf


def test_the_fudge_rule_fits_fully_first_then_the_fewest_values_then_the_smallest_move():
    # a body past one edge by a quarter width needs one value moved; a class past two edges needs two, so one value beats two however small the two are
    earth = body()
    full = cf.classify(earth)
    assert (full.letter, full.fudge_values, full.fudge_amount) == ("M", (), 0)
    width = cf.WIDTHS["temperature_k"]["edges"][275.0]
    near = cf.classify(body(surface_temperature_k=275.0 - 0.1 * width))
    assert near.fudge_values and 0.0 < near.match < 1.0 and near.fudge_amount <= cf.CAP_TIERS[0]
    assert cf.classify(body(surface_temperature_k=275.0 - 0.6 * width)).match < near.match      # still above the 273.15 K line of W, so M with a larger move


def test_the_smaller_cap_tier_is_tried_before_the_larger():
    within_first = {"gravity_g": 0.2}
    within_second = {"gravity_g": 0.9}
    past = {"gravity_g": 1.0}
    assert cf.tier_of(within_first) == 0 and cf.tier_of(within_second) == 1 and cf.tier_of(past) is None
    assert cf.match_of({}) == 1.0 and cf.match_of(past) == 0.0
    assert cf.match_of({"gravity_g": 0.5, "pressure_pa": 0.5}) == pytest.approx(0.5)


def test_a_body_inside_a_class_box_keeps_the_trees_letter_at_full_match():
    for letter, values in CASES.items():
        best, fit, runner, runner_fit = cf.best_class(values)
        assert pc.classify(values) == letter
        assert best == letter, f"{letter}: best fit {best} ({fit:.2f}), runner-up {runner} ({runner_fit:.2f})"
        assert fit == pytest.approx(1.0, abs=1e-9) or letter in ("K", "P", "E", "B")  # these cases sit on a shoulder of one quantity


def test_best_class_is_total_over_gated_bodies():
    for values in CASES.values():
        letter, fit, _runner, _runner_fit = cf.best_class(values)
        assert letter in cf.LETTERS and 0.0 < fit <= 1.0


def test_a_hot_thin_air_world_is_class_z_and_a_body_no_class_matches_has_no_best_letter():
    hot_thin = body(surface_temperature_k=414.0, atmospheric_pressure_pa=1000.0, ocean_fraction=0.0, hydrosphere="vapor", gravity_g=1.4)
    assert cf.best_class(hot_thin)[:2] == ("Z", 1.0)      # GEN.1221: the hot thin-air world has its own class
    nothing = body(surface_temperature_k=300.0, atmospheric_pressure_pa=10000.0, ocean_fraction=0.0, hydrosphere="vapor", gravity_g=0.9)
    assert cf.best_class(nothing)[:2] == (None, 0.0)


def test_a_cold_hydrogen_envelope_is_a_gas_dwarf_and_not_a_deep_frozen_world():
    values = body(surface_temperature_k=72.0, mass_earth=0.7, share_h2=0.98, share_o2=0.0, hydrosphere="ice", ocean_fraction=1.0,
                  p_o2_kpa=0.0, p_n2_kpa=0.0, p_ar_kpa=0.0, p_co2_kpa=0.0)
    put = dict(values)
    assert pc.classify(values) == "T" and cf.best_class(values)[0] == "T"


def test_one_condition_missed_outright_zeroes_the_match_and_a_soft_miss_only_lowers_it():
    earth = cf.scores(EARTH_VALUES := body())
    assert earth["M"] == pytest.approx(1.0)
    warm = cf.scores(body(surface_temperature_k=274.7))      # water still freezes at 273.15 K: the phase gate (GEN.1221) is a hard condition
    assert 0.0 < warm["M"] < 1.0
    cold = cf.scores(body(surface_temperature_k=250.0))
    assert cold["M"] == 0.0
    assert cf.scores(body(surface_temperature_k=271.0))["M"] == 0.0      # below freezing no water family applies


def test_a_reducing_wet_world_is_f_and_not_m_by_its_own_fit():
    from tests.test_planet_classes import no_oxygen
    values = no_oxygen(surface_temperature_k=300.0, ocean_fraction=0.6, share_ch4=0.06)
    letter, fit, runner, _ = cf.best_class(values)
    assert letter == "F" and fit == pytest.approx(1.0)
    assert cf.scores(values).get("M", 0.0) == 0.0


def test_gates_are_kinds_not_degrees():
    assert "J" not in cf.scores(body())
    assert "M" not in cf.scores(airless())
    assert "C" not in cf.scores(body())
    assert set(cf.scores(body(body_type="g", mass_earth=20.0))) == {"T", "I", "J"}


def test_the_ice_gas_giant_overlap_is_a_shoulder_not_a_cut():
    scored = cf.scores(body(body_type="g", mass_earth=40.0))
    assert scored["J"] == 1.0 and 0.0 < scored["I"] < 1.0


def test_a_plain_nitrogen_air_with_no_oxygen_is_f_and_a_carbon_dioxide_rich_one_is_g():
    """GEN.1281: heavy worlds with an inert nitrogen air, no oxygen, no methane and no hydrogen fit no class before F was widened from the reducing world to the anoxic world."""
    from tests.test_planet_classes import no_oxygen
    inert = no_oxygen(surface_temperature_k=300.0, gravity_g=1.2, share_ch4=0.0, share_h2=0.0, share_co2=0.01, ocean_fraction=0.5)
    detail = pc.classify_detail(inert)
    assert (detail.letter, detail.match) == ("F", pytest.approx(1.0)) and detail.digit == 0
    assert pc.class_of("F", "t").name == "Anoxic world"
    rich = no_oxygen(surface_temperature_k=300.0, gravity_g=1.2, share_ch4=0.0, share_h2=0.0, share_co2=0.5, ocean_fraction=0.5)
    assert pc.classify(rich) == "G"
    heavy_hot_thin = no_oxygen(surface_temperature_k=380.0, gravity_g=1.8, atmospheric_pressure_pa=9000.0, ocean_fraction=0.0, hydrosphere="vapor")
    assert pc.classify(heavy_hot_thin) == "Z"      # F stops at 1.4 g so the heavy thin-air class Z keeps its bodies


def test_a_heavy_oxygen_ocean_world_of_1_6_g_is_a_temperate_class_and_not_a_hole():
    # GEN.1306 part 2 (registry AW3): M, H, O, P and R reach V's 2 g line
    heavy = body(gravity_g=1.6)
    letter, fit, _runner, _fit2 = cf.best_class(heavy)
    assert letter == "M" and fit > 0.9
    assert cf.best_class(body(gravity_g=2.1))[0] == "V"


def test_a_carbon_dioxide_air_is_g_whatever_its_oxygen_and_a_thin_air_is_k_whatever_its_oxygen():
    # GEN.1306 parts 2 and 3 (registry AW4, AW5)
    rich = body(surface_temperature_k=290.0, atmospheric_pressure_pa=15000.0, share_co2=0.5, share_o2=0.35, ocean_fraction=1.0, gravity_g=1.5)
    assert cf.best_class(rich)[:2] == ("G", 1.0)
    thin = body(surface_temperature_k=290.0, atmospheric_pressure_pa=1000.0, share_co2=0.01, share_o2=0.63, ocean_fraction=0.0, hydrosphere="vapor", gravity_g=1.3)
    assert cf.best_class(thin)[:2] == ("K", 1.0)


def test_a_hot_thin_air_is_z_under_1_4_g_too_and_at_the_clamp_edges_of_1_and_10_kpa():
    # GEN.1306 parts 1 and 5 (registry AW6, AW7): Z reaches down to 0.5 g and up to 1500 K; the pressure edges hold through the shoulder
    for pressure in (989.0, 1000.0, 10000.0):
        hot = body(surface_temperature_k=626.0, atmospheric_pressure_pa=pressure, ocean_fraction=0.0, hydrosphere="vapor", gravity_g=1.14, share_o2=0.0)
        letter, fit, _runner, _fit2 = cf.best_class(hot)
        assert letter == "Z" and fit > (0.95 if pressure < 1000.0 else 0.99)      # 989 Pa is 1 percent past the 1 kPa edge: a tiny fudge
    very_hot = body(surface_temperature_k=1298.0, atmospheric_pressure_pa=10000.0, ocean_fraction=0.0, hydrosphere="vapor", gravity_g=1.59, share_o2=0.0, share_co2=0.88)
    assert cf.best_class(very_hot)[:2] == ("Z", 1.0)
    warm = body(surface_temperature_k=350.0, atmospheric_pressure_pa=5000.0, ocean_fraction=0.0, hydrosphere="vapor", gravity_g=1.0, share_o2=0.0)
    assert cf.best_class(warm)[0] == "F"       # F holds more conditions than Z, so the warm anoxic air stays F


def test_a_cold_air_world_between_minus_40_c_and_freezing_is_class_w_unless_class_p_fits_and_a_colder_one_is_e():
    # GEN.1306 part 4 (registry AW8, Boss 11-06): W is the frozen world, P the Earth-like cold world, E the deep-frozen one
    cold = body(surface_temperature_k=250.0, atmospheric_pressure_pa=10000.0, share_o2=0.0, ocean_fraction=0.1, hydrosphere="ice", gravity_g=0.83)
    assert cf.best_class(cold)[:2] == ("W", 1.0)
    assert cf.best_class(dict(cold, surface_temperature_k=220.0))[:2] == ("E", 1.0)
    assert cf.best_class(dict(cold, surface_temperature_k=274.0))[0] != "W"
    earthlike = body(surface_temperature_k=250.0, hydrosphere="ice", ocean_fraction=0.4, band1_mean_k=235.0, band4_mean_k=276.0)
    letter, fit, runner, _fit2 = cf.best_class(earthlike)
    assert letter == "P" and fit == pytest.approx(1.0) and runner == "W"      # both fit fully; P holds more conditions, so P wins the tie


def test_a_fifth_carbon_dioxide_is_g_whatever_the_oxygen():
    # GEN.1306 (registry AW5), with L (Boss's low-oxygen mask world) unchanged: 1 to 20 percent carbon dioxide (a share, class G's edge) with 3 to 16 kPa of oxygen (GEN.1321)
    assert cf.best_class(body(share_co2=0.3, share_o2=0.2))[0] == "G"
    assert cf.best_class(body(share_co2=0.05, share_o2=0.1))[0] == "L"


def test_oxygen_is_judged_by_pressure_and_not_by_share():
    """GEN.1321 (registry OP1): a thick air whose oxygen is a few percent of the molecules but 25 kPa is a temperate world; a thin air of 90 percent oxygen
    that is only 9 kPa is a marginal one; the same oxygen share on a thin air is not enough."""
    thick = body(atmospheric_pressure_pa=1.0e6, p_o2_kpa=25.0, p_n2_kpa=950.0, p_ar_kpa=9.0, p_co2_kpa=0.4)
    detail = cf.classify(thick)
    assert (detail.letter, detail.digit) == ("M", 9)                          # 2.5 percent of the air, 25 kPa of oxygen: M, limited by its air amount (M's air is 50 to 250 kPa)
    thin = body(atmospheric_pressure_pa=10000.0, p_o2_kpa=9.0, p_n2_kpa=0.9, p_ar_kpa=0.05, p_co2_kpa=0.05)
    assert cf.best_class(thin)[0] != "M"                                      # 90 percent of the air, only 9 kPa
    assert cf.scores(body(p_o2_kpa=40.0))["M"] < 1.0 and cf.scores(body(p_o2_kpa=35.4))["M"] == pytest.approx(1.0)    # the fire limit is the upper edge


def test_carbon_dioxide_of_the_oxygen_water_classes_is_judged_by_pressure():
    """GEN.1321 (registry OP1): under 1 kPa of carbon dioxide M fits fully, over it L takes over from M, and a thick air whose carbon dioxide is under 1 percent of the molecules
    but over 1 kPa is not M; class G keeps its share."""
    assert cf.scores(body(p_co2_kpa=0.5))["M"] == pytest.approx(1.0)
    assert cf.scores(body(p_co2_kpa=1.5))["M"] < 1.0
    deep = body(atmospheric_pressure_pa=2.0e5, p_co2_kpa=1.8, p_o2_kpa=40.0, p_n2_kpa=150.0, p_ar_kpa=2.0)       # 0.9 percent of the molecules, 1.8 kPa
    assert cf.scores(deep)["M"] < 1.0
    assert cf.best_class(body(share_co2=0.3, share_o2=0.2))[0] == "G"


def test_the_o_edge_is_95_percent_of_liquid_cover_and_m_runs_up_to_it():
    # GEN.1570 (Boss's edge, registry CL6): M is 20 to 95 percent liquid cover, O is 95 percent or more. Each is 1 on its own side of the edge and falls across the
    # 0.05 cover shoulder (CF5) on the other, so a body AT 0.95 scores 1.0 in both and the narrower description (M has more conditions) takes it; a full ocean is O at 1.0 and M at 0.
    def at(cover):
        found = cf.scores(body(surface_temperature_k=285.0, ocean_fraction=cover))
        return found["M"], found["O"]
    assert at(0.8) == (1.0, 0.0)                                              # the old O edge: a body of 80 percent is now plain M
    assert at(0.949)[0] == 1.0 and at(0.949)[1] < 1.0
    assert at(0.95) == (1.0, 1.0) and cf.best_class(body(surface_temperature_k=285.0, ocean_fraction=0.95))[0] == "M"
    assert at(0.951)[1] == 1.0 and cf.best_class(body(surface_temperature_k=285.0, ocean_fraction=0.951))[0] == "O"
    assert at(1.0) == (0.0, 1.0)


def test_a_cool_global_ocean_with_oxygen_is_o_and_a_cool_damp_desert_is_h_and_a_world_with_no_water_is_neither():
    # GEN.1342 (Boss): O and H take their mean temperature from M's floor (275 K) up to their old ceilings (312 and 315 K)
    for kelvin in (276.0, 285.0, 300.0):
        ocean = body(surface_temperature_k=kelvin, ocean_fraction=0.97)
        assert cf.best_class(ocean)[:2] == ("O", 1.0), kelvin
        desert = body(surface_temperature_k=kelvin, ocean_fraction=0.1, groundwater="aquifer")
        assert cf.best_class(desert)[0] == "H" and cf.best_class(desert)[1] > (0.9 if kelvin < 300.0 else 0.5), kelvin      # at 300 K the wet bulb (25.7 C) is 0.4 width past the 25 C line
    assert cf.best_class(body(surface_temperature_k=285.0, ocean_fraction=0.5))[0] == "M"
    assert cf.best_class(body(surface_temperature_k=305.0, ocean_fraction=0.97))[0] == "O"      # 309 K has a wet bulb of 33.7 C, past the 31 C line by 1.8 widths: a hole, nearest O at match 0
    dry = body(surface_temperature_k=285.0, ocean_fraction=0.0, groundwater="none")
    assert cf.best_class(dry)[0] != "H"


def test_a_nan_quantity_is_no_match_and_classification_refuses_it():
    """GEN.1593: every comparison with NaN is false, which read as a perfect pass; a NaN in a class quantity now scores 0 and the generation path stops on it."""
    from tests.test_planet_classes import body
    import pytest
    from planetgen.physics import class_fit as cf, planet_classes as pc

    assert cf.membership({"quantity": "gravity_g", "min": 0.75}, {"gravity_g": float("nan")}) == 0.0
    assert cf.membership({"quantity": "gravity_g", "lt": 2.0}, {"gravity_g": float("nan")}) == 0.0
    assert cf.membership({"quantity": "groundwater_present", "equals": 1.0}, {"groundwater_present": float("nan")}) == 0.0
    sound = cf.scores(body())
    nan_gravity = cf.scores(body(gravity_g=float("nan")))
    assert nan_gravity.get("M", 0.0) == 0.0 < sound["M"]
    assert cf.scores(body(atmospheric_pressure_pa=float("nan"))) == {}      # no side of the air gate
    with pytest.raises(pc.ClassError, match="not a number"):
        cf.classify(body(gravity_g=float("nan")))
    with pytest.raises(pc.ClassError, match="gravity_g"):
        pc.classify_detail(body(gravity_g=float("nan")))
