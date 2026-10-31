# tests/test_atmosphere_recipes.py

"""
GEN.478: the atmosphere possibility tables (`planetgen.physics.atmosphere_recipes`): the 16 tables and 104 rows of part 09, each
with a source, a state and a gate; the worked examples a to p of the part reproduced from the rows; the corrections of the later
research; the open rows listed with their research item; the draws, deterministic at any worker count; and the loader's refusals.
"""

import copy
import json
import math
import concurrent.futures
import multiprocessing
import os
import sys

import pytest

from planetgen.generation import possibility_table as pt
from planetgen.generation.possibility_table import Subject
from planetgen.physics import envelope_escape
from planetgen.physics import atmosphere_recipes as ar
from planetgen.physics.atmosphere_recipes import INPUTS, RECIPES, Recipes
from planetgen.physics.elements import DataError
from planetgen.util import draw

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
BOOLEANS = {"is_moon", "long_lived_geology", "ocean_present", "water_lost", "co2_rich", "h2_poor", "life_chemical",
            "sep_reducing", "lightning", "gas_has_co_h2s", "has_solid_ground", "h2o_present", "carbon_bearing", "h2_air",
            "thin_disc_fgk", "planetesimals_early"}

SOLAR_HOST = {"n_over_c": 1.0, "ar_over_c": 1.0, "s_over_c": 1.0}
"""dict: The host ratios a comet's species rows read (GEN.1005, k7), at the proto-solar values."""


def document():
    with open(ar.DATA_FILE, encoding="utf-8") as handle:
        return json.load(handle)


def broken(change):
    """A loader over a copy of the data file edited by `change(document)`."""
    doc = copy.deepcopy(document())
    change(doc)
    return Recipes(doc, register=False)


def row_of(doc, row_id):
    return next(e for e in doc["entries"] if e["id"] == row_id)


def main(table, floor=1e-5):
    """The rows of a probability table with more than a rounding error of probability."""
    return {name: p for name, p in dict(table).items() if p > floor}


def probs(group, body_class, **quantities):
    return main(RECIPES.probabilities(group, Subject(body_class, quantities)))


def full(body_class="terrestrial", **overrides):
    """A subject that carries every stored quantity a gate or weight reads, set to something plain."""
    q = {name: (False if name in BOOLEANS else 1.0) for name in INPUTS}
    q.update(zone="e", host_type="G", family="stony", redox_diw=0.5, mass_earth=1.0, mass_jupiter=0.003, radius_km=6371.0,
            v_esc_kms=11.186, teq_k=255.0, ts_k=288.0, tday_k=300.0, tp_k=300.0, t_layer_k=250.0, teff_k=300.0,
            pressure_bar=1.0, insolation=1.0, luminosity_solar=1.0, xuv_index=1.0, snow_ratio=0.5, host_mass_solar=1.0,
            age_gy=4.5, core_radius_earth=1.0, perihelion_au=5.0, profile_t_min_k=200.0, profile_t_max_k=300.0,
            thin_disc_fgk=True, long_lived_geology=True, ocean_present=True, has_solid_ground=True, h2o_present=True,
            o2_pal=1.0, carbon_bearing=True, u_star=0.5, u_star_threshold=0.3, sulfur_flux_earth=1.0, grain_um=1.0)
    q.update(overrides)
    return Subject(body_class, q)


# --- the tables --------------------------------------------------------------------------------------------------------

def test_the_part_has_16_tables_and_104_rows_with_the_row_counts_of_each_table():
    assert len(RECIPES) == 104
    counts = {letter: len(RECIPES.in_category(letter)) for letter in "abcdefghijklmnop"}
    assert counts == {"a": 7, "b": 6, "c": 6, "d": 7, "e": 6, "f": 8, "g": 6, "h": 5, "i": 6, "j": 5, "k": 7, "l": 3,
                      "m": 11, "n": 9, "o": 6, "p": 6}
    assert [t.name for t in RECIPES.tables] == [f"atmosphere-{c}" for c in "abcdefghijklmnop"]
    assert all(t.name in pt.REGISTRY for t in RECIPES.tables)


def test_the_status_of_each_row_is_the_parts_except_the_rows_the_later_decisions_moved():
    statuses = {}
    for entry in RECIPES.entries:
        statuses[entry.status] = statuses.get(entry.status, 0) + 1
    # the part: SOURCED 8, DERIVED 65, DESIGN 31 (b1 to b4 moved by GEN.802; a1 and a2 by GEN.482; f7, h2 to h4, i6 and k1 to k7 by GEN.484 to GEN.488); Boss's fork 7 (2026-10-18) makes j2 to j4 outputs of the retention physics
    assert statuses == {"SOURCED": 8, "DERIVED": 70, "DESIGN": 26}
    assert {e.id for e in RECIPES.entries if e.status == "DERIVED"} >= {"j2", "j3", "j4"}


def test_every_row_has_a_source_a_state_and_a_gate_and_the_states_follow_the_rule():
    for entry in RECIPES.entries:
        assert entry.source.strip(), entry.id
        assert entry.state in pt.STATES, entry.id
        assert entry.classes, entry.id
        if entry.state == pt.PARTLY_VERIFIED:
            assert entry.limit, entry.id
        if entry.state == pt.STAND_IN:
            assert entry.replacement.startswith("GEN."), entry.id
        if entry.state == pt.VERIFIED:
            assert set(entry.labels) <= {"READ", "COMPUTED", "DATA"}, (entry.id, entry.labels)
        if entry.state == pt.OPEN_STATE:
            assert entry.status == "DESIGN" and entry.items, entry.id
    assert RECIPES.state_counts() == {"verified": 37, "partly verified": 54, "open": 0, "stand-in": 13}


def test_every_research_item_gen_480_to_494_is_answered_by_a_row_state():
    """The research found no measurement for some rows: they ship as labelled stand-ins that name their builder."""
    assert RECIPES.research == {f"R{n}": f"GEN.{479 + n}" for n in range(1, 16)}
    assert RECIPES.open_rows() == [], "no row waits on research any more"
    researched = {r for entry in RECIPES.entries for r in entry.research}
    assert researched == set(RECIPES.research), "every research item concerns at least one row"
    builders = {f"GEN.{n}" for n in range(509, 526)} | {"GEN.595", "GEN.1021"}      # GEN.1021: the 67P calibration of row k7 (GEN.1005 replaced the rest with models)
    for entry in RECIPES.entries:
        if entry.state == pt.STAND_IN:
            assert entry.replacement in builders, (entry.id, entry.replacement)


def test_the_stand_ins_are_rows_with_no_measurement_or_built_from_solar_system_bodies():
    stand_ins = {e.id: e.replacement for e in RECIPES.entries if e.state == pt.STAND_IN}
    assert {"b5", "n7", "o3"} <= set(stand_ins)
    # GEN.1002: a row whose weights rest on Solar System bodies is a flagged stand-in with its replacing item; GEN.1005 replaced a1, a2 and b1 with models
    # (the bodies are their tests) and left k7 resting on the one comet measured, until GEN.1021
    assert stand_ins["k7"] == "GEN.1021" and not {"a1", "a2", "b1"} & set(stand_ins)
    assert not {"k1", "k2", "k3", "k4", "k5"} & set(stand_ins), "the comet rows follow the Robinson database (GEN.520)"
    assert not {"m3", "m4", "m11"} & set(stand_ins), "the life-linked gases are design flags (GEN.525)"
    assert not {"a3", "a4", "f7", "p6"} & set(stand_ins), "rows with a published statistic are not stand-ins"


def test_the_registered_tables_report_the_stand_ins_and_limits():
    report = pt.REGISTRY.report()
    assert "atmosphere-a (version 1)" in report
    assert "atmosphere-b: b1 " in report                      # a stand-in row of table b
    assert "atmosphere-a: a2 " in report and "replaced by GEN.513" in report
    assert "limit:" in report


def test_the_digest_of_a_table_changes_with_a_weight_expression():
    other = broken(lambda d: row_of(d, "a5").update(weight="phi((0 - redox_diw)/0.4)"))
    assert other.tables[0].digest() != RECIPES.tables[0].digest()
    same = broken(lambda d: None)
    assert same.tables[0].digest() == RECIPES.tables[0].digest()


def test_the_data_file_is_what_the_builder_writes():
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    try:
        import build_atmosphere_recipes as builder
    finally:
        sys.path.pop(0)
    with open(builder.OUT, encoding="utf-8") as handle:
        assert handle.read() == builder.render(builder.build())


# --- the worked examples of part 09, reproduced from the rows ---------------------------------------------------------

def test_example_a_the_redox_classes_at_one_drawn_value():
    p = probs("a", "terrestrial", redox_diw=0.5)
    assert p["a5"] == pytest.approx(0.159, abs=5e-4) and p["a6"] == pytest.approx(0.841, abs=5e-4)
    assert p.get("a7", 0) == pytest.approx(0.0, abs=5e-4)
    # GEN.1005: the magma-ocean mean is the model of physics.magma_ocean_redox, called from the row
    assert RECIPES.entry("a1").spec["dist"]["mu"] == "magma_ocean_diw(mass_earth)"
    d = RECIPES._dist_of(RECIPES.entry("a1"), Subject("terrestrial", {"mass_earth": 0.05}))
    assert d["mu"] == pytest.approx(-1.0743, abs=1e-3)
    d = RECIPES._dist_of(RECIPES.entry("a1"), Subject("terrestrial", {"mass_earth": 1.0}))
    assert d["mu"] == pytest.approx(1.0724, abs=1e-3)
    later = RECIPES._dist_of(RECIPES.entry("a2"), Subject("terrestrial", {"mass_earth": 1.0}))
    assert later["mu"] == pytest.approx(3.4724, abs=1e-3) and later["sigma"] == 1.0


def test_the_magma_ocean_redox_follows_the_model_and_the_first_hand_ranges_it_is_tested_against_gen_1005():
    """The row calls physics.magma_ocean_redox; the Moon, Mars and Earth ranges of Hirschmann 2022 section 6.3 and Deng 2020 are checks (tests/test_magma_ocean_redox.py)."""
    from planetgen.physics import magma_ocean_redox as mor

    def a1(mass):
        return RECIPES._dist_of(RECIPES.entry("a1"), Subject("terrestrial", {"mass_earth": mass}))

    for mass in (0.0123, 0.107, 1.0, 5.4):
        assert a1(mass)["mu"] == pytest.approx(mor.surface_diw(mass)) and a1(mass)["sigma"] == pytest.approx(mor.diw_sigma(mass))
    assert -2.2 - a1(0.0123)["sigma"] <= a1(0.0123)["mu"] <= -1.5 + a1(0.0123)["sigma"]
    assert a1(30.0)["mu"] == pytest.approx(a1(6.0)["mu"]), "held where the base pressure passes the 140 GPa of Deng 2020"
    later = RECIPES._dist_of(RECIPES.entry("a2"), Subject("terrestrial", {"mass_earth": 30.0}))
    assert later["mu"] == pytest.approx(mor.later_mantle_diw(1.0)), "no source for the later mantle of a super-Earth: the value at one Earth mass is held"



def _delivery_law(x, mass):
    """The class shares of the formation-distance law of GEN.802 from the coefficients of its data file, written out again."""
    import json
    import os
    from planetgen.physics import elements
    with open(os.path.join(elements.DATA_DIRECTORY, "ngpps_formation_water.json"), encoding="utf-8") as handle:
        model = json.load(handle)["model"]
    lx = math.log(min(max(x, model["x_range"][0]), model["x_range"][1]))
    lm = math.log10(min(max(mass, model["mass_range_earth"][0]), model["mass_range_earth"][1]))
    hinge = max(0.0, lx - math.log(model["knot"]))
    events = [1 / (1 + math.exp(-(a + b * lx + c * lm + h * hinge))) for a, b, c, h in model["coefficients"]]
    events[1], events[2] = min(events[1], events[0]), min(events[2], min(events[1], events[0]))
    return 1 - events[0], events[0] - events[1], events[1] - events[2], events[2]


def test_example_b_the_delivery_split_by_the_formation_distance_over_the_snow_line_and_the_mass():
    for x, mass in ((0.447, 1.0), (0.75, 5.0), (1.5, 5.0), (0.1, 0.3), (6.0, 20.0)):
        p = probs("b", "terrestrial", snow_ratio=x, mass_earth=mass, **{"const:b:c_b4": 1.0})
        assert [p.get(row, 0.0) for row in ("b1", "b2", "b3", "b4")] == pytest.approx(list(_delivery_law(x, mass)), abs=2e-6), (x, mass)
        assert sum(p.values()) == pytest.approx(1.0, abs=2e-5)           # `probs` drops the shares under 1e-5
    # no gate at the snow line any more: 60 percent of the planets that formed at three quarters of it are water-rich (the population)
    assert probs("b", "terrestrial", snow_ratio=0.75, mass_earth=5.0, **{"const:b:c_b4": 1.0})["b4"] == pytest.approx(0.564, abs=2e-3)
    # the plateau moves the water-rich planets it removes into the 'high' class
    full_law = probs("b", "terrestrial", snow_ratio=0.75, mass_earth=5.0, **{"const:b:c_b4": 1.0})
    half = probs("b", "terrestrial", snow_ratio=0.75, mass_earth=5.0, **{"const:b:c_b4": 0.5})
    assert half["b4"] == pytest.approx(0.5 * full_law["b4"]) and half["b3"] == pytest.approx(full_law["b3"] + 0.5 * full_law["b4"])
    assert half["b1"] == pytest.approx(full_law["b1"]) and half["b2"] == pytest.approx(full_law["b2"])


def test_example_c_the_outgassed_regime_under_the_redox_class_and_the_joint_probability():
    s = full("terrestrial", redox_diw=0.5, ts_k=1000.0, pressure_bar=100.0, **{"drawn.a": "a5"})
    p = dict(RECIPES.probabilities("c", s))
    assert p == pytest.approx({"c1": 0.463, "c2": 0.537}, abs=1e-3)
    s3 = full("terrestrial", ts_k=1000.0, **{"drawn.a": "a6"})
    assert dict(RECIPES.probabilities("c", s3)) == {"c3": 1.0}
    after = full("terrestrial", ts_k=400.0, ocean_present=True, **{"drawn.c": "c3"})
    assert dict(RECIPES.probabilities("c-after", after)) == {"c5": 1.0}
    steam = full("terrestrial", ts_k=700.0, **{"drawn.c": "c3"})
    assert dict(RECIPES.probabilities("c-after", steam)) == {"c6": 1.0}
    dry = full("terrestrial", ts_k=400.0, ocean_present=False, **{"drawn.c": "c3"})
    with draw.bound(1):
        assert RECIPES.draw_group("c-after", dry).name == "c6"          # no row passes the weights: the fallback wins


def test_example_d_the_shoreline_is_the_published_logistic_and_the_regimes_share_its_probability():
    def at(flux, v=11.186, lum=1.0):
        return dict(RECIPES.probabilities("d", Subject("terrestrial", dict(insolation=flux, v_esc_kms=v, luminosity_solar=lum))))
    # at the shoreline flux of this escape speed and luminosity (log10 f0 = 2.62) half the planets keep an atmosphere
    shore = 10 ** 2.62
    assert at(shore)["d1"] == pytest.approx(0.5, abs=1e-9)
    assert at(shore * 10 ** 0.247)["d1"] == pytest.approx(1 - 1 / 11, abs=1e-9)         # one width above: 1 / (1 + 10^-1)... 
    assert at(1.0)["d5"] > 0.99
    far = at(5000.0)
    assert far["d1"] > 0.999
    for flux in (1, 10, 100, 300, 1000, 5000):
        p = at(flux)
        assert sum(p.values()) == pytest.approx(1.0, abs=1e-12)
        # the four atmosphere windows add up to P(atmosphere): airless is 1 - P(atmosphere)
        delta = math.log10(flux) - 2.62
        p_atm = 1 / (1 + 10 ** (delta / 0.247))
        assert p.get("d1", 0.0) == pytest.approx(1 - p_atm, abs=1e-9)
    # the slope in the escape speed is 5.9 and in the luminosity 1.17 (Berta-Thompson 2025)
    assert at(shore * 2 ** 5.9, v=11.186 * 2)["d1"] == pytest.approx(0.5, abs=1e-9)
    assert at(shore * 10 ** 1.17, lum=10.0)["d1"] == pytest.approx(0.5, abs=1e-9)
    assert RECIPES.derive("d6", Subject("terrestrial", {"drawn.d": "d3"})) == {"p_max_bar": 0.01}
    assert RECIPES.derive("d6", Subject("terrestrial", {"drawn.d": "d5"})) == {"p_max_bar": None}


def test_the_retired_k_window_is_gone_and_the_shoreline_is_marked_as_corrected():
    d1 = RECIPES.entry("d1")
    assert "2507.02136" in d1.source and "Corrected 2026-10-18" in d1.notes
    assert "K window" in d1.notes and "13 Solar System bodies" in d1.notes
    assert "3.0e-4" not in d1.weight.text


def test_the_hydrogen_flag_of_d7_uses_the_jeans_parameter_of_the_part():
    gate = dict(mass_earth=5.0, teq_k=300.0, xuv_index=1.5)
    assert RECIPES.flag("d7", Subject("giant", dict(gate, v_esc_kms=25.0)))
    assert not RECIPES.flag("d7", Subject("giant", dict(gate, v_esc_kms=5.0)))


def test_example_e_the_envelope_the_water_world_and_the_hycean_share_the_water_rich_probability():
    s = Subject("giant", dict(mass_earth=5.0, core_radius_earth=5 ** 0.27, insolation=20.0, snow_ratio=1.5, host_mass_solar=1.0, luminosity_solar=1.0,
                          age_gy=5.0, teq_k=255.0 * 20.0 ** 0.25, **{"const:e:f_env0": 0.05, "const:e:eta": 0.224, "const:b:c_b4": 0.1454}))
    assert RECIPES.probability_of("b4", s) == pytest.approx(0.1454 * _delivery_law(1.5, 5.0)[3], abs=1e-9)
    p = dict(RECIPES.probabilities("e", s))
    lost = envelope_escape.lost_fraction(5.0, 0.05, 20.0, 1.0, 1.0, 5.0, 255.0 * 20.0 ** 0.25)
    assert 0.2 < lost < 0.9, "a 5 Earth-mass planet at 20 times the Earth's insolation has lost part of its envelope"
    assert p["e1"] / (p["e1"] + p["e2"]) == pytest.approx(lost, abs=1e-9)
    assert p["e3"] == pytest.approx(0.111, abs=1e-3) and p["e4"] == pytest.approx(0.032, abs=1e-3)
    assert sum(p.values()) == pytest.approx(1.0)
    bins = RECIPES.bins("e5", Subject("giant", {"mass_earth": 5.0}))
    assert [round(b[2], 3) for b in bins] == [0.268, 0.294, 0.268, 0.169]
    z = [0.02, 0.0548, 0.173, 0.424]
    mu = [RECIPES.derive("e5", Subject("giant", {"mass_earth": 5.0, "drawn.e": "e2"}), value=v)["mean_molecular_weight"] for v in z]
    assert mu[0] == pytest.approx(2.4 + 3.1 * math.log10(2) / math.log10(58))
    assert RECIPES.probability("e6", Subject("giant", {"mass_earth": 5.0, "teq_k": 500.0, "drawn.e": "e2"})) == pytest.approx(0.5)


def test_example_f_the_giant_regimes_the_metallicity_bins_and_the_carbon_to_oxygen_bins():
    s = Subject("giant", dict(mass_jupiter=0.3, teq_k=950.0, zone="h"))
    p = main(RECIPES.probabilities("f", s))
    assert p == pytest.approx({"f2": 0.0062, "f3": 0.9936, "f4": 0.0002}, abs=1e-4)
    assert [round(b[2], 3) for b in RECIPES.bins("f6", s)] == [0.0, 0.013, 0.369, 0.566, 0.052]
    carbon = RECIPES.bins("f7", s)
    assert sum(b[2] for b in carbon) == pytest.approx(1.0, abs=1e-12)
    assert [(b[0], b[1]) for b in carbon] == [(0.2, 0.35), (0.35, 0.7), (0.7, 0.8), (0.8, 1.0), (1.0, 1.5)]


def test_the_giants_carbon_to_oxygen_correction_the_tail_is_explicit_and_not_a_log_normal_spill():
    s = Subject("giant", dict(mass_jupiter=0.3, teq_k=950.0, zone="h"))
    tail = math.sqrt(0.043 * 0.087)
    above = RECIPES.probability_in("f7", s, 0.8, None)
    assert above == pytest.approx(tail, abs=1e-9)                  # the core stops at 0.8, so only the tail is above it
    assert 0.043 <= above <= 0.087
    # part 05's recipe, 0.55 * 10^N(0, 0.15) truncated to [0.2, 1.5], put 13.9 percent above 0.8
    lo, hi, x = math.log10(0.2 / 0.55) / 0.15, math.log10(1.5 / 0.55) / 0.15, math.log10(0.8 / 0.55) / 0.15
    old = (ar.phi(hi) - ar.phi(x)) / (ar.phi(hi) - ar.phi(lo))
    assert old == pytest.approx(0.139, abs=3e-3)
    assert "13.9 percent" in RECIPES.entry("f7").notes
    # the tail is uniform on 0.8 to 1.5, the core centre is 0.57
    assert RECIPES.probability_in("f7", s, 1.0, None) == pytest.approx(tail * 0.5 / 0.7, abs=1e-9)
    d = RECIPES._dist_of(RECIPES.entry("f7"), s)
    assert d["core"]["mu"] == pytest.approx(math.log10(0.57))
    # the tail share draws inside the data's bracket
    with draw.bound(3):
        shares = [RECIPES.draw_group("e", full("giant", mass_earth=5.0, core_radius_earth=1.5, snow_ratio=0.5)).constants[("e", "f_env0")]
                  for _ in range(5)]
    assert all(0.05 <= v <= 0.2 for v in shares)


def test_example_g_h_and_i_the_brown_dwarf_the_hot_jupiter_and_the_lava_windows():
    g = Subject("brown-dwarf", {"teff_k": 848.4})
    assert main(RECIPES.probabilities("g", g)) == pytest.approx({"g3": 0.9801, "g4": 0.0199}, abs=1e-3)
    assert RECIPES.probability("g6", g) == pytest.approx(0.994, abs=1e-3)
    assert main(RECIPES.probabilities("h", Subject("giant", {"tday_k": 2500.0, "zone": "h", "photosphere_bar": 0.1}))) == pytest.approx(
        {"h2": 0.9983, "h3": 0.0017}, abs=1e-4)
    assert main(RECIPES.probabilities("h", Subject("giant", {"tday_k": 900.0, "zone": "h", "photosphere_bar": 0.1}))) == pytest.approx(
        {"h1": 0.7454, "h2": 0.2546}, abs=1e-4), "the methane line at 0.1 bar is 966 K, not 1500 K"
    s = Subject("giant", {"tday_k": 2150.0, "zone": "h", "teq_k": 2350.0, "drawn.h": "h2"})
    assert RECIPES.probability("h4", s) == 0.0, "TiO is not seen at 2350 K and below"
    assert RECIPES.probability("h4", Subject("giant", dict(s.quantities, teq_k=2641.0))) == 1.0, "and is seen at 2641 K"
    lava = lambda t: dict(RECIPES.probabilities("i", Subject("terrestrial", {"tp_k": t, "zone": "h"})))  # noqa: E731
    assert main(lava(2060.0)) == pytest.approx({"i1": 0.326, "i2": 0.674}, abs=1e-3)
    assert main(lava(2500.0)) == pytest.approx({"i2": 0.929, "i3": 0.071}, abs=1e-3)
    s = Subject("terrestrial", {"tp_k": 2060.0, "zone": "h", "mass_earth": 3.0})
    assert RECIPES.derive("i6", s)["p_eq_bar"] == pytest.approx(3.5e-4, rel=0.02)


def test_the_lava_world_boundaries_are_the_text_edges_of_miguel_2011_and_the_difference_is_63_k():
    edges = [2033, 2588, 2890, 3168]
    for row, (lo, hi) in zip(("i2", "i3", "i4"), zip(edges, edges[1:])):
        assert f"window(tp_k, {lo}, {hi}, 60)" == RECIPES.entry(row).weight.text
    assert 2096 - 2033 == 63                                  # the caption list differs by 63 K and not by 130 K


def test_example_j_the_frost_and_evolved_atmospheres_are_outputs_of_the_retention_physics():
    titan = full("terrestrial", is_moon=True, radius_km=2500.0, ts_k=61.0, insolation=0.012, **{"drawn.d": "d5"})
    frost = full("terrestrial", ts_k=61.0, **{"drawn.d": "d5"})
    w3 = RECIPES.raw_weight(RECIPES.entry("j3"), titan)
    w4 = RECIPES.raw_weight(RECIPES.entry("j4"), frost)
    assert w3 / (w3 + w4) == pytest.approx(0.471, abs=1e-3)   # the two windows overlap from 57 to 67 K
    mars = full("terrestrial", mass_earth=0.1, ts_k=220.0, **{"drawn.d": "d3"})
    assert dict(RECIPES.probabilities("j", mars)) == {"j2": 1.0}
    none = full("terrestrial", mass_earth=0.1, ts_k=300.0, **{"drawn.d": "d3"})
    with draw.bound(1):
        assert RECIPES.draw_group("j", none).name == "none"   # the explicit outcome: none of these applies
    for entry in (RECIPES.entry("j2"), RECIPES.entry("j3"), RECIPES.entry("j4")):
        assert entry.status == "DERIVED" and not entry.research


def test_example_k_the_coma_classes_and_the_carbon_dioxide_to_water_bins():
    comet = Subject("comet", {"perihelion_au": 1.0, "family": "x", **SOLAR_HOST})
    assert dict(RECIPES.probabilities("k", comet)) == pytest.approx({"k1": 0.52, "k2": 0.36, "k3": 0.12}, abs=1e-3)
    assert [round(b[2], 3) for b in RECIPES.bins("k4", comet)] == [0.136, 0.335, 0.335, 0.194]
    quiet = Subject("comet", {"perihelion_au": 5.0, "family": "x", **SOLAR_HOST})
    with draw.bound(2):
        assert RECIPES.draw_group("k", quiet).name == "none"
    icy = Subject("asteroid", {"perihelion_au": 2.0, "family": "icy", **SOLAR_HOST})
    assert set(dict(RECIPES.probabilities("k", icy))) == {"k1", "k2", "k3"}


def test_example_l_m_n_and_o_the_white_dwarf_atmosphere_the_trace_gases_the_clouds_and_the_radioisotopes():
    wd = Subject("terrestrial", {"host_type": "WD", "ts_k": 1800.0})
    p = dict(RECIPES.probabilities("l", wd))
    assert p["l1"] == pytest.approx(0.997, abs=1e-3) and RECIPES.derive("l1", wd)["p_eq_bar"] == pytest.approx(2.39e-5, rel=0.01)
    giant = Subject("giant", {"t_layer_k": 900.0, "gas_has_co_h2s": True})
    assert RECIPES.probability("m7", giant) == pytest.approx(0.030, abs=1e-3)
    assert RECIPES.probability("m8", giant) == pytest.approx(1.0, abs=1e-3)
    assert RECIPES.derive("m9", giant)["k_eq"] == pytest.approx(0.034, abs=1e-3)
    cloud = Subject("giant", {"profile_t_min_k": 300.0, "profile_t_max_k": 2000.0, "teq_k": 800.0})
    assert RECIPES.crossed_windows("n1", cloud) == ["Al-Ti oxides", "Fe", "forsterite", "enstatite", "Cr", "MnS"]
    assert RECIPES.crossed_windows("n2", cloud) == ["Na2S", "KCl/ZnS"]
    assert RECIPES.crossed_windows("n3", cloud) == []
    assert RECIPES.probability("n6", cloud) == pytest.approx(0.067, abs=1e-3)
    sun = Subject("terrestrial", {"thin_disc_fgk": True, "age_gy": 3.0})
    assert [round(b[2], 4) for b in RECIPES.bins("o1", sun)] == [0.0, 0.8729, 0.1271]
    assert RECIPES.derive("o5", sun)["heating_ratio"] == pytest.approx(1.432, abs=1e-3)
    assert RECIPES.derive("o5", Subject("terrestrial", {"age_gy": 8.0}))["heating_ratio"] == pytest.approx(0.589, abs=1e-3)


def test_example_p_the_uv_regime_depends_on_the_stellar_mass_through_the_saturation_time():
    early_m = Subject("terrestrial", {"host_type": "M", "host_mass_solar": 0.55, "age_gy": 1.0})
    v = RECIPES.derive("p4", early_m)
    assert v["tau_sat_gy"] == pytest.approx(0.162) and v["lx_over_lbol"] == pytest.approx(8.5e-5, rel=0.02)
    for mass, tau in ((0.25, 0.6), (0.15, 2.3)):
        s = Subject("terrestrial", {"host_type": "M", "host_mass_solar": mass, "age_gy": 1.0})
        assert RECIPES.derive("p4", s)["tau_sat_gy"] == pytest.approx(tau)
    young_m = Subject("terrestrial", {"host_type": "M", "host_mass_solar": 0.15, "age_gy": 0.5})
    p = dict(RECIPES.probabilities("p", young_m))
    assert p["p4"] > 0.99 and p["p3"] < 0.01
    sun = Subject("terrestrial", {"host_type": "G", "host_mass_solar": 1.0, "age_gy": 4.5})
    assert dict(RECIPES.probabilities("p", sun)) == {"p1": 1.0}
    infant = Subject("terrestrial", {"host_type": "G", "host_mass_solar": 1.0, "age_gy": 0.05})
    assert dict(RECIPES.probabilities("p", infant)) == {"p5": 1.0}
    assert "replaces LU[0.162, 2]" in RECIPES.entry("p4").notes


def test_the_corrected_white_dwarf_redox_ensemble_replaces_the_six_star_centre():
    d4 = RECIPES._dist_of(RECIPES.entry("a4"), Subject("terrestrial", {}))
    d3 = RECIPES._dist_of(RECIPES.entry("a3"), Subject("terrestrial", {}))
    assert d4["mu"] == -1.25 and d4["sigma"] == 0.22 and d3["mu"] == -1.25 and d3["sigma"] == 0.22
    assert "Corrected 2026-10-18" in RECIPES.entry("a3").notes and "0.55 dex" in RECIPES.entry("a3").notes


# --- the sum rule: for every class the weights of the rows that pass add up to 100 percent -------------------------------

@pytest.mark.parametrize("group", sorted(RECIPES.groups))
def test_for_every_class_a_group_s_probabilities_add_up_to_one_or_the_group_does_not_apply(group):
    spec = RECIPES.groups[group]
    classes = sorted({c for r in spec.rows for c in RECIPES.entry(r).classes if c != "*"} or {"E"})
    for body_class in classes:
        for ratio in (0.2, 1.5):
            subject = full(body_class, snow_ratio=ratio, ts_k=700.0, **{f"drawn.{g}": None for g in RECIPES.groups})
            table = RECIPES.group_probabilities(group, subject)
            if table:
                assert math.fsum(p for _r, p in table) == pytest.approx(1.0, abs=1e-12)
                assert all(0.0 < p <= 1.0 for _r, p in table)


def test_a_group_with_no_supporting_row_for_the_class_is_an_error_that_names_the_class_and_the_table():
    with pytest.raises(pt.PossibilityError, match="atmosphere-g.*'terrestrial'"):
        RECIPES.draw_group("g", full("terrestrial"))


# --- draws -------------------------------------------------------------------------------------------------------------

def test_a_group_draw_follows_its_probabilities_and_stores_the_design_constants_inside_their_brackets():
    s = full("giant", mass_earth=5.0, core_radius_earth=1.54, insolation=30.0, v_esc_kms=20.0, snow_ratio=1.5)
    counts = {}
    with draw.bound(11):
        for _ in range(3000):
            picked = RECIPES.draw_group("e", s)
            counts[picked.name] = counts.get(picked.name, 0) + 1
            assert 0.05 <= picked.constants[("e", "f_env0")] <= 0.2 and 0.1 <= picked.constants[("e", "eta")] <= 0.5
    assert counts.get("e1", 0) > 0 and counts.get("e2", 0) > 0 and set(counts) <= {"e1", "e2", "e3", "e4"}
    given = full("giant", mass_earth=5.0, core_radius_earth=1.54, insolation=100.0, snow_ratio=1.5)
    with draw.bound(1):
        assert 0 < RECIPES.draw_group("e", given).probability <= 1.0


def test_a_continuous_row_draws_inside_its_range_and_its_bins_describe_the_draw():
    s = Subject("giant", dict(mass_jupiter=0.3, teq_k=950.0, zone="h"))
    with draw.bound(5):
        z = [RECIPES.draw_value("f5", s) for _ in range(4000)]
        c = [RECIPES.draw_value("f7", s) for _ in range(4000)]
    assert all(10 ** -0.5 <= v <= 100 for v in z) and all(0.2 <= v <= 1.5 for v in c)
    bins = RECIPES.bins("f6", s)
    for lo, hi, p in bins:
        got = sum(1 for v in z if lo <= v < hi) / len(z)
        assert got == pytest.approx(p, abs=0.03)
    above = sum(1 for v in c if v > 0.8) / len(c)
    assert above == pytest.approx(RECIPES.probability_in("f7", s, 0.8, None), abs=0.02)
    assert RECIPES.draw_value("a4", Subject("terrestrial", {"host_type": "G"})) is None          # its gate does not pass


def test_a_row_with_several_values_draws_a_dict_and_the_log_uniform_ranges_hold():
    s = Subject("comet", {})
    with draw.bound(8):
        got = [RECIPES.draw_value("k6", s) for _ in range(200)]
    assert set(got[0]) == {"CH3OH", "NH3", "HCN"}
    assert all(0.002 <= g["CH3OH"] <= 0.026 and 0.001 <= g["NH3"] <= 0.0073 and 0.0009 <= g["HCN"] <= 0.0029 for g in got)
    species = RECIPES.draw_value("k7", Subject("comet", dict(SOLAR_HOST)))
    assert set(species) == {"CH4", "C2H6", "H2CO", "H2S", "N2", "O2", "OCS", "SO2", "Ar", "Kr", "Xe"}


def test_a_binary_row_is_present_with_its_probability_and_carries_a_value_when_it_has_one():
    s = full("giant", t_layer_k=962.0)
    with draw.bound(4):
        hits = sum(RECIPES.draw_binary("m7", s)[0] for _ in range(4000))
    assert hits / 4000 == pytest.approx(0.5, abs=0.03)
    rogue = full("giant", zone="r")
    with draw.bound(4):
        present, fraction = RECIPES.draw_binary("b6", rogue)
        counts = sum(RECIPES.draw_binary("b6", rogue)[0] for _ in range(2000))
    assert counts / 2000 == pytest.approx(0.8, abs=0.04)
    assert RECIPES.draw_binary("b6", full("terrestrial", zone="e")) == (False, None)


def test_the_rotation_percentile_is_uniform_and_only_drawn_for_a_host_in_the_saturated_or_young_branch():
    s = full("terrestrial", **{"drawn.p": "p4"})
    with draw.bound(6):
        picks = [RECIPES.draw_value("p6", s) for _ in range(4000)]
    assert min(picks) >= 0.02 and max(picks) <= 0.98
    assert sum(picks) / 4000 == pytest.approx(0.5, abs=0.01)
    assert RECIPES.draw_value("p6", full("terrestrial", **{"drawn.p": "p1"})) is None


def test_draw_all_draws_the_tables_in_order_and_the_joint_probability_is_parent_times_conditional():
    s = full("terrestrial", redox_diw=0.5, ts_k=1000.0, pressure_bar=100.0, snow_ratio=1.5, insolation=1.0)
    with draw.bound(21):
        result = RECIPES.draw_all(s)
    assert list(result.groups)[:3] == ["a", "b", "d"]
    a, c = result.groups["a"], result.groups["c"]
    parent = result.subject.quantities["drawn.a"]
    assert parent == a.name and RECIPES.entry(c.name).conditions[0].quantity == "drawn.a"
    assert result.joint["c"] == pytest.approx(a.probability * c.probability)
    assert result.subject.quantities["drawn.d"] == result.groups["d"].name
    assert "a1" in result.values and "redox_diw" in result.subject.quantities
    assert "m10" in result.binaries and "o1" in result.values
    assert result.windows["n4"] == ["water cloud or ice"]


def _worker_draw(seed):
    s = full("terrestrial", redox_diw=0.5, ts_k=1000.0, pressure_bar=100.0, snow_ratio=1.5, insolation=3.0)
    with draw.bound(seed):
        r = RECIPES.draw_all(s)
    return ([(g, d.name, d.probability) for g, d in sorted(r.groups.items())],
            sorted((k, v if not isinstance(v, dict) else sorted(v.items())) for k, v in r.values.items()),
            sorted((k, v) for k, v in r.binaries.items()), sorted(r.choices.items()))


def test_a_draw_is_deterministic_from_its_seed_at_one_and_four_workers():
    seeds = list(range(8))
    serial = [_worker_draw(seed) for seed in seeds]
    assert serial == [_worker_draw(seed) for seed in seeds]
    assert serial[0] != serial[1]
    with concurrent.futures.ProcessPoolExecutor(max_workers=4, mp_context=multiprocessing.get_context("spawn")) as pool:
        assert list(pool.map(_worker_draw, seeds)) == serial


# --- the loader's refusals ----------------------------------------------------------------------------------------------

def test_a_row_without_a_source_or_with_a_state_it_does_not_know_is_refused():
    with pytest.raises(DataError, match="a1.*source"):
        broken(lambda d: row_of(d, "a1").update(source=""))
    with pytest.raises(DataError, match="state"):
        broken(lambda d: row_of(d, "a1").update(state="guessed"))
    with pytest.raises(DataError, match="a5.*partly verified and names no limit"):
        broken(lambda d: row_of(d, "a5").update(state="partly verified", limit=""))
    with pytest.raises(DataError, match="b5.*stand-in and names no replacement"):
        broken(lambda d: row_of(d, "b5").update(replacement=""))
    with pytest.raises(DataError, match="b1.*open and names no research item"):
        broken(lambda d: row_of(d, "b1").update(state="open", items=[]))


def test_a_class_a_kind_or_a_quantity_it_does_not_know_is_refused():
    with pytest.raises(DataError, match="unknown kinds Q9"):
        broken(lambda d: row_of(d, "a5").update(classes=["Q9"]))
    with pytest.raises(DataError, match="kind"):
        broken(lambda d: row_of(d, "a5").update(kind="guess"))
    with pytest.raises(DataError, match="no stored quantity"):
        broken(lambda d: row_of(d, "a1").update(conditions=[{"quantity": "mood", "minimum": 1}]))
    with pytest.raises(DataError, match="mood"):
        broken(lambda d: row_of(d, "a5").update(weight="phi(mood)"))


def test_an_expression_may_only_do_arithmetic_on_the_names_it_is_given():
    for text in ("__import__('os').system('true')", "(lambda: 1)()", "teq_k.real", "[x for x in range(3)]", "open('f')"):
        with pytest.raises(DataError):
            broken(lambda d, t=text: row_of(d, "a5").update(weight=t))
    with pytest.raises(DataError, match="cannot be read"):
        broken(lambda d: row_of(d, "a5").update(weight="phi((0 - redox_diw"))
    with pytest.raises(DataError, match="does not exist"):
        broken(lambda d: row_of(d, "a5").update(weight="w('zz9')"))


def test_weights_that_refer_to_each_other_in_a_circle_are_an_error_naming_the_circle():
    def change(d):
        row_of(d, "a5").update(weight="w('a6')")
        row_of(d, "a6").update(weight="w('a5')")
    other = broken(change)
    with pytest.raises(pt.TableError, match="circle"):
        other.probabilities("a", full("terrestrial"))


def test_a_weight_that_comes_out_negative_or_a_probability_above_one_is_an_error():
    other = broken(lambda d: row_of(d, "a5").update(weight="0 - 1"))
    with pytest.raises(pt.TableError, match="not a finite number of 0 or more"):
        other.probabilities("a", full("terrestrial"))
    other = broken(lambda d: row_of(d, "e6").update(weight="1.5"))
    with pytest.raises(pt.TableError, match="above 1"):
        other.probability("e6", full("giant", **{"drawn.e": "e2"}))
    other = broken(lambda d: row_of(d, "b1").update(lets={}, weight="log(snow_ratio)"))
    with pytest.raises(pt.TableError, match="cannot be evaluated"):
        other.probabilities("b", full("terrestrial", snow_ratio=0.0))


def test_groups_distributions_and_windows_that_are_broken_are_refused():
    with pytest.raises(DataError, match="fallback"):
        broken(lambda d: d["groups"]["d"].update(fallback="e1"))
    with pytest.raises(DataError, match="missing or belongs elsewhere"):
        broken(lambda d: d["groups"]["a"]["rows"].append("zz9"))
    with pytest.raises(DataError, match="distribution kind"):
        broken(lambda d: row_of(d, "a1")["dist"].update(kind="cauchy"))
    with pytest.raises(DataError, match="lacks sigma"):
        broken(lambda d: row_of(d, "a1")["dist"].pop("sigma"))
    with pytest.raises(DataError, match="low < high"):
        broken(lambda d: row_of(d, "n1")["windows"].append(["bad", 5, 4]))
    with pytest.raises(DataError, match="not a bracket"):
        broken(lambda d: d["groups"]["b"]["constants"].update(c_b4=[0.4, 0.05, "log"]))
    with pytest.raises(DataError, match="every row of group d"):
        broken(lambda d: row_of(d, "d6")["table"].pop("d5"))


def test_a_condition_on_a_drawn_parent_that_was_not_drawn_does_not_pass_and_is_no_error():
    entry = RECIPES.entry("c1")
    assert not entry.supports(Subject("terrestrial", {"ts_k": 1000.0, "pressure_bar": 100.0, "drawn.a": None}))
    assert not entry.supports(Subject("terrestrial", {"ts_k": 1000.0, "pressure_bar": 100.0}))
    assert entry.supports(Subject("terrestrial", {"ts_k": 1000.0, "pressure_bar": 100.0, "drawn.a": "a5"}))


def test_every_quantity_a_gate_or_weight_reads_is_documented_in_inputs():
    names = set()
    for entry in RECIPES.entries:
        names |= {c.quantity for c in entry.conditions if not c.quantity.startswith("drawn.")}
    assert names <= set(INPUTS)
    assert all(INPUTS[n].strip() for n in INPUTS)
    for entry in RECIPES.entries:
        if entry.weight is not None:
            allowed = set(INPUTS) | set(RECIPES.physical) | {"body_class"} | set(entry.constants) | set(entry.lets)
            if entry.group:
                allowed |= set(RECIPES.groups[entry.group].constants)
            assert entry.weight.names <= allowed, entry.id


def test_the_condition_one_of_of_the_possibility_table_passes_a_member_and_refuses_the_rest():
    c = pt.Condition("host_type", one_of=("F", "G"))
    assert c.holds(Subject("terrestrial", {"host_type": "G"})) and not c.holds(Subject("terrestrial", {"host_type": "M"}))
    assert pt.Table("t", 1, [pt.Row("a", ("E",), 1.0, "S 2020", conditions=(c,))]).digest() != pt.Table(
        "t", 1, [pt.Row("a", ("E",), 1.0, "S 2020", conditions=(pt.Condition("host_type", one_of=("F",)),))]).digest()
