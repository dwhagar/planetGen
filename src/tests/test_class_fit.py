"""
The best-fit class table and scorer (GEN.1205, GEN.1206): `physics/data/planet_class_fit.json` gives every class its own positive
conditions and `physics/class_fit.py` scores a body against them. Nothing stored changes: `planet_classes.classify` still decides.
"""

import json
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


def test_every_quantity_a_condition_reads_has_a_shoulder_and_a_registry_basis():
    for cls in cf.CLASSES:
        for condition in cls.when:
            if "equals" not in condition:
                assert condition["quantity"] in cf.SHOULDERS
    for name, entry in cf.SHOULDERS.items():
        assert entry["basis"].startswith("CF")
        assert ("absolute" in entry) != ("relative" in entry)


def test_the_shoulder_basis_rows_are_registered():
    text = open(os.path.join(os.path.dirname(__file__), "..", "..", "docs", "design", "research", "gap-filling-plan.md"), encoding="utf-8").read()
    for entry in cf.SHOULDERS.values():
        assert f"| {entry['basis']} |" in text


def test_membership_is_a_straight_line_across_the_shoulder():
    edge_min = {"quantity": "temperature_k", "min": 273.15}
    assert cf.membership(edge_min, {"temperature_k": 273.15}) == 1.0
    assert cf.membership(edge_min, {"temperature_k": 268.15}) == pytest.approx(0.5)
    assert cf.membership(edge_min, {"temperature_k": 263.15}) == 0.0
    edge_lt = {"quantity": "temperature_k", "lt": 320.0}
    assert cf.membership(edge_lt, {"temperature_k": 319.0}) == 1.0
    assert cf.membership(edge_lt, {"temperature_k": 325.0}) == pytest.approx(0.5)
    assert cf.membership(edge_lt, {"temperature_k": 330.0}) == 0.0
    assert cf.membership(edge_lt, {"temperature_k": None}) == 0.0


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


def test_a_body_no_class_matches_has_no_best_letter():
    hot_thin = body(surface_temperature_k=414.0, atmospheric_pressure_pa=1000.0, ocean_fraction=0.0, hydrosphere="vapour", gravity_g=1.4)
    assert cf.best_class(hot_thin)[:2] == (None, 0.0)


def test_a_cold_hydrogen_envelope_is_a_gas_dwarf_and_not_a_deep_frozen_world():
    values = body(surface_temperature_k=72.0, mass_earth=0.7, share_h2=0.98, share_o2=0.0, hydrosphere="ice", ocean_fraction=1.0,
                  p_o2_kpa=0.0, p_n2_kpa=0.0, p_ar_kpa=0.0, p_co2_kpa=0.0)
    put = dict(values)
    assert pc.classify(values) == "T" and cf.best_class(values)[0] == "T"


def test_one_condition_missed_outright_zeroes_the_match_and_a_soft_miss_only_lowers_it():
    earth = cf.scores(EARTH_VALUES := body())
    assert earth["M"] == pytest.approx(1.0)
    warm = cf.scores(body(surface_temperature_k=271.0))
    assert 0.0 < warm["M"] < 1.0
    cold = cf.scores(body(surface_temperature_k=250.0))
    assert cold["M"] == 0.0


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
