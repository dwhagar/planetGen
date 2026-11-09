"""
The stand-in sensitivity harness (GEN.1311, scripts/stand_in_sensitivity.py): reads a range from a registry cell, finds the constant a row names, overrides it for the length of a
block and puts it back, and sorts every registry row into overridable or not with a reason. No system is generated here.
"""

import math
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import stand_in_sensitivity as harness  # noqa: E402


@pytest.mark.parametrize("knobs, value, expected", [
    ("0.05 to 0.2", 0.1, (0.05, 0.2)),
    ("0.05 to 0.2 [R12]; the form", 0.1, (0.05, 0.2)),
    ("1.5 K to 3 K", 2.0, (1.5, 3.0)),
    ("a factor of 2", 10.0, (5.0, 20.0)),
    ("none", 1.0, None),
    ("1.0 Gyr, fixed", 1.0, None),
])
def test_parse_range(knobs, value, expected):
    assert harness.parse_range(knobs, value) == expected


def test_override_sets_and_restores_a_module_constant():
    from planetgen.physics import class_fit
    original = class_fit.HOLE_FLOOR
    with harness.override("planetgen.physics.class_fit", "HOLE_FLOOR", 0.5):
        assert class_fit.HOLE_FLOOR == 0.5
    assert class_fit.HOLE_FLOOR == original


def test_override_keeps_an_integer_constant_an_integer():
    from planetgen.physics import orbit_elements
    name = next(key for key, value in vars(orbit_elements).items() if key.isupper() and type(value) is int)
    original = getattr(orbit_elements, name)
    with harness.override("planetgen.physics.orbit_elements", name, float(original) * 2):
        assert type(getattr(orbit_elements, name)) is int
    assert getattr(orbit_elements, name) == original


def test_every_registry_row_is_sorted_with_a_reason():
    rows = harness.classify_rows(harness.read_rows())
    assert rows
    assert all(("target" in row) != ("reason" in row) for row in rows)
    overridable = [row for row in rows if "target" in row]
    assert overridable
    for row in overridable:
        target = row["target"]
        assert target["low"] < target["high"]
        assert isinstance(target["value"], (int, float))


def test_movement_is_a_log_ratio_for_a_pressure_and_a_difference_for_a_share():
    assert harness.movement("pressure_pa", 100.0, 1000.0) == pytest.approx(2.302585, rel=1e-6)
    assert harness.movement("share_n2", 0.2, 0.25) == pytest.approx(0.05)


def test_compare_pairs_bodies_and_counts_a_changed_structure():
    base = [{(0, -1): {"temperature_k": 280.0, "gravity_g": 1.0}}, {(0, -1): {"temperature_k": 300.0}}]
    moved = [{(0, -1): {"temperature_k": 290.0, "gravity_g": 1.0}}, {(0, -1): {"temperature_k": 300.0}, (1, -1): {"temperature_k": 100.0}}]
    totals, structure = harness.compare(base, moved)
    assert structure == 1
    summary = harness.summarize(totals)
    assert summary["temperature_k"]["rms"] == pytest.approx(10.0)
    assert summary["gravity_g"]["moved"] == 0.0


def test_combine_pools_the_passes_and_counts_a_resumed_row_once():
    def row(fill_id, rms, moved, structure=0):
        return {"fill_id": fill_id, "low": {"quantities": {"temperature_k": {"rms": rms, "moved": moved}}, "structure": structure},
                "high": {"quantities": {"temperature_k": {"rms": rms / 2, "moved": moved}}, "structure": 0}}
    first = {"results": [row("A1", 3.0, 0.5), row("A1", 3.0, 0.5), row("B1", 1.0, 0.1, structure=2)]}
    second = {"results": [row("A1", 4.0, 0.7)]}
    matrix = harness.combine([first, second])
    assert set(matrix) == {"A1", "B1"}
    assert matrix["A1"]["quantities"]["temperature_k"]["rms"] == pytest.approx(math.sqrt((9.0 + 16.0) / 2))
    assert matrix["A1"]["quantities"]["temperature_k"]["passes"] == 2 and matrix["A1"]["passes"] == 2
    assert matrix["B1"]["structure"] == 2
