"""
The fudge table (GEN.1313, scripts/build_fudge_table.py): the committed file is what the script builds from the matrix, the class table and the registry; every quantity the class fit gives a
shoulder has a width; no width passes its cap; a width is the root sum of squares of the movements its rows make, each widened by the factor of its tier.
"""

import json
import math
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_fudge_table as builder  # noqa: E402


@pytest.fixture(scope="module")
def table():
    with open(builder.OUT, encoding="utf-8") as handle:
        return json.load(handle)


def test_the_committed_table_is_what_the_script_builds():
    assert open(builder.OUT, encoding="utf-8").read() == builder.make(), "run python scripts/build_fudge_table.py"


def test_the_build_is_reproducible():
    assert builder.make() == builder.make()


def test_every_shoulder_quantity_has_a_positive_width(table):
    from planetgen.physics import class_fit

    assert set(table["quantities"]) == set(class_fit.SHOULDERS)
    for quantity, entry in table["quantities"].items():
        assert entry["width"] > 0.0 or not entry["edges"], quantity        # a quantity no class edge reads has no width to speak of
        assert entry["basis"] in ("measured", "fallback"), quantity
        assert (entry["tier"], entry["tweak"]) == (builder.TIER, builder.TWEAK)
        assert entry["unit"] == ("ln" if quantity in builder.LOG_QUANTITIES else "plain")


def test_a_measured_width_is_two_sigma_of_the_widened_movements(table):
    for quantity, entry in table["quantities"].items():
        if entry["basis"] != "measured":
            assert entry["rows"] == 0, quantity
            assert entry["sigma"] == 0.0, quantity
            continue
        assert entry["rows"] > 0, quantity
        assert entry["width"] == pytest.approx(builder.SIGMA_MULTIPLE * entry["sigma"], rel=1e-4, abs=1e-5), quantity
        assert sum(item["share"] for item in entry["contributors"]) <= 1.0 + 1e-3, quantity
        assert [item["share"] for item in entry["contributors"]] == sorted((item["share"] for item in entry["contributors"]), reverse=True), quantity


def test_no_width_in_use_passes_its_cap(table):
    for quantity, entry in table["quantities"].items():
        for edge in entry["edges"]:
            if edge["cap"] is None:
                assert edge["width"] == entry["width"]
                continue
            assert edge["width"] <= edge["cap"] + 1e-9, (quantity, edge)
            assert (edge["width"] == edge["cap"]) or not edge["binds"], (quantity, edge)        # a width that was cut is the cap


def test_tier_factor_widens_an_argued_range():
    assert builder.tier_factor("T1") == 1.0
    assert builder.tier_factor("T3") == 1.0
    assert builder.tier_factor("T4") == 1.5
    assert builder.tier_factor("T8") == 1.5
    assert builder.tier_factor("T9") == 2.0
    assert builder.tier_factor(None) == 2.0


def test_build_sums_independent_movements_in_quadrature():
    class Cls:
        def __init__(self, when):
            self.when = when

    classes = [Cls([{"quantity": "temperature_k", "min": 250.0}]), Cls([{"quantity": "temperature_k", "lt": 330.0}])]
    shoulders = {"temperature_k": {"absolute": 10.0}}
    matrix = {"sample": {}, "rows": {"A1": {"quantities": {"temperature_k": {"rms": 3.0}}}, "B1": {"quantities": {"temperature_k": {"rms": 4.0}}}}}
    rows = [{"fill_id": "A1", "tier": "T2", "tweak": "L0"}, {"fill_id": "B1", "tier": "T9", "tweak": "L3"}]
    entry = builder.build(matrix, rows, classes, shoulders)["quantities"]["temperature_k"]
    assert entry["sigma"] == pytest.approx(math.hypot(3.0, 8.0), rel=1e-5)
    assert entry["width"] == pytest.approx(2.0 * math.hypot(3.0, 8.0), rel=1e-5)
    assert entry["contributors"][0]["fill_id"] == "B1"
    assert entry["tuned_share"] == pytest.approx(64.0 / 73.0, abs=1e-3)
    assert all(edge["width"] <= 40.0 for edge in entry["edges"])        # half the 80 K gap between the two edges


def test_a_movement_past_the_half_gap_counts_as_the_half_gap():
    cell = {"rms": 50.0, "bodies": 10, "steps": [100.0, 2.0]}
    assert builder.clipped_rms(cell, None) == 50.0
    assert builder.clipped_rms(cell, 4.0) == pytest.approx(math.sqrt((16.0 + 4.0) / 10.0))
    assert builder.clipped_rms({"rms": 7.0}, 4.0) == 7.0
