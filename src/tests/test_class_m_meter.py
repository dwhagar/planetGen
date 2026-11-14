# tests/test_class_m_meter.py

"""
`scripts/class_m_meter.py`: the funnel and the per-1,000 table on a small fixed file give a fixed report, the `sunlike` sample list is fixed, and one real
system draws under a tier (Class M abundance tiers, thread 5).
"""

import importlib.util
import json
import pathlib
import sys

import pytest

SCRIPTS = pathlib.Path(__file__).resolve().parents[2] / "scripts"


@pytest.fixture(scope="module")
def meter():
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location("class_m_meter", SCRIPTS / "class_m_meter.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _body(**fields):
    body = {"body_type": "t", "is_moon": False, "zone": "e", "mass": 5.972e24, "hydrosphere": "surface ocean", "gravity": 1.0, "ocean_fraction": 0.7,
            "planet_class": "T", "overall_tier": "Red"}
    body.update(fields)
    return body


def _records():
    return [
        {"tag": "sunlike", "seed": 30000, "bodies": [_body(planet_class="M", overall_tier="Blue"), _body(zone="h"), _body(mass=1.0e22)]},
        {"tag": "sunlike", "seed": 30001, "bodies": [_body(hydrosphere=None), _body(planet_class="O", overall_tier="Green"), _body(is_moon=True, planet_class="M")]},
        {"tag": "sunlike", "seed": 30002, "bodies": []},
        {"tag": "sunlike", "seed": 30003, "error": "boom", "bodies": [_body(planet_class="M")]},
    ]


def test_the_funnel_counts_a_candidate_down_to_the_first_row_it_fails(meter):
    records = [r for r in _records() if "error" not in r]
    numbers = meter.measure(records)
    assert numbers["systems"] == 3
    # candidates: seed 30000 has 2 (the zone h body is not one), seed 30001 has 2 (the moon is not one); the light body fails mass, the bare one fails ground
    assert numbers["funnel"] == [4, 3, 2, 2, 2, 2, 1, 1, 1]
    assert numbers["class_m_worlds"] == 2 and numbers["class_m_systems"] == 2
    assert numbers["blue_worlds"] == 1 and numbers["blue_systems"] == 1


def test_the_report_of_a_small_file_is_fixed(meter, tmp_path, capsys):
    with open(tmp_path / "science.jsonl", "w") as handle:
        for record in _records():
            handle.write(json.dumps(record) + "\n")
    meter.report(type("Args", (), {"directory": str(tmp_path)})())
    out = capsys.readouterr().out
    assert "| science | 3 | 4 | 666.7 (" in out and " | none (what the laws give) |" in out
    assert "### science: 3 Sun-like systems" in out
    assert "| Class M" in out and "| Blue" in out
    assert out == _again(meter, tmp_path, capsys)


def _again(meter, tmp_path, capsys):
    meter.report(type("Args", (), {"directory": str(tmp_path)})())
    return capsys.readouterr().out


def test_the_wilson_interval_brackets_the_rate(meter):
    low, high = meter.wilson(10, 100)
    assert low < 0.1 < high
    assert meter.wilson(0, 0) == (0.0, 0.0)
    assert meter.wilson(0, 400)[0] == 0.0


def test_the_sunlike_sample_is_fixed(meter):
    census = sys.modules["class_census"]
    jobs = census.jobs_for("sunlike", 25)
    assert len(jobs) == 25
    assert jobs[0] == ("sunlike", 30000, {"STAR_TYPE": "F7V"})
    assert jobs[10] == ("sunlike", 30010, {"STAR_TYPE": "F7V"})
    assert jobs[9][2]["STAR_TYPE"] == "K2V"
    assert len(census.jobs_for("sunlike")) == 400


def test_a_system_draws_under_a_tier(meter):
    census = sys.modules["class_census"]
    census.CLASS_M = "abundant"
    try:
        record = census.draw_one(("sunlike", 30003, {"STAR_TYPE": "G2V"}))
    finally:
        census.CLASS_M = None
    assert "error" not in record
    assert record["class_m"] == "abundant"
    assert all("overall_tier" in body for body in record["bodies"])
