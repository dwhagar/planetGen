# tests/test_class_census_script.py

"""
TEST.719: `scripts/class_census.py` draws the class census sample and prints its tables in one markdown format. The table maths (the funnel, the crosstab, the age bins)
is held on a small stored sample, the draw on one real system, and the sample lists on their fixed seeds, so a run of one train compares line by line with the next.
"""

import importlib.util
import json
import pathlib

import pytest

SCRIPT = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "class_census.py"


@pytest.fixture(scope="module")
def census():
    spec = importlib.util.spec_from_file_location("class_census", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _body(klass, **found):
    quantities = {"body_type": "t", "pressure_pa": 1.0e5, "gravity_g": 1.0, "temperature_k": 288.0, "phase_water": 1, "ocean_fraction": 0.7,
                  "o2_partial_pa": 2.1e4, "co2_partial_pa": 40.0, "reducing_share": 0.0, "share_h2": 0.0, "share_so2": 0.0, "hottest_region_k": 300.0,
                  "coldest_region_k": 260.0, "night_temperature_k": 270.0, "wet_bulb_c": 20.0, **found}
    return {"planet_class": klass, "is_moon": False, "host": 0, "zone": "h", "name": klass, "distance": 1.0, "planet_type": "terrestrial", "found": quantities,
            "family_conditions": {"M": [["pressure_pa", "min", 5.0e4, 1.0], ["temperature_k", "lt", 305.0, 1.0 if quantities["temperature_k"] < 305 else 0.0]]},
            "fit": [klass, "T", 0.9, 0.5, 1], "life_stage": "photosynthesis"}


def _write(directory, records):
    with open(directory / "sample.jsonl", "w") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")


def _record(tag, seed, bodies, age=4.6, disk=None):
    star = {"type": "G2V", "luminosity": 3.828e26, "age": age, "disk_kind": disk, "mass": 2.0e30}
    return {"tag": tag, "seed": seed, "stars": [star], "belts": [{"mass_kg": 1.0e23, "largest_body_km": 500.0}], "bodies": bodies}


def test_the_report_prints_every_table_in_one_markdown_format(census, tmp_path, capsys):
    _write(tmp_path, [_record("natural", 1, [_body("M"), _body("T", temperature_k=700.0, pressure_pa=5.0)]),
                      _record("natural", 2, [_body("E", gravity_g=0.1)], age=0.02),
                      _record("hab-G2V", 3, [_body("M")], age=0.001, disk="class-ii")])
    census.report(type("Args", (), {"directory": str(tmp_path), "caps_off_dir": None})())
    out = capsys.readouterr().out
    for title in ("A. Samples", "B. Class counts", "C. Natural sample", "C2.", "C3.", "C4.", "D. M family", "E. Funnel", "E2.", "E3.", "E4.", "F. Rocky",
                  "G. Pressure pile-ups", "H. Oxygen chain", "I. Stars", "J. Other value checks"):
        assert f"### {title}" in out
    assert "| set | systems | planets | moons | belts | errors |" in out
    assert "| natural | 2 | 3 | 0 | 2 | 0 |" in out
    assert "total systems 3 bodies 4" in out
    assert "| (0, 3] Myr | 1 | 1 | 1 |" in out


def test_the_funnel_counts_a_body_at_the_first_stage_it_fails(census):
    good, thin, hot = ({"q_" + k: v for k, v in _body(c, **extra)["found"].items()} | {"planet_class": c}
                       for c, extra in (("M", {}), ("T", {"pressure_pa": 1.0}), ("E", {"temperature_k": 400.0})))
    counts = census.funnel([good, thin, hot])
    stages = [stage for stage, _ in census.FUNNEL_STAGES]
    assert counts[stages.index("all planets and moons")] == 3
    assert counts[stages.index("air at 10 Pa or more")] == 2
    assert counts[stages.index("mean temperature 258 to 343 K (liquid-water span)")] == 1
    assert counts[-1] == 1


def test_the_crosstab_drops_a_body_with_no_class_and_gives_row_shares(census):
    items = [("a", "T"), ("a", "T"), ("a", "C"), ("b", None), ("b", "C")]
    table = census.crosstab(items, lambda i: i[0], lambda i: i[1])
    assert table.splitlines()[0] == "|  | C | T |" and "| a | 1 | 2 |" in table and "| b | 1 | 0 |" in table
    shares = census.crosstab(items, lambda i: i[0], lambda i: i[1], percent=True)
    assert "| a | 33.3 | 66.7 |" in shares


def test_the_spectral_letter_and_the_age_bins(census):
    assert [census.spectral(t) for t in ("G2V", "M5.5V", "DA white dwarf", None)] == ["G", "M", "WD", "?"]
    assert census.cut(3.0, [0, 3, 10]) == 0 and census.cut(3.1, [0, 3, 10]) == 1 and census.cut(0.0, [0, 3, 10]) is None and census.cut(11.0, [0, 3, 10]) is None


def test_the_samples_have_fixed_seeds_and_the_names_the_command_takes(census):
    assert [len(census.jobs_for(name)) for name in census.SETS] == [600, 160, 120, 40, 100, 190]
    assert census.jobs_for("natural")[0] == ("natural", 1000, {}) and census.jobs_for("youngforced")[-1][1] == 20099
    assert census.jobs_for("hab")[0][2] == {"STAR_TYPE": "G2V", "HABITABLE_WORLD": True}
    with pytest.raises(SystemExit):
        census.jobs_for("nothing")


def test_one_system_is_drawn_into_a_record_the_report_reads(census, tmp_path, capsys):
    record = census.draw_one(("natural", 1001, {"STAR_TYPE": "G2V"}))
    assert "error" not in record and record["stars"] and record["bodies"]
    body = record["bodies"][0]
    assert body["fit"][0] == body["planet_class"] and "found" in body and "family_conditions" in body
    _write(tmp_path, [json.loads(json.dumps(record))])
    census.report(type("Args", (), {"directory": str(tmp_path), "caps_off_dir": None})())
    assert "total systems 1 bodies" in capsys.readouterr().out
