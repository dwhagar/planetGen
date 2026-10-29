# tests/test_prevalence_occurrence.py

"""
TEST.395: planet occurrence against observation (part E of the prevalence suite): the figures the generator's planets give per host
group (giants, hot Jupiters, cold giants, 1 to 20 Earth-mass planets inside 0.5 AU, Earth-mass planets in the ecosphere, the
super-Earth to Neptune ratio) and the observed ranges they are checked against. The live check shares the system census's draw
(`test_prevalence_systems_census.py`).
"""

import re

import pytest

from tests.audit import prevalence_occurrence as occurrence

ITEM = re.compile(r"^[A-Z]+\.\d+$")


def _system(star, planets, hz=(1.0, 1.7)):
    return {"star": star, "orbits": [list(p) for p in planets], "hz": list(hz), "planets": len(planets)}


def test_the_figures_count_planets_per_star_by_mass_and_orbit():
    records = [
        _system("G2V", [(300.0, 0.05, "J"), (400.0, 3.0, "J"), (5.0, 0.3, "S"), (1.0, 1.2, "M")]),
        _system("K5V", [(4.0, 0.2, "S"), (30.0, 0.4, "I"), (12.0, 20.0, "I")]),
        _system("F5V", []),
        _system("M2V", [(150.0, 0.04, "J")]),
        {"star": "M2V", "abandoned": True},
    ]
    kgf = occurrence.summarize(records)["KGF"]["values"]
    assert kgf["giants"] == pytest.approx(2 / 3)
    assert kgf["hot_jupiters"] == pytest.approx(1 / 3)
    assert kgf["cold_giants"] == pytest.approx(1 / 3)
    assert kgf["inner_planets"] == pytest.approx(2 / 3)
    assert kgf["earth_mass_in_hz"] == pytest.approx(1 / 3)
    assert kgf["super_earths_per_neptune"] == pytest.approx(2 / 2)
    m = occurrence.summarize(records)["M"]
    assert (m["n"], m["abandoned"], m["values"]["giants"], m["values"]["hot_jupiters"]) == (1, 1, 1.0, 1.0)


def test_a_rate_the_generator_misses_is_flagged_and_a_rate_it_meets_is_not():
    expectations = occurrence.build_expectations(100)
    records = [_system("M2V", [(300.0, 0.05, "J")] * 3)] * 10
    findings = occurrence.compare(occurrence.summarize(records), expectations)
    assert {(f.cell, f.measure, f.kind) for f in findings} >= {("M", "giants", "out-of-range"), ("M", "hot_jupiters", "out-of-range")}
    sun_like = ([_system("G2V", [(5.0, 0.3, "S"), (1.0, 1.2, "M")])] * 30       # an inner super-Earth, an Earth in the ecosphere
                + [_system("G2V", [(5.0, 0.3, "S"), (200.0, 3.0, "J")])] * 20  # an inner super-Earth and a cold giant
                + [_system("G2V", [(30.0, 0.8, "I")])] * 10                     # a Neptune
                + [_system("G2V", [])] * 40)
    quiet = [_system("M2V", [(5.0, 1.0, "S")])] * 100 + sun_like
    assert occurrence.compare(occurrence.summarize(quiet), expectations) == []


def test_every_observed_range_has_a_kind_a_source_and_a_known_measure():
    expectations = occurrence.build_expectations(100)
    for group, spec in expectations["types"].items():
        assert group in occurrence.GROUPS
        for name, bounds in spec["measures"].items():
            assert name in occurrence.MEASURES and bounds["kind"] == "observed" and bounds["source"], (group, name)
            assert 0.0 <= bounds["min"] < bounds["max"], (group, name)


def test_the_stored_file_holds_the_same_observed_ranges_and_names_its_gaps():
    stored = occurrence.load_expectations()
    built = occurrence.build_expectations(stored["reference"]["systems_per_type"])
    for group, spec in built["types"].items():
        for name, bounds in spec["measures"].items():
            kept = stored["types"][group]["measures"][name]
            assert {k: kept[k] for k in ("min", "max", "kind", "source")} == bounds, (group, name)
            if "anomaly" in kept:
                assert ITEM.match(kept["anomaly"]) and kept["anomaly"] in stored["anomalies"], (group, name)
