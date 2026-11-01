# tests/test_prevalence_occurrence.py

"""
TEST.395: planet occurrence against observation (part E of the prevalence suite): the figures the generator's planets give per host
group (Kepler-like systems, giants, hot Jupiters, habitable-zone Earths, spacing) and the observed ranges they are checked against,
which come from the sourced bands of `planet_occurrence`. The live check shares the system census's draw
(`test_prevalence_systems_census.py`).
"""

import re

import pytest

from planetgen.physics import planet_occurrence as po
from tests.audit import prevalence_occurrence as occurrence

ITEM = re.compile(r"^[A-Z]+\.\d+$")


def _system(star, planets, hz=(0.9, 1.6), star_mass=1.0):
    """`planets` are `(mass in Earth masses, distance in AU, radius in Earth radii)`."""
    return {"star": star, "star_mass": star_mass, "orbits": [[m, a, "S", r] for m, a, r in planets], "hz": list(hz), "planets": len(planets)}


def test_the_figures_count_planets_per_star_by_radius_and_period():
    records = [
        _system("G2V", [(5.0, 0.1, 1.8), (3.0, 0.2, 1.5), (300.0, 0.04, 11.0), (300.0, 4.0, 11.0), (1.0, 1.2, 1.0)]),
        _system("G2V", [(2.0, 0.3, 1.3)]),
        _system("G2V", []),
        _system("G2V", [(1.0, 20.0, 1.0)]),
        {"star": "G2V", "abandoned": True},
        _system("M2V", [(2.0, 0.05, 1.3), (3.0, 0.07, 1.6)], hz=(0.1, 0.2), star_mass=0.4),
    ]
    summary = occurrence.summarize(records)
    g = summary["G"]["values"]
    assert summary["G"]["n"] == 4 and summary["G"]["abandoned"] == 1
    assert g["inner_fraction"] == pytest.approx(2 / 4)          # the 20 AU Earth has a period of 89 years
    assert g["inner_multiplicity"] == pytest.approx((3 + 1) / 2)       # the hot Jupiter counts, the 4 AU giant and the 1.2 AU Earth do not
    assert g["planets_per_star"] == pytest.approx(4 / 4)
    assert g["small_planets_per_star"] == pytest.approx(3 / 4)
    assert g["giants_2_8_au"] == pytest.approx(1 / 4)
    assert g["hot_jupiter_rate"] == pytest.approx(1 / 4)
    assert g["earth_in_hz"] == pytest.approx(1 / 4)
    m = summary["M"]["values"]
    assert m["inner_fraction"] == 1.0 and m["small_planets_per_star"] == 2.0 and m["hot_jupiter_rate"] == 0.0


def test_the_spacing_is_counted_in_mutual_hill_radii_between_neighbours():
    wide = _system("G2V", [(1.0, 1.0, 1.0), (1.0, 2.0, 1.0)])
    tight = _system("G2V", [(1.0, 1.0, 1.0), (1.0, 1.02, 1.0)])
    assert occurrence.summarize([wide])["G"]["values"]["pairs_over_10_hill_radii"] == 1.0
    assert occurrence.summarize([tight])["G"]["values"]["pairs_over_10_hill_radii"] == 0.0
    assert "pairs_over_10_hill_radii" not in occurrence.summarize([_system("G2V", [(1.0, 1.0, 1.0)])])["G"]["values"]


def _sun_like_sky(inner_hosts, planets_each, giants=0.0):
    """100 Sun-like systems: `inner_hosts` of them with `planets_each` planets of 1 to 4 Earth radii in 10 to 100 days, none else."""
    records = []
    for i in range(100):
        planets = []
        if i < inner_hosts:
            planets = [(3.0, 0.1 * 1.8 ** j, 1.5) for j in range(planets_each)]
        if i < giants * 100:
            planets.append((300.0, 4.0, 11.0))
        records.append(_system("G2V", planets))
    return records


def test_a_rate_the_generator_misses_is_flagged_and_a_rate_it_meets_is_not():
    expectations = occurrence.build_expectations(100)
    sparse = _sun_like_sky(inner_hosts=60, planets_each=1)
    flags = {(f.cell, f.measure, f.kind) for f in occurrence.compare(occurrence.summarize(sparse), expectations)}
    assert ("G", "inner_fraction", "out-of-range") in flags
    assert ("G", "inner_multiplicity", "out-of-range") in flags
    met = _sun_like_sky(inner_hosts=31, planets_each=3, giants=0.14)
    found = {f.measure for f in occurrence.compare(occurrence.summarize(met), expectations) if f.cell == "G"}
    assert not {"inner_fraction", "inner_multiplicity", "giants_2_8_au"} & found


def test_every_observed_range_has_a_kind_a_source_and_a_known_measure():
    expectations = occurrence.build_expectations(100)
    for group, spec in expectations["types"].items():
        assert group in occurrence.GROUPS
        for name, bounds in spec["measures"].items():
            assert name in occurrence.MEASURES and bounds["kind"] == "observed" and "[R" in bounds["source"], (group, name)
            assert 0.0 <= bounds["min"] < bounds["max"], (group, name)
            if name in occurrence.SHARE_MEASURES:
                assert bounds["max"] <= 1.0


def test_a_larger_sample_narrows_the_range_to_the_published_interval():
    small, large = occurrence.build_expectations(100), occurrence.build_expectations(10000)
    for group, spec in small["types"].items():
        for name, bounds in spec["measures"].items():
            wide = bounds["max"] - bounds["min"]
            narrow = large["types"][group]["measures"][name]["max"] - large["types"][group]["measures"][name]["min"]
            assert narrow < wide, (group, name)


def test_every_band_a_measure_reads_exists_and_the_gaps_name_items():
    for measure in occurrence.MEASURES:
        assert po.bands_of(occurrence.BAND_MEASURES[measure]), measure
        item, note = occurrence.GAP_ITEMS[measure]
        assert ITEM.match(item) and note


def test_the_stored_file_holds_the_same_observed_ranges_and_names_its_gaps():
    stored = occurrence.load_expectations()
    built = occurrence.build_expectations(stored["reference"]["systems_per_type"])
    assert set(stored["types"]) == set(built["types"])
    for group, spec in built["types"].items():
        assert set(stored["types"][group]["measures"]) == set(spec["measures"])
        for name, bounds in spec["measures"].items():
            kept = stored["types"][group]["measures"][name]
            assert {k: kept[k] for k in ("min", "max", "kind", "source")} == bounds, (group, name)
            if "anomaly" in kept:
                assert ITEM.match(kept["anomaly"]) and kept["anomaly"] in stored["anomalies"], (group, name)
