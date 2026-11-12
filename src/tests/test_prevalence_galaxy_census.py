# tests/test_prevalence_galaxy_census.py

"""
TEST.395, the galaxy-placed repeat: the planets of the stars a seeded galaxy places round the reference solar position, against the same
observed bands as the system draw's occurrence census (`test_prevalence_systems_census.py`), with the host mix the galaxy's own. The
live census draws sectors and fails on a flag no known gap explains; the rest of the file holds the harness to its rules on small
hand-made records.

Classes O (ocean world), P (glaciated world) and Q are given no target rate here (Boss, 2026-10-30: physics drives, no game rates, and no
survey measures a class): they are only tested for being possible when conditions allow them (the class sweep), and the report lists how
many of each the sample drew.
"""

import re

import pytest

from tests.audit import prevalence as census
from tests.audit import prevalence_galaxy as galaxy
from tests.audit import prevalence_occurrence as occurrence

ITEM = re.compile(r"^[A-Z]+\.\d+$")


def _system(group, planets, hz=(0.9, 1.6), star_mass=None, classes=None):
    """`planets` are `(mass in Earth masses, distance in AU, radius in Earth radii)`."""
    mass = star_mass or {"M": 0.4, "K": 0.7, "G": 1.0, "F": 1.3}[group]
    return {"star": galaxy.REPRESENTATIVE[group], "group": group, "star_mass": mass, "hz": list(hz), "planets": len(planets),
            "orbits": [[m, a, "S", r] for m, a, r in planets], "classes": classes or {"S": len(planets)}, "life_bodies": 0}


@pytest.fixture(scope="module")
def stored():
    return occurrence.load_expectations()


def _sky(group, n, inner_hosts, planets_each):
    """`n` systems of `group`, `inner_hosts` of them with `planets_each` planets of 1 to 4 Earth radii in 10 to 100 days."""
    return [_system(group, [(3.0, 0.1 * 1.8 ** j, 1.5) for j in range(planets_each)] if i < inner_hosts else []) for i in range(n)]


def test_a_star_is_filed_under_the_group_of_its_mass():
    assert galaxy.REPRESENTATIVE == {"M": "M2V", "K": "K5V", "G": "G2V", "F": "F5V"}
    assert [occurrence.po.host_group(m) for m in (0.4, 0.7, 1.0, 1.3, 2.5)] == ["M", "K", "G", "F", None]


def test_a_group_under_the_minimum_is_a_cell_only_in_the_pool(stored):
    records = _sky("M", galaxy.MIN_GROUP, 20, 3) + _sky("G", galaxy.MIN_GROUP - 1, 10, 3)
    summary = galaxy.summarize(records)
    expectations = galaxy.expectations_for(summary, stored)
    assert set(expectations["types"]) == {"M", galaxy.POOL}
    assert summary["G"]["n"] == galaxy.MIN_GROUP - 1 and summary[galaxy.POOL]["n"] == 2 * galaxy.MIN_GROUP - 1


def test_the_pool_is_the_stars_of_the_groups_that_have_a_band_for_the_measure(stored):
    records = _sky("M", 40, 30, 3) + _sky("G", 10, 3, 1)
    summary = galaxy.summarize(records)
    pooled = summary[galaxy.POOL]["values"]
    weights = summary[galaxy.POOL]["weights"]
    # the M dwarf band for 1 to 4 Earth radii under 200 days is the only one of its kind, so only M stars are pooled for it
    assert weights["small_planets_per_star"] == {"M": 40}
    assert pooled["small_planets_per_star"] == pytest.approx(summary["M"]["values"]["small_planets_per_star"])
    # the Kepler-like share has a band for both, so both are pooled by their stars
    assert set(weights["inner_fraction"]) == {"M", "G"}
    mixed = (summary["M"]["values"]["inner_fraction"] * 40 + summary["G"]["values"]["inner_fraction"] * 10) / 50
    assert pooled["inner_fraction"] == pytest.approx(mixed)


def test_the_pooled_range_is_the_groups_ranges_averaged_by_their_stars(stored):
    summary = galaxy.summarize(_sky("M", 60, 40, 2) + _sky("G", 20, 8, 2))
    expectations = galaxy.expectations_for(summary, stored)
    pooled = expectations["types"][galaxy.POOL]["measures"]["inner_fraction"]
    m = occurrence.build_expectations(60)["types"]["M"]["measures"]["inner_fraction"]
    g = occurrence.build_expectations(20)["types"]["G"]["measures"]["inner_fraction"]
    assert pooled["min"] == pytest.approx((m["min"] * 60 + g["min"] * 20) / 80, abs=1e-3)
    assert pooled["max"] == pytest.approx((m["max"] * 60 + g["max"] * 20) / 80, abs=1e-3)


def test_a_flag_names_the_cell_and_the_measure_that_left_its_band(stored):
    records = _sky("M", 60, 3, 1)          # 5 percent of M dwarfs with a Kepler-like planet, against 30 to 100 observed
    summary = galaxy.summarize(records)
    flags = galaxy.compare(summary, galaxy.expectations_for(summary, stored))
    named = {(f.cell, f.measure, f.kind) for f in flags}
    assert ("M", "inner_fraction", "out-of-range") in named
    assert (galaxy.POOL, "inner_fraction", "out-of-range") in named
    assert all(f.cell in ("M", galaxy.POOL) for f in flags)
    assert any(f.line().startswith("OUT-OF-RANGE M inner_fraction") for f in flags)


def test_a_flag_whose_group_has_a_known_gap_carries_the_gap_item(stored):
    tagged = {(g, m): b["anomaly"] for g, spec in stored["types"].items() for m, b in spec["measures"].items() if b.get("anomaly")}
    assert tagged, "the stored occurrence file names at least one known gap"
    (group, measure), item = next(iter(sorted(tagged.items())))
    summary = galaxy.summarize(_sky(group, galaxy.MIN_GROUP, 0, 1))
    expectations = galaxy.expectations_for(summary, stored)
    assert expectations["types"][group]["measures"][measure]["anomaly"] == item
    assert ITEM.match(item) and item in expectations["anomalies"]
    pooled = expectations["types"].get(galaxy.POOL, {"measures": {}})["measures"]
    if measure in pooled:
        assert pooled[measure]["anomaly"] == item


def test_a_system_that_raised_or_ran_over_is_flagged_and_not_averaged(stored):
    records = _sky("M", 30, 20, 3) + [{"star": "M2V", "group": "M", "error": "ValueError: x"},
                                      {"star": "other", "group": None, "abandoned": True}]
    summary = galaxy.summarize(records)
    assert summary["M"]["n"] == 30 and summary["M"]["error_count"] == 1
    flags = galaxy.compare(summary, galaxy.expectations_for(summary, stored))
    assert any(f.kind == "error" for f in flags)


def test_a_star_outside_every_group_is_counted_and_never_scored(stored):
    records = _sky("M", 30, 20, 3) + [dict(_system("G", [], star_mass=3.0), star="other", group=None)]
    summary = galaxy.summarize(records)
    assert summary["other"]["n"] == 1 and summary["M"]["n"] == 30 and summary[galaxy.POOL]["n"] == 30


def test_the_class_table_counts_final_classes_and_the_report_names_the_life_bodies():
    records = [_system("M", [(1.0, 0.1, 1.0)], classes={"S": 2, "asteroid": 1}), _system("M", [], classes={"C": 1})]
    records[0]["life_bodies"] = 1
    table, counted = galaxy.class_table(records)
    assert counted == 2 and dict((letter, count) for letter, count, _share in table) == {"S": 2, "asteroid": 1, "C": 1}
    assert sum(share for *_x, share in table) == pytest.approx(1.0)
    summary = galaxy.summarize(records)
    text = galaxy.format_report(summary, records, [], "header")
    assert "Bodies with life, from the stored life class: 1" in text and "| S | 2 |" in text


def test_a_sector_of_the_solar_position_is_the_same_on_every_draw():
    first = galaxy.draw_sector(0)
    again = galaxy.draw_sector(0)
    assert first == again and first


@pytest.mark.slow
@pytest.mark.timeout(3600)
def test_the_galaxy_placed_stars_raise_no_new_occurrence_flag(stored):
    """Draws `SECTORS` seeded sectors at the solar position on every core, writes the latest-run report and fails on a flag no known gap
    explains, naming its cell (a host group or `galaxy`) and the measure that left its band."""
    records = galaxy.draw_sectors()
    summary = galaxy.summarize(records)
    expectations = galaxy.expectations_for(summary, stored)
    findings = galaxy.compare(summary, expectations)
    header = (f"{galaxy.SECTORS} seeded sectors at the reference solar position; host groups of {galaxy.MIN_GROUP} stars or more are cells "
              f"of their own, every group counts in the `{galaxy.POOL}` cell.")
    path = census.write_report(galaxy.format_report(summary, records, findings, header), name=galaxy.REPORT_NAME)
    new = [f for f in findings if not f.known]
    assert not new, f"new galaxy occurrence flags (report: {path}):\n" + "\n".join(f.line() for f in new)
