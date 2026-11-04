# tests/test_redraw_rates.py

"""TEST.476 (GEN.1006): the redraw-rate test. `draw.redrawn` counters beside every retry loop of a system build, the registry of their
kinds (`generation/redraw_kinds.py`), the census harness (`tests/audit/redraw_rates.py`) and its two gates: a ceiling per kind and
group (a rate that doubles trips it) and the exact-zero invariant for the overlap kinds (expected to fail until GEN.1015)."""

import json
import pathlib
import re

import pytest

from planetgen.generation import redraw_kinds
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from tests.audit import prevalence as census
from tests.audit import redraw_rates as rr
from tests.fuzz_support import deterministic_entropy

SOURCE = pathlib.Path(__file__).resolve().parent.parent / "planetgen"
SAMPLE_PER_GROUP = 12
"""int: Systems per group in the sample the gates read (the baseline file was measured on a far larger one; a run's own sampling error is
part of the gate)."""

CALL = re.compile(r"draw\.redrawn\(([^)]*)\)")
LITERAL = re.compile(r"\"([a-z_]+\.[a-z_]+)\"")


def _call_sites():
    """`{kind: [file:line]}` for every `draw.redrawn(...)` in the package, from the source text."""
    found = {}
    for path in sorted(SOURCE.rglob("*.py")):
        text = path.read_text(encoding="utf-8")
        for match in CALL.finditer(text):
            if path.name == "draw.py" and "kind" in match.group(1) and "count" in match.group(1):
                continue
            line = text.count("\n", 0, match.start()) + 1
            for kind in LITERAL.findall(match.group(1)):
                found.setdefault(kind, []).append(f"{path.relative_to(SOURCE)}:{line}")
    return found


# --- the counter -------------------------------------------------------------------------------------------------------------

def test_a_redraw_outside_a_counting_block_is_not_counted():
    draw.redrawn("pair.kozai")
    with draw.counting() as tally:
        pass
    assert tally == {}


def test_counting_counts_by_kind_and_nests():
    with draw.counting() as outer:
        draw.redrawn("pair.kozai")
        with draw.counting() as inner:
            draw.redrawn("pair.kozai", 2)
            draw.redrawn("tree.member_dropped")
        draw.redrawn("pair.kozai")
    assert inner == {"pair.kozai": 2, "tree.member_dropped": 1}
    assert outer == {"pair.kozai": 4, "tree.member_dropped": 1}


def test_counting_reads_no_stream():
    with deterministic_entropy(7):
        draw.set_run_seed(7)
        plain = [draw.random() for _ in range(5)]
    with deterministic_entropy(7):
        draw.set_run_seed(7)
        with draw.counting():
            draw.redrawn("pair.kozai", 3)
            counted = [draw.random() for _ in range(5)]
    assert plain == counted


@pytest.mark.parametrize("stars", [None, 3])
def test_a_counted_system_is_the_same_system(stars):
    """The counters change nothing that is generated: a system built with a tally on equals the one built without."""
    def build(counted):
        config = SystemConfig()
        if stars:
            config.STARS = stars
        with deterministic_entropy(11):
            draw.set_run_seed(11)
            if counted:
                with draw.counting():
                    return StarSystem(system_config=config).to_dict()
            return StarSystem(system_config=config).to_dict()

    assert json.dumps(build(True), sort_keys=True, default=str) == json.dumps(build(False), sort_keys=True, default=str)


# --- the registry ------------------------------------------------------------------------------------------------------------

def test_every_call_names_a_registered_kind():
    unknown = {kind: sites for kind, sites in _call_sites().items() if kind not in redraw_kinds.KINDS}
    assert not unknown, f"add these to redraw_kinds.KINDS: {unknown}"


def test_every_registered_kind_has_a_call_site():
    assert not set(redraw_kinds.KINDS) - set(_call_sites()), "a kind with no counter beside any loop"


def test_every_kind_has_a_scope_and_a_summary_and_the_overlap_kinds_are_the_stars():
    for name, kind in redraw_kinds.KINDS.items():
        assert re.fullmatch(r"[a-z_]+\.[a-z_]+", name)
        assert kind.scope and len(kind.summary) > 20, name
    assert redraw_kinds.OVERLAP_KINDS == ("companions.overlap", "pair.overlap", "tree.restart_overlap")


# --- the harness -------------------------------------------------------------------------------------------------------------

def _summary(events_by_group, n=100):
    return {"groups": {group: {"n": n, "events": events, "systems": {k: min(v, n) for k, v in events.items()}, "abandoned": 0,
                               "errors": [], "error_count": 0} for group, events in events_by_group.items()},
            "by_stars": {}, "by_class": {}}


def _expectations(rate, ceiling):
    return {"groups": {"stars-3": {"n_ref": 1000, "kinds": {"pair.kozai": {"baseline": rate, "events": 0, "ceiling": ceiling}}}}}


def test_a_rate_that_doubles_trips_its_ceiling_and_a_repeat_does_not():
    expectations = _expectations(0.5, 2 * 0.5 + 0.001)
    assert not rr.ceilings_broken(_summary({"stars-3": {"pair.kozai": 50}}, 100), expectations)
    assert not rr.ceilings_broken(_summary({"stars-3": {"pair.kozai": 100}}, 100), expectations)
    broken = rr.ceilings_broken(_summary({"stars-3": {"pair.kozai": 180}}, 100), expectations)
    assert [(group, kind) for group, kind, *_ in broken] == [("stars-3", "pair.kozai")]


def test_a_few_events_of_a_rare_kind_are_not_a_doubling():
    expectations = _expectations(0.0, 0.001)
    assert not rr.ceilings_broken(_summary({"stars-3": {"pair.kozai": 2}}, 100), expectations)
    assert rr.ceilings_broken(_summary({"stars-3": {"pair.kozai": 60}}, 100), expectations)


def test_the_overlap_kinds_are_never_held_under_a_ceiling():
    expectations = {"groups": {"stars-3": {"n_ref": 10, "kinds": {"pair.overlap": {"baseline": 0.0, "events": 0, "ceiling": 0.0}}}}}
    assert not rr.ceilings_broken(_summary({"stars-3": {"pair.overlap": 90}}, 100), expectations)
    assert rr.overlap_events(_summary({"stars-3": {"pair.overlap": 90, "pair.kozai": 4}}, 100)) == {"pair.overlap": 90}


def test_wilson_and_the_rate_bound():
    low, high = rr.wilson(0, 100)
    assert low == 0.0 and 0.0 < high < 0.05
    low, high = rr.wilson(50, 100)
    assert 0.39 < low < 0.41 and 0.59 < high < 0.61
    assert rr.rate_lower_bound(0, 100) == 0.0
    assert rr.rate_lower_bound(100, 100) == pytest.approx(0.8)
    assert rr.rate_lower_bound(100, 100, squares=100) == pytest.approx(1.0 - 2.0 * 0.0)  # one event in every system: no spread
    assert rr.rate_lower_bound(100, 100, squares=10000) < 0.01  # all in one system: the mean of 1.0 is within its own spread of zero


def test_build_expectations_sets_the_ceiling_to_twice_the_baseline_plus_one_event():
    records = [{"group": "stars-3", "stars": 3, "class": "GV", "counts": {"pair.kozai": 2 if i % 2 else 0}} for i in range(100)]
    built = rr.build_expectations(records, "abc", "test")
    bounds = built["groups"]["stars-3"]["kinds"]["pair.kozai"]
    assert bounds["baseline"] == pytest.approx(1.0) and bounds["ceiling"] == pytest.approx(2.0 + 0.01)
    assert built["reference"]["tier"] == "T9" and built["reference"]["tweak"] == "L3"


def test_the_report_lists_every_kind_by_group_stars_and_class():
    records = [{"group": "stars-3", "stars": 3, "class": "GV", "counts": {"pair.kozai": 1}},
               {"group": "ordinary", "stars": 1, "class": "KV", "counts": {}}]
    text = rr.format_report(rr.summarize(records))
    for kind in redraw_kinds.KINDS:
        assert f"| {kind} |" in text
    assert "by number of stars" in text and "by class of the placed star" in text
    assert "| pair.kozai | 1.0000 | 0.0000 |" in text.split("by class of the placed star")[1]


def test_bucket_of_and_group_names():
    assert [rr.bucket_of(n) for n in (1, 2, 3, 4, 5, 6, 9)] == ["1", "2", "3", "4", "5+", "5+", "5+"]
    assert set(rr.GROUPS) == {"ordinary", "stars-1", "stars-2", "stars-3", "stars-4", "stars-5", "stars-6", "world-requested", "world-barred"}


# --- the registered ceilings -------------------------------------------------------------------------------------------------

def test_the_expectations_register_every_ceiling_with_a_tier_and_a_tweak_level():
    expectations = rr.load_expectations()
    reference = expectations["reference"]
    assert reference["tier"] == "T9" and reference["tweak"] == "L3" and reference["commit"] and reference["ceiling_factor"] == 2.0
    assert set(expectations["groups"]) == set(rr.GROUPS)
    for group, spec in expectations["groups"].items():
        assert spec["n_ref"] > 0, group
        for kind, bounds in spec["kinds"].items():
            assert kind in redraw_kinds.KINDS, kind
            assert bounds["ceiling"] >= rr.CEILING_FACTOR * bounds["baseline"], (group, kind)


@pytest.fixture(scope="module")
def sample():
    """The seeded sample the gates read: `SAMPLE_PER_GROUP` systems of every group."""
    summary = rr.summarize(rr.draw_groups(SAMPLE_PER_GROUP, chunk=SAMPLE_PER_GROUP))
    census.write_report(rr.format_report(summary, f"{SAMPLE_PER_GROUP} systems per group; the baseline file holds the reference run."),
                        name=rr.REPORT_NAME)
    return summary


def test_the_sample_raises_no_error_and_hangs_nowhere(sample):
    for group, mine in sample["groups"].items():
        assert mine["n"] == SAMPLE_PER_GROUP, (group, mine["errors"])
        assert not mine["abandoned"] and not mine["error_count"], (group, mine["errors"])


def test_no_unregistered_kind_fires_in_the_sample(sample):
    assert rr.unregistered_kinds(sample) == []


def test_the_forced_star_groups_come_out_with_the_stars_they_were_forced_to():
    records = rr.draw_group("stars-3", 3)
    assert all(r["stars"] == 3 for r in records), [r.get("stars") for r in records]


def test_no_redraw_rate_is_above_its_registered_ceiling(sample):
    """TEST.476: a change that suddenly redraws more is caught. Ceilings are twice the baseline measured before the overlap-free build;
    the run's own sampling error is held against them. A change that moves a rate on purpose re-measures the baseline in the same change:
    `python -m tests.audit.redraw_rates measure --count 600 --out records.json` (per group, see `--groups`), then `build records.json --commit <sha>
    --label <why>`, and says so in its change note."""
    broken = rr.ceilings_broken(sample, rr.load_expectations())
    assert not broken, "\n".join(f"{g} {k}: {e} events in {n} systems = {rate:.4f} per system, ceiling {c:.4f}" for g, k, e, n, rate, c in broken)


OVERLAP_SYSTEMS = (("stars-3", 59, "companions.overlap"), ("stars-3", 6, "tree.restart_overlap"), ("stars-3", 1213, "pair.overlap"),
                   ("stars-4", 58, "tree.restart_overlap"), ("stars-3", 1906, "pair.overlap"))
"""tuple: `(group, number, kind)` of seeded systems of the baseline run that redraw because a star overlaps what it holds, one of each overlap
kind or more. A random sample of a few systems seldom meets an overlap (about one system in forty does), so the invariant reads these."""


@pytest.fixture(scope="module")
def overlap_systems():
    """The overlap-kind counts of every system in `OVERLAP_SYSTEMS`: `{(group, number): {kind: events}}`."""
    found = {}
    for group, number, _ in OVERLAP_SYSTEMS:
        (record,) = rr.draw_group(group, 1, first=number)
        assert "counts" in record, record
        found[(group, number)] = {kind: n for kind, n in record["counts"].items() if kind in redraw_kinds.OVERLAP_KINDS}
    return found


def test_the_overlap_redraw_count_is_exactly_zero(sample, overlap_systems):
    """GEN.1006's Done clause: no pair, tree or companion set is drawn again because a star overlaps what it holds. Holds from the member
    record (GEN.1009), the keep-out rule (GEN.1010) and the truncated draws (GEN.1011); the systems that used to redraw are the sample
    it reads."""
    assert rr.overlap_events(sample) == {}
    assert not {key: counts for key, counts in overlap_systems.items() if counts}
