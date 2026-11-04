"""
PERF.92: the time each generated object takes, per object type: exclusive of the objects inside it, kept as per-kind
aggregates with a repeatable 1-in-1,000 sample (a row for every object only on request), and saved under the run.
"""

import types

import pytest

from planetgen.db import store
from planetgen.generation import object_timing, run_sector

from tests import first_fill
from tests.test_galaxy_gen import REFERENCE_SHAPE_ARGV
from tests.test_sector_gen import _real_sector_args


_CLOCK = [0.0]    # a fake clock the toy classes move, so no timing assertion depends on how fast the machine is


def _spend(seconds):
    _CLOCK[0] += seconds


class Outer:
    def __init__(self, inner=None):
        _spend(0.004)
        self.inner = inner() if inner else None


class Inner:
    def __init__(self):
        _spend(0.004)


class SubInner(Inner):
    def __init__(self):
        super().__init__()


@pytest.fixture
def toy_classes(monkeypatch):
    """`Outer`, `Inner` and `SubInner` registered as a system, a star and a binary star, for one test."""
    monkeypatch.setattr(object_timing, "_KIND_OF", {})
    _CLOCK[0] = 0.0
    monkeypatch.setattr(object_timing, "time", types.SimpleNamespace(perf_counter=lambda: _CLOCK[0]))
    saved = {cls: cls.__init__ for cls in (Outer, Inner, SubInner)}
    object_timing.track(Outer, "system")
    object_timing.track(Inner, "star")
    object_timing.track(SubInner, "binary-star")
    yield
    for cls, init in saved.items():
        cls.__init__ = init


def test_every_kind_has_a_registered_class_and_the_other_way_round():
    import planetgen.generation.system  # noqa: F401 -- importing registers every generator
    import planetgen.generation.run_sector  # noqa: F401
    kinds = set(object_timing._KIND_OF.values()) | {"moon"}    # a moon is a Planet that knows it is one
    assert kinds == set(object_timing.KINDS)


def test_a_sector_reports_every_kind_it_built_and_nothing_else():
    args = _real_sector_args(6, "TimedSector")
    with object_timing.unit("1,2,3") as timings:
        _name, sector = run_sector.generate_sector(args)
    assert set(timings.aggregates) <= set(object_timing.KINDS)
    assert timings.aggregates["system"].count == len(sector.entries)
    assert {"system", "star", "planet"} <= set(timings.aggregates)
    for aggregate in timings.aggregates.values():
        assert aggregate.count == sum(aggregate.histogram) and aggregate.minimum <= aggregate.mean <= aggregate.maximum


def test_an_object_outside_a_unit_costs_nothing_and_records_nothing(toy_classes):
    Outer(Inner)
    assert object_timing.active() is None


def test_time_spent_in_an_inner_object_is_not_counted_twice(toy_classes):
    with object_timing.unit("a") as timings:
        Outer(Inner)
    outer, inner = timings.aggregates["system"], timings.aggregates["star"]
    assert outer.total == pytest.approx(0.004) and inner.total == pytest.approx(0.004)    # 4 ms each, not 8 ms for the outer


def test_a_subclass_building_through_its_base_is_one_object_of_its_own_kind(toy_classes):
    with object_timing.unit("a") as timings:
        SubInner()
    assert set(timings.aggregates) == {"binary-star"} and timings.aggregates["binary-star"].count == 1


def test_units_do_not_nest(toy_classes):
    with object_timing.unit("outer") as outer, object_timing.unit("inner") as inner:
        Inner()
    assert inner is outer and outer.aggregates["star"].count == 1


def test_the_sample_is_about_one_in_a_thousand_and_repeats():
    chosen = [n for n in range(200_000) if object_timing.sampled("0,0,0", "star", n)]
    assert 150 < len(chosen) < 250
    assert chosen == [n for n in range(200_000) if object_timing.sampled("0,0,0", "star", n)]
    assert chosen != [n for n in range(200_000) if object_timing.sampled("0,0,1", "star", n)]


def test_only_the_sample_is_kept_unless_every_object_is_asked_for(toy_classes):
    with object_timing.unit("a", full=False) as sample_only:
        for _ in range(30):
            Inner()
    with object_timing.unit("b", full=True) as every:
        for _ in range(30):
            Inner()
    assert len(sample_only.samples) == sum(object_timing.sampled("a", "star", n) for n in range(30))
    assert [sample[1] for sample in every.samples] == list(range(30))


def test_the_histogram_buckets_double_from_a_microsecond():
    assert [object_timing.bucket(s) for s in (5e-7, 1e-6, 1.9e-6, 2e-6, 1e-3)] == [0, 1, 1, 2, 10]
    assert object_timing.bucket(1000.0) == object_timing.HISTOGRAM_BUCKETS - 1


@pytest.mark.db
def test_a_unit_saves_under_its_run_and_the_run_adds_its_sectors_up(mysql_config, monkeypatch, toy_classes):
    monkeypatch.setenv("PLANETGEN_GENERATION_STATS", "1")
    conn = store.get_connection(mysql_config)
    try:
        run_id = store.start_generation_run(conn, "galaxy", ["galaxy"])
        for label in ("a", "b"):
            with object_timing.unit(label, full=True) as timings:
                Outer(Inner)
                Inner()
            assert object_timing.save(conn, run_id, timings) == 2 + 3
        conn.commit()
        rows = {row["kind"]: row for row in object_timing.kind_summaries(conn)}
        assert rows["system"]["objects"] == 2 and rows["star"]["objects"] == 4
        assert rows["star"]["total_seconds"] > 0.012 and 0.003 < rows["star"]["mean_seconds"] < 0.008
        assert rows["star"]["fastest_seconds"] <= rows["star"]["mean_seconds"] <= rows["star"]["slowest_seconds"]
        # the kind that took the most time first; a unit counts each object's own time: the four stars take 0.016 s, the two systems 0.008 s
        ordered = object_timing.kind_summaries(conn)
        assert [row["kind"] for row in ordered] == ["star", "system"]
        assert ordered[0]["total_seconds"] > ordered[1]["total_seconds"]
        assert conn.execute("SELECT COUNT(*) AS n FROM object_timing_samples WHERE run_id = ?",
                            (run_id,)).fetchone()["n"] == 6
    finally:
        conn.close()


@pytest.mark.db
def test_nothing_is_saved_without_a_run_or_with_the_stats_off(mysql_config, monkeypatch, toy_classes):
    conn = store.get_connection(mysql_config)
    try:
        with object_timing.unit("a") as timings:
            Inner()
        assert object_timing.save(conn, None, timings) == 0
        monkeypatch.setenv("PLANETGEN_GENERATION_STATS", "0")
        assert object_timing.save(conn, 1, timings) == 0
    finally:
        conn.close()


@pytest.mark.db
def test_a_galaxy_run_saves_each_sectors_timings_under_its_run(mysql_config, monkeypatch):
    from planetgen.galaxy import seed as galaxySeed
    from tests.test_version_key import SEED, _generate
    monkeypatch.setenv("PLANETGEN_GENERATION_STATS", "1")
    first_fill.without_the_deepest_floor(monkeypatch)   # the counts below are of the systems asked for (GEN.633)
    _generate(mysql_config, monkeypatch, "plan", "--nuclear-amplitude", "0", "--edge=160pc",
              "--seed", galaxySeed.format_seed(SEED), *REFERENCE_SHAPE_ARGV)
    _generate(mysql_config, monkeypatch, "sectors", "--shell", "1", "--limit", "2", "--num-systems", "2")
    conn = store.get_connection(mysql_config)
    try:
        run_id = conn.execute("SELECT id FROM generation_runs WHERE command = 'sectors'").fetchone()["id"]
        sectors = {row["id"] for row in conn.execute("SELECT id FROM sectors").fetchall()}
        rows = conn.execute("SELECT sector_id, kind, object_count FROM object_timing_stats WHERE run_id = ?",
                            (run_id,)).fetchall()
        assert {row["sector_id"] for row in rows} == sectors
        counts = {}
        for row in rows:
            counts[row["kind"]] = counts.get(row["kind"], 0) + row["object_count"]
        assert counts["system"] == 4 and counts["star"] >= 4
    finally:
        conn.close()
