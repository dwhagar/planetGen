"""OPS.161: the queue engine chooses a job's split from its declaration, the profile and the free resources, and says why (`planetgen.queue.splitter`)."""

import json

import pytest

from planetgen.queue import breakpoints as bp, commands, runrecord, splitprofile, splitter
from planetgen.queue.splitprofile import Reading

ROOMY = splitter.Resources(available_mb=64000, cores=8, workers=8)
FAST = Reading(seconds_per_item=0.5, mb_base=1000, mb_per_item=2.0, samples=10)


def dynamic(items=10000, unit=100, whole=(), **options):
    cuts = tuple(c for c in range(unit, items, unit) if not any(a < c < b for a, b in whole))
    return bp.Declaration(items, cuts, split=bp.BYTE_IDENTICAL, whole=whole, **options)


def test_a_job_that_declares_no_flag_keeps_its_suggested_split_whatever_the_profile_says():
    got = splitter.choose("k", bp.even(1000, 100), FAST, ROOMY)
    assert got.by == splitter.BY_SUGGESTED and got.cuts == bp.even(1000, 100).suggested
    assert "output depends on the cut" in got.reason
    assert (got.old_units, got.new_units, got.old_items, got.new_items) == (10, 10, 100, 100)


def test_with_no_usable_profile_the_suggested_split_stands_and_the_decision_says_so():
    for reading in (None, Reading(0.5, samples=2), Reading(0.0, samples=50), Reading(0.5, samples=5, trust=0.0)):
        got = splitter.choose("k", dynamic(), reading, ROOMY)
        assert got.by == splitter.BY_SUGGESTED and got.cuts == dynamic().suggested, reading
        assert "nothing usable" in got.reason
    assert splitter.choose("k", dynamic(), FAST, None).by == splitter.BY_SUGGESTED, "no resources given: nothing to choose from"


def test_a_profile_cuts_the_stage_into_units_of_the_target_seconds():
    got = splitter.choose("k", dynamic(), FAST, ROOMY, seconds=60)
    assert got.by == splitter.BY_PROFILE and got.limit == "time"
    assert got.new_items == 120 and got.new_units == 84 and got.old_units == 100
    assert got.cuts[:3] == (120, 240, 360)
    assert all(b - a <= 120 for a, b in zip((0, *got.cuts), (*got.cuts, got.items)))


def test_every_cut_is_legal_and_no_span_that_stays_whole_is_cut():
    decl = dynamic(whole=((1000, 1500), (4000, 4010)))
    for seconds in (5, 20, 60, 400):
        got = splitter.choose("k", decl, FAST, ROOMY, seconds=seconds)
        assert set(got.cuts) <= set(decl.legal_cuts), seconds
        assert not any(1000 < c < 1500 or 4000 < c < 4010 for c in got.cuts)
    long_span = dynamic(whole=((100, 9000),))
    got = splitter.choose("k", long_span, FAST, ROOMY, seconds=10)
    assert 100 in got.cuts and 9000 in got.cuts and not any(100 < c < 9000 for c in got.cuts), "the long span is one unit, however long it runs"


def test_memory_limits_the_unit_when_a_unit_would_not_fit_its_share():
    tight = splitter.Resources(available_mb=12000, cores=8, workers=4)
    got = splitter.choose("k", dynamic(), FAST, tight, seconds=600)
    assert got.limit == "memory"
    assert got.new_items == pytest.approx(250, abs=1), "(12000 x 0.5 / 4 workers - 1000 MB base) / 2 MB an item"
    assert splitter.choose("k", dynamic(), FAST, ROOMY, seconds=600).new_items > got.new_items


def test_the_typed_split_is_the_default_and_the_bound():
    slow = Reading(seconds_per_item=60.0, mb_per_item=0.0, samples=9)
    got = splitter.choose("k", dynamic(), slow, ROOMY, seconds=60)
    assert got.limit == "floor" and got.new_items == 25, "one item a unit would be a quarter of the suggested unit at the smallest"
    quick = Reading(seconds_per_item=0.0001, samples=9)
    got = splitter.choose("k", dynamic(items=100000), quick, ROOMY, seconds=60)
    assert got.limit == "ceiling" and got.new_items == 400
    got = splitter.choose("k", dynamic(floor_items=7, ceiling_items=9), slow, ROOMY)
    assert got.new_items == 7


def test_enough_units_are_made_to_keep_the_workers_busy():
    got = splitter.choose("k", dynamic(items=800, unit=100, ceiling_items=800), Reading(0.001, samples=9), splitter.Resources(64000, 8, 8), seconds=60)
    assert got.limit == "workers" and got.new_units == 8 and got.new_items == 100


def test_never_more_units_than_the_queue_holds():
    got = splitter.choose("k", dynamic(items=100000, unit=10, floor_items=1, ceiling_items=1), Reading(1.0, samples=9), ROOMY, seconds=1, max_units=1000)
    assert got.new_units <= 1000 and got.items == 100000


def test_a_job_of_one_unit_by_nature_or_no_item_is_never_cut():
    for decl in (bp.whole_job(split=bp.BYTE_IDENTICAL), bp.Declaration(0, split=bp.BYTE_IDENTICAL)):
        got = splitter.choose("k", decl, FAST, ROOMY)
        assert got.cuts == () and got.by == splitter.BY_SUGGESTED


def test_the_same_inputs_give_the_same_cut():
    one = splitter.choose("k", dynamic(), FAST, ROOMY, seconds=60)
    two = splitter.choose("k", dynamic(), FAST, ROOMY, seconds=60)
    assert one == two


def test_the_decision_is_plain_json_that_carries_the_cuts_only_when_they_were_chosen():
    chosen = splitter.choose("k", dynamic(), FAST, ROOMY, seconds=60).as_dict()
    assert json.loads(json.dumps(chosen))["cuts"][:2] == [120, 240] and chosen["by"] == "profile" and chosen["samples"] == 10
    suggested = splitter.choose("k", bp.even(100, 10), FAST, ROOMY).as_dict()
    assert "cuts" not in suggested and suggested["by"] == "suggested"


def test_the_message_names_the_old_and_new_figures_and_the_reason():
    text = splitter.choose("k", dynamic(), FAST, ROOMY, seconds=60).message()
    assert "changed from 100 unit(s) of about 100 item(s) to 84 of about 120" in text and "10 sample(s)" in text and "limited by time" in text
    assert "suggested split" in splitter.choose("k", bp.even(100, 10), None, None).message()


def test_spans_of_a_decision_cover_every_item_once():
    got = splitter.choose("k", dynamic(), FAST, ROOMY, seconds=60)
    spans = splitter.spans(got.items, got)
    assert spans[0][0] == 0 and spans[-1][1] == got.items and all(a[1] == b[0] for a, b in zip(spans, spans[1:]))


def test_decide_reads_the_stub_profile_and_the_host(monkeypatch):
    monkeypatch.setattr(splitprofile, "STUB", {"k": FAST})
    monkeypatch.setattr(splitter, "resources_now", lambda kind=None: ROOMY)
    monkeypatch.setattr(splitter, "target_seconds", lambda: 60.0)
    assert splitter.decide("k", dynamic()).by == splitter.BY_PROFILE
    assert splitter.decide("other", dynamic()).by == splitter.BY_SUGGESTED, "no reading for the kind"
    assert splitter.decide("k", bp.even(1000, 100)).by == splitter.BY_SUGGESTED, "flagged suggested"
    monkeypatch.setattr(splitprofile, "lookup", lambda kind: 1 / 0)
    assert splitter.decide("k", dynamic()).by == splitter.BY_SUGGESTED, "learning never blocks a job"


def test_the_stub_profile_is_empty_and_a_test_run_reads_no_table():
    assert splitprofile.STUB == {} and splitprofile.lookup("scatter-mass") is None


def _row(items, spi, peak, memory=32000, cores=16, epoch=146):
    return {"items": items, "size_bucket": splitprofile_bucket(items), "seconds_per_item": spi, "peak_mb": peak, "host_memory_mb": memory, "host_cores": cores, "generator_epoch": epoch}


def splitprofile_bucket(items):
    from planetgen.queue import jobprofiles
    return jobprofiles.size_bucket(items)


def test_the_reading_of_this_machines_rows_is_in_full_trust_with_a_memory_line():
    rows = [_row(100, 0.02, 500 + 100 * 2.0), _row(110, 0.03, 500 + 110 * 2.0), _row(120, 0.01, 500 + 120 * 2.0), _row(1000, 0.02, 500 + 1000 * 2.0), _row(1100, 0.02, 500 + 1100 * 2.0)]
    got = splitprofile.from_rows(rows, 32000, 16, 146)
    assert got.trust == 1.0 and got.samples == 5 and got.seconds_per_item == pytest.approx(0.02)
    assert got.mb_per_item == pytest.approx(2.0, rel=0.05) and got.mb_base == pytest.approx(500, rel=0.1) and got.usable()


def test_one_size_bucket_gives_the_raw_peak_over_items_and_no_base():
    got = splitprofile.from_rows([_row(100, 0.02, 600), _row(105, 0.02, 650), _row(110, 0.02, 700)], 32000, 16, 146)
    assert got.mb_base == 0.0 and got.mb_per_item == pytest.approx(650 / 105)


def test_rows_from_another_machine_or_epoch_count_at_half_trust_and_too_few_rows_give_nothing():
    rows = [_row(100, 0.02, 700, memory=64000, cores=32) for _ in range(3)]
    assert splitprofile.from_rows(rows, 32000, 16, 146).trust == splitprofile.OTHER_TRUST
    assert splitprofile.from_rows(rows[:2], 32000, 16, 146) is None
    assert splitprofile.from_rows([_row(100, 0.02, 700, epoch=145) for _ in range(3)], 32000, 16, 146).trust == 0.5
    assert splitprofile.from_rows([_row(100, None, 700) for _ in range(5)], 32000, 16, 146) is None


def test_lookup_reads_the_unit_rows_of_the_kind_and_fails_open(monkeypatch):
    from planetgen.queue import jobprofiles, memorygate, profileclaim
    monkeypatch.setattr(jobprofiles, "enabled", lambda: True)
    monkeypatch.setattr(profileclaim, "enabled", lambda: True)
    monkeypatch.setattr(memorygate, "read_host", lambda: memorygate.Host(32000, 20000, 16, 0.1))
    seen = []
    monkeypatch.setattr(jobprofiles, "unit_rows", lambda kind, *a, **k: seen.append(kind) or [_row(100, 0.02, 700)] * 4)
    got = splitprofile.lookup("scatter-mass")
    assert seen == ["unit:scatter-mass"] and got is not None and got.samples == 4
    monkeypatch.setattr(jobprofiles, "unit_rows", lambda *a, **k: 1 / 0)
    assert splitprofile.lookup("scatter-mass") is None, "a table that cannot answer leaves the suggested split"
    monkeypatch.setattr(profileclaim, "enabled", lambda: False)
    monkeypatch.setattr(jobprofiles, "unit_rows", lambda *a, **k: seen.append("read") or [])
    assert splitprofile.lookup("scatter-mass") is None and "read" not in seen, "the self-tuning switch turns the reading off"


def test_the_resources_come_from_the_gate_and_the_stage_cap(monkeypatch):
    from planetgen.queue import memorygate, workers
    monkeypatch.setattr(memorygate, "read_host", lambda: memorygate.Host(32000, 20000, 16, 0.1))
    monkeypatch.setattr(workers, "stage_worker_cap", lambda kind: 5)
    assert splitter.resources_now("scatter-mass") == splitter.Resources(20000, 16, 5)
    monkeypatch.setattr(workers, "stage_worker_cap", lambda kind: None)
    assert splitter.resources_now("plan-skeleton") == splitter.Resources(20000, 16, 16)
    monkeypatch.setattr(memorygate, "read_host", lambda: 1 / 0)
    assert splitter.resources_now() == splitter.Resources(0.0, 1, 1)


def test_the_target_is_a_share_of_the_unit_ceiling_setting():
    assert splitter.target_seconds() == pytest.approx(0.2 * 300)


@pytest.fixture
def conn(redis_server, key_prefix, monkeypatch):
    from planetgen.queue import redisqueue
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    return redisqueue.connect(redis_server)


def test_a_decision_is_written_to_the_runs_log_with_the_old_and_new_values(conn):
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("s", "S", units=1)])
    splitter.announce(conn, run_id, 0, splitter.choose("scatter-mass", dynamic(), FAST, ROOMY, seconds=60))
    events = [e for e in runrecord.log_records(conn, run_id) if "split" in str(e)]
    text = runrecord.log_text(conn, run_id)
    assert "Split of scatter-mass changed from 100 unit(s)" in text
    assert events


def test_a_stage_opened_with_a_decision_stores_it_in_its_plan_and_logs_it(conn, monkeypatch):
    from planetgen.queue import units, wrapper
    decision = splitter.choose("probe", dynamic(items=500, unit=50), FAST, ROOMY, seconds=60)

    _EXPANSION[:] = [commands.Expansion(units=[({"i": i}, 1.0) for i in range(3)], plan={"x": 1}, split=decision)]
    monkeypatch.setitem(commands.REGISTRY, "probe", commands.CommandSpec("probe", "A probe", unit=("tests.test_splitter", "unit"), expand=("tests.test_splitter", "expand")))
    queued = []
    monkeypatch.setattr(units, "enqueue_units", lambda *a, **k: queued.append(a))
    step = commands.Step("probe", "probe", "A probe", {})
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("probe", "A probe", units=1)], steps=[step.as_dict()])
    assert wrapper._open(conn, run_id, 0) == 3
    plan = runrecord.read_plan(conn, run_id, 0)
    assert plan["split"]["by"] == "profile" and plan["split"]["cuts"][:2] == [63, 126] and plan["detail"] == {"x": 1}
    assert "Split of probe changed" in runrecord.log_text(conn, run_id)
    assert queued


_EXPANSION = []


def expand(params):
    """The probe command's cut: the expansion the test put in `_EXPANSION`."""
    return _EXPANSION[0]


def unit(args, report):
    return {}
