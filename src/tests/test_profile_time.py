"""OPS.159: the expected seconds of a unit from the job profile (queue/profiletime.py), read by the untimed estimate and the heartbeat judge."""

import logging

import pytest

from planetgen.queue import beatjudge, commandgaps, memorygate, profileclaim, profiletime, redisqueue, roughtime, runrecord
from planetgen.queue.memorygate import Host

HOST = (16000.0, 8, 146)
OTHER = (32000.0, 16, 146)


@pytest.fixture
def conn(redis_server, key_prefix, monkeypatch):
    monkeypatch.setattr(memorygate, "PREFIX", f"{key_prefix}:gate")
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    monkeypatch.setenv(profileclaim.ENV_VAR, "1")
    monkeypatch.setattr(memorygate, "read_host", lambda: Host(HOST[0], 8000.0, HOST[1], 0.1))
    from planetgen.util import generator_epoch
    monkeypatch.setattr(generator_epoch, "generator_epoch", lambda: HOST[2])
    return redisqueue.connect(redis_server)


def rows(*seconds, shape=HOST):
    return lambda kind: [{"seconds": s, "outcome": "done", "host_memory_mb": shape[0], "host_cores": shape[1], "generator_epoch": shape[2]} for s in seconds]


def test_it_is_the_median_of_this_machines_units_and_nothing_with_too_few(conn):
    assert profiletime.unit_seconds(conn, "link", token="r1", read_rows=rows(10.0, 12.0)) is None
    assert profiletime.unit_seconds(conn, "link", token="r1", read_rows=rows(10.0, 12.0, 90.0)) == 12.0


def test_it_moves_a_quarter_a_run_and_is_held_for_the_run(conn):
    first = profiletime.unit_seconds(conn, "link", token="r1", typed_seconds=20.0, read_rows=rows(*[40.0] * 5))
    assert first == pytest.approx(25.0)                                    # 20 s typed, 40 s measured, limited to a quarter
    assert profiletime.unit_seconds(conn, "link", token="r1", read_rows=rows(*[400.0] * 5)) == first
    assert profiletime.unit_seconds(conn, "link", token="r2", read_rows=rows(*[40.0] * 5)) == pytest.approx(31.25)


def test_another_machines_samples_count_half(conn):
    assert profiletime.unit_seconds(conn, "link", token="r1", typed_seconds=20.0, read_rows=rows(*[30.0] * 5, shape=OTHER)) == pytest.approx(25.0)
    assert profiletime.target(rows(*[30.0] * 5, shape=OTHER)("x"), None, *HOST)[0] == 30.0


def test_off_and_a_broken_database_leave_the_readers_what_they_had(conn, monkeypatch):
    def broken(kind):
        raise RuntimeError("no control database")
    assert profiletime.unit_seconds(conn, "link", token="r1", read_rows=broken) is None
    monkeypatch.setenv(profileclaim.ENV_VAR, "0")
    assert profiletime.unit_seconds(conn, "link", token="r1", read_rows=rows(*[5.0] * 5)) is None


def test_a_change_is_logged_with_the_kind_old_new_samples_and_reason(conn, caplog):
    caplog.set_level(logging.INFO, logger="planetgen")
    profiletime.unit_seconds(conn, "link", token="r1", typed_seconds=20.0, read_rows=rows(*[40.0] * 5))
    assert "expected unit seconds of unit:link was 20.0, now 25.0" in caplog.text and "5 samples" in caplog.text
    (change,) = profileclaim.changes(conn)
    assert (change["what"], change["kind"], change["old"], change["samples"]) == ("expected unit seconds", "unit:link", 20.0, 5)


def test_the_untimed_estimate_uses_the_profiles_unit_time_before_the_typed_one():
    assert roughtime.stage_seconds(8, 4) == 2 * roughtime.UNIT_SECONDS
    assert roughtime.stage_seconds(8, 4, unit_seconds=5.0) == 10.0


def test_a_stage_carries_the_profiles_unit_time_to_the_eta():
    from planetgen.queue import runeta
    stage = {"units": 8, "workers": 4, "started": None, "unit_seconds": 5.0}
    assert runeta._rough_remaining(stage, 0.0)[0] == 10.0
    assert runeta._rough_remaining({**stage, "unit_seconds": None}, 0.0)[0] == 2 * roughtime.UNIT_SECONDS


def test_the_judge_flags_a_unit_far_slower_than_the_profiles_item_time():
    beat = beatjudge.Beat(beat=1000.0, started=100.0, done=100.0, estimate=10000.0, moved=990.0)
    slow = beatjudge.judge(beat, beatjudge.History(expected_item_seconds=0.01), 1000.0)          # 100 items a second expected, 0.11 shown
    assert slow.state == beatjudge.SLOW
    assert beatjudge.judge(beat, beatjudge.History(expected_item_seconds=5.0), 1000.0).state == beatjudge.WORKING


def test_the_judges_history_reads_the_profile_per_item(conn):
    profiletime.item_seconds(conn, "link", token="r1", read_rows=rows(*[0.5] * 5))
    history = commandgaps.history_of(conn, "r1")("link")
    assert history.expected_item_seconds == 0.5
