"""OPS.156: the memory claim read from the job profile (queue/profileclaim.py). No history gives today's claims, history lowers a claim, a memory kill raises it."""

import logging

import pytest

from planetgen.queue import memorygate, profileclaim, redisqueue, runrecord
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


def row(peak, shape=HOST, outcome="done"):
    return {"peak_mb": peak, "outcome": outcome, "host_memory_mb": shape[0], "host_cores": shape[1], "generator_epoch": shape[2]}


def history(*peaks, **kw):
    return lambda kind: [row(p, **kw) for p in peaks]


def test_the_margin_falls_from_125_at_three_samples_to_110_at_thirty():
    assert profileclaim.margin(3) == pytest.approx(1.25) and profileclaim.margin(30) == pytest.approx(1.10) and profileclaim.margin(100) == pytest.approx(1.10)
    assert profileclaim.margin(16) < 1.25 and profileclaim.margin(16) > 1.10


def test_the_95th_percentile_is_taken_between_ranks():
    assert profileclaim.percentile([100.0]) == 100.0
    assert profileclaim.percentile([float(n) for n in range(1, 101)]) == pytest.approx(95.05)


def test_no_history_gives_todays_claim(conn):
    assert profileclaim.claim(conn, "unit:x", 1000.0, token="r1", read_rows=history()) is None
    assert profileclaim.claim(conn, "unit:x", 1000.0, token="r1", read_rows=history(500.0, 500.0)) is None         # two samples are not enough
    assert memorygate.claim_mb(conn, "unit:x", 1000.0, known={}) == 1000.0


def test_history_lowers_a_claim_below_the_typed_figure_but_not_by_more_than_a_quarter_a_run(conn):
    first = profileclaim.claim(conn, "unit:x", 2000.0, token="r1", read_rows=history(*[500.0] * 5))
    assert first == pytest.approx(1500.0)                                    # wanted 625, limited to 75 percent of the typed 2000
    assert profileclaim.claim(conn, "unit:x", 2000.0, token="r1", read_rows=history(*[100.0] * 5)) == first       # frozen for the run
    second = profileclaim.claim(conn, "unit:x", 2000.0, token="r2", read_rows=history(*[500.0] * 5))
    assert second == pytest.approx(1125.0)
    third = profileclaim.claim(conn, "unit:x", 2000.0, token="r3", read_rows=history(*[500.0] * 5))
    assert third == pytest.approx(843.75)


def test_it_settles_on_the_profile_and_can_rise_by_a_quarter_a_run(conn):
    assert profileclaim.claim(conn, "unit:x", 1000.0, token="r1", read_rows=history(*[800.0] * 30)) == pytest.approx(880.0)       # 800 x 1.10
    assert profileclaim.claim(conn, "unit:x", 1000.0, token="r2", read_rows=history(*[2000.0] * 30)) == pytest.approx(880.0 * 1.25)


def test_a_memory_kill_raises_the_claim_at_once(conn):
    rows = lambda kind: [row(1500.0, outcome="killed_memory")] + [row(400.0)] * 5       # noqa: E731
    assert profileclaim.claim(conn, "unit:x", 1000.0, token="r1", read_rows=rows) == pytest.approx(1500.0 * memorygate.LEARN_MARGIN)


def test_a_profile_from_another_machine_counts_half(conn):
    rows = history(*[1000.0] * 5, shape=OTHER)
    assert profileclaim.claim(conn, "unit:x", 2000.0, token="r1", read_rows=rows) == pytest.approx(1625.0)       # (1250 + 2000) / 2, inside the quarter
    assert profileclaim.target([row(1000.0, OTHER)] * 5, None, *HOST)[0] == pytest.approx(1250.0)
    assert profileclaim.target([row(1000.0, OTHER)] * 5, 2000.0, *HOST)[0] == pytest.approx(1625.0)


def test_the_claim_goes_through_the_gates_claim_and_a_retry_still_grows_it(conn):
    assert memorygate.claim_mb(conn, "unit:x", 1000.0, known={}) == 1000.0
    profileclaim.claim(conn, "unit:x", 1000.0, token="r1", read_rows=history(*[800.0] * 30))
    assert memorygate.claim_mb(conn, "unit:x", 1000.0, known={"unit:x": 5000.0}, token="r1") == pytest.approx(880.0)    # the old largest peak no longer holds it up
    memorygate.grow(conn, "unit:x", "job-1", 880.0, 0.0)
    assert memorygate.claim_mb(conn, "unit:x", 1000.0, job_id="job-1", token="r1") == pytest.approx(880.0 * memorygate.GROW_FACTOR)


def test_off_gives_todays_claim(conn, monkeypatch):
    monkeypatch.setenv(profileclaim.ENV_VAR, "0")
    assert profileclaim.claim(conn, "unit:x", 1000.0, token="r1", read_rows=history(*[200.0] * 30)) is None


def test_a_profile_that_cannot_be_read_gives_the_typed_claim(conn):
    def broken(kind):
        raise RuntimeError("no control database")
    assert profileclaim.claim(conn, "unit:x", 1000.0, token="r1", read_rows=broken) is None


def test_every_change_is_logged_and_kept_for_the_queue_page(conn, caplog):
    caplog.set_level(logging.INFO, logger="planetgen")
    profileclaim.claim(conn, "unit:x", 1000.0, token="r1", read_rows=history(*[800.0] * 30))
    assert "claim of unit:x was 1,000 MB, now 880 MB" in caplog.text and "30 samples" in caplog.text
    (change,) = profileclaim.changes(conn)
    assert (change["kind"], change["old"], change["new"], change["source"], change["samples"]) == ("unit:x", 1000.0, pytest.approx(880.0), "own", 30)
