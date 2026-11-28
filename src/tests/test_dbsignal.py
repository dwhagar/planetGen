"""OPS.248: MySQL as a signal (`planetgen/queue/dbsignal.py`): the calm figure the database learns of itself, the reading of a probe against it, a probe still out reading as slow, and the hold that is the memory
controller's own machine with this signal. The probe statement is run against the real control schema once."""

import pytest

from planetgen.queue import dbsignal, pressure

from planetgen.db import store

from tests.test_launcher import control  # noqa: F401 -- the fixture


class Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


class Conn:
    """A connection whose statement takes `took` seconds of the clock."""

    def __init__(self, clock, took, fail=False):
        self.clock, self.took, self.fail, self.sent, self.ended, self.closed = clock, took, fail, [], 0, False

    def execute(self, sql):
        self.sent.append(sql)
        self.clock.now += self.took
        if self.fail:
            raise OSError("the database is gone")
        return self

    def fetchone(self):
        return None

    def rollback(self):
        self.ended += 1

    def close(self):
        self.closed = True


def _learned(baseline, ms, times=dbsignal.BASELINE_PROBES):
    for _ in range(times):
        baseline.level(ms)
    return baseline


def test_nothing_is_pressed_until_the_calm_figure_has_been_learned_but_a_probe_over_the_starved_floor_still_counts():
    baseline = dbsignal.Baseline()
    assert baseline.level(200.0) == pressure.CALM, "no calm figure yet to judge a slow one by, and it is under the starved floor"
    assert baseline.level(dbsignal.STARVED_FLOOR_MS) == pressure.STARVED
    assert baseline.count == 1 and baseline.level(1.0) == pressure.CALM


def test_a_probe_is_pressed_and_starved_at_a_factor_of_the_calm_figure_and_never_under_the_floors():
    fast = _learned(dbsignal.Baseline(), 1.0)
    assert fast.level(20.0) == pressure.CALM, "twenty times calm, but under the pressed floor"
    assert fast.level(dbsignal.PRESSED_FLOOR_MS) == pressure.PRESSED
    assert fast.level(dbsignal.STARVED_FLOOR_MS) == pressure.STARVED
    slow = _learned(dbsignal.Baseline(), 40.0)
    assert slow.level(100.0, learn=False) == pressure.CALM and slow.level(dbsignal.PRESSED_FACTOR * 40.0) == pressure.PRESSED
    assert slow.level(dbsignal.STARVED_FACTOR * 40.0) == pressure.STARVED


def test_a_pressed_probe_teaches_nothing_and_a_probe_counted_twice_teaches_once():
    baseline = _learned(dbsignal.Baseline(), 2.0)
    before = (baseline.mean_ms, baseline.count)
    assert baseline.level(500.0) == pressure.STARVED and (baseline.mean_ms, baseline.count) == before
    assert baseline.level(2.0, learn=False) == pressure.CALM and (baseline.mean_ms, baseline.count) == before
    assert baseline.level(2.0) == pressure.CALM and baseline.count == before[1] + 1


def test_the_prober_times_the_statement_and_a_probe_still_out_reads_as_the_time_it_has_been_out():
    clock = Clock()
    prober = dbsignal.Prober(lambda: None, clock=clock)
    assert prober.latest() is None
    conn = Conn(clock, 0.004)
    prober.probe(conn)
    assert conn.sent == [dbsignal.PROBE_SQL] and prober.latest() == pytest.approx(4.0) and prober.completed == 1
    prober.sent_at = clock.now          # a probe is sent and does not come back
    clock.now += 2.0
    assert prober.latest() == pytest.approx(2000.0), "the figure is the time it has been out"
    clock.now += 1.0
    assert prober.latest() == pytest.approx(3000.0)


def test_the_look_learns_each_probe_once_however_often_the_controller_looks():
    clock = Clock()
    prober = dbsignal.Prober(lambda: None, clock=clock)
    look = dbsignal.Look(prober)
    assert look().level == pressure.CALM and look.baseline.count == 0, "no probe yet"
    prober.probe(Conn(clock, 0.003))
    for _ in range(4):
        look()
    assert look.baseline.count == 1
    prober.probe(Conn(clock, 0.003))
    look()
    assert look.baseline.count == 2


def test_the_look_tells_the_record_each_probe_it_judged_and_nothing_before_the_first():
    clock, told = Clock(), []
    prober = dbsignal.Prober(lambda: None, clock=clock)
    look = dbsignal.Look(prober, record=lambda latency, calm: told.append((round(latency), round(calm))))
    look()
    assert told == [], "no probe yet"
    prober.probe(Conn(clock, 0.004))
    look()
    prober.probe(Conn(clock, 0.040))
    look()
    assert told == [(4, 0), (40, 4)], "each probe with the calm figure it was judged against, before it taught it anything"


def test_a_database_that_stops_answering_holds_admission_and_lets_go_when_it_answers_again():
    clock, wall = Clock(), Clock()
    prober = dbsignal.Prober(lambda: None, clock=clock)
    controller = dbsignal.controller(prober, clock=wall)
    said = []

    class Boss:
        def act(self, what, reason=""):
            said.append((what, reason))
    controller.attach(Boss())
    for _ in range(dbsignal.BASELINE_PROBES + 1):
        prober.probe(Conn(clock, 0.002))
        controller.turn()
        wall.now += dbsignal.PROBE_SECONDS
    assert controller.holding() is None
    prober.sent_at = clock.now          # the next probe goes out and does not come back
    for _ in range(int(pressure.HOLD_AFTER_SECONDS) + 2):
        clock.now += 1.0
        wall.now += 1.0
        controller.turn()
    assert "the database is the limit" in controller.holding() and said[-1][0] == "stopped admitting bulk and upkeep jobs"
    assert "a ledger read took" in said[-1][1]
    prober.sent_at = None
    prober.probe(Conn(clock, 0.002))
    for _ in range(int(pressure.RELEASE_AFTER_SECONDS) + 2):
        wall.now += 1.0
        controller.turn()
    assert controller.holding() is None and said[-1][0] == "admitting again"


def test_the_memory_controllers_own_reason_is_unchanged():
    one = pressure.Controller(look=lambda: pressure.Reading("psi", 20.0, 0.0, 1000.0, 2000.0, pressure.PRESSED), clock=Clock())
    one.turn()
    one.clock.now += pressure.HOLD_AFTER_SECONDS
    one.turn()
    assert "the host is short of memory" in one.holding()


def test_the_probe_statement_runs_against_the_ledger_table(control):
    prober = dbsignal.Prober(lambda: control)
    prober.probe(control)
    assert prober.latest() is not None and prober.completed == 1


def test_a_probe_ends_its_transaction_whether_it_was_answered_or_not():
    """OPS.323: the probe's connection lives as long as the supervisor and does not autocommit, so a read left its transaction open for good."""
    clock = Clock()
    prober = dbsignal.Prober(lambda: None, clock=clock)
    answered = Conn(clock, 0.001)
    prober.probe(answered)
    assert answered.ended == 1
    gone = Conn(clock, 0.001, fail=True)
    with pytest.raises(OSError):
        prober.probe(gone)
    assert gone.ended == 1


def test_a_connection_that_failed_is_closed_and_the_one_at_the_end_is_closed_too(monkeypatch):
    """OPS.323: a dropped connection kept its transaction until the garbage collector found it."""
    clock = Clock()
    made = [Conn(clock, 0.001, fail=True), Conn(clock, 0.001)]
    prober = dbsignal.Prober(lambda: made[len([c for c in made if c.sent])], clock=clock)
    monkeypatch.setattr(dbsignal, "PROBE_SECONDS", 0.0)
    seen = []
    original = prober.probe

    def probe(conn):
        original(conn)
        seen.append(conn)
        if conn is made[1]:
            prober.stop.set()

    prober.probe = probe
    prober.run()
    assert made[0].closed and made[1].closed and made[1].ended == 1


def test_a_second_connection_can_alter_the_ledger_table_while_the_probe_has_run(control, mysql_config):
    """OPS.323: after a probe the connection holds no metadata lock on `engine_jobs` (it held a shared one for the life of the supervisor, which stopped a test's DROP for 81 minutes)."""
    prober = dbsignal.Prober(lambda: control)
    prober.probe(control)
    other = store.get_control_connection(mysql_config)
    try:
        other.execute("SET SESSION lock_wait_timeout = 3")
        other.execute("ALTER TABLE engine_jobs COMMENT = 'ops323'")
        other.execute("ALTER TABLE engine_jobs COMMENT = ''")
    finally:
        other.close()
