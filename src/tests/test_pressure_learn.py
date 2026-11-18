"""OPS.249: the memory controller's edges are learned from what it recorded (`planetgen/queue/pressurelearn.py`, control schema v25)."""

import pytest

from planetgen.db import store
from planetgen.queue import pressure, pressurelearn
from planetgen.queue.pressurelearn import Minute, learn

HOUR = 3600.0


def minute(start, some=0.0, level=0, held=0.0, yields=0, full=0.0, psi=True, available=8000.0, reserve=2000.0):
    return Minute(start, some, full, available, reserve, level, held, yields, psi)


def calm_hours(hours):
    return [minute(60.0 * n) for n in range(int(hours * 60))]


def test_no_record_gives_the_typed_edges():
    edges = learn([], [])
    assert edges == pressurelearn.typed() and edges.sensitivity == 1.0
    assert (edges.hold_some, edges.yield_some, edges.yield_full) == (pressure.HOLD_SOME, pressure.YIELD_SOME, pressure.YIELD_FULL)
    assert (edges.pressure_share, edges.starved_share) == (pressure.PRESSURE_SHARE, pressure.STARVED_SHARE)


def test_kills_that_the_hold_did_not_see_coming_make_the_edges_more_sensitive():
    minutes = calm_hours(6)
    kills = []
    for at in (3600.0, 7200.0, 10800.0):           # pressure of 5 percent before each kill, under the typed edge of 10
        minutes[int(at // 60) - 1] = minute(at - 60.0, some=5.0)
        kills.append(at)
    edges = learn(minutes, kills)
    assert edges.sensitivity > 1.0 and edges.hold_some < pressure.HOLD_SOME
    assert edges.hold_some <= 5.0 / pressurelearn.MARGIN + 1e-9, "the reading before a kill is now pressed, with the margin to spare"
    assert edges.pressure_share > pressure.PRESSURE_SHARE and edges.yield_some < pressure.YIELD_SOME
    assert "3 memory kills" in edges.words()


def test_too_few_kills_or_kills_the_hold_saw_say_nothing():
    minutes = calm_hours(4)
    minutes[119] = minute(7140.0, some=5.0)
    assert learn(minutes, [7200.0, 7260.0]).sensitivity == 1.0                         # two, under MIN_KILLS
    seen = calm_hours(4)
    for at in (3600.0, 7200.0, 10800.0):
        seen[int(at // 60) - 1] = minute(at - 60.0, some=40.0, level=2, held=60.0)       # the hold had them in sight
    assert learn(seen, [3600.0, 7200.0, 10800.0]).sensitivity == 1.0
    quiet = calm_hours(4)
    assert learn(quiet, [3600.0, 7200.0, 10800.0]).sensitivity == 1.0                    # a kill on a calm host: no edge could have caught it


def test_the_sensitivity_never_passes_its_clamp():
    minutes = calm_hours(6)
    kills = []
    for at in (3600.0, 7200.0, 10800.0):
        minutes[int(at // 60) - 1] = minute(at - 60.0, some=0.1)          # 1 percent of the edge
        kills.append(at)
    assert learn(minutes, kills).sensitivity == pressurelearn.MAX_SENSITIVITY


def _holds(count, starved_every=None, spacing=1800.0, length=5):
    minutes = calm_hours(count * spacing / HOUR + 1)
    for n in range(count):
        for step in range(length):
            index = int((n * spacing) // 60) + step
            starved = starved_every is not None and n % starved_every == 0
            minutes[index] = minute(60.0 * index, some=14.0, level=2 if starved and step == 2 else 1, held=60.0, yields=1 if starved and step == 2 else 0)
    return minutes


def test_holds_that_were_never_needed_and_no_kills_raise_the_edges():
    edges = learn(_holds(10), [])
    assert edges.sensitivity < 1.0 and edges.sensitivity >= pressurelearn.MIN_SENSITIVITY
    assert edges.hold_some > pressure.HOLD_SOME and edges.hold_some <= pressure.HOLD_SOME * pressurelearn.MAX_RAISE + 1e-9
    assert edges.pressure_share < pressure.PRESSURE_SHARE
    assert "10 of 10 holds" in edges.words()


def test_holds_that_ended_in_a_yield_or_a_starved_host_keep_the_typed_edges():
    assert learn(_holds(10, starved_every=2), []).sensitivity == 1.0       # half were needed: under IDLE_SHARE idle
    assert learn(_holds(5), []).sensitivity == 1.0                         # under MIN_EPISODES
    assert learn(_holds(10), [3600.0]).sensitivity == 1.0                  # a kill in the window: the edge is not raised


def test_without_psi_the_reserve_is_the_signal():
    minutes = calm_hours(6)
    kills = []
    for at in (3600.0, 7200.0, 10800.0):
        # available 1900 MB of a 2000 MB reserve: pressed would be under 1500, so the strength is 0.79
        minutes[int(at // 60) - 1] = minute(at - 60.0, psi=False, available=1900.0, reserve=2000.0)
        kills.append(at)
    edges = learn(minutes, kills)
    assert edges.sensitivity > 1.0
    assert edges.pressure_share * 2000.0 > 1900.0 and edges.pressure_share <= 1.0, "1,900 MB of the reserve is now pressed"


def test_a_reading_judged_by_learned_edges_is_pressed_where_the_typed_ones_are_calm():
    host = pressure.memorygate.Host(16000.0, 8000.0, 8, 0.1)
    edges = pressurelearn.scaled(2.0, "test")
    assert (edges.hold_some, edges.yield_some) == (pressure.HOLD_SOME / 2.0, pressure.YIELD_SOME / 2.0)
    reading = pressure.read(read_host=lambda: host, path="/nonexistent", edges=None)
    assert reading.level == pressure.CALM

    class Psi:
        pass
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".psi", delete=False) as handle:
        handle.write("some avg10=6.00 avg60=0.00 avg300=0.00 total=0\nfull avg10=0.00 avg60=0.00 avg300=0.00 total=0\n")
    assert pressure.read(read_host=lambda: host, path=handle.name).level == pressure.CALM
    assert pressure.read(read_host=lambda: host, path=handle.name, edges=edges).level == pressure.PRESSED


class FakeConn:
    def __init__(self):
        self.statements = []

    def execute(self, sql, params=()):
        self.statements.append((sql, params))
        return self

    def commit(self):
        pass

    def close(self):
        pass


def _reading(some=0.0, level=0, source="psi"):
    return pressure.Reading(source, some, 0.0, 8000.0, 2000.0, level)


def test_the_recorder_sums_a_minute_up_in_one_row_and_writes_it_when_the_minute_ends():
    conn, clock = FakeConn(), [600.0]
    recorder = pressurelearn.Recorder(lambda: conn, clock=lambda: clock[0])
    for at, some, level in ((600.0, 2.0, 0), (610.0, 12.0, 1), (620.0, 30.0, 2), (630.0, 4.0, 0)):
        clock[0] = at
        recorder.note(_reading(some, level), held=at >= 620.0)
    assert conn.statements == []                                            # the minute is still open
    clock[0] = 660.0
    recorder.note(_reading(1.0, 0), held=False)
    (sql, params), = conn.statements
    assert "INSERT INTO pressure_minutes" in sql
    assert params[0] == 600.0 and params[1] == 30.0 and params[5] == 2      # the minute, the highest some, the highest level
    assert params[6] == pytest.approx(20.0) and params[8] == 1              # held for the 20 s from 610 to 630, read by PSI


def test_the_recorder_never_raises_and_backs_off_after_a_failed_write():
    calls, clock, mono = [], [0.0], [0.0]

    def broken():
        calls.append(1)
        raise OSError("down")
    recorder = pressurelearn.Recorder(broken, clock=lambda: clock[0], monotonic=lambda: mono[0])
    recorder.note(_reading(), False)
    clock[0] = 60.0
    recorder.note(_reading(), False)                      # writes the first minute: fails, and the writer leaves the database alone for BACKOFF_SECONDS
    clock[0], mono[0] = 120.0, 30.0
    recorder.note(_reading(), False)                      # the second minute is dropped while backed off
    assert len(calls) == 1
    mono[0] = 61.0
    clock[0] = 180.0
    recorder.note(_reading(), False)                      # after the back-off it tries again
    assert len(calls) == 2
    assert recorder.edges() == pressurelearn.typed()      # and the typed edges stand while the record cannot be read


def test_off_records_nothing_and_gives_the_typed_edges(monkeypatch):
    monkeypatch.setenv(pressurelearn.ENV_VAR, "0")
    conn = FakeConn()
    recorder = pressurelearn.Recorder(lambda: conn, clock=lambda: 600.0)
    recorder.note(_reading(), False)
    recorder.current = None
    assert conn.statements == [] and recorder.edges() == pressurelearn.typed()


def test_the_controller_judges_by_the_learned_edges_and_says_when_they_move(monkeypatch):
    said, seen = [], []

    class Supervisor:
        def act(self, what, reason=""):
            said.append((what, reason))

    class Recorder:
        edges_now = pressurelearn.typed()

        def edges(self):
            return self.edges_now

        def note(self, *args, **kwargs):
            seen.append(args)

    def look(edges=None, **kwargs):
        seen.append(("edges", edges))
        return pressure.Reading("psi", 0.0, 0.0, 8000.0, 2000.0, 0)

    monkeypatch.setattr(pressure, "read", look)
    recorder = Recorder()
    controller = pressure.Controller(look=look, recorder=recorder, clock=lambda: 0.0)
    controller.attach(Supervisor())
    controller.turn()
    assert said == [] and ("edges", recorder.edges_now) in seen
    recorder.edges_now = pressurelearn.scaled(1.5, "test edges")
    controller.turn()
    assert said == [("moved the memory edges", "test edges")] and ("edges", recorder.edges_now) in seen
    own = pressure.Controller(look=lambda: pressure.Reading("db", 0.0, 0.0, 0.0, 0.0, 0), recorder=recorder, clock=lambda: 0.0)
    own.turn()
    assert not any(args[0] == "db" for args in seen if args and isinstance(args[0], str)), "a look of the caller's own is not recorded as the host's"


@pytest.fixture
def control(mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    conn.execute("DELETE FROM pressure_minutes")
    conn.commit()
    try:
        yield conn
    finally:
        conn.close()


def test_the_record_is_written_and_read_back_from_the_control_database(control):
    recorder = pressurelearn.Recorder(lambda: _Same(control))
    recorder.clock = lambda: 1_800_000_000.0
    recorder.note(_reading(12.0, 1), True)
    recorder.clock = lambda: 1_800_000_060.0
    recorder.note(_reading(3.0, 0), False)
    row = control.execute("SELECT some_max, level_max, psi FROM pressure_minutes").fetchone()
    assert (row["some_max"], row["level_max"], row["psi"]) == (12.0, 1, 1)
    minutes, kills = pressurelearn.read_record(control)
    assert kills == [] and [(m.some, m.level, m.psi) for m in minutes] == [(12.0, 1, True)]


class _Same:
    """A connection whose `close` leaves the test's own open."""

    def __init__(self, conn):
        self.conn = conn

    def execute(self, *args):
        return self.conn.execute(*args)

    def commit(self):
        self.conn.commit()

    def close(self):
        pass
