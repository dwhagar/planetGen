"""OPS.231: the one memory controller (`planetgen/queue/pressure.py`): the reading (PSI, or the host's available memory without it), the hold and the ask to the youngest bulk job, and what the launcher and
the light jobs' duty do with them. The launcher's rows are real; the templates are a fake forker."""

import os

import pytest

from planetgen.db import store
from planetgen.queue import jobhost, launcher, ledger, lightjobs, memorygate, pressure, redisqueue

from tests.test_launcher import FakeForker, HOST, _launcher, _ready, control  # noqa: F401 -- the fixture and the fakes


PSI = "some avg10={some} avg60=0.00 avg300=0.00 total=1\nfull avg10={full} avg60=0.00 avg300=0.00 total=1\n"


def _file(tmp_path, some, full=0.0):
    path = tmp_path / "memory"
    path.write_text(PSI.format(some=some, full=full), encoding="utf-8")
    return str(path)


def test_the_kernels_pressure_file_is_read_and_a_text_that_is_not_one_is_refused():
    assert pressure.parse_psi(PSI.format(some="12.50", full="3.00")) == (12.5, 3.0)
    assert pressure.parse_psi("some avg10=1.00 avg60=0 avg300=0 total=1\n") == (1.0, 0.0)
    assert pressure.parse_psi("") is None and pressure.parse_psi("some avg10=x\n") is None


def test_the_level_follows_psi_where_there_is_one_and_the_available_memory_against_the_reserve_where_there_is_not(tmp_path):
    host = memorygate.Host(16000.0, 9000.0, 8, 0.2)
    reserve = memorygate.reserve_mb(host)
    levels = [pressure.read(lambda: host, _file(tmp_path, some, full)).level for some, full in ((0.5, 0.0), (11.0, 0.0), (30.0, 0.0), (2.0, 6.0))]
    assert levels == [pressure.CALM, pressure.PRESSED, pressure.STARVED, pressure.STARVED]
    absent = str(tmp_path / "nothing")
    scarce = lambda mb: (lambda: memorygate.Host(16000.0, mb, 8, 0.2))      # noqa: E731
    assert pressure.read(scarce(reserve * 2), absent).level == pressure.CALM
    assert pressure.read(scarce(reserve * 0.7), absent).level == pressure.PRESSED
    assert pressure.read(scarce(reserve * 0.4), absent).level == pressure.STARVED
    assert pressure.read(scarce(reserve * 0.4), absent).source == "memory" and pressure.read(lambda: host, _file(tmp_path, 1)).source == "psi"


class Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


class Bulk:
    """A launcher that only says whom it would ask."""

    def __init__(self):
        self.asked = []

    def yield_youngest(self):
        self.asked.append(True)
        return 7, "scatter-mass"


def _controller(level, launcher_=None):
    clock, readings = Clock(), [level]
    done = []

    class Boss:
        def act(self, what, reason=""):
            done.append((what, reason))
    one = pressure.Controller(launcher=lambda: launcher_, look=lambda: pressure.Reading("psi", 20.0, 0.0, 1000.0, 2000.0, readings[0]), clock=clock)
    one.attach(Boss())
    return one, clock, readings, done


def test_lasting_pressure_holds_admission_and_lasting_calm_lets_it_go():
    one, clock, readings, done = _controller(pressure.PRESSED)
    one.turn()
    assert one.holding() is None, "one look is not lasting"
    clock.now += pressure.HOLD_AFTER_SECONDS
    one.turn()
    assert "short of memory" in one.holding() and done[-1][0] == "stopped admitting bulk and upkeep jobs"
    readings[0] = pressure.CALM
    clock.now += 1
    one.turn()
    clock.now += pressure.RELEASE_AFTER_SECONDS - 1
    one.turn()
    assert one.holding() is not None, "calm for less than the release time"
    clock.now += 2
    one.turn()
    assert one.holding() is None and done[-1][0] == "admitting again"


def test_a_pressure_that_eases_before_the_hold_never_holds():
    one, clock, readings, _done = _controller(pressure.PRESSED)
    one.turn()
    clock.now += pressure.HOLD_AFTER_SECONDS - 1
    readings[0] = pressure.CALM
    one.turn()
    readings[0] = pressure.PRESSED
    clock.now += 2
    one.turn()
    clock.now += pressure.HOLD_AFTER_SECONDS - 1
    one.turn()
    assert one.holding() is None, "the clock began again when it eased"


def test_lasting_starvation_under_a_hold_asks_the_youngest_bulk_job_once_a_cooldown():
    bulk = Bulk()
    one, clock, _readings, done = _controller(pressure.STARVED, bulk)
    one.turn()
    clock.now += pressure.HOLD_AFTER_SECONDS
    one.turn()
    assert bulk.asked == [] and one.holding()
    clock.now += pressure.YIELD_AFTER_SECONDS
    assert one.turn() == 7 and len(bulk.asked) == 1
    assert "asked job 7 (scatter-mass)" in done[-1][0]
    clock.now += pressure.YIELD_COOLDOWN_SECONDS - 1
    assert one.turn() is None and len(bulk.asked) == 1, "the first has not had its cooldown"
    clock.now += 2
    assert one.turn() == 7 and len(bulk.asked) == 2


def test_pressure_that_is_not_starvation_holds_but_asks_nobody():
    bulk = Bulk()
    one, clock, _readings, _done = _controller(pressure.PRESSED, bulk)
    one.turn()
    clock.now += pressure.HOLD_AFTER_SECONDS + pressure.YIELD_AFTER_SECONDS + 5
    one.turn()
    assert one.holding() and bulk.asked == []


def test_a_gate_that_is_off_holds_nothing(monkeypatch):
    monkeypatch.setenv(memorygate.DISABLE_ENV, "0")
    one, clock, _readings, _done = _controller(pressure.STARVED, Bulk())
    one.turn()
    clock.now += 1000
    one.turn()
    assert one.holding() is None


# -- the launcher ---------------------------------------------------------------------------------------------------------------

def test_a_held_launcher_hands_out_nothing_and_every_waiting_row_says_why(control):
    run = _ready(control)
    forker = FakeForker()
    one = launcher.Launcher(control, forker, read_host=lambda: HOST, hold=lambda: "admission is held: the host is short of memory")
    assert one.turn() == [] and forker.forked == []
    waited = {row["waited"] for row in control.execute("SELECT waited FROM engine_jobs WHERE run_id = ? AND item_key <> ''", (run,)).fetchall()}
    assert waited == {"admission is held: the host is short of memory"}
    one.hold = lambda: None
    assert len(one.turn()) == 3, "the hold lifted: the rows go"


def _running(control, forker, count=3):
    run = _ready(control, count=count)
    one = _launcher(control, forker)
    clock = iter(range(100, 200))
    one.clock = lambda: float(next(clock))
    assert len(one.turn()) == count
    return run, one


def test_the_youngest_bulk_unit_is_asked_to_stop_once_and_never_an_opening(control):
    forker = FakeForker()
    run, one = _running(control, forker)
    sent = []
    for held in (item[3] for item in forker.forked):
        held.send = lambda message, held=held: sent.append((held.id, message))
    asked = one.yield_youngest()
    assert asked[1] == "fill-unit" and sent == [(3, {"op": "stop"})], "the last one started"
    assert one.yield_youngest()[0] != asked[0] and len(sent) == 2, "then the one before it"
    one.yield_youngest()
    assert one.yield_youngest() is None and len(sent) == 3, "each is asked once"


def test_a_unit_that_stopped_at_the_ask_goes_back_to_ready_without_using_a_retry(control):
    forker = FakeForker()
    run, one = _running(control, forker, count=1)
    kind, runner, job, held = forker.forked[0]
    held.send = lambda message: None
    row_id = one.yield_youngest()[0]
    ledger.begin_attempt(control, row_id, pid=held.pid)
    control.execute("UPDATE engine_jobs SET state = 'waiting' WHERE id = ?", (row_id,))        # a stopped unit's own end leaves its row waiting (`attempts.end`)
    ledger.end_attempt(control, row_id, 1, "lost", error_kind=ledger.STOPPED)
    control.commit()
    held.said = [{"event": "unit", "result": {"stopped": True}}]
    forker.exits = [("plain", jobhost.Exit(held.id, held.pid, 0, None))]
    assert one.ended() == [row_id]
    row = control.execute("SELECT state, waited FROM engine_jobs WHERE id = ?", (row_id,)).fetchone()
    assert (row["state"], row["waited"]) == ("ready", launcher.YIELDED_TEXT)
    assert ledger.stopped_attempts(control, row_id) == 1
    again = one.turn()
    assert again == [row_id]
    assert forker.forked[1][2]["retries_left"] == 1, "the stop did not use the job's retry"


def test_a_stop_that_nobody_asked_for_is_settled_as_before(control):
    forker = FakeForker()
    _run, one = _running(control, forker, count=1)
    held = forker.forked[0][3]
    held.said = [{"event": "unit", "result": {"stopped": True}}]
    forker.exits = [("plain", jobhost.Exit(held.id, held.pid, 0, None))]
    one.ended()
    row = control.execute("SELECT state FROM engine_jobs WHERE item_key <> '' LIMIT 1").fetchone()
    assert row["state"] == "ready", "the old rule: a row left running goes once more"
