"""OPS.191: the engine owns the loop of a command with `pieces` and `do` (planetgen/queue/engine.py): counts into the beat, a cancel and a pause honoured between two pieces,
`ctx` with the database and the one commit, and a command without `pieces` run as it always was."""

import threading
import time

import pytest

from planetgen.db import store
from planetgen.queue import commands, engine, pause, redisqueue, runrecord, units
from planetgen.queue.runrecord import Stage

DONE = []


def listed(args, cursor):
    DONE.append(("pieces", cursor))
    return list(range(args["count"]))


def one(args, piece, ctx):
    DONE.append(("do", piece, ctx.cursor))
    return piece * 10


def together(args, results):
    return {"sum": sum(results)}


def committing(args, piece, ctx):
    ctx.db.execute("SELECT 1")
    ctx.commit()
    return None


STATE = {}


def fails(args, piece, ctx):
    ctx.db
    raise RuntimeError("boom")


def sql_piece(args, piece, ctx):
    return ctx.sql("DELETE FROM t", "Emptying t", ("x",))


def cancel_after_first(args, piece, ctx):
    DONE.append(("do", piece))
    runrecord.request_cancel(STATE["conn"], STATE["run_id"])
    STATE["now"] += 2.0
    return piece


def tick(args, piece, ctx):
    STATE["now"] += 0.5
    return piece


def own_unit(args, report):
    report(1, 1)
    return {"own": True}


def never(*_):
    raise AssertionError("not called")


@pytest.fixture(autouse=True)
def clear():
    DONE.clear()
    STATE.clear()


@pytest.fixture
def command(monkeypatch):
    def make(kind="engine-test", **fields):
        spec = commands.CommandSpec(kind, "A test", expand=(__name__, "never"), **{"pieces": (__name__, "listed"), "do": (__name__, "one"), **fields})
        monkeypatch.setitem(commands.REGISTRY, kind, spec)
        return spec
    return make


class Report:
    def __init__(self):
        self.seen = []

    def __call__(self, done, estimate):
        self.seen.append((done, estimate))


def test_the_engine_counts_the_pieces_done_of_the_pieces_in_all_into_the_beat(command):
    command()
    report = Report()
    assert engine.run("engine-test", {"count": 3}, report) == [0, 10, 20]
    assert report.seen == [(0, 3), (1, 3), (2, 3), (3, 3)]


def test_the_pieces_are_asked_for_after_the_cursor_and_every_piece_is_done_in_order(command):
    command()
    engine.run("engine-test", {"count": 2}, Report())
    assert DONE == [("pieces", 0), ("do", 0, 0), ("do", 1, 0)]


def test_a_command_may_make_its_result_of_the_results_of_its_pieces(command):
    command(gather=(__name__, "together"))
    assert engine.run("engine-test", {"count": 4}, Report()) == {"sum": 60}


def test_a_unit_with_no_piece_reports_nothing_done_of_nothing_and_ends(command):
    command()
    report = Report()
    assert engine.run("engine-test", {"count": 0}, report) == []
    assert report.seen == [(0, 0)]


def test_a_command_with_its_own_unit_function_runs_as_it_always_did(monkeypatch):
    monkeypatch.setitem(commands.REGISTRY, "own-test", commands.CommandSpec("own-test", "A test", unit=(__name__, "own_unit"), expand=(__name__, "never")))
    report = Report()
    assert engine.run("own-test", {}, report) == {"own": True}
    assert report.seen == [(1, 1)]


def test_a_command_run_alone_goes_through_the_engine(command):
    command()
    assert commands.run_alone("engine-test", {"count": 2}) == [0, 10]


def test_pieces_and_do_go_together_and_replace_the_unit_function():
    with pytest.raises(ValueError, match="go together"):
        commands.CommandSpec("bad", "x", pieces=(__name__, "listed"))
    with pytest.raises(ValueError, match="go together"):
        commands.CommandSpec("bad", "x", do=(__name__, "one"))
    with pytest.raises(ValueError, match="not both"):
        commands.CommandSpec("bad", "x", unit=(__name__, "own_unit"), pieces=(__name__, "listed"), do=(__name__, "one"))
    with pytest.raises(ValueError, match="needs `pieces`"):
        commands.CommandSpec("bad", "x", unit=(__name__, "own_unit"), gather=(__name__, "together"))


def test_a_command_with_pieces_is_not_inline_and_has_no_unit_function(command):
    spec = command()
    assert not spec.inline
    with pytest.raises(ValueError, match="engine"):
        commands.unit_function("engine-test")


def test_ctx_gives_the_database_once_commits_through_ctx_and_closes_it_at_the_end(command, monkeypatch):
    command(do=(__name__, "committing"))
    seen = []

    class Connection:
        def execute(self, statement, params=()):
            seen.append(("execute", statement))

        def commit(self):
            seen.append(("commit",))

        def close(self):
            seen.append(("close",))

    opened = []
    monkeypatch.setattr(store, "get_connection", lambda config, **kw: opened.append(config) or Connection())
    args = {"count": 2, "database": {"host": "h", "port": 3306, "user": "u", "database": "d"}}
    engine.run("engine-test", args, Report())
    assert len(opened) == 1 and opened[0].database == "d"
    assert seen == [("execute", "SELECT 1"), ("commit",), ("execute", "SELECT 1"), ("commit",), ("close",)]


def test_the_connection_is_closed_when_a_piece_fails_and_a_unit_that_never_asked_opens_none(command, monkeypatch):
    command(do=(__name__, "fails"))
    closed = []

    class Connection:
        def close(self):
            closed.append(True)

    monkeypatch.setattr(store, "get_connection", lambda config, **kw: Connection())
    with pytest.raises(RuntimeError, match="boom"):
        engine.run("engine-test", {"count": 1, "database": {"host": "h", "port": 1, "user": "u", "database": "d"}}, Report())
    assert closed == [True]
    command()
    monkeypatch.setattr(store, "get_connection", never)
    engine.run("engine-test", {"count": 1}, Report())          # no `ctx.db` asked for: no connection opened, and `ctx.commit()` has nothing to do


def test_ctx_sql_names_the_statement_in_the_beat(command, monkeypatch):
    from planetgen.queue import longsql

    command(do=(__name__, "sql_piece"))
    seen = []
    monkeypatch.setattr(longsql, "run", lambda conn, statement, label, params=(), config=None: seen.append((conn, statement, label, params, config.database)) or "cursor")
    monkeypatch.setattr(store, "get_connection", lambda config, **kw: type("C", (), {"close": lambda self: None})())
    assert engine.run("engine-test", {"count": 1, "database": {"host": "h", "port": 1, "user": "u", "database": "d"}}, Report()) == ["cursor"]
    assert seen[0][1:] == ("DELETE FROM t", "Emptying t", ("x",), "d")


@pytest.fixture
def started(redis_server):
    conn = redisqueue.connect(redis_server)
    run_id = runrecord.start_run(conn, "Engine test", [Stage("s", "Stage", units=1)], who="boss", now=100.0)
    runrecord.begin_stage(conn, run_id, 0, now=101.0)
    return conn, run_id


def test_a_run_that_is_canceled_ends_the_unit_between_two_pieces_with_the_pieces_done_kept(command, started):
    conn, run_id = started
    STATE.update(conn=conn, run_id=run_id, now=0.0)
    command(do=(__name__, "cancel_after_first"))
    report = Report()
    with pytest.raises(engine.UnitStopped, match="canceling"):
        engine.run("engine-test", {"count": 5}, report, redis=conn, run_id=run_id, stage=0, unit=0, job=units.job_id(run_id, 0, 0), clock=lambda: STATE["now"])
    assert DONE == [("pieces", 0), ("do", 0)]
    assert report.seen == [(0, 5), (1, 5)]


def test_the_run_is_read_at_most_once_a_second_not_at_every_piece(command, started, monkeypatch):
    conn, run_id = started
    command()
    reads = []
    real = conn.hget
    monkeypatch.setattr(conn, "hget", lambda *a, **k: reads.append(a) or real(*a, **k))
    now = [0.0]
    engine.run("engine-test", {"count": 50}, Report(), redis=conn, run_id=run_id, stage=0, unit=0, job="j", clock=lambda: now[0])
    assert len(reads) == 1            # the clock never moved
    reads.clear()

    command(do=(__name__, "tick"))
    STATE["now"] = 0.0
    engine.run("engine-test", {"count": 10}, Report(), redis=conn, run_id=run_id, stage=0, unit=0, job="j", clock=lambda: STATE["now"])
    assert 4 <= len(reads) <= 6       # ten pieces of half a second: one look a second


def test_a_pause_holds_the_unit_between_two_pieces_until_the_queues_are_let_go(command, started, monkeypatch):
    conn, run_id = started
    command()
    monkeypatch.setattr(engine, "PAUSE_POLL_SECONDS", 0.05)
    pause.pause("boss", conn)
    done = []
    thread = threading.Thread(target=lambda: done.append(engine.run("engine-test", {"count": 2}, Report(), redis=conn, run_id=run_id, stage=0, unit=0, job="j")))
    thread.start()
    time.sleep(0.4)
    assert thread.is_alive() and DONE == [("pieces", 0)]          # held before the first piece
    pause.resume(conn)
    thread.join(timeout=10)
    assert done == [[0, 10]]


def test_a_pause_is_shown_in_the_beat_as_waiting_on_pause(command, started, monkeypatch):
    from planetgen.queue import unitbeat
    conn, run_id = started
    command()
    monkeypatch.setattr(engine, "PAUSE_POLL_SECONDS", 0.05)
    pause.pause("boss", conn)
    beat = unitbeat.UnitBeat(conn, run_id, 0, 0)
    with beat:
        thread = threading.Thread(target=lambda: engine.run("engine-test", {"count": 1}, beat, redis=conn, run_id=run_id, stage=0, unit=0, job="j"))
        thread.start()
        time.sleep(0.3)
        assert beat.record()["waiting_on"] == engine.PAUSE_TAG
        pause.resume(conn)
        thread.join(timeout=10)
    assert beat.record()["waiting_on"] == ""


def test_a_run_canceled_while_the_unit_is_held_by_a_pause_ends_the_hold(command, started, monkeypatch):
    conn, run_id = started
    command()
    monkeypatch.setattr(engine, "PAUSE_POLL_SECONDS", 0.05)
    pause.pause("boss", conn)
    caught = []

    def hold():
        try:
            engine.run("engine-test", {"count": 2}, Report(), redis=conn, run_id=run_id, stage=0, unit=0, job="j")
        except engine.UnitStopped as stopped:
            caught.append(str(stopped))

    thread = threading.Thread(target=hold)
    thread.start()
    time.sleep(0.2)
    runrecord.request_cancel(conn, run_id)
    thread.join(timeout=10)
    assert caught == ["the run is canceling"]
    pause.resume(conn)


def test_a_solo_hold_lets_its_own_unit_through(command, started):
    conn, run_id = started
    command()
    pause.hold_solo(run_id, 0, 0, units.job_id(run_id, 0, 0), conn)
    assert engine.run("engine-test", {"count": 2}, Report(), redis=conn, run_id=run_id, stage=0, unit=0, job=units.job_id(run_id, 0, 0)) == [0, 10]
    pause.resume(conn)


def test_a_unit_of_a_run_that_ended_leaves_the_stage_unfinished_and_says_so(command, started):
    conn, run_id = started
    command()
    runrecord.request_cancel(conn, run_id)
    assert units._do_unit(conn, run_id, 0, 0, "engine-test", {"count": 3}, 1.0) == {"stopped": True}
    assert not conn.sismember(runrecord.done_key(run_id, 0), 0)


def test_a_unit_of_the_engine_is_counted_and_done_like_any_other(command, started, monkeypatch):
    conn, run_id = started
    command(gather=(__name__, "together"))
    monkeypatch.setattr(units, "peak_resident_mb", lambda: 100.0)
    units._do_unit(conn, run_id, 0, 0, "engine-test", {"count": 3}, 1.0)
    assert conn.sismember(runrecord.done_key(run_id, 0), 0)
    import json
    assert json.loads(conn.hget(runrecord.results_key(run_id, 0), 0)) == {"sum": 30}
