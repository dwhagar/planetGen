"""
OPS.264: `ctx.commit_piece_in(conn, result)`, the piece's cursor and result written inside the transaction of a save that opens a connection of its own (`store.save_sector`, `save_system`). The rows, the cursor
past the piece and the result land together or not at all, so a unit stopped after the save starts at the piece that follows and one stopped before it does the piece again; and the two commands that save so
(`generate-unplaced-sector`, `generate-system`) make their sector or system once.
"""

import types

import pytest

from planetgen.db import store
from planetgen.queue import attempts, commands, engine, ledger, wrapper

pytestmark = pytest.mark.fresh_database

PIECES = 4
DONE = []
STOP = {}


def listed(args, cursor):
    return list(range(args["count"]))[cursor:]


def saving(args, piece, ctx):
    """A save that opens its own connection and writes the piece's row and cursor in one transaction, killed just after it when asked."""
    DONE.append(piece)
    conn = store.get_connection(ctx.config)
    try:
        conn.execute("INSERT INTO commit_in_probe (piece) VALUES (?)", (piece,))
        result = {"piece": piece}
        ctx.commit_piece_in(conn, result)
        if STOP.get("before_commit") == piece:
            conn.rollback()
            raise RuntimeError("the save is lost before its commit")
        conn.commit()
    finally:
        conn.close()
    if STOP.get("after") == piece:
        raise RuntimeError("the unit is killed after the save")
    return result


def together(args, results, ctx):
    return {"results": results}


@pytest.fixture(autouse=True)
def clear():
    DONE.clear()
    STOP.clear()


@pytest.fixture
def job(mysql_config, monkeypatch):
    spec = commands.CommandSpec("commit-in-test", "A test", expand=(__name__, "listed"), pieces=(__name__, "listed"), do=(__name__, "saving"), gather=(__name__, "together"), resumes=True)
    monkeypatch.setitem(commands.REGISTRY, "commit-in-test", spec)
    control, galaxy = store.get_control_connection(store.control_mysql_config(mysql_config)), store.get_connection(mysql_config)
    galaxy.execute("CREATE TABLE IF NOT EXISTS commit_in_probe (piece INT NOT NULL)")
    galaxy.execute("DELETE FROM commit_in_probe")
    galaxy.commit()
    run = ledger.add_run(control, "generate", "g")
    ledger.add_jobs(control, run, [{"stage": "s", "command": "commit-in-test", "item_key": "0"}])
    control.commit()
    row = control.execute("SELECT id FROM engine_jobs WHERE run_id = ?", (run,)).fetchone()["id"]
    args = {"count": PIECES, "database": wrapper.database_args(mysql_config)}

    def take():
        return attempts.Attempt(control, row, ledger.begin_attempt(control, row, galaxy), mysql_config)

    try:
        yield control, galaxy, row, args, take
    finally:
        galaxy.execute("DROP TABLE IF EXISTS commit_in_probe")
        galaxy.commit()
        control.close()
        galaxy.close()


def _rows(galaxy):
    galaxy.commit()
    return [row["piece"] for row in galaxy.execute("SELECT piece FROM commit_in_probe ORDER BY piece").fetchall()]


def test_a_unit_killed_after_a_piece_saved_resumes_past_it_and_the_row_is_there_once(job):
    control, galaxy, row, args, take = job
    STOP["after"] = 1
    with pytest.raises(RuntimeError):
        engine.run("commit-in-test", args, lambda done, total: None, attempt=take())
    galaxy.commit()
    assert ledger.item_state(galaxy, row)["cursor_key"] == "2", "the cursor is past the piece whose save committed"
    STOP.clear()
    DONE.clear()
    result = engine.run("commit-in-test", args, lambda done, total: None, attempt=take())
    assert DONE == [2, 3], "the piece that saved is not done again"
    assert _rows(galaxy) == [0, 1, 2, 3]
    assert result["results"] == [{"piece": number} for number in range(PIECES)], "its stored result is the one gather reads"


def test_a_save_lost_before_its_commit_leaves_no_row_and_no_cursor_so_the_piece_is_done_again(job):
    control, galaxy, row, args, take = job
    STOP["before_commit"] = 2
    with pytest.raises(RuntimeError):
        engine.run("commit-in-test", args, lambda done, total: None, attempt=take())
    galaxy.commit()
    assert ledger.item_state(galaxy, row)["cursor_key"] == "2" and _rows(galaxy) == [0, 1]
    STOP.clear()
    DONE.clear()
    engine.run("commit-in-test", args, lambda done, total: None, attempt=take())
    assert DONE == [2, 3] and _rows(galaxy) == [0, 1, 2, 3]


def test_a_unit_another_attempt_took_over_saves_nothing(job):
    control, galaxy, row, args, take = job
    first = take()
    take()          # a second attempt replaces the first
    with pytest.raises(attempts.Replaced):
        engine.run("commit-in-test", args, lambda done, total: None, attempt=first)
    assert _rows(galaxy) == [], "the replaced attempt's save was rolled back, not committed"


def test_commit_piece_in_does_nothing_for_a_unit_with_no_attempt(job):
    control, galaxy, row, args, take = job
    ctx = engine.Context("commit-in-test", args, lambda done, total: None)
    conn = store.get_connection(store.MySQLConfig.from_args(args["database"]))
    try:
        assert ctx.commit_piece_in(conn, {"piece": 0}) is None
        assert conn.execute("SELECT COUNT(*) AS n FROM engine_piece_results").fetchone()["n"] == 0
    finally:
        conn.close()


# --- the two commands that save through it ---------------------------------------------------------------------------------------

def _unit(mysql_config, args, **more):
    return {"database": wrapper.database_args(mysql_config), "args": args, "run_seed": 7, **more}


@pytest.fixture
def named(mysql_config, monkeypatch):
    """`run_galaxy._sector_args_of` that names the test database and asks for nothing else."""
    from planetgen.generation import run_galaxy
    args = types.SimpleNamespace(density=1.0, num_systems=None, output=False, format="text", strict=False, mysql_host=mysql_config.host, mysql_port=mysql_config.port,
                                 mysql_user=mysql_config.user, mysql_password=mysql_config.password, mysql_database=mysql_config.database)
    monkeypatch.setattr(run_galaxy, "_sector_args_of", lambda unit: args)
    return args


def _kill_after(monkeypatch, name):
    """`store.<name>` that saves and then dies, once: the unit is killed between the save's commit and the engine's own checkpoint."""
    real = getattr(store, name)
    state = {"killed": False}

    def killing(*args, **kwargs):
        saved = real(*args, **kwargs)
        if not state["killed"]:
            state["killed"] = True
            raise RuntimeError("the unit is killed after the save")
        return saved

    monkeypatch.setattr(store, name, killing)


def _count(galaxy, table):
    galaxy.commit()
    return galaxy.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]


def test_a_sector_run_unit_killed_after_a_sector_saved_does_not_make_it_twice(job, mysql_config, named, monkeypatch):
    from planetgen.galaxy.sector import SpaceSector
    from planetgen.generation import run_sector
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem

    def made(args, **kwargs):
        sector = SpaceSector("Piece sector", edge_ly=10.0)
        config = SystemConfig()
        config.STAR_TYPE, config.PLANETS, config.BINARY_SYSTEM = "G2V", False, False
        sector.add_system(StarSystem(system_config=config), position=(1.0, 1.0, 1.0), system_config=config)
        return sector.name, sector

    monkeypatch.setattr(run_sector, "generate_directed_sector", made)
    monkeypatch.setattr(run_sector, "sector_generation_summary_lines", lambda sector, args: ["a line"])
    control, galaxy, row, args, take = job
    unit = _unit(mysql_config, {}, numbers=[0, 1, 2])
    # the second sector saves and dies: its row and the cursor past it were committed together
    real = store.save_sector
    seen = {"saved": 0}

    def killing_second(*a, **k):
        saved = real(*a, **k)
        seen["saved"] += 1
        if seen["saved"] == 2:
            raise RuntimeError("the unit is killed after the second sector")
        return saved

    monkeypatch.setattr(store, "save_sector", killing_second)
    with pytest.raises(RuntimeError):
        engine.run("generate-unplaced-sector", unit, lambda done, total: None, attempt=take())
    assert _count(galaxy, "sectors") == 2 and ledger.item_state(galaxy, row)["cursor_key"] == "2"
    monkeypatch.setattr(store, "save_sector", real)
    result = engine.run("generate-unplaced-sector", unit, lambda done, total: None, attempt=take())
    assert _count(galaxy, "sectors") == 3, "the sector saved before the kill is not made again"
    assert [saved["number"] for saved in result["sectors"]] == [0, 1, 2] and all(saved["sector_id"] for saved in result["sectors"])


def test_a_system_unit_killed_after_its_system_saved_does_not_save_it_twice(job, mysql_config, named, monkeypatch):
    from planetgen.generation import run_system
    from planetgen.generation.config import SystemConfig

    def config(args):
        made = SystemConfig()
        made.STAR_TYPE, made.PLANETS, made.BINARY_SYSTEM = "G2V", False, False
        return made

    monkeypatch.setattr(run_system, "build_system_config", config)
    control, galaxy, row, args, take = job
    unit = _unit(mysql_config, {})
    _kill_after(monkeypatch, "save_system")
    with pytest.raises(RuntimeError):
        engine.run("generate-system", unit, lambda done, total: None, attempt=take())
    assert _count(galaxy, "star_systems") == 1 and ledger.item_state(galaxy, row)["cursor_key"] == "2"
    result = engine.run("generate-system", unit, lambda done, total: None, attempt=take())
    assert _count(galaxy, "star_systems") == 1, "the system is not saved a second time"
    assert result["star_system_id"] and result["name"]


def test_a_system_unit_that_starts_at_its_finish_piece_draws_the_system_it_drew_before(job, mysql_config, named, monkeypatch):
    from planetgen.generation import run_system
    from planetgen.generation.config import SystemConfig

    def config(args):
        made = SystemConfig()
        made.STAR_TYPE, made.PLANETS, made.BINARY_SYSTEM = "G2V", False, False
        return made

    monkeypatch.setattr(run_system, "build_system_config", config)
    control, galaxy, row, args, take = job
    unit = _unit(mysql_config, {})
    drawn = []
    real = run_system.class_request_report

    def failing_once(system, system_config):
        drawn.append((system.name, system.stars[0].mass))
        if len(drawn) == 1:
            raise RuntimeError("the unit is killed before the finish piece saves")
        return real(system, system_config)

    monkeypatch.setattr(run_system, "class_request_report", failing_once)
    with pytest.raises(RuntimeError):
        engine.run("generate-system", unit, lambda done, total: None, attempt=take())
    assert ledger.item_state(galaxy, row)["cursor_key"] == "1", "the generate piece's cursor is stored, its state is gone"
    engine.run("generate-system", unit, lambda done, total: None, attempt=take())
    assert len(drawn) == 2 and drawn[0] == drawn[1], "the second try drew the system the first drew"
    assert _count(galaxy, "star_systems") == 1
