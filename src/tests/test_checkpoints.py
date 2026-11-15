"""OPS.219: checkpoints. `ctx.commit()` stores the cursor of the last finished piece and what the finished pieces returned in the transaction that holds the data, and a unit of a command that `resumes`
that is run again starts after them and gives the result an undisturbed unit gives (`planetgen/queue/engine.py`, `ledger.resume_point`)."""

import pytest

from planetgen.db import store
from planetgen.queue import attempts, commands, engine, ledger, wrapper

pytestmark = pytest.mark.fresh_database

PIECES = 6
DONE = []
STOP_AT = {}


def listed(args, cursor):
    return list(range(args["count"]))[cursor:]


def probe(args, piece, ctx):
    """Writes a row for the piece, in the galaxy, and returns something JSON-able (nothing for an odd piece, a tuple for an even one)."""
    DONE.append(piece)
    if STOP_AT.get("piece") == piece:
        raise RuntimeError("the unit is killed here")
    ctx.db.execute("INSERT INTO checkpoint_probe (piece) VALUES (?)", (piece,))
    return None if piece % 2 else {"piece": piece, "pair": (piece, piece * 2)}


def counted_probe(args, piece, ctx):
    """Writes a row for the piece and commits it as the piece (`ctx.commit_piece`), then is killed there when asked: the death between the commit and the engine's own checkpoint."""
    DONE.append(piece)
    ctx.db.execute("INSERT INTO checkpoint_probe (piece) VALUES (?)", (piece,))
    result = {"piece": piece}
    ctx.commit_piece(result)
    if STOP_AT.get("after") == piece:
        raise RuntimeError("the unit is killed after the piece's commit")
    return result


def together(args, results, ctx):
    return {"results": results, "count": len(results)}


@pytest.fixture(autouse=True)
def clear():
    DONE.clear()
    STOP_AT.clear()


@pytest.fixture
def command(monkeypatch):
    def make(resumes=True, kind="resume-test"):
        spec = commands.CommandSpec(kind, "A test", expand=(__name__, "listed"), pieces=(__name__, "listed"), do=(__name__, "probe"), gather=(__name__, "together"), resumes=resumes)
        monkeypatch.setitem(commands.REGISTRY, kind, spec)
        return spec
    return make


@pytest.fixture
def job(mysql_config):
    """A job row in the ledger and the galaxy's probe table; yields `(control, galaxy, row, args, attempt)` where `attempt()` takes the next attempt."""
    control, galaxy = store.get_control_connection(store.control_mysql_config(mysql_config)), store.get_connection(mysql_config)
    galaxy.execute("CREATE TABLE IF NOT EXISTS checkpoint_probe (piece INT NOT NULL)")
    galaxy.execute("DELETE FROM checkpoint_probe")
    galaxy.commit()
    run = ledger.add_run(control, "generate", "g")
    ledger.add_jobs(control, run, [{"stage": "s", "command": "resume-test", "item_key": "0"}])
    control.commit()
    row = control.execute("SELECT id FROM engine_jobs WHERE run_id = ?", (run,)).fetchone()["id"]
    args = {"count": PIECES, "database": wrapper.database_args(mysql_config)}

    def take():
        return attempts.Attempt(control, row, ledger.begin_attempt(control, row, galaxy), mysql_config)

    try:
        yield control, galaxy, row, args, take
    finally:
        galaxy.execute("DROP TABLE IF EXISTS checkpoint_probe")
        galaxy.commit()
        control.close()
        galaxy.close()


def _rows(galaxy):
    galaxy.commit()
    return [row["piece"] for row in galaxy.execute("SELECT piece FROM checkpoint_probe ORDER BY piece").fetchall()]


def _undisturbed(command, job):
    command()
    control, galaxy, row, args, take = job
    result = engine.run("resume-test", args, lambda done, total: None, attempt=take())
    assert DONE == list(range(PIECES))
    return result


@pytest.mark.parametrize("killed", range(PIECES))
def test_a_unit_killed_at_any_piece_starts_after_the_pieces_it_finished_and_gives_the_undisturbed_result(command, job, killed):
    control, galaxy, row, args, take = job
    undisturbed = _undisturbed(command, job)
    galaxy.execute("DELETE FROM checkpoint_probe")
    galaxy.execute("UPDATE engine_items SET cursor_key = NULL WHERE job_id = ?", (row,))
    galaxy.execute("DELETE FROM engine_piece_results WHERE job_id = ?", (row,))
    galaxy.commit()
    DONE.clear()
    STOP_AT["piece"] = killed
    with pytest.raises(RuntimeError):
        engine.run("resume-test", args, lambda done, total: None, attempt=take())
    assert DONE == list(range(killed + 1))
    STOP_AT.clear()
    DONE.clear()
    reports = []
    again = engine.run("resume-test", args, lambda done, total: reports.append((done, total)), attempt=take())
    assert DONE == list(range(killed, PIECES)), "the second attempt did only the pieces the first did not finish"
    assert again == undisturbed, "the results of the pieces done before the kill are the ones the unit gives, including a piece that returned nothing"
    assert _rows(galaxy) == list(range(PIECES)), "no piece wrote its data twice and none was lost"
    assert reports[0] == (killed, PIECES), "the bar starts where the unit stands"


def test_the_cursor_and_the_results_are_stored_with_the_data_and_a_piece_that_returned_nothing_has_no_row(command, job):
    command()
    control, galaxy, row, args, take = job
    STOP_AT["piece"] = 4
    with pytest.raises(RuntimeError):
        engine.run("resume-test", args, lambda done, total: None, attempt=take())
    galaxy.commit()
    assert ledger.item_state(galaxy, row)["cursor_key"] == "4"
    stored = galaxy.execute("SELECT piece_index FROM engine_piece_results WHERE job_id = ? ORDER BY piece_index", (row,)).fetchall()
    assert [r["piece_index"] for r in stored] == [0, 2], "pieces 1 and 3 returned nothing and have no row"
    assert ledger.resume_point(galaxy, row, 1) == (4, [{"pair": [0, 0], "piece": 0}, None, {"pair": [2, 4], "piece": 2}, None])


def test_a_piece_that_commits_itself_as_the_piece_is_not_done_again_by_a_unit_killed_just_after(command, job, monkeypatch):
    """OPS.222: `ctx.commit_piece()` lands the piece's data, the cursor past it and its result in one transaction, so a unit killed before the engine's checkpoint starts at the next piece, not the same one."""
    monkeypatch.setitem(commands.REGISTRY, "resume-test", commands.CommandSpec("resume-test", "A test", expand=(__name__, "listed"), pieces=(__name__, "listed"), do=(__name__, "counted_probe"),
                                                                            gather=(__name__, "together"), resumes=True))
    control, galaxy, row, args, take = job
    STOP_AT["after"] = 2
    with pytest.raises(RuntimeError):
        engine.run("resume-test", args, lambda done, total: None, attempt=take())
    galaxy.commit()
    assert ledger.item_state(galaxy, row)["cursor_key"] == "3", "the cursor is past the piece that committed itself"
    STOP_AT.clear()
    DONE.clear()
    result = engine.run("resume-test", args, lambda done, total: None, attempt=take())
    assert DONE == [3, 4, 5], "the piece that committed itself is not done again"
    assert _rows(galaxy) == list(range(PIECES)), "its row is there once"
    assert result["results"] == [{"piece": number} for number in range(PIECES)], "its result is the one it stored"


def test_a_piece_that_commits_itself_commits_plainly_with_no_attempt(command, job, monkeypatch):
    monkeypatch.setitem(commands.REGISTRY, "resume-test", commands.CommandSpec("resume-test", "A test", expand=(__name__, "listed"), pieces=(__name__, "listed"), do=(__name__, "counted_probe"),
                                                                            gather=(__name__, "together"), resumes=True))
    control, galaxy, row, args, take = job
    result = engine.run("resume-test", args, lambda done, total: None)
    assert result["count"] == PIECES and _rows(galaxy) == list(range(PIECES))


def test_a_piece_result_goes_through_json_so_a_resumed_unit_and_an_undisturbed_one_hand_gather_the_same_values(command, job):
    command()
    control, galaxy, row, args, take = job
    result = engine.run("resume-test", args, lambda done, total: None, attempt=take())
    assert result["results"][0] == {"pair": [0, 0], "piece": 0}, "a tuple the piece returned is a list for gather, resumed or not"


def test_a_command_that_does_not_resume_starts_at_its_first_piece_and_stores_nothing(command, job):
    command(resumes=False)
    control, galaxy, row, args, take = job
    STOP_AT["piece"] = 3
    with pytest.raises(RuntimeError):
        engine.run("resume-test", args, lambda done, total: None, attempt=take())
    galaxy.commit()
    assert ledger.item_state(galaxy, row)["cursor_key"] is None
    assert galaxy.execute("SELECT COUNT(*) AS n FROM engine_piece_results WHERE job_id = ?", (row,)).fetchone()["n"] == 0
    STOP_AT.clear()
    DONE.clear()
    engine.run("resume-test", args, lambda done, total: None, attempt=take())
    assert DONE == list(range(PIECES)), "every piece again"


def test_a_unit_with_no_attempt_runs_from_its_first_piece_and_commits_plainly(command, job):
    command()
    control, galaxy, row, args, take = job
    result = engine.run("resume-test", args, lambda done, total: None)
    assert result["count"] == PIECES and DONE == list(range(PIECES))


def test_a_replaced_attempt_stops_at_its_next_checkpoint_and_stores_nothing(command, job):
    command()
    control, galaxy, row, args, take = job
    first = take()
    take()                                  # a retry started while the first was alive
    with pytest.raises(attempts.Replaced):
        engine.run("resume-test", args, lambda done, total: None, attempt=first)
    galaxy.commit()
    assert DONE == [0], "the first piece ran and its commit found the unit taken over"
    assert _rows(galaxy) == [], "its data rolled back with the fence"
    assert galaxy.execute("SELECT COUNT(*) AS n FROM engine_piece_results WHERE job_id = ?", (row,)).fetchone()["n"] == 0


def test_a_unit_that_finished_and_is_run_again_does_nothing_and_gives_its_result(command, job):
    command()
    control, galaxy, row, args, take = job
    first = engine.run("resume-test", args, lambda done, total: None, attempt=take())
    DONE.clear()
    again = engine.run("resume-test", args, lambda done, total: None, attempt=take())
    assert DONE == [] and again == first and _rows(galaxy) == list(range(PIECES))


def test_a_new_attempt_keeps_the_cursor_a_replaced_one_may_not_read(command, job):
    command()
    control, galaxy, row, args, take = job
    STOP_AT["piece"] = 2
    with pytest.raises(RuntimeError):
        engine.run("resume-test", args, lambda done, total: None, attempt=take())
    galaxy.commit()
    second = take()
    galaxy.commit()
    assert ledger.item_state(galaxy, row)["cursor_key"] == "2", "the next attempt's number is written, the cursor stays"
    assert second.resume_point()[0] == 2
    assert ledger.resume_point(galaxy, row, second.number - 1) == (0, []), "a number that is not the unit's current attempt reads nothing"


# --- the commands that say they resume: each one proved, none left to a guess --------------------------------

RESUMES = {
    "scatter-mass": "pieces: one for each shell of the unit's run, then the far shells and the components, from the unit's own arguments; every piece draws on its own random stream and writes its rows",
    "scatter-luminosity": "the same pieces as scatter-mass",
    "scatter-band": "the same pieces as scatter-mass",
    "scatter-phenomena": "the same shell, far and component pieces from the unit's arguments, each piece writing its own rows",
    "scatter-backfill": "a sky unit is one piece; a unit of sectors is the sorted items in chunks of the unit's own arguments; what the pieces keep in `ctx.kept` is read again or is a cache",
    "refresh-run-blocks": "two fixed pieces, the nebulae and then the sectors saved since the run began, each a catch-up that does nothing a second time",
    "fill-map-luminosity": "ranges of ids from the unit's low and high, each piece touching only rows still empty",
    "draw-nebula-field": "one slab of the field for each index from the unit's first, each written to `nebula_field`",
    "link-neighbors": "the unit's chunk of sectors in `store.link_order`, in batches, each batch one transaction under the neighbour lock",
    "settle-sector-paths": "the unit's sectors that still exist, in batches that sit together",
    "finish-phenomena": "four fixed pieces; the drawn special rows are the first piece's stored result, the write and the record commit through `ctx.commit_piece()` so their data and the cursor land together, the indexes are built only where missing",
}
"""dict: The commands that `resumes`, each with why its pieces may be started after a cursor (the three conditions of `commands.CommandSpec`)."""

RESTARTS = {
    "plan-skeleton": "the pieces hand each other the skeleton through ctx.state",
    "refresh-block-lights": "one piece",
    "store-scatter-limits": "one piece",
    "rescore-habitability": "the list is the systems whose stamp is not current: a piece done changes the list, so a cursor into it skips the wrong ones",
    "refill-effect-ledger": "the list is the systems whose stamp is not current, as rescore-habitability",
    "generate-start-sector": "a unit run again (`again`) lists only the sectors not yet saved, and its result counts the process that began it (`pid`)",
    "generate-neighborhood": "as generate-start-sector",
    "generate-sectors": "as generate-start-sector",
    "generate-unplaced-sector": "its result names the process that ran it (`pid`), which the starter counts its sectors by",
    "run-population": "the first piece forgets every species for a rescan and the list is read from what the scan has done",
    "run-territories": "the pieces are read from the polities the earlier pieces changed",
    "generate-system": "two pieces that hand each other the system through ctx.state",
    **dict.fromkeys(("clear-bright-stars", "empty-phenomenon-scatter", "empty-nebula-field"), "the look piece decides whether the pieces after it do anything, and a clear stopped before its last statement is done again in full (the TRUNCATE comes last)"),
    "drop-phenomenon-indexes": "each piece drops an index only where the table still has it, and the whole is a few statements: run again from the first costs nothing",
    "clear-scatter": "the first piece chooses the clear and hands it to the next through ctx.state",
    "clean-nebula-field": "the pieces hand each other the clouds read and checked through ctx.state",
    "reset-galaxy": "the facts of the tables and the foreign-key setting are read once for the unit, and a wipe stopped half way is run again whole",
}
"""dict: The commands with pieces that start from their first piece when run again, each with why (the warm-ups and checks, which are short or read-only and restart cheaply, are listed by prefix below)."""

RESTART_PREFIXES = ("check-", "warm-", "end-map-rebuild")
"""tuple: Kinds that restart because a run again costs next to nothing (read-only checks) or only rebuilds a cache (warm-ups)."""


def test_every_command_with_pieces_has_been_audited_for_resuming():
    pieced = {kind for kind, spec in commands.REGISTRY.items() if spec.pieces is not None}
    unaudited = sorted(kind for kind in pieced if kind not in RESUMES and kind not in RESTARTS and not kind.startswith(RESTART_PREFIXES))
    assert not unaudited, f"say whether these resume (RESUMES, with why it is safe) or start again (RESTARTS, with why not): {unaudited}"
    assert {kind for kind in pieced if commands.REGISTRY[kind].resumes} == set(RESUMES), "the commands that say `resumes` are exactly the audited list"
    assert not sorted(kind for kind in RESTARTS if kind in commands.REGISTRY and commands.REGISTRY[kind].resumes)
    assert not sorted(kind for kind in {*RESUMES, *RESTARTS} if kind not in commands.REGISTRY), "listed and not registered"


def test_a_command_that_resumes_needs_pieces():
    with pytest.raises(ValueError, match="needs `pieces`"):
        commands.CommandSpec("never", "A test", unit=(__name__, "listed"), resumes=True)


SHELL_UNIT = {"shells": [(0, 1), (1, 2), (2, 3)], "far": (3, 6), "component_parts": [1, 2]}
SECTOR_ITEMS = [[x, 0, 0, 1.0] for x in range(5 * 200 + 1)]

PIECE_ARGS = {
    "scatter-mass": SHELL_UNIT, "scatter-luminosity": SHELL_UNIT, "scatter-band": SHELL_UNIT, "scatter-phenomena": SHELL_UNIT,
    "scatter-backfill": {"items": SECTOR_ITEMS},
    "refresh-run-blocks": {"since": 0},
    "fill-map-luminosity": {"low": 1, "high": 100_000},
    "draw-nebula-field": {"i": 3, "count": 7},
}


@pytest.mark.parametrize("kind", sorted(PIECE_ARGS))
def test_the_pieces_of_a_command_after_a_cursor_are_its_pieces_less_the_first_cursor_of_them_whatever_the_cursor(kind):
    """Condition (1) of `CommandSpec`, for each command whose list comes from its arguments alone."""
    listing = commands.pieces_function(kind)
    whole = listing(PIECE_ARGS[kind], 0)
    assert len(whole) > 1
    for cursor in range(len(whole) + 2):
        assert listing(PIECE_ARGS[kind], cursor) == whole[cursor:]


def test_a_stopped_map_luminosity_fill_run_again_does_the_pieces_left_and_fills_what_an_undisturbed_unit_fills(mysql_config, monkeypatch, job):
    """A real resumable command, killed at every piece boundary: the second attempt does only the pieces left and the unit's result and the galaxy are those of one that was never stopped."""
    from planetgen.generation import run_plan
    from tests.test_phenomenon_scatter import _scatter_args, _seed_galaxy
    control, galaxy, row, args, take = job
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_scatter_args(mysql_config))
    galaxy.commit()
    top = int(galaxy.execute("SELECT MAX(id) AS n FROM phenomenon_scatter").fetchone()["n"])
    monkeypatch.setattr(run_plan, "MAP_LUMINOSITY_PIECE_IDS", max(top // 5, 1))
    params = {"low": 1, "high": top, "database": args["database"]}
    pieces = len(run_plan.fill_map_luminosity_pieces(params, 0))
    assert pieces >= 4

    def values():
        galaxy.commit()
        return [(r["id"], r["map_luminosity_sol"]) for r in galaxy.execute("SELECT id, map_luminosity_sol FROM phenomenon_scatter ORDER BY id").fetchall()]

    def empty():
        galaxy.execute("UPDATE phenomenon_scatter SET map_luminosity_sol = NULL")
        galaxy.commit()

    empty()
    undisturbed = engine.run("fill-map-luminosity", params, lambda done, total: None, attempt=take())
    expected = values()
    assert undisturbed["filled"] > 0
    real = run_plan.fill_map_luminosity_do
    for killed in range(pieces):
        empty()
        galaxy.execute("UPDATE engine_items SET cursor_key = NULL WHERE job_id = ?", (row,))
        galaxy.execute("DELETE FROM engine_piece_results WHERE job_id = ?", (row,))
        galaxy.commit()
        calls = []

        def stopping(unit_params, piece, ctx, killed=killed):
            if len(calls) == killed:
                raise RuntimeError("killed")
            calls.append(piece)
            return real(unit_params, piece, ctx)

        monkeypatch.setattr(run_plan, "fill_map_luminosity_do", stopping)
        with pytest.raises(RuntimeError):
            engine.run("fill-map-luminosity", params, lambda done, total: None, attempt=take())
        assert len(calls) == killed
        monkeypatch.setattr(run_plan, "fill_map_luminosity_do", lambda unit_params, piece, ctx: calls.append(piece) or real(unit_params, piece, ctx))
        again = engine.run("fill-map-luminosity", params, lambda done, total: None, attempt=take())
        assert len(calls) == pieces, "every piece was done once, the ones before the kill by the first attempt and the rest by the second"
        assert again == undisturbed and values() == expected
