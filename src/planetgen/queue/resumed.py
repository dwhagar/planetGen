# planetgen/queue/resumed.py

"""
Resume from rows (OPS.222 part B, L8 of `docs/design/operations/job-engine.md`): a run that carries on from a stopped one subtracts what the stopped run's rows say is done.

`generation/finish_run.py` plans the stopped run again: a stage the run finished is carried over whole (`wrapper.FINISHED_BY`), and every other stage is cut again and marked with the id of the run it
carries on from (`RESUMES_FROM`). When such a stage is cut (`wrapper._open`) `carry` reads the earlier run's rows of the same stage and, for each unit the new cut has the same as the earlier run's
(same item key, same arguments, same cost):

* **A unit the rows say is done** is marked done in the new run's rows, its stored result is copied into the new run's results (the stage's `conclude` reads them), and it joins the stage's done count.
  It is not queued: it did its work once, and the galaxy holds it.
* **A unit the rows say is not done** is queued as any unit is, and a command that `resumes` (`commands.CommandSpec`) starts it after the pieces the earlier run's attempts finished: the cursor and
  the stored piece results of the earlier unit are copied to the new unit's rows in the galaxy (`engine_items`, `engine_piece_results`), where `ledger.begin_attempt` leaves them for the attempt
  that follows.

* **A unit of a command that cannot resume** (`commands.CommandSpec.resumes` false: it starts from its first piece) is queued with `AGAIN` among its arguments, so a command whose pieces skip what an earlier attempt
  saved (the sector stages) knows it is run again. The row it was cut with is not changed (a unit row's arguments are `{args, cost}`, as the plan wrote them).

A unit whose cut differs (a stage that plans from what the galaxy holds now cuts fewer units, numbered from 0 again) matches nothing and runs whole, as it would in a run begun from nothing. Nothing
here is a second record of what is done: the earlier run's rows are, and they are only read.
"""

import json

from planetgen.queue import commands, ledger, planrows, runrecord
from planetgen.util import log

RESUMES_FROM = "resumes_from"
"""str: The step parameter that marks a stage a run carries on from a stopped one (`generation/finish_run.py`): the id of the run whose rows are subtracted when the stage is cut. A stage is handed its
parameters without it, as a stage that finished (`wrapper.FINISHED_BY`) is."""


AGAIN = "again"
"""str: The argument a unit is queued with (and its row is not cut with) when it belongs to a stage that carries on from a stopped run and its command cannot resume: the unit starts from its first piece, and
a command that lists what it has not done yet (`run_galaxy.generate_sectors_pieces`: the sectors not yet saved) says so by it. A command that `resumes` is given its cursor instead."""


class Carried:
    """
    What `carry` found for stage `index`.

    Attributes:
        done (dict[int, tuple]): Unit number -> `(stored result as text, cost)` of each unit the earlier run's rows say is done and whose result is stored.
        resumed (int): How many units that are not done start after pieces the earlier run's attempts finished.
    """

    def __init__(self):
        self.done = {}
        self.resumed = 0


def _same(row, args, cost):
    """Whether the earlier row's arguments are the ones the new cut gives the unit (compared as the JSON both are stored as, without the `AGAIN` mark)."""
    try:
        earlier, new = json.loads(row["arguments"]), json.loads(json.dumps({"args": args, "cost": cost}, sort_keys=True))
    except (TypeError, ValueError):
        return False
    for body in (earlier, new):
        if isinstance(body.get("args"), dict):
            body["args"].pop(AGAIN, None)       # what a run that carries on tells a unit that cannot resume (`planrows.add_to_unit_args`) is no part of what the unit is: a unit of a run that carried on from another is the same unit
    return earlier == new


def _earlier_rows(control, run, stage, keys=None):
    """`{item key: row}` of the units of `stage` in the earlier run's ledger run; only the units whose item keys are in `keys` when it is given (a window of a stage written in windows, OPS.240)."""
    where, values = "", ()
    if keys is not None:
        keys = [str(key) for key in keys]
        where, values = " AND item_key IN (" + ", ".join(["?"] * len(keys)) + ")", tuple(keys)
    return {row["item_key"]: row for row in control.execute(
        "SELECT id, item_key, state, arguments FROM engine_jobs WHERE run_id = ? AND stage = ? AND item_key <> ''" + where, (run, stage, *values)).fetchall()}


def carry(conn, run_id, index, prior, kind, rows, units, first=0):
    """
    Subtracts the earlier run `prior`'s done units from the cut of stage `index` of `run_id` and gives the rest the earlier run's checkpoints (see the module).

    A unit is done when its row says so, or when its done mark in the galaxy is set (`engine_items.done_at`, OPS.217, true before the result is written), or when it is in the earlier run's done set (the
    barrier Redis still holds until OPS.218's follow-ups remove it): a unit that finished its work and was then failed by the end of its stage (the conclusion that raised) has the second and third and not the
    first. Whichever says so, its stored result must be there.

    Args:
        rows (list[int]): The new run's unit rows, in unit order (`planrows.record_cut`).
        units (list[tuple]): `(args, cost)` per unit, as the cut gave them.
        first (int): The number of the first of `units`, for a window of a stage written in windows (OPS.240): `rows` and `units` are that window's and the earlier run's rows are read for its numbers only.

    Returns:
        Carried: The units done (their rows are marked done here), and how many others resume. Empty when the earlier run has no rows for the stage.
    """
    from planetgen.db import store
    found = Carried()
    stage = planrows._stage_key(conn, run_id, index)
    stages = runrecord.declared_stages(conn, prior)
    engine = conn.hget(runrecord.run_key(prior), planrows.ENGINE_RUN_FIELD)
    if not engine or index >= len(stages) or runrecord._s(stages[index]["key"]) != stage:
        return found
    control = planrows.control_connection(conn, run_id)
    galaxies = {}
    try:
        earlier = _earlier_rows(control, int(runrecord._s(engine)), stage, range(first, first + len(units)))
        resumes = commands.spec(kind).resumes
        for place, (args, cost) in enumerate(units):
            unit = first + place
            row = earlier.get(str(unit))
            if row is None or not _same(row, args, cost):
                continue
            named = (args or {}).get("database")
            item = _item(galaxies, named, row["id"], store) if named else None
            raw = conn.hget(runrecord.results_key(prior, index), unit)
            if raw is not None and (row["state"] == "done" or (item is not None and item["done_at"] is not None) or conn.sismember(runrecord.done_key(prior, index), unit)):
                found.done[unit] = (runrecord._s(raw), cost)
                ledger.set_job_state(control, rows[place], "done")
            elif resumes and item is not None and _cursor(item) > 0:
                _copy_checkpoint(galaxies[json.dumps(named, sort_keys=True)], row["id"], rows[place], _cursor(item))
                found.resumed += 1
        control.commit()
    finally:
        control.close()
        for galaxy in galaxies.values():
            try:
                galaxy.close()
            except Exception as exc:  # noqa: BLE001 -- a connection that will not close is not the run's failure
                log.swallowed("queue/resumed", "close", exc)
    return found


def _item(galaxies, named, job_id, store):
    """The earlier unit's `engine_items` row (`cursor_key`, `done_at`) in the galaxy `named`, reading it on a connection kept for the stage; `None` when it has none."""
    key = json.dumps(named, sort_keys=True)
    if key not in galaxies:
        galaxies[key] = store.get_connection(store.MySQLConfig.from_args(named))
    return galaxies[key].execute("SELECT cursor_key, done_at FROM engine_items WHERE job_id = ?", (job_id,)).fetchone()


def _cursor(item):
    """The number of pieces the earlier unit finished (0 when it stored none)."""
    text = str(item["cursor_key"] or "")
    return int(text) if text.isdigit() else 0


def _copy_checkpoint(galaxy, earlier, new, cursor):
    """
    Copies the cursor and the piece results of the earlier unit `earlier` to the new unit's rows `new` in the galaxy. Attempt 0 is what the new unit's row says until `ledger.begin_attempt` raises
    it, which leaves the cursor and the results as they are.
    """
    galaxy.execute("INSERT INTO engine_items (job_id, attempt, cursor_key, updated_at) VALUES (?, 0, ?, NOW(6))"
                   " ON DUPLICATE KEY UPDATE cursor_key = ?, updated_at = NOW(6)", (new, str(cursor), str(cursor)))
    galaxy.execute("INSERT INTO engine_piece_results (job_id, piece_index, result) SELECT job_id, piece_index, copied FROM"
                   " (SELECT ? AS job_id, piece_index, result AS copied FROM engine_piece_results WHERE job_id = ? AND piece_index < ?) AS earlier"
                   " ON DUPLICATE KEY UPDATE result = copied", (new, earlier, cursor))
    galaxy.commit()


def settle(conn, run_id, index, found):
    """
    Counts the units `carry` found done in the stage's totals and results (`runrecord.unit_done`, the stored result first) so the stage's bar, its ETA and its conclusion read them as they would
    a unit that just finished. A stage with no unit left to run is ended by the opener (`wrapper._open`) as the unit that completes a stage ends it.
    """
    for unit, (raw, cost) in sorted(found.done.items()):
        conn.hset(runrecord.results_key(run_id, index), unit, raw)
        runrecord.unit_done(conn, run_id, index, unit, cost=cost, seconds=0.0)
