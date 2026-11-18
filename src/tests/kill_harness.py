# tests/kill_harness.py

"""
Kills in the middle of a run's units, for the proof of the job engine's stage 3 (OPS.225): the units of a run are done by the inline launcher on a real MySQL, and a seeded plan kills some of them at a
commit, as a SIGKILL would. A kill before the commit loses what the piece wrote since the last one; a kill after it leaves the data, the cursor and the results committed and the unit's end unwritten. Either way
the launcher finds the job process gone, closes its attempt `lost` and hands the unit out again, and the galaxy the run ends with must be the galaxy of an undisturbed run.

A unit is killed at most once (the launcher retries a lost unit once, `launcher.RETRIES`), so a plan names global commit numbers and the first commit of each unit that falls on one of them is the kill.
"""

import hashlib
import os
import signal

from planetgen.queue import attempts, engine, forkedjob, jobhost
from tests import inline_queue


class Killed(BaseException):
    """What a kill looks like to the code it interrupts: nothing in the engine may catch it by accident."""


class Plan:
    """Which commits of a run are killed: `kills` maps the number of a commit (counted from 0, over every unit of the run) to `"before"` or `"after"`."""

    def __init__(self, kills=None):
        self.kills = dict(kills or {})
        self.seen = 0
        self.units = set()          # the units already killed once
        self.fired = []             # (commit number, mode, unit) of every kill that happened
        self._depth = 0
        self.dead = False           # a killed job process writes nothing more: its handlers run in this process, and are made to do nothing

    def hook(self, real):
        plan = self

        def hooked(ctx, *args, **kwargs):
            if plan._depth:                 # `commit_piece` calls `commit`: one commit is one number
                return real(ctx, *args, **kwargs)
            number, plan.seen = plan.seen, plan.seen + 1
            mode = plan.kills.get(number)
            unit = (ctx.run_id, ctx.stage, ctx.unit)
            if mode is None or unit in plan.units or ctx.attempt is None:
                plan._depth += 1
                try:
                    return real(ctx, *args, **kwargs)
                finally:
                    plan._depth -= 1
            plan.units.add(unit)
            plan.fired.append((number, mode, unit))
            if mode == "after":
                plan._depth += 1
                try:
                    real(ctx, *args, **kwargs)
                finally:
                    plan._depth -= 1
            plan.dead = True
            raise Killed(f"commit {number} {mode}")
        return hooked


class KillingForker(inline_queue.InlineForker):
    """The inline forker, where a unit that is killed ends as a job process that died of SIGKILL: no end written, no message."""

    plan = None

    def fork(self, kind, runner, job=None, template=None):
        self.made += 1
        link = inline_queue._Link()
        paused = signal.alarm(0)
        try:
            try:
                exit_ = jobhost.Exit(self.made, os.getpid(), code=forkedjob.run(job, link))
            except Killed:
                self.plan.dead = False
                exit_ = jobhost.Exit(self.made, os.getpid(), signal=signal.SIGKILL)
        finally:
            if paused:
                signal.alarm(paused)
        self.exits.append(("inline", exit_))
        if os.environ.get("KILL_DEBUG"):
            from planetgen.db import store
            control = store.get_control_connection()
            print("ROWS", [(r["id"], r["command"], r["item_key"], r["state"], r["attempt"]) for r in control.execute("SELECT id, command, item_key, state, attempt FROM engine_jobs WHERE run_id = (SELECT MAX(id) FROM engine_runs) ORDER BY id").fetchall()], flush=True)
            control.close()
        return "inline", inline_queue._InlineJob(self.made, link)


def _forker(plan):
    forker = KillingForker()
    forker.plan = plan
    return forker


def install(monkeypatch, plan):
    """Puts `plan` on every commit a unit makes (`engine.Context.commit`, `commit_piece`, `commit_piece_in`) and the killing forker under the inline launcher."""
    for name in ("commit", "commit_piece", "commit_piece_in"):
        monkeypatch.setattr(engine.Context, name, plan.hook(getattr(engine.Context, name)))
    real_end = attempts.Attempt.end

    def end(attempt, *args, **kwargs):
        return None if plan.dead else real_end(attempt, *args, **kwargs)

    monkeypatch.setattr(attempts.Attempt, "end", end)
    monkeypatch.setattr(inline_queue, "FORKER", lambda: _forker(plan))


def digest_of_galaxy(config):
    """What the galaxy is, without the numbers the database gave it: the charted sectors, the systems, stars and planets by their own columns, in an order of their own, and how many paths join the sectors (a rogue planet in a cell nobody charted gets a sector row of its own, GEN.1111, which is not a sector the run made)."""
    from planetgen.db import store
    db = store.get_connection(config)
    try:
        queries = {
            "sectors": "SELECT s.shell_index, s.band_index, s.wedge_index, s.name, s.center_x_pc, s.center_y_pc, s.center_z_pc FROM sectors s JOIN chart_ledger c ON c.shell_index = s.shell_index "
                       "AND c.band_index = s.band_index AND c.wedge_index = s.wedge_index ORDER BY s.shell_index, s.band_index, s.wedge_index",
            "systems": "SELECT cell_key, name, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc, is_binary FROM star_systems ORDER BY cell_key, name",
            "stars": "SELECT y.cell_key, y.name AS system_name, s.role, s.mass_kg, s.temperature_k, s.luminosity_w FROM stars s JOIN star_systems y ON y.id = s.star_system_id ORDER BY y.cell_key, y.name, s.role, s.mass_kg",
            "planets": "SELECT y.cell_key, y.name AS system_name, p.orbital_index, p.body_type, p.planet_class, p.mass_kg, p.radius_km, p.distance_km, p.atmospheric_pressure_pa "
                       "FROM planets p JOIN star_systems y ON y.id = p.star_system_id ORDER BY y.cell_key, y.name, p.orbital_index, p.body_type, p.mass_kg",
        }
        found = {}
        for name, query in queries.items():
            sha, count = hashlib.sha256(), 0
            for row in db.execute(query).fetchall():
                sha.update(repr(tuple(row.values())).encode())
                count += 1
            found[name] = (count, sha.hexdigest())
        found["paths"] = db.execute("SELECT COUNT(*) AS n FROM sector_paths").fetchone()["n"]
        return found
    finally:
        db.close()
