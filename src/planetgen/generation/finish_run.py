# planetgen/generation/finish_run.py

"""
Finish run (PERF.182)
=====================

`planetgen finish-run`, the Generate page's "Finish the last run" and the typed command "finish run" carry on a generation
run that stopped partway (a killed process, an out-of-memory MySQL, a cancel), from the first stage that did not finish,
instead of wiping the galaxy and starting again.

**The record is the run's own.** A run is data (`queue/runrecord.py`): its steps, each a command and its complete
parameters, and for each stage the clock it began, the clock it ended and the conclusion its command stored. A run also
records the galaxy's seed and the build (the version key) it was started on. Nothing else is needed to carry it on: the
steps are queued again from the record, so the arguments are not worked out again and a random start cannot move.

**What it does.** `plan` reads the record of the newest generation run that stopped (or the one named), refuses when the
galaxy's seed or the build is not the one the run was started on or a stage after the first unfinished one has no way to
carry on (`FINISHABLE`), and gives the steps of a new run (`Plan.steps`):

* a stage the run finished is *carried over*: it opens and ends at once with no unit, and the conclusion the earlier run
  stored (`wrapper.FINISHED_BY`);
* the stage that stopped runs again from its start, its `prepare` included, which empties what the stage rewrites
  (the scatters, the nebula field), and a sector stage skips the sectors already saved (`again`);
* a stage that cannot carry on from the middle of the one before it runs that one again too (`REWIND_WITH`): the
  phenomena indexes and the special phenomena are written into the table the scatter's shells fill, and only the scatter's
  `prepare` empties it;
* every stage after it runs as it would have.

The new run records the run it carries on from (`finishes`), so a finish that is itself stopped can be finished in turn.
"""

import contextlib
import dataclasses
import sys

from planetgen.db import store
from planetgen.queue import commands, redisqueue, runrecord, wrapper

FINISHABLE = frozenset({
    "generate-start-sector", "generate-neighborhood", "generate-sectors", "link-neighbours", "scatter-mass",
    "scatter-luminosity", "scatter-phenomena", "finish-phenomena", "draw-nebula-field", "refresh-block-lights",
    "scatter-backfill", "run-population", "run-territories", "settle-sector-paths", "refresh-run-blocks"})
"""frozenset[str]: The commands whose stage can be run again from its start after it stopped partway, and give the
galaxy an uninterrupted run would have: each either empties what it rewrites when it opens (`prepare`) or does work that
saves nothing twice (a sector already saved is skipped, a total made again is the same total)."""

REWIND_WITH = {"finish-phenomena": "scatter-phenomena", "scatter-luminosity": "scatter-mass"}
"""dict: The command whose stage stopped -> the command of the stage before it that runs again with it. The special
phenomena are added to the rows the shells wrote, and the bright stars of a run that draws both passes are split from the
massive ones at a limit the run stores when it ends, so neither can be emptied alone."""

AGAIN = "again"
"""str: The step parameter that marks the stage that stopped: it is run again from its start (`generate_sectors_unit`
leaves the sectors already saved alone)."""


class CannotFinish(RuntimeError):
    """There is no run to carry on, or the one found cannot be carried on; the message says why, in words for the page."""


@dataclasses.dataclass
class Plan:
    """
    A run to carry on.

    Attributes:
        run (dict): The stopped run (`runrecord.read_run`).
        steps (list[commands.Step]): The steps of the new run, every stage of the stopped one: those it finished carry
            `wrapper.FINISHED_BY`, the one that stopped `AGAIN`.
        first (int): The index of the first stage that runs again (stages before it are carried over).
        stopped (int): The index of the first stage that did not finish.
    """

    run: dict
    steps: list
    first: int
    stopped: int

    @property
    def label(self):
        """The first stage that runs again, as the pages say it."""
        return self.steps[self.first].label


def stage_finished(conn, run_id, index, step):
    """
    Whether stage `index` of the run finished: every unit done (the counter marks the stage `ended`) and, for a command
    with a `conclude`, the conclusion stored. A stage cut to no unit never gets an `ended` mark; the run having
    begun the stage after it says it finished.
    """
    fields = {runrecord._s(name): runrecord._s(value)
              for name, value in conn.hgetall(runrecord.stage_key(run_id, index)).items()}
    if "started" not in fields:
        return False
    if conn.hget(runrecord.stage_key(run_id, index + 1), "started") is not None:
        return True
    if "ended" not in fields:
        return False
    return commands.conclude_function(step.kind) is None or "concluded" in fields


def first_unfinished(conn, run_id, steps):
    """The index of the first stage of the run that did not finish, `len(steps)` when every one did."""
    for index, step in enumerate(steps):
        if not stage_finished(conn, run_id, index, step):
            return index
    return len(steps)


def last_incomplete(conn):
    """
    The id of the newest generation run, when it stopped before it was done (failed, cancelled, or live with nobody
    following it), `None` when the newest generation run was done or there is none.

    Raises:
        CannotFinish: The newest generation run is still running.
    """
    newest = None
    for run_id in runrecord.recent_runs(conn, limit=200):      # ids sort by the second a run began in: two can tie
        fields = conn.hmget(runrecord.run_key(run_id), "state", "kind", "created")
        state, kind = runrecord._s(fields[0]), runrecord._s(fields[1]) or runrecord.GENERATION
        if kind == runrecord.GENERATION and state:
            created = float(runrecord._s(fields[2]) or 0)
            if newest is None or created > newest[0]:
                newest = (created, run_id, state)
    if newest is not None:
        _created, run_id, state = newest
        if state == "done":
            return None
        if state in runrecord.LIVE_STATES:
            if not runrecord.holder_is_dead(conn, run_id):
                raise CannotFinish(f"run {run_id} is still running.")
            runrecord.finish_run(conn, run_id, "failed",
                                 error="its process was gone (killed or restarted) before the run ended")
            conn.register_script(runrecord._RELEASE)(keys=[runrecord.active_key(runrecord.GENERATION)], args=[run_id])
        return run_id
    return None


def plan(conn, run_id, galaxy_seed, build, database):
    """
    The run `run_id` as a run to carry on.

    Args:
        galaxy_seed (str): The seed of the galaxy in the database now, as hex, and `build` the running code's version key:
            the run is carried on only on the galaxy and the build it was started on. `build=None` leaves the build
            unchecked (the Generate page asks before the command, which checks it).
        database (str): The database this process works on: the run's.

    Raises:
        CannotFinish: There is no such run, it is not a stopped generation run, it was started on another galaxy or
            build or database, or a stage it would run has no way to carry on.
    """
    run = runrecord.read_run(conn, run_id)
    if run is None:
        raise CannotFinish(f"There is no record of run {run_id}.")
    if run["kind"] != runrecord.GENERATION:
        raise CannotFinish(f"Run {run_id} ({run['title']}) does not change the galaxy.")
    if run["state"] in runrecord.LIVE_STATES:
        if not runrecord.holder_is_dead(conn, run_id):
            raise CannotFinish(f"Run {run_id} ({run['title']}) is still running.")
        runrecord.finish_run(conn, run_id, "failed", error="its process was gone (killed or restarted) before the run ended")
        run = runrecord.read_run(conn, run_id)
    if run["state"] == "done":
        raise CannotFinish(f"Run {run_id} ({run['title']}) was done; there is nothing to finish.")
    carried_on = next((other for other in runrecord.recent_runs(conn, limit=200)
                       if runrecord._s(conn.hget(runrecord.run_key(other), "finishes")) == run_id), None)
    if carried_on:
        raise CannotFinish(f"Run {run_id} ({run['title']}) was already carried on by run {carried_on}.")
    if not run["galaxy_seed"] or not run["build"]:
        raise CannotFinish(f"Run {run_id} ({run['title']}) did not record the galaxy and the build it worked on, so it "
                           f"cannot be carried on; run it again.")
    if run["database"] != database:
        raise CannotFinish(f"Run {run_id} worked on the database {run['database']}, not {database}.")
    if run["galaxy_seed"] != galaxy_seed:
        raise CannotFinish(f"The galaxy was planned again since run {run_id} started (its seed is not the one the run "
                           f"worked on), so the run cannot be carried on.")
    if build is not None and run["build"] != build:
        raise CannotFinish(f"PlanetGen changed since run {run_id} started (version key {run['build']} then, {build} now): a "
                           f"galaxy made by two builds is not the one either makes, so the run cannot be carried on. "
                           f"Run it again.")
    steps = [commands.Step.from_dict(item) for item in runrecord.read_steps(conn, run_id)]
    queued = wrapper.queued(steps)
    if not queued or len(queued) != len(run["stages"]):
        raise CannotFinish(f"Run {run_id} ({run['title']}) holds no steps to carry on.")
    stopped = first_unfinished(conn, run_id, queued)
    if stopped >= len(queued):
        raise CannotFinish(f"Every stage of run {run_id} ({run['title']}) finished; there is nothing to carry on.")
    first = stopped
    while True:
        earlier = REWIND_WITH.get(queued[first].kind)
        if earlier is None or (queued[first].kind == "scatter-luminosity" and queued[first].params.get("clear")):
            break
        found = next((index for index in range(first - 1, -1, -1) if queued[index].kind == earlier), None)
        if found is None:
            break
        first = found
    unable = [step for step in queued[first:] if step.kind not in FINISHABLE]
    if unable:
        raise CannotFinish(f"The stage \"{unable[0].label}\" of run {run_id} cannot be carried on from the middle: run it "
                           f"again.")
    carried = []
    for index, step in enumerate(queued):
        # a stage that opened in the run being carried on keeps the clock it recorded (PERF.213); one that is run again opens anew
        params = {name: value for name, value in wrapper.resolve_params(conn, run_id, step.params, strict=False).items()
                  if name not in (wrapper.FINISHED_BY, AGAIN)}
        if index < first:
            params[wrapper.FINISHED_BY] = run_id
        elif index == stopped:
            params[AGAIN] = True
        carried.append(commands.Step(step.key, step.kind, step.label, params))
    return Plan(run, carried, first, stopped)


def run_description(run):
    """A run as one line for the pages: its title and the stage it stopped in."""
    stage = run["stages"][min(max(run["current"], 0), len(run["stages"]) - 1)]["label"] if run["stages"] else ""
    return f"{run['title']}, stopped in \"{stage}\""


def find(conn, galaxy_seed, build, database, run_id=None):
    """
    The `Plan` to carry on `run_id`, else the last run that did not finish (`last_incomplete`).

    Raises:
        CannotFinish: See `plan`; or there is no run that stopped.
    """
    if run_id is None:
        run_id = last_incomplete(conn)
        if run_id is None:
            raise CannotFinish("No run stopped before it was done; there is nothing to finish.")
    return plan(conn, run_id, galaxy_seed, build, database)


def add_arguments(parser):
    """Adds the options of `planetgen finish-run`: `--run`, then the common ones (`cli.generate.add_fill_arguments`)."""
    parser.add_argument("--run", metavar="ID", default=None,
                        help="The run to carry on (the Queue page lists the ids). Default: the last run that stopped "
                             "before it was done.")
    from planetgen.cli import generate
    generate.add_fill_arguments(parser)


def current_galaxy(args):
    """`(galaxy seed as hex, build)` of the database `args` name and the running code, `(None, build)` when it is not planned."""
    from planetgen.galaxy import version_key
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        seed = store.get_galaxy_seed(conn)
    finally:
        conn.close()
    return (None if seed is None else seed.hex(), version_key.version_key())


def run_finish(args):
    """
    `planetgen finish-run`: finds the run (`find`), says what it carries over and what it runs again, and follows the new
    run to its end.

    Raises:
        SystemExit: There is nothing to finish, or the run cannot be carried on (the reason is logged).
    """
    from planetgen.generation import run_common, stages
    from planetgen.util import log
    mysql_config = store.mysql_config_from_args(args)
    seed, build = current_galaxy(args)
    conn = redisqueue.connect()
    try:
        found = find(conn, seed, build, mysql_config.database, args.run)
    except CannotFinish as exc:
        log.error(str(exc))
        raise SystemExit(1) from exc
    steps, old = found.steps, found.run
    log.normal(f"Carrying on run {old['id']} ({run_description(old)}): {found.first} stage(s) finished before, "
               f"\"{found.label}\" and the {len(steps) - found.first - 1} after it run now.")
    stages.begin([stages.Stage(step.key, step.label) for step in steps], args, "finish-run")
    followers = [RecordedStage(args, step) for step in steps]
    workers = run_common._worker_count(args)

    def on_plan(index, plan_):
        followers[index].on_plan(plan_)

    def on_result(index, unit, result):
        followers[index].on_result(unit, result)

    def on_stage_done(index, conclusion):
        followers[index].end(conclusion)

    with run_common._generation_progress() as progress:
        log.set_console(progress.console)
        for follower in followers:
            follower.progress, follower.workers = progress, workers
        try:
            wrapper.run(f"Finish: {old['title']}", steps, on_plan, on_result, workers_count=workers,
                        who=getattr(args, "who", None), database=mysql_config.database, argv=sys.argv[1:],
                        listing=steps, on_stage_done=on_stage_done, galaxy_seed=seed, build=build, finishes=old["id"])
        except BaseException as exc:
            for follower in followers:
                follower.abort(type(exc), exc, exc.__traceback__)
            raise
        finally:
            log.reset_console()
            run_common._finish_stats(args)
    stages.finish()


class RecordedStage:
    """
    The starting process's side of one stage of a run carried on from a record: it knows only the stage's key and label, as
    the record does. A stage the run before finished is announced as skipped; any other opens a bar of its units.
    """

    def __init__(self, args, step):
        self.args, self.step = args, step
        self.progress = None
        self.workers = 1
        self.bar = None
        self._held = contextlib.ExitStack()

    def on_plan(self, plan_):
        from planetgen.generation import stages, steps as stepBars
        from planetgen.util import log
        carried = (plan_["detail"] or {}).get(wrapper.FINISHED_BY)
        if carried:
            stages.skip(self.step.key, f"finished by run {carried}")
            return
        stages.enter(self.step.key, queued=True)
        if plan_["units"]:
            log.normal(f"{self.step.label}: {plan_['units']:,} unit(s).")
            self.bar = self._held.enter_context(stepBars.Step(
                self.step.label, None, max(plan_["units"], 1), args=self.args, progress=self.progress,
                workers=self.workers, record=False))

    def on_result(self, unit, result):
        if self.bar is not None:
            self.bar.update(advance=1)

    def end(self, conclusion):
        self.bar = None
        self._held.close()

    def abort(self, exc_type, exc, traceback):
        self.bar = None
        self._held.__exit__(exc_type, exc, traceback)


__all__ = ["AGAIN", "CannotFinish", "FINISHABLE", "Plan", "REWIND_WITH", "add_arguments", "find", "last_incomplete", "plan",
           "run_finish", "stage_finished"]
