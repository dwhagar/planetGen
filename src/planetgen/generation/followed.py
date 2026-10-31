# planetgen/generation/followed.py

"""
A run followed from the process that started it (PERF.172).

A run is a list of stages the wrapper queues (`queue/wrapper.py`); the process that started it only follows. A
`Follower` is that process's side of one stage: the parameters the stage's step needs, and what to show and say as the
stage's plan arrives, its units come back and the stage ends (bars, lines of the log, the summary the caller wants).
`run_followed` puts any list of followers in one run: `planetgen link` is the run of the link follower, and
`planetgen galaxy` the run of the followers of every stage it takes, so a galaxy is one run in the queue, with one
estimate over all of it.

A follower does no work of the stage's: that is the command's (`queue/commands.py`), in a unit or in the command's own
`prepare` and `conclude`, and runs in a worker. Everything a follower needs is in the step's parameters, which are
complete before the run starts, and in the stage's plan and results.
"""

import contextlib

from planetgen.db import store
from planetgen.generation import finish_run, run_common, stages
from planetgen.queue import runrecord, wrapper
from planetgen.util import log


class Follower:
    """
    The side of stage `key` the starting process takes. A subclass names its `key` and `label`, gives `params` (the
    step's parameters besides `database` and `workers`) and overrides what it shows or counts: `on_plan`, `on_result`,
    `on_progress` and `on_done`. A bar or a block it opens for the stage is held with `hold`, and ends when the stage
    does (or when the run stops, and then without recording its speed).

    Attributes:
        args (argparse.Namespace): The run's arguments.
        progress (rich.progress.Progress): The run's one display, set when the run starts.
        workers (int): The workers the run uses.
    """

    key = None
    label = ""
    own_row = True
    """bool: Whether the unit that completes the stage stores its `generation_stage_runs` row (false for a stage the
    process that started the run keeps open past its units, and stores for all of it)."""
    repeatable = False
    """bool: Whether a run that several commands put the stage in does it once, at the last place (PERF.203): the stage
    redoes a total that each of them changed (`dedupe`)."""

    def __init__(self, args):
        self.args = args
        self.progress = None
        self.workers = 1
        self._held = contextlib.ExitStack()

    def params(self):
        """The step's parameters besides `database` and `workers`: plain JSON, complete before the run starts."""
        return {}

    def hold(self, context):
        """Enters a context manager (a bar, a job node) that lasts until the stage ends; returns what it gives."""
        return self._held.enter_context(context)

    def on_plan(self, plan):
        """The stage has been cut: `plan` is `{"kind", "units", "cost", "detail"}`."""

    def on_result(self, unit, result):
        """Unit `unit` is done and came back with `result`."""

    def on_progress(self, unit, done, estimate):
        """Unit `unit`, still working, has done `done` of `estimate`."""

    def on_done(self, conclusion):
        """Every unit is done and the command's `conclude` has run; `conclusion` is what it returned (or `None`)."""

    def end(self, conclusion):
        """The stage ends: `on_done`, then whatever it held is closed."""
        try:
            self.on_done(conclusion)
        finally:
            self._held.close()

    def abort(self, exc_type, exc, traceback):
        """The run stops with the stage open: what it held is closed as a failure."""
        self._held.__exit__(exc_type, exc, traceback)


class Planned:
    """
    A run a verb is going to make, before it is started (PERF.203): the followers of its stages in order, the run's title, the
    options `run_followed` takes besides them (`seed`, `kind`, `workers`, `galaxy`) and `result`, the callable that gives what the
    verb returns once the run has ended (it reads the followers). `planetgen <verb>` starts its own with `run_planned`; the one
    run of a whole Generate job (`planetgen run-job`) takes the followers of every selected verb from the same functions, so the
    two cannot differ.

    Attributes:
        followers (list[Follower]): One for each stage, in order.
        title (str): What the run is called.
        options (dict): The keyword arguments of `run_followed`.
        result (callable or None): `result()` is what the verb returns.
    """

    __slots__ = ("followers", "title", "options", "result")

    def __init__(self, followers, title, result=None, **options):
        self.followers, self.title, self.options, self.result = list(followers), title, options, result


def run_planned(args, planned):
    """Runs a `Planned` as one run (`run_followed`) and returns what its `result` gives, `None` without one."""
    run_followed(args, planned.followers, planned.title, **planned.options)
    return planned.result() if planned.result is not None else None


def dedupe(followers):
    """
    The followers of a run made of several commands' (PERF.203): a `repeatable` stage that more than one of them ends with
    (the block light totals) is kept at its last place only, since the last redoes what all of them changed.

    Raises:
        ValueError: A stage that is not `repeatable` is there twice: the run's steps are keyed by their stage.
    """
    last = {f.key: place for place, f in enumerate(followers) if f.repeatable}
    kept = [f for place, f in enumerate(followers) if not f.repeatable or last[f.key] == place]
    keys = [f.key for f in kept]
    for key in keys:
        if keys.count(key) > 1:
            raise ValueError(f"the stage {key!r} comes twice in one run")
    return kept


def _enter(follower):
    """Enters the follower's stage for the run's numbering, as one whose row the queue stores unless it keeps its own."""
    key = follower.key
    if not follower.own_row:
        stages.enter(key)
    elif stages.open_key() == key:
        stages.hand_over_open_stage()
    else:
        stages.enter(key, queued=True)


def run_followed(args, followers, title, seed=None, kind=runrecord.GENERATION, workers=None, galaxy=True, settings=None):
    """
    Runs the stages of `followers` as one run in the queue and follows it to its end: the steps are the followers'
    (`wrapper.plan_steps`, each with its parameters), the stage rows are handed to the queue up front
    (`stages.row_for`), and each stage is entered (`stages.enter`) as its plan arrives.

    Args:
        args (argparse.Namespace): The run's arguments.
        followers (list[Follower]): One for each stage, in order.
        title (str): What the run is called.
        seed (int, optional): The seed to record with the run.
        kind (str): The kind of run (`runrecord.GENERATION`, or `runrecord.MAINTENANCE` for one that saves nothing).
        workers (int, optional): The workers to use, the run's `--workers` when not given.
        galaxy (bool): Whether to record the stored galaxy's seed and build with the run, so `planetgen finish-run` can
            carry it on (false for a run that belongs to no galaxy: one sector or one system).
        settings (dict, optional): `{key: settings}` of the stages of a run of several commands (`wrapper.plan_steps`).

    Raises:
        redisqueue.Unavailable: No Redis server answers.
        runrecord.Busy: Another run is running and this one may not wait.
        units.RunStopped: A unit failed or the run was cancelled.
    """
    followers = dedupe(followers)
    mysql_config = store.mysql_config_from_args(args)
    workers = workers or run_common._worker_count(args)
    database = wrapper.database_args(mysql_config)
    galaxy_seed, build = finish_run.current_galaxy(args) if galaxy else (None, None)
    parts = [(f.key, f.label, {"database": database, "workers": workers, **f.params()}) for f in followers]
    steps, listing = wrapper.plan_steps(args, parts) if settings is None else wrapper.plan_steps(args, parts, settings)
    rows = {f.key: row for f in followers if f.own_row and (row := stages.row_for(f.key)) is not None}
    open_ones = []

    def on_plan(index, plan):
        _enter(followers[index])
        open_ones.append(followers[index])
        followers[index].on_plan(plan)

    def on_result(index, unit, result):
        followers[index].on_result(unit, result)

    def on_progress(index, unit, done, estimate):
        followers[index].on_progress(unit, done, estimate)

    def on_stage_done(index, conclusion):
        followers[index].end(conclusion)
        open_ones.remove(followers[index])

    # a document written to stdout (`--output -`) must start at its first line: the bar's display, even a finished one, writes a blank line
    with run_common._generation_progress(disable=getattr(args, "output", None) == "-") as progress:
        log.set_console(progress.console)
        for follower in followers:
            follower.progress, follower.workers = progress, workers
        try:
            wrapper.run(title, steps, on_plan, on_result, on_progress, workers_count=workers,
                        who=getattr(args, "who", None), database=mysql_config.database,
                        argv=getattr(args, "argv", None), seed=seed, listing=listing, rows=rows or None, kind=kind,
                        on_stage_done=on_stage_done, galaxy_seed=galaxy_seed, build=build)
        except BaseException as exc:
            for follower in open_ones:
                follower.abort(type(exc), exc, exc.__traceback__)
            raise
        finally:
            log.reset_console()
            run_common._finish_stats(args)
