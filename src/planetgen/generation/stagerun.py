# planetgen/generation/stagerun.py

"""
The stages of a generating verb as parts of ONE wrapper run (PERF.202).

Every generating verb used to carry its own starter around a one-step `wrapper.run`: the step's parameters, the
stage row it hands over, the progress bar it draws and the totals it folds out of the units' results were all
written into the starter, so a verb could only be run as a run of its own. A `Part` holds what the starter knew
about its stage, apart from the starting, and `run_parts` follows any number of them, in order, in one run: one
run record, one progress display, a bar for the stage that is working, and the same refusals (a stage that
cannot go on says so from its `expand` or its unit, and the run stops with that message, as `refusal_of` hands
back).

What a stage does is in its command's `expand` and unit (`queue.commands`); what a part adds is only how the
process that follows the run shows it. The Generate page and the API jobs (PERF.203) build the parts of the
selected commands and give them to `run_parts` once.
"""

import contextlib
import sys
from dataclasses import dataclass, field
from typing import Callable, Optional

from planetgen.db import store
from planetgen.generation import run_common, steps
from planetgen.queue import commands, units, wrapper
from planetgen.util import log

REFUSAL_MARKER = "RuntimeError: "
"""str: What precedes a unit's or an expansion's refusal in the text of the run's error (`units.RunStopped`)."""


@dataclass
class Part:
    """
    One stage of a run, as the process that follows the run shows it.

    Args:
        step (commands.Step): The stage's command and parameters (`wrapper.step_for`).
        listing (list[commands.Step]): The whole run's steps with this one in place (`wrapper.step_for`).
        bar (str, optional): The progress bar's description; none draws no bar.
        kind (str, optional): The `generation_stats` kind the bar's speed is predicted from (`steps.Step`).
        total (callable, optional): `total(plan)`, the bar's total from the stage's plan (default: the plan's
            `detail["sectors"]`, at least 1).
        advance (callable, optional): `advance(result)`, how far a unit's result moves the bar (default: none).
        sums (tuple[str]): Names of integer fields of the units' results that are added up into `totals`.
        fold (callable, optional): `fold(result)`, anything more the starter takes out of a unit's result.
        row (dict, optional): The stage's `generation_stage_runs` row (`stages.hand_over_open_stage`).
        totals (dict): What `sums` add up to, by field name, once the stage is done.
    """

    step: commands.Step
    listing: list = field(default_factory=list)
    bar: Optional[str] = None
    kind: Optional[str] = None
    total: Optional[Callable] = None
    advance: Optional[Callable] = None
    sums: tuple = ()
    fold: Optional[Callable] = None
    row: Optional[dict] = None
    totals: dict = field(default_factory=dict)

    def __post_init__(self):
        self.totals = {name: 0 for name in self.sums}

    def take(self, result):
        """What this part does with one unit's result: adds up its `sums` and calls its `fold`."""
        for name in self.sums:
            self.totals[name] += result[name]
        if self.fold is not None:
            self.fold(result)

    def bar_total(self, plan):
        return max(self.total(plan) if self.total is not None else plan["detail"]["sectors"], 1)


def listing_of(parts):
    """The run's whole step list for `parts`: the first part's `listing`, with each later part's step in place of the
    step of its key (a part whose key the list lacks is added at the end)."""
    listing = list(parts[0].listing) if parts[0].listing else [parts[0].step]
    for part in parts[1:]:
        keys = [item.key for item in listing]
        if part.step.key in keys:
            listing[keys.index(part.step.key)] = part.step
        else:
            listing.append(part.step)
    return listing


def refusal_of(exc):
    """
    The refusal a stage made, as a `RuntimeError` of its own words, from the `units.RunStopped` a run raises when its
    expansion or unit raised one (the error text is "RuntimeError: ..." after the worker's own); `exc` itself when
    it is not one.
    """
    text = str(exc)
    if isinstance(exc, units.RunStopped) and REFUSAL_MARKER in text:
        return RuntimeError(text.split(REFUSAL_MARKER, 1)[1])
    return exc


def run_parts(args, title, parts, *, workers=None, on_progress=None, show_bars=True, finally_=None):
    """
    Runs `parts` as the stages of one wrapper run, in order, and follows them: the bar of the stage that is working is
    drawn on the run's one progress display, each unit's result goes to its part (`Part.take`), and the part's bar
    advances. A refusal a stage makes is raised as the `RuntimeError` of its own words (`refusal_of`).

    Args:
        args (argparse.Namespace): The run's arguments (the database, `who`, `argv`).
        title (str): The run's title.
        parts (list[Part]): The stages, in order.
        workers (int, optional): The workers the run uses (default: the run's `--workers`).
        on_progress (callable, optional): `on_progress(index, unit, done, estimate)` for the units still working.
        show_bars (bool): Whether to draw the bars (a caller that draws its own passes `False`).
        finally_ (callable, optional): Called once the run has ended, however it ended.

    Returns:
        str: The run id.
    """
    mysql_config = store.mysql_config_from_args(args)
    workers = run_common._worker_count(args) if workers is None else workers
    bars = {}
    rows = {part.step.key: part.row for part in parts if part.row is not None}

    def on_plan(index, plan):
        part = parts[index]
        if part.bar is not None and show_bars:
            bar = steps.Step(part.bar, part.kind, 1, args=args, progress=progress, workers=workers)
            bar.__enter__()
            bars[index] = bar
            bar.update(total=part.bar_total(plan))

    def on_result(index, _unit, result):
        part = parts[index]
        part.take(result)
        if part.advance is not None and index in bars:
            bars[index].update(advance=part.advance(result))

    def on_stage_done(index):
        bar = bars.pop(index, None)
        if bar is not None:
            bar.__exit__(None, None, None)

    with run_common._generation_progress() as progress:
        log.set_console(progress.console)
        try:
            try:
                return wrapper.run(title, [part.step for part in parts], on_plan, on_result, on_progress,
                                   workers_count=workers, who=getattr(args, "who", None),
                                   database=mysql_config.database, argv=getattr(args, "argv", None), seed=None,
                                   listing=listing_of(parts), rows=rows or None, on_stage_done=on_stage_done)
            except units.RunStopped as exc:
                refusal = refusal_of(exc)
                if refusal is exc:
                    raise
                raise refusal from exc
            finally:
                for bar in list(bars.values()):
                    with contextlib.suppress(Exception):
                        bar.__exit__(*sys.exc_info())
                bars.clear()
                if finally_ is not None:
                    finally_()
        finally:
            log.reset_console()
