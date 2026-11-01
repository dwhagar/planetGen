# planetgen.cli.warm_progress

"""
What `planetgen warm --progress` shows while the warm-up runs (update.sh stage 1 and 2, OPS.88).

A stage that is cut into several units (chunks: the Galaxy Map's nodes `NODE_CHUNK` at a time, the object pages
`OBJECT_CHUNK` at a time) reports its own position, `galaxy-deep-nodes chunk 12 of 180`, the total known as soon as the
stage is queued, and the overall bar advances with the chunks, not only when a stage ends.

* With a terminal attached the `rich` bars are drawn (as the maintenance commands and the wipe step do): one for the
  whole warm-up and one for the stage in progress, each with elapsed time, time left and, for the stage, the seconds a
  chunk takes. Lines the warm-up logs are printed above the bars.
* Without one (a log file, cron, a pipe) each chunk leaves a compact line with its number, the nodes or pages it held,
  the built and already stored counts so far, the elapsed time and the time left, one line per chunk up to
  `LOG_CHUNK_LINES` chunks and one per several chunks beyond that (always the last), then a summary line for the stage.
  A stage of one unit prints its own line, as before.
"""

import math
import sys
import time

LOG_CHUNK_LINES = 20
"""int: Most chunk lines a stage leaves in a log; a stage of more chunks prints every Nth, and the last."""


def clock_text(seconds):
    """`seconds` as `1:05:12` (hours when there are any), or `?` when it is not known."""
    if seconds is None or seconds != seconds or seconds < 0:
        return "?"
    seconds = int(round(seconds))
    hours, rest = divmod(seconds, 3600)
    minutes, secs = divmod(rest, 60)
    return f"{hours}:{minutes:02d}:{secs:02d}" if hours else f"{minutes}:{secs:02d}"


def text_bar(done, total, label, width=30):
    """One line of the plain bar: `[#####-----] 3/10 label`."""
    filled = width * done // total if total else width
    return f"[{'#' * filled}{'-' * (width - filled)}] {done}/{total} {label}"


class WarmProgress:
    """
    Follows the stages of one warm-up: tell it `plan` when a stage is queued and `result` as each of its units
    finishes; it prints (`say`) the lines for a log, or draws the bars when `rich_bars` is set. Use it as a context
    manager so the bars are cleared at the end.

    Args:
        keys (list[str]): The stages' keys, in order.
        say (callable): Prints one line (the caller's log; replaced by the bars' console while they are drawn).
        show (bool): Show progress at all (`--progress`).
        rich_bars (bool): Draw the `rich` bars (a terminal is attached) instead of the log lines.
        clock (callable): Seconds, monotonic; replaced in tests.
    """

    def __init__(self, keys, say, show=True, rich_bars=False, clock=time.monotonic):
        self.keys = list(keys)
        self.say = say
        self.show = show
        self.rich_bars = show and rich_bars
        self.clock = clock
        self.progress = None
        self.overall = None
        self.task = None
        self.index = -1
        self.units = 0
        self.done = 0
        self.sums = {}
        self.noun = "items"
        self.biggest = 0
        self.started = 0.0
        self.last_printed = 0

    def __enter__(self):
        if self.rich_bars:
            from rich.progress import (BarColumn, Progress, ProgressColumn, TextColumn, TimeElapsedColumn,
                                       TimeRemainingColumn)
            from rich.text import Text

            class PerChunk(ProgressColumn):
                def render(self, task):
                    if not task.completed or not task.elapsed or task.fields.get("single"):
                        return Text("")
                    return Text(f"{task.elapsed / task.completed:.1f} s/chunk")

            self.progress = Progress(TextColumn("{task.description}", markup=False), BarColumn(),
                                     TimeElapsedColumn(), TimeRemainingColumn(), PerChunk(), transient=True)
            self.overall = self.progress.add_task("warm-up", total=len(self.keys), single=True)
            self.progress.start()
            self.say = self._console_say
        return self

    def __exit__(self, exc_type, exc, tb):
        if self.progress is not None:
            self.progress.stop()
        return False

    def _console_say(self, line):
        self.progress.console.print(line, markup=False, highlight=False)

    def plan(self, index, plan):
        """Stage `index` has been queued as `plan["units"]` units."""
        self.index = index
        self.units = int(plan.get("units") or 0)
        self.done = 0
        self.sums = {"warmed": 0, "built": 0, "cached": 0, "left_out": 0, "seconds": 0.0}
        self.noun = "items"
        self.biggest = 0
        self.last_printed = 0
        self.started = self.clock()
        if self.progress is not None:
            if self.task is not None:
                self.progress.remove_task(self.task)
            self.task = self.progress.add_task(self._description(), total=max(self.units, 1),
                                               single=self.units <= 1)
            self.progress.update(self.overall, completed=index, description=self._overall_text())

    def _description(self):
        key = self.keys[self.index]
        return f"{key} chunk {self.done} of {self.units}" if self.units > 1 else key

    def _overall_text(self):
        return f"warm-up, step {self.index + 1} of {len(self.keys)}"

    def result(self, result):
        """One unit of the current stage finished with `result`; prints what that calls for."""
        key = self.keys[self.index]
        self.done += 1
        for name in ("warmed", "built", "cached", "left_out"):
            self.sums[name] += int(result.get(name) or 0)
        self.sums["seconds"] += float(result.get("seconds") or 0.0)
        self.noun = result.get("noun") or self.noun
        self.biggest = max(self.biggest, int(result.get("warmed") or 0) + int(result.get("left_out") or 0))
        if self.progress is not None:
            self.progress.update(self.task, completed=self.done, description=self._description())
            self.progress.update(self.overall, completed=self.index + self.done / max(self.units, 1))
        if self.units <= 1:
            self.say(result["line"])
            if self.show and not self.rich_bars:
                self.say(text_bar(self.index + 1, len(self.keys), key))
            return
        elapsed = self.clock() - self.started
        if self.done == self.units:
            self.say(self.summary(key, elapsed))
            if self.show and not self.rich_bars:
                self.say(text_bar(self.index + 1, len(self.keys), key))
        elif self.show and not self.rich_bars and self._due():
            self.say(self.chunk_line(key, result, elapsed))

    def _due(self):
        step = max(1, math.ceil(self.units / LOG_CHUNK_LINES))
        return self.done % step == 0

    def chunk_line(self, key, result, elapsed):
        """`galaxy-deep-nodes chunk 12 of 180: this chunk 64 nodes in 38.2 s; so far ... elapsed 7:31, about 1:05:00 left`."""
        left = elapsed / self.done * (self.units - self.done)
        return (f"{key} chunk {self.done} of {self.units}: this chunk {int(result.get('warmed') or 0):,} {self.noun} in "
                f"{float(result.get('seconds') or 0.0):.1f} s; so far {self.sums['warmed']:,} {self.noun} "
                f"({self.sums['built']:,} built, {self.sums['cached']:,} already stored), {clock_text(elapsed)} "
                f"elapsed, about {clock_text(left)} left")

    def summary(self, key, elapsed):
        """The one line a stage of several chunks ends with."""
        out = f", {self.sums['left_out']:,} left out" if self.sums["left_out"] else ""
        return (f"{key}: {self.sums['warmed']:,} {self.noun} warmed ({self.sums['built']:,} built, "
                f"{self.sums['cached']:,} already stored{out}) in {clock_text(elapsed)} over {self.units} chunks of up "
                f"to {self.biggest:,} {self.noun}")


def stdout_is_terminal():
    """Whether lines go to a terminal (the rich bars can be drawn)."""
    return bool(getattr(sys.stdout, "isatty", lambda: False)())
