# planetgen.cli.warm_progress

"""
What `planetgen warm --progress` shows while the warm-up runs (update.sh stage 1 and 2, OPS.88).

A stage that is cut into several units (chunks: the Galaxy Map's nodes `NODE_CHUNK` at a time, the object pages
`OBJECT_CHUNK` at a time) reports its own position, `galaxy-deep-nodes chunk 12 of 180`, the total known as soon as the
stage is queued, and the overall bar advances with the chunks, not only when a stage ends.

The overall bar weighs a step by what it is expected to cost (PERF.234): a step already queued by the cost its units
declare (`warmup.ITEM_SECONDS` an item), a step not yet queued by `expected` (`warmup.EXPECTED_SECONDS`, measured on a
test galaxy), so the long object steps move it as far as they take, not one step in eleven. Its time left is a rough guess
while a later step has not been queued yet, and says so. No time left is shown for the stage or the whole warm-up from
fewer than `MIN_CHUNKS` chunks (PERF.237): the first results are the first of the workers that run side by side.

A unit that reports how far it has got (`unit_progress`: the node chunks and the object slices say so as they work) moves the
bars between the chunk ends too: its fraction counts as that share of a chunk done, for the stage's position, the cost done
and the time left, so a stage of a few big chunks (a big galaxy) is not still at zero until the first one ends.

The seconds a chunk reports are its own worker's. Several chunks run at once, so the line also says how many (the
chunks' own seconds added up over the elapsed time) and how long a chunk takes overall (the elapsed time over the chunks
done); the time left comes from that overall rate only.

* With a terminal attached the `rich` bars are drawn (as the maintenance commands and the wipe step do): one for the
  whole warm-up and one for the stage in progress, each with elapsed time, time left and, for the stage, the seconds a
  chunk takes overall and how many run at once. Lines the warm-up logs are printed above the bars.
* Without one (a log file, cron, a pipe) each chunk leaves a compact line with its number, the nodes or pages it held,
  its own worker's seconds, the chunks at once, the built and already stored counts so far, the elapsed time and the
  time left, one line per chunk up to `LOG_CHUNK_LINES` chunks and one per several chunks beyond that (always the last),
  then a summary line for the stage. A stage of one unit prints its own line, as before.
"""

import math
import sys
import time

LOG_CHUNK_LINES = 20
"""int: Most chunk lines a stage leaves in a log; a stage of more chunks prints every Nth, and the last."""

MIN_CHUNKS = 5
"""int: Chunks that must be done before a time left is shown (PERF.237): fewer say nothing of how long the rest takes."""

MIN_SECONDS = 0.5
"""float: Elapsed time below which no rate is worked out."""

FALLBACK_SECONDS = 1.0
"""float: What a step with no expected cost and no queued cost weighs in the overall bar."""


def clock_text(seconds):
    """`seconds` as `1:05:12` (hours when there are any), or `?` when it is not known."""
    if seconds is None or seconds != seconds or seconds < 0:
        return "?"
    seconds = int(round(seconds))
    hours, rest = divmod(seconds, 3600)
    minutes, secs = divmod(rest, 60)
    return f"{hours}:{minutes:02d}:{secs:02d}" if hours else f"{minutes}:{secs:02d}"


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
        expected (dict): `{stage key: cost}` the steps are expected to cost, in one worker's seconds, before each is queued
            (`warmup.EXPECTED_SECONDS`); a step queued weighs the cost its plan declares instead.
        mirror (callable): Also given each line the drawn bars print (they bypass the logging layer), for the warm-up log (OPS.96).
    """

    def __init__(self, keys, say, show=True, rich_bars=False, clock=time.monotonic, expected=None, mirror=None):
        self.keys = list(keys)
        self.say = say
        self.mirror = mirror
        self.show = show
        self.rich_bars = show and rich_bars
        self.clock = clock
        self.expected = dict(expected or {})
        self.weights = [max(float(self.expected.get(key) or 0.0), FALLBACK_SECONDS) for key in self.keys]
        """list[float]: Each step's cost: its `expected` until the step is queued, then its queued cost."""
        self.queued = 0
        """int: Steps queued so far."""
        self.progress = None
        self.overall = None
        self.task = None
        self.index = -1
        self.units = 0
        self.done = 0
        self.sums = {}
        self.inflight = {}
        """dict[int, float]: Unit number to the fraction (0 to 1) of it done, for the units of the current stage that have reported."""
        self.noun = "items"
        self.biggest = 0
        self.started = 0.0
        self.run_started = None
        self.last_printed = 0

    def __enter__(self):
        if self.rich_bars:
            from rich.progress import BarColumn, Progress, ProgressColumn, TextColumn, TimeElapsedColumn
            from rich.text import Text

            class Field(ProgressColumn):
                """Shows a text the progress keeps in the task's fields."""

                def __init__(self, name):
                    super().__init__()
                    self.name = name

                def render(self, task):
                    return Text(str(task.fields.get(self.name) or ""))

            self.progress = Progress(TextColumn("{task.description}", markup=False), BarColumn(),
                                     TimeElapsedColumn(), Field("left"), Field("rate"), transient=True)
            self.overall = self.progress.add_task("warm-up", total=sum(self.weights), single=True, left="", rate="")
            self.progress.start()
            self.say = self._console_say
        return self

    def __exit__(self, exc_type, exc, tb):
        if self.progress is not None:
            self.progress.stop()
        return False

    def _console_say(self, line):
        self.progress.console.print(line, markup=False, highlight=False)
        if self.mirror is not None:      # the bars' console is not the logging layer: the warm-up log gets the line here (OPS.96)
            self.mirror(line)

    def plan(self, index, plan):
        """Stage `index` has been queued as `plan["units"]` units."""
        self.index = index
        self.queued = max(self.queued, index + 1)
        self.units = int(plan.get("units") or 0)
        self.done = 0
        self.inflight = {}
        self.sums = {"warmed": 0, "built": 0, "cached": 0, "left_out": 0, "seconds": 0.0}
        self.noun = "items"
        self.biggest = 0
        self.last_printed = 0
        self.started = self.clock()
        if self.run_started is None:
            self.run_started = self.started
        cost = float(plan.get("cost") or 0.0)
        if cost > 0:
            self.weights[index] = cost
        if self.progress is not None:
            if self.task is not None:
                self.progress.remove_task(self.task)
            self.task = self.progress.add_task(self._description(), total=max(self.units, 1),
                                               single=self.units <= 1, left="", rate="")
            self.progress.update(self.overall, total=sum(self.weights), completed=self._weight_done(),
                                 description=self._overall_text(), left=self._overall_left_text())

    def _description(self):
        key = self.keys[self.index]
        return f"{key} chunk {self.done} of {self.units}" if self.units > 1 else key

    def _overall_text(self):
        return f"warm-up, step {self.index + 1} of {len(self.keys)}"

    def _weight_done(self):
        """The cost of the steps before this one and the share of this one's that its chunks have done."""
        share = min(1.0, self._position() / self.units) if self.units > 0 else 0.0
        return sum(self.weights[:self.index]) + self.weights[self.index] * share

    def _position(self):
        """Chunks done, with the share of each chunk that reported a part of itself (`unit_progress`) added."""
        return self.done + sum(self.inflight.values())

    def _enough(self):
        return self._position() >= MIN_CHUNKS

    def unit_progress(self, index, unit, done, estimate):
        """Unit `unit` of stage `index` has done `done` of `estimate` of its work (`wrapper.run`'s `on_progress`); moves the bars."""
        if index != self.index or estimate is None or estimate <= 0:
            return
        fraction = max(0.0, min(1.0, float(done) / float(estimate)))
        if fraction <= self.inflight.get(unit, 0.0):
            return
        self.inflight[unit] = fraction
        if self.progress is not None:
            elapsed = self.clock() - self.started
            stage_left = self._stage_left(elapsed)
            self.progress.update(self.task, completed=min(self._position(), self.units),
                                 left="" if stage_left is None else f"about {clock_text(stage_left)} left")
            self.progress.update(self.overall, completed=self._weight_done(), left=self._overall_left_text())

    def _stage_left(self, elapsed):
        """Seconds the stage still takes at the rate it has run at (`None` while too few chunks are done to say)."""
        if not self._enough() or elapsed < MIN_SECONDS:
            return None
        return elapsed / self._position() * (self.units - self._position())

    def _at_once(self, elapsed):
        """How many chunks run side by side: the chunks' own seconds added up over the elapsed time (`None` if unknown)."""
        if not self._enough() or elapsed < MIN_SECONDS or self.sums["seconds"] <= 0:
            return None
        return max(1, round(self.sums["seconds"] / elapsed))

    def overall_left(self):
        """Seconds the whole warm-up still takes at the rate it has run at, by cost (`None` while too few chunks are done)."""
        if self.run_started is None or not self._enough():
            return None
        elapsed = self.clock() - self.run_started
        done = self._weight_done()
        if elapsed < MIN_SECONDS or done <= 0:
            return None
        return elapsed * (sum(self.weights) - done) / done

    def rough(self):
        """Whether a step after this one has not been queued yet, so its cost is a guess."""
        return self.queued < len(self.keys)

    def _overall_left_text(self):
        left = self.overall_left()
        if left is None:
            return ""
        return f"about {clock_text(left)} left" + (" (rough guess)" if self.rough() else "")

    def result(self, result, unit=None):
        """One unit (number `unit`) of the current stage finished with `result`; prints what that calls for."""
        key = self.keys[self.index]
        self.done += 1
        self.inflight.pop(unit, None)
        for name in ("warmed", "built", "cached", "left_out"):
            self.sums[name] += int(result.get(name) or 0)
        self.sums["seconds"] += float(result.get("seconds") or 0.0)
        self.noun = result.get("noun") or self.noun
        self.biggest = max(self.biggest, int(result.get("warmed") or 0) + int(result.get("left_out") or 0))
        elapsed = self.clock() - self.started
        if self.progress is not None:
            stage_left = self._stage_left(elapsed)
            at_once = self._at_once(elapsed)
            self.progress.update(
                self.task, completed=self.done, description=self._description(),
                left="" if stage_left is None else f"about {clock_text(stage_left)} left",
                rate="" if at_once is None else f"{elapsed / self.done:.1f} s/chunk overall, {at_once} at once")
            self.progress.update(self.overall, completed=self._weight_done(),
                                 left=self._overall_left_text())
        if self.units <= 1:
            self.say(result["line"])
            if self.show and not self.rich_bars:
                self.say(self.bar_line(key, finished=True))
            return
        if self.done == self.units:
            self.say(self.summary(key, elapsed))
            if self.show and not self.rich_bars:
                self.say(self.bar_line(key, finished=True))
        elif self.show and not self.rich_bars and self._due():
            self.say(self.chunk_line(key, result, elapsed))

    def _due(self):
        step = max(1, math.ceil(self.units / LOG_CHUNK_LINES))
        return self.done % step == 0

    def bar_line(self, key, finished=False):
        """The plain bar at a stage's end: weighted by the steps' cost, with the step count and the whole warm-up's time left."""
        done = sum(self.weights[:self.index + 1]) if finished else self._weight_done()
        filled = int(30 * done / sum(self.weights)) if sum(self.weights) else 30
        left = self._overall_left_text()
        return (f"[{'#' * filled}{'-' * (30 - filled)}] {self.index + 1}/{len(self.keys)} {key}"
                + (f", warm-up {left}" if left else ""))

    def chunk_line(self, key, result, elapsed):
        """
        `galaxy-deep-nodes chunk 12 of 180: this chunk 64 nodes took 38.2 s on its worker; about 10 chunks at once, 3.8 s a
        chunk overall; so far ..., 7:31 elapsed, about 1:05:00 left`. The time left, from the overall rate, and the chunks at
        once are left out until `MIN_CHUNKS` chunks are done.
        """
        left = self._stage_left(elapsed)
        at_once = self._at_once(elapsed)
        together = ("" if at_once is None else
                    f"; about {at_once} chunk{'' if at_once == 1 else 's'} at once, {elapsed / self.done:.1f} s a chunk overall")
        wait = "time left not known yet" if left is None else f"about {clock_text(left)} left"
        return (f"{key} chunk {self.done} of {self.units}: this chunk {int(result.get('warmed') or 0):,} {self.noun} took "
                f"{float(result.get('seconds') or 0.0):.1f} s on its worker{together}; so far {self.sums['warmed']:,} "
                f"{self.noun} ({self.sums['built']:,} built, {self.sums['cached']:,} already stored), "
                f"{clock_text(elapsed)} elapsed, {wait}")

    def summary(self, key, elapsed):
        """The one line a stage of several chunks ends with."""
        out = f", {self.sums['left_out']:,} left out" if self.sums["left_out"] else ""
        return (f"{key}: {self.sums['warmed']:,} {self.noun} warmed ({self.sums['built']:,} built, "
                f"{self.sums['cached']:,} already stored{out}) in {clock_text(elapsed)} over {self.units} chunks of up "
                f"to {self.biggest:,} {self.noun}")


def stdout_is_terminal():
    """Whether lines go to a terminal (the rich bars can be drawn)."""
    return bool(getattr(sys.stdout, "isatty", lambda: False)())
