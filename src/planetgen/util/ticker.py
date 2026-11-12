# planetgen.util.ticker

"""
The one progress ticker (OPS.169): how a long loop says where it is
====================================================================

A step that works for minutes says where it is every `tuning.PROGRESS_LOG_SECONDS` and no more often, so a log read without the bar (a job's
stage log, `tail -f`, the run log) never goes quiet and never floods. Before this module four places kept their own clock and their own
constant for it (the block light units, the uncharted index, the download lines and the warm-up lines) and one built its own elapsed stamp.

* `Ticker(label, total, say)`: `tick(done, detail)` says `Label: 1,240 of 4,000 done (31%), 2.1 per s, 4:07 elapsed.` when the cadence has
  passed and returns whether it did; `finish(done)` says the last line whatever the clock. `percent=True` is for a counter on a bar's scale (0 to
  100), where a count and a rate mean nothing: `Label: 37% done, 4:07 elapsed.`
* `Pulse(say, text)`: a thread that says `text()` every cadence for as long as it runs (a context manager), for a stage whose unit may be silent
  between two counts, such as one long query.
* `stamp(text, began)`: `text` with the time since `began`, `Text (4:07 elapsed).`, the one form of the hand stamp the block light lines had.

The clock is `time.monotonic` and can be given, so a test needs no wait. The text is one sentence with a capital and a period (the logging
standard's message text). Time is `logrecord.duration_text`.
"""

import threading
import time

from planetgen import tuning
from planetgen.util import log, logrecord


def _seconds(every):
    return tuning.PROGRESS_LOG_SECONDS if every is None else every


def stamp(text, began, clock=time.monotonic):
    """`text` (a trailing period or ellipsis is dropped) with how long it has been since `began`: `Adding up the levels (4:07 elapsed).`"""
    return f"{text.rstrip('. ')} ({logrecord.duration_text(clock() - began)} elapsed)."


class Ticker:
    """
    Says where a loop is at most every `every` seconds (`tuning.PROGRESS_LOG_SECONDS` by default).

    Args:
        label (str): What is being counted, `Block light totals, units`.
        total (int, optional): How many there will be; left out, a line says only how many are done.
        say (callable, optional): `say(text)`; `log.normal` by default.
        every (float, optional): Seconds between two lines.
        clock (callable): `time.monotonic` by default.
        percent (bool): `done` is on a bar's scale (`total` is its end): say the share, not the count and the rate.
    """

    def __init__(self, label, total=None, say=None, every=None, clock=time.monotonic, percent=False):
        self.label, self.total, self.every, self.clock, self.percent = label, total, _seconds(every), clock, percent
        self._say = say if say is not None else log.normal
        self.began = self._last = self.clock()

    def text(self, done, detail=None):
        """The line for `done`: the count and share, the rate, the time since the start and `detail` (what else the loop knows)."""
        elapsed = self.clock() - self.began
        if self.percent:
            body = f"{done / self.total:.0%} done" if self.total else f"{done:,.0f} done"
        else:
            body = f"{done:,.0f} of {self.total:,.0f} done ({done / self.total:.0%})" if self.total else f"{done:,.0f} done"
            if done and elapsed >= 1.0:
                body += f", {done / elapsed:,.1f} per s"
        tail = f", {detail}" if detail else ""
        return f"{self.label}: {body}{tail}, {logrecord.duration_text(elapsed)} elapsed."

    def tick(self, done, detail=None):
        """Says the line when `every` seconds have passed since the last one (or the start). Returns whether it did."""
        now = self.clock()
        if now - self._last < self.every:
            return False
        self._last = now
        self._say(self.text(done, detail))
        return True

    def finish(self, done, detail=None):
        """Says the closing line, whatever the clock."""
        self._last = self.clock()
        self._say(self.text(done, detail))


class Pulse:
    """
    A thread that says `text()` every `every` seconds while the block runs: `with Pulse(say, lambda: "Still working.") as pulse:`. For work that
    can go quiet between two counts (one long query), where a `Ticker` that waits for the next count would say nothing. A `text` that raises is
    skipped for that round and never stops the work.
    """

    def __init__(self, say, text, every=None):
        self._say, self._text, self.every = say, text, _seconds(every)
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, name="progress-pulse", daemon=True)

    def _run(self):
        while not self._stop.wait(self.every):
            try:
                self._say(self._text())
            except Exception as exc:  # noqa: BLE001 -- a progress line is never worth the work it reports on
                log.swallowed("util/ticker", "Pulse._run", exc)
                continue

    def start(self):
        self._thread.start()
        return self

    def stop(self):
        self._stop.set()

    def __enter__(self):
        return self.start()

    def __exit__(self, *exc_info):
        self.stop()
        return False
