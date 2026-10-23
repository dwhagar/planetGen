# planetgen/util/step_timer.py

"""
A list of how long each step of a command took, for the commands that print a timing table (the plan, PERF.163).
"""

import contextlib
import time


class StepTimer:
    """The steps of one piece of work, in the order they ran, with the seconds each took."""

    def __init__(self):
        self.steps = []

    @contextlib.contextmanager
    def step(self, label):
        """Times the block under `label`; a block that raises is not recorded."""
        started = time.perf_counter()
        yield
        self.steps.append((label, time.perf_counter() - started))

    def total(self):
        """The seconds of all the steps."""
        return sum(seconds for _label, seconds in self.steps)

    def lines(self, title):
        """The table as lines of text: the `title`, then each step with its seconds and share of the total, the
        costliest first."""
        total = self.total()
        shown = sorted(self.steps, key=lambda item: -item[1])
        return [f"{title} ({total:.2f}s):"] + [
            f"  {seconds:8.2f}s {100.0 * seconds / total if total else 0.0:5.1f}%  {label}" for label, seconds in shown]
