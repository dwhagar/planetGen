# planetgen.cli.stage_progress

"""
The progress bar of the maintenance commands (`planetgen.cli.reset`,
`planetgen.cli.orbits`): one bar over the command's steps and, under it, a
second bar over the items inside the current step -- the same two-bar shape
a galaxy run draws (`run_common._generation_progress`), so the terminal and
the Generate page's job status (`PLANETGEN_PROGRESS_FILE`) show these runs
the way they show a generation run.
"""

from planetgen.generation import run_common


class StageProgress:
    """
    A step bar with an optional item bar for the step in progress. Use it as
    a context manager; call `stage` as each step starts and `detail` as the
    items inside it finish.

    Args:
        total (int): How many steps the run has.
        disable (bool): Draw nothing (the progress file is still written).
    """

    def __init__(self, total, disable=False):
        self.total = total
        self.progress = run_common._generation_progress(disable=disable)
        self.task = None
        self.detail_task = None
        self.done = 0

    def __enter__(self):
        self.progress.start()
        self.task = self.progress.add_task("", total=self.total)
        self.progress.main_task = self.task
        return self

    def __exit__(self, exc_type, exc, tb):
        if exc_type is None:
            self._clear_detail()
            self.progress.update(self.task, completed=self.total, total=self.total)
        self.progress.stop()
        return False

    def stage(self, description):
        """Starts the next step: the one before it counts as done and its
        item bar goes away."""
        self._clear_detail()
        # A `total` makes the progress file write at once, so a slow step
        # shows its name before its first item finishes.
        self.progress.update(self.task, completed=self.done, total=self.total, description=description)
        self.done += 1

    def detail(self, label, done, total):
        """Shows `done` of `total` items of the current step under the step
        bar (the signature the store's `on_progress` hooks call)."""
        label = f"  {label}"
        if self.detail_task is None:
            self.detail_task = self.progress.add_task(label, total=max(total, 1), completed=done)
            self.progress.detail_task = self.detail_task
        self.progress.update(self.detail_task, completed=done, total=max(total, 1), description=label)

    def _clear_detail(self):
        if self.detail_task is not None:
            self.progress.remove_task(self.detail_task)
            self.progress.detail_task = None
            self.detail_task = None
