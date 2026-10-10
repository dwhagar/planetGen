# planetgen.cli.stage_progress

"""
The progress bar of the maintenance commands (`planetgen.cli.reset`,
`planetgen.cli.orbits`, the population pass): one bar over the command's
steps and, under it, a second bar over the items inside the current step --
the same two-bar shape a galaxy run draws (`run_common._generation_progress`),
so the terminal and the Generate page's job status
(`PLANETGEN_PROGRESS_FILE`) show these runs the way they show a generation
run. Both are `planetgen.generation.steps` (PERF.51): each draws when it is
predicted to take over 15 seconds (or has no recorded speed yet).
"""

from planetgen.generation import run_common, steps


class StageProgress:
    """
    A step bar with an optional item bar for the step in progress. Use it as
    a context manager; call `stage` as each step starts and `detail` as the
    items inside it finish.

    Args:
        total (int): How many steps the run has.
        disable (bool): Draw nothing (the progress file is still written).
        kind (str): What the steps' speed is recorded under (`generation_stats`).
        args (argparse.Namespace, optional): The run's arguments, for its stored speeds.
        stats (GenerationStats, optional): The stored speeds, when there are no `args` to read them from.
    """

    def __init__(self, total, disable=False, kind="stages", args=None, stats=None):
        self.total = total
        self.kind = kind
        self.args = args
        self.stats = stats
        self.progress = run_common._generation_progress(disable=disable)
        self.bar = None
        self.detail_bar = None
        self.detail_label = None
        self.done = 0

    def __enter__(self):
        self.progress.start()
        self.bar = steps.Step("", self.kind, self.total, args=self.args, stats=self.stats,
                              progress=self.progress).__enter__()
        return self

    def __exit__(self, exc_type, exc, tb):
        self._clear_detail()
        self.bar.close(success=exc_type is None)
        self.progress.stop()
        return False

    def stage(self, description):
        """Starts the next step: the one before it counts as done and its
        item bar goes away."""
        self._clear_detail()
        # A `total` makes the progress file write at once, so a slow step
        # shows its name before its first item finishes.
        self.bar.update(completed=self.done, total=self.total, description=description)
        self.done += 1

    def detail(self, label, done, total):
        """Shows `done` of `total` items of the current step under the step
        bar (the signature the store's `on_progress` hooks call)."""
        text = f"  {label}"
        if self.detail_bar is None or self.detail_label != label:
            self._clear_detail()
            self.detail_bar = steps.Step(text, f"{self.kind}:{label}", max(total, 1), args=self.args,
                                         stats=self.stats, progress=self.progress).__enter__()
            self.detail_label = label
        self.detail_bar.update(completed=done, total=max(total, 1), description=text)

    def _clear_detail(self):
        if self.detail_bar is not None:
            self.detail_bar.close(success=self.detail_bar.done >= (self.detail_bar.total or 0))
            self.detail_bar = None
            self.detail_label = None
