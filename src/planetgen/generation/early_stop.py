# planetgen/generation/early_stop.py

"""
Early stop for the layer-walking scatter passes: the layers are
walked from the galactic plane outwards, and once `tuning.SCATTER_DRY_LAYERS`
layers in a row, in walk order, have produced nothing, the rest are not
walked. The count is of walk order, not of finishing order, so any number of
workers stops at the same layer and a layer that finishes late cannot hide a
productive one behind a run of empty ones.
"""

from planetgen import tuning


class DryStreak:
    """Counts the empty layers in a row, in walk order, as layers report.

    Args:
        total (int): How many layers the pass walks.
        limit (int): How many empty layers in a row end the walk
            (`None` for `tuning.SCATTER_DRY_LAYERS`, 0 for never).
    """

    def __init__(self, total, limit=None):
        self.total = total
        self.limit = tuning.SCATTER_DRY_LAYERS if limit is None else limit
        self._results = {}
        self._next = 0
        self.dry = 0
        self.stopped_at = None

    def record(self, position, produced):
        """A layer (its place in the walk) finished having produced `produced` (a count)."""
        self._results[position] = bool(produced)
        while self.stopped_at is None and self._next in self._results:
            self.dry = 0 if self._results[self._next] else self.dry + 1
            self._next += 1
            if self.limit and self.dry >= self.limit and self._next < self.total:
                self.stopped_at = self._next

    @property
    def stopped(self):
        return self.stopped_at is not None

    def left_out(self, walked):
        """How many layers were never walked, once `walked` were queued."""
        return self.total - walked

    def message(self, label, walked):
        """The line for the log when the walk ended early."""
        left = self.left_out(walked)
        return (f"{label}: stopped after {self.limit:,} layers in a row with nothing; "
                f"{left:,} of {self.total:,} layers not walked.")
