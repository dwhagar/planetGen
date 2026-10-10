# planetgen/generation/layer_groups.py

"""
Grouping empty layers in the layer-walking scatter passes (PERF.57).

The passes walk the layers from the galactic plane outwards, each layer a
task of its own. Once `tuning.SCATTER_EMPTY_LAYERS_BEFORE_GROUP` layers in a
row took nothing, the next `tuning.SCATTER_GROUP_LAYERS` layers are drawn as
one group: the expected count over all of the group's sectors, the objects
placed in sectors by their density. A group that places nothing is followed
by one twice the size, and so on until something is placed or the layers run
out; after a placement the walk goes back to single layers.

The walk is decided by the layers' results in walk order, so the same seed
gives the same objects with any number of workers: a layer is only started
early when no result still to come could have put the walk into a group by
then. A layer expected to hold `tuning.SCATTER_CERTAIN_OBJECTS` or more is
counted as certain to take some, so the dense middle of the galaxy still runs
side by side.
"""

from planetgen import tuning


class Walk:
    """
    The order in which a pass draws its layers, singly and in groups.

    Args:
        expected (list): The expected objects of each layer, in walk order.
        empty_before (int): Layers in a row that took nothing before the next
            are grouped (`None` for the tuning constant).
        group_layers (int): Size of the first group (`None` for the constant).
        certain (float): Expected objects from which a layer counts as certain
            to take some (`None` for the constant).
    """

    def __init__(self, expected, empty_before=None, group_layers=None, certain=None):
        self.expected = list(expected)
        self.total = len(self.expected)
        self.empty_before = tuning.SCATTER_EMPTY_LAYERS_BEFORE_GROUP if empty_before is None else empty_before
        self.first_group = tuning.SCATTER_GROUP_LAYERS if group_layers is None else group_layers
        self.certain = tuning.SCATTER_CERTAIN_OBJECTS if certain is None else certain
        self.group_layers = self.first_group
        self.status = [None] * self.total   # None until it reports, then whether it took anything
        self.group_sizes = []               # the size of each group drawn, in order
        self.singles = 0

    # -- results -----------------------------------------------------------

    def record(self, position, produced):
        """Position `position` (a single layer) took `produced` objects."""
        self.status[position] = bool(produced)

    def record_group(self, start, size, produced):
        """The group of `size` layers from `start` took `produced` objects in all."""
        for position in range(start, start + size):
            self.status[position] = bool(produced)
        self.group_layers = self.first_group if produced else self.group_layers * 2

    # -- decisions ---------------------------------------------------------

    def _trailing(self, position):
        """The layers just before `position`, back to the last one that took (or is certain to take) something:
        `(empty layers in a row, whether any of them has yet to report)`."""
        run, waiting = 0, False
        for earlier in range(position - 1, -1, -1):
            status = self.status[earlier]
            if status is None:
                if self.expected[earlier] >= self.certain:
                    break
                waiting = True
            elif status:
                break
            run += 1
        return run, waiting

    def must_wait(self, position):
        """Whether the layer at `position` cannot be chosen yet: results to come could still put the walk in a group."""
        if not self.empty_before:
            return False
        run, waiting = self._trailing(position)
        return waiting and run >= self.empty_before

    def grouped(self, position):
        """Whether the walk is in a group at `position` (call after `must_wait` says no)."""
        if not self.empty_before:
            return False
        run, waiting = self._trailing(position)
        return run >= self.empty_before and not waiting

    def group_at(self, position):
        """How many layers the group starting at `position` takes."""
        return min(self.group_layers, self.total - position)

    def settled(self, start, size):
        """Whether every layer from `start` for `size` has reported."""
        return all(status is not None for status in self.status[start:start + size])

    def drive(self, queue, submit_single, submit_group):
        """
        Walks every position: `submit_single(position)` for a layer on its own,
        `submit_group(start, size)` for a group; the work queue's `wait_until`
        is what lets earlier results arrive when the next choice depends on them.
        """
        position = 0
        while position < self.total:
            if self.must_wait(position):
                queue.wait_until(lambda: not self.must_wait(position))
                if self.must_wait(position):
                    raise RuntimeError(f"the scatter walk is waiting on layers that were never started (position {position})")
                continue
            if self.grouped(position):
                size = self.group_at(position)
                self.group_sizes.append(size)
                submit_group(position, size)
                queue.wait_until(lambda: self.settled(position, size))
                position += size
            else:
                self.singles += 1
                submit_single(position)
                position += 1

    def left_ungrouped(self):
        return self.total - sum(self.group_sizes)
