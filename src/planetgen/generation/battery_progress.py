# planetgen/generation/battery_progress.py

"""
What the generator battery's digest process says of its progress (OPS.184)
=========================================================================

`generator_battery.start_running_digest` works the digest out in a process of its own while a plan does other work, and the plan then waits for it. On a code key not seen
before the digest takes minutes (the bright-star sampling tables are built cold), so the wait has to show real counts: the process counts the battery's cases and every sampling
table it builds in two shared numbers (made by `create`, handed to the process when it starts, `adopt`), and the waiting side reads them (`read`) as often as it likes. Nothing is
written or counted in any other process: `case` and `table_built` do nothing until `adopt` has been called.
"""

_shared = None
"""tuple or None: `(cases done, cases in all, tables built)`, shared numbers of the digest process; `None` anywhere else."""


def create(context):
    """The shared numbers for a digest process of `context` (a `multiprocessing` context). One process writes them and one reads, so they take no lock."""
    return (context.Value("i", 0, lock=False), context.Value("i", 0, lock=False), context.Value("q", 0, lock=False))


def adopt(shared):
    """Called first in the digest process: the numbers it writes."""
    global _shared
    _shared = shared


def case(done, total):
    """The battery's `progress(done, total)`: `done` of its `total` cases are finished."""
    if _shared is not None:
        _shared[0].value, _shared[1].value = done, total


def table_built():
    """A sampling table was worked out (not read back from a kept copy)."""
    if _shared is not None:
        _shared[2].value += 1


def read(shared):
    """`(cases done, cases in all, tables built)` as the digest process has them now."""
    return shared[0].value, shared[1].value, shared[2].value
