# planetgen/queue/breakpoints.py

"""
What a job declares when it enters the queue (OPS.161, docs/design/operations/dynamic-queue-chain.md): where it may be cut, whether the engine may choose the cut, and what it
depends on. The queue engine (`queue/splitter.py`) reads this with the job profile (`queue/splitprofile.py`) and the machine's free resources and picks the split; a job
declares nothing about the machine.

A declaration is plain data and has no side effects. It says, for the `items` the job will work on, numbered 0 to `items - 1` in the order the job works them:

* **the suggested split** (`suggested`): the cut the job would make today, as the positions it is cut at (a unit holds the items between two cuts). It is the cut when the flag
  says so and the cut whenever the engine has nothing to go on. Boss, 2026-11-12, 3:02 PM EST: "use the suggested split in the job information".
* **the split flag** (`split`; for a registered command the flag of its `CommandSpec.split`, which `commands.declaration` puts here): `BYTE_IDENTICAL` says the job's output is the same at any legal cut, so the engine may choose the cut from the profile and the free resources;
  `SUGGESTED` says the output depends on the cut, so the engine uses the suggested cut exactly. A job that declares nothing is `SUGGESTED`. PERF.260 holds a job flagged
  `BYTE_IDENTICAL` to a test that gives the same output at two cuts.
* **where a cut may fall** (`boundaries`, `whole`): `boundaries` are the item positions a cut may be at (`None` is any position), `whole` are spans of items `(start, end)` that
  must stay in one unit, so no cut falls strictly inside. `legal_cuts` is what is left.
* **its dependencies** (`after`, `reads`, `writes`): `after` are the kinds of command that must be done first (the same edges as `commands.REQUIRES`, which a command that
  declares none takes); `reads` and `writes` are named facts, tables or files the job uses and changes. From them `conflicts` says whether two jobs touch the same thing, and so
  whether they may overlap (OPS.131/OPS.132: a pair overlaps only when both declared and nothing conflicts; the byte-for-byte proof of a pair stays the rule of those items).
* **the clamp** (`floor_items`, `ceiling_items`): the fewest and most items the engine may put in a unit. Unset, a unit may hold a quarter to four times the suggested unit
  (`CLAMP_FACTOR`), so the typed split stays the default and the bound.
"""

import statistics
from dataclasses import dataclass

from planetgen.queue.commands import SPLIT_BYTE_IDENTICAL as BYTE_IDENTICAL, SPLIT_SUGGESTED as SUGGESTED

FLAGS = (BYTE_IDENTICAL, SUGGESTED)

CLAMP_FACTOR = 4
"""int: How far from the suggested unit size the engine may go when the job gives no floor and no ceiling."""


class DeclarationError(ValueError):
    """A declaration that contradicts itself."""


def _cuts(values, items, what):
    out = tuple(int(value) for value in values)
    if list(out) != sorted(set(out)) or any(not 0 < value < items for value in out):
        raise DeclarationError(f"{what} must be increasing positions between 1 and {items - 1}, not {list(out)[:8]}")
    return out


@dataclass(frozen=True)
class Declaration:
    """
    One job's declaration (the module docstring). All positions count items from 0; a cut at `c` ends a unit after item `c - 1`. A declaration of no items has no cut.
    """

    items: int
    suggested: tuple = ()
    split: str = SUGGESTED
    boundaries: tuple = None
    whole: tuple = ()
    after: tuple = ()
    reads: tuple = ()
    writes: tuple = ()
    floor_items: int = None
    ceiling_items: int = None

    def __post_init__(self):
        if self.items < 0:
            raise DeclarationError("a job has no fewer than 0 items")
        if self.split not in FLAGS:
            raise DeclarationError(f"the split flag is one of {FLAGS}, not {self.split!r}")
        object.__setattr__(self, "suggested", _cuts(self.suggested, self.items, "the suggested split"))
        if self.boundaries is not None:
            object.__setattr__(self, "boundaries", _cuts(self.boundaries, self.items, "the boundaries"))
        spans = tuple((int(start), int(end)) for start, end in self.whole)
        if any(not 0 <= start < end <= self.items for start, end in spans):
            raise DeclarationError(f"a span that stays whole lies inside the {self.items} items, not {list(spans)[:4]}")
        object.__setattr__(self, "whole", spans)
        for name in ("after", "reads", "writes"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        for name in ("floor_items", "ceiling_items"):
            value = getattr(self, name)
            if value is not None and value < 1:
                raise DeclarationError(f"{name} is at least 1 item")
        if self.floor_items and self.ceiling_items and self.floor_items > self.ceiling_items:
            raise DeclarationError("the floor is above the ceiling")
        illegal = [cut for cut in self.suggested if cut not in self.legal_cuts]
        if illegal:
            raise DeclarationError(f"the suggested split cuts at {illegal[:4]}, where a cut may not fall")

    @property
    def legal_cuts(self):
        """The positions a cut may be at, in order: the boundaries (every position when there are none), less any inside a span that stays whole."""
        return _legal(self)

    @property
    def dynamic(self):
        """Whether the engine may choose the cut."""
        return self.split == BYTE_IDENTICAL

    def sizes(self, cuts=None):
        """The number of items of each unit of `cuts` (the suggested split by default)."""
        edges = [0, *(self.suggested if cuts is None else cuts), self.items]
        return [end - start for start, end in zip(edges, edges[1:])] if self.items else []

    def spans(self, cuts=None):
        """The `(start, end)` item range of each unit of `cuts` (the suggested split by default)."""
        edges = [0, *(self.suggested if cuts is None else cuts), self.items]
        return list(zip(edges, edges[1:])) if self.items else []

    def suggested_unit_items(self):
        """The usual number of items of a suggested unit (the median), or 0 for a job with no item."""
        sizes = self.sizes()
        return int(statistics.median(sizes)) if sizes else 0

    def clamp(self):
        """`(floor, ceiling)` in items the engine keeps a unit between (the module docstring)."""
        usual = max(self.suggested_unit_items(), 1)
        floor = self.floor_items or max(1, usual // CLAMP_FACTOR)
        ceiling = self.ceiling_items or usual * CLAMP_FACTOR
        return floor, max(floor, ceiling)


def _legal(decl):
    allowed = range(1, decl.items) if decl.boundaries is None else decl.boundaries
    return tuple(cut for cut in allowed if not any(start < cut < end for start, end in decl.whole))


def even(items, unit_items, **options):
    """A declaration whose suggested split is units of `unit_items` items (the last takes what is left); `options` are the other fields."""
    if unit_items < 1:
        raise DeclarationError("a unit holds at least 1 item")
    return Declaration(items, tuple(range(unit_items, items, unit_items)), **options)


def from_chunks(chunks, **options):
    """A declaration whose suggested split is `chunks`, the lists of items a stage already cuts into (`units.static_chunks`), in order."""
    sizes = [len(chunk) for chunk in chunks if len(chunk)]
    cuts, total = [], 0
    for size in sizes[:-1]:
        total += size
        cuts.append(total)
    return Declaration(sum(sizes), tuple(cuts), **options)


def whole_job(**options):
    """The declaration of a job that is one unit by nature (`CommandSpec.whole`): one item, no cut, so the suggested split is the only split."""
    return Declaration(1, (), **options)


def conflicts(first, second):
    """
    Whether two declarations touch the same thing so that the two jobs may not run together: one writes what the other reads or writes. A job that declares neither reads nor
    writes is unknown, not free: it conflicts with everything, so a job that has not declared yet runs in plan order, as it does today.
    """
    if not (first.reads or first.writes) or not (second.reads or second.writes):
        return True
    one, other = set(first.writes), set(second.writes)
    return bool(one & (set(second.reads) | other) or other & set(first.reads))


def may_overlap(first, second):
    """Whether two jobs may run together on the strength of what they declare (no conflict). Whether a pair does is for OPS.131/OPS.132 to prove byte for byte."""
    return not conflicts(first, second)
