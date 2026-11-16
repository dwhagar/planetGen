# tests/pair_overlap_support.py

"""
The proof of one pair (OPS.242): two stages overlapped store the same bytes as the same two stages one after the other.

A pair goes in `queue/independence.json` only when `prove(recipe, context)` says `passed`. What it does, from the one prepared database each time (`cut_neutral_support.dump` and `restore`):

* **Capture.** The recipe's `capture` runs the real command line once (the old queue, inline), and a hook on `commands.expand_function` and `commands.end_function` records what each stage was
  cut with and what each of its end commands was given (`Captured`). The database is then put back. Every run below replays those, so the params are the ones a run really gives, settings that
  are worked out from the stored galaxy included.
* **The runs.** Plan order (the first stage, then the second, each cut, its units one after another, then its end commands), plan order again (it must reproduce itself), the second stage first, and
  the stages overlapped: both cut up front, before either has stored anything (an overlapped stage opens before the other ends), their units mixed one at a time (`workers` 1) or on a pool
  of 15 threads, in both orders of who starts first, the end commands in the order of the start. A unit is `engine.run`, the same as a worker's.
* **The verdict.** Every run is dumped and compared with plan order, table by table, row by row (`cut_neutral_support._comparable`: the clock columns left out, the recipe's `ignore_columns`
  and `views` for a row id that no output depends on). Besides: each stage must store something in plan order (a pair that stores nothing proves nothing), and what each command read
  and wrote, seen by the SQL tap (`queue/sqltap.py`) on the plan-order run, must be what `queue/touches.py` declares for it, on data that has the stage doing its work.

A pair that fails names the rows in `Verdict.differences`; its entry is taken out of `independence.json` until it is fixed.
"""

import concurrent.futures
import contextlib
import json
from dataclasses import dataclass, field

from planetgen.queue import commands, engine, planrows, sqltap, touches
from tests import cut_neutral_support as support

WORKERS = (1, 15)
"""tuple: How many workers the overlapped runs mix the units of the two stages on: one (the units one at a time, alternating) and fifteen (a pool of threads, many units at once)."""


@dataclass(frozen=True)
class PairRecipe:
    """
    How to prove one pair.

    Attributes:
        first, second (str): The two commands' kinds (`commands.REGISTRY`), in plan order.
        prepare (callable): `prepare(context)` puts the database in the state both stages work on.
        capture (callable): `capture(context)` runs the real command lines that make both stages (`Capture`); the database is put back after it.
        ignore_columns (tuple): `(table, column)` pairs left out of the comparison (a surrogate key nothing refers to), each with a reason in `note`.
        views (tuple or callable): `(name, sql)` pairs compared as well (a row's parent seen through its join), or `views(context)` that gives them.
        note (str): Why the recipe is shaped as it is.
    """

    first: str
    second: str
    prepare: object
    capture: object
    ignore_columns: tuple = ()
    views: tuple = ()
    note: str = ""


@dataclass
class Captured:
    """What the real run gave each stage: `expanded[kind]` the params its `expand` was called with, `ended[(kind, role)]` the params and the detail its end command was given."""

    expanded: dict = field(default_factory=dict)
    ended: dict = field(default_factory=dict)


def _without_run(value):
    """`value` with every `run_id` set to `None`: the run record of the capture is gone when the database is put back, and a row that refers to it would break its foreign key."""
    if isinstance(value, dict):
        return {key: None if key == "run_id" else _without_run(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_without_run(item) for item in value]
    return value


@contextlib.contextmanager
def capturing(monkeypatch):
    """Records, for the stages that run in the block, the params their `expand` and their end commands get (`Captured`)."""
    found = Captured()
    real_expand, real_end = commands.expand_function, commands.end_function

    def expand(kind):
        function = real_expand(kind)

        def cut(params):
            found.expanded[kind] = _without_run(params)
            return function(params)
        return cut

    def end(kind, role):
        function = real_end(kind, role)
        if function is None:
            return function

        def hook(params, detail, results, seconds):
            found.ended[(kind, role)] = (_without_run(params), detail)
            return function(params, detail, results, seconds)
        return hook

    monkeypatch.setattr(commands, "expand_function", expand)
    monkeypatch.setattr(commands, "end_function", end)
    try:
        yield found
    finally:
        monkeypatch.setattr(commands, "expand_function", real_expand)
        monkeypatch.setattr(commands, "end_function", real_end)


class _Stage:
    """One stage being replayed: cut, its units run, its end commands called."""

    def __init__(self, kind, captured, tap):
        self.kind, self.captured, self.tap = kind, captured, tap
        self.params = captured.expanded[kind]
        self.expansion = None
        self.results = []

    def cut(self):
        with self.tap.running(self.kind):
            self.expansion = commands.expand_function(self.kind)(self.params)
        self.results = [None] * len(self.expansion.units)

    def units(self):
        return [(self, index, params) for index, (params, _cost) in enumerate(self.expansion.units)]

    def run_unit(self, index, params):
        with self.tap.running(self.kind):
            self.results[index] = json.loads(json.dumps(engine.run(self.kind, params, lambda *args, **kwargs: None)))

    def end(self):
        for role in (planrows.CONCLUDE, planrows.FINISH):
            hook = commands.end_function(self.kind, role)
            if hook is None:
                continue
            params, _detail = self.captured.ended.get((self.kind, role), (self.params, None))
            with self.tap.running(commands.end_name(self.kind, role)):
                hook(params, self.expansion.plan, self.results, 0.0)


def _run(captured, kinds, tap, *, overlapped, workers=1):
    """
    One run of the two stages in `kinds` (the one that starts first first): one after another, or overlapped (both cut before either runs a unit, the units mixed, the end commands
    in the order of the start).
    """
    stages = [_Stage(kind, captured, tap) for kind in kinds]
    if not overlapped:
        for stage in stages:
            stage.cut()
            for _stage, index, params in stage.units():
                stage.run_unit(index, params)
            stage.end()
        return
    for stage in stages:
        stage.cut()
    mixed = []
    lists = [stage.units() for stage in stages]
    for position in range(max(len(units) for units in lists)):
        for units in lists:
            if position < len(units):
                mixed.append(units[position])
    if workers == 1:
        for stage, index, params in mixed:
            stage.run_unit(index, params)
    else:
        with concurrent.futures.ThreadPoolExecutor(workers) as pool:
            for future in [pool.submit(stage.run_unit, index, params) for stage, index, params in mixed]:
                future.result()
    for stage in stages:
        stage.end()


@dataclass
class Verdict:
    """What `prove` found for one pair."""

    first: str
    second: str
    runs: list
    differences: list
    undeclared: list
    stored: dict

    @property
    def passed(self):
        return not self.differences and not self.undeclared and all(self.stored.values())

    def line(self):
        words = "INDEPENDENT" if self.passed else "FAILS"
        why = [text for text, bad in (("the runs differ in the stored bytes", bool(self.differences)),
                                      ("a command touches what it does not declare", bool(self.undeclared)),
                                      ("a stage stored nothing, so it proves nothing: " + ", ".join(kind for kind, rows in self.stored.items() if not rows), not all(self.stored.values()))) if bad]
        return f"{words} {self.first} and {self.second}: {len(self.runs)} runs compared with plan order" + (f" ({'; '.join(why)})" if why else "")


def _comparable(config, recipe, context):
    return {**support._comparable(support.dump(config), set(recipe.ignore_columns)), **_views(config, recipe, context)}


def _views(config, recipe, context):
    class Holder:
        views = recipe.views(context) if callable(recipe.views) else recipe.views
    return support._view_rows(config, Holder)


def _undeclared(tap, kinds):
    """Lines for what the tap saw a command touch that `touches` does not declare for it (reads are covered by the writes, the engine's own tables are left out)."""
    lines = []
    for kind in kinds:
        slot = tap.touched.get(kind)
        if slot is None:
            continue
        declared = {name: set(names) for name, names in (("reads", touches.reads_of(kind)), ("writes", touches.writes_of(kind)))}
        extra_reads = slot["reads"] - declared["reads"] - declared["writes"] - touches.ENGINE_TABLES if touches.ANY not in touches.reads_of(kind) else set()
        extra_writes = slot["writes"] - declared["writes"] - touches.ENGINE_TABLES if touches.ANY not in touches.writes_of(kind) else set()
        if extra_reads:
            lines.append(f"{kind} reads {sorted(extra_reads)}, not declared")
        if extra_writes:
            lines.append(f"{kind} writes {sorted(extra_writes)}, not declared")
    return lines


def prove(recipe, context):
    """The verdict on the pair of `recipe` (module note)."""
    recipe.prepare(context)
    base = support.dump(context.mysql_config)
    with capturing(context.monkeypatch) as captured:
        recipe.capture(context)
    support.restore(context.mysql_config, base)
    first, second = recipe.first, recipe.second
    ends = [commands.end_name(kind, role) for kind in (first, second) for role in (planrows.CONCLUDE, planrows.FINISH) if commands.end_function(kind, role)]
    kinds = [first, second, *ends]
    tap = sqltap.Tap()
    runs = [("plan order", (first, second), False, 1), ("plan order again", (first, second), False, 1), ("second stage first", (second, first), False, 1)]
    runs += [(f"overlapped, {order[0]} first, {workers} worker(s)", order, True, workers)
             for workers in WORKERS for order in ((first, second), (second, first))]
    found = {}
    with tap:
        for name, order, overlapped, workers in runs:
            support.restore(context.mysql_config, base)
            _run(captured, order, tap, overlapped=overlapped, workers=workers)
            found[name] = _comparable(context.mysql_config, recipe, context)
    differences = []
    for name, rows in found.items():
        if name == "plan order":
            continue
        differences += [f"{name}: {line}" for line in support._differences(found["plan order"], rows)]
    stored = _stored_by_stage(context, recipe, captured, base, first, second)
    return Verdict(first, second, [name for name, *_ in runs], differences, _undeclared(tap, kinds), stored)


def _comparable_of(dumped, recipe):
    return support._comparable(dumped, set(recipe.ignore_columns))


def _stored_by_stage(context, recipe, captured, base, first, second):
    """What each stage alone stores, from the prepared database: `{kind: [tables that changed]}` (a stage that stores nothing proves nothing)."""
    out = {}
    before = _comparable_of(base, recipe)
    for kind in (first, second):
        support.restore(context.mysql_config, base)
        _run(captured, (kind,), sqltap.Tap(), overlapped=False)
        after = support._comparable(support.dump(context.mysql_config), set(recipe.ignore_columns))
        out[kind] = [table for table, rows in after.items() if rows != before.get(table)]
    return out
