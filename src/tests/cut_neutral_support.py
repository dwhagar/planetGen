# tests/cut_neutral_support.py

"""
The cut-neutral check (PERF.260): a stage whose split is flagged "dynamic split is byte identical" must store the same bytes whatever size the queue cuts its units to.

A command's stage is cut into units by a knob (a module constant such as `run_galaxy.LINK_UNIT_SECTORS`, or a size the planner will one day choose from the job profile).
The queue engine may pick any cut for a stage flagged `commands.SPLIT_BYTE_IDENTICAL`; this module shows whether the flag is true:

* a `Recipe` says, for one command kind, how to put a database in the state the stage works on (`prepare`), what params its stage is expanded with, which knob cuts it and two
  values of it (`cuts`). What the stage stored is everything it changed anywhere in the database (the whole database is dumped and put back, `dump` and `restore`);
* `check(recipe, context)` runs the stage at the first cut, at the second and at the first again, from the same stored state each time, one unit after another in unit order, and
  compares the dumps with `repr` (the exact bits of a float, the clock columns left out): a `Verdict` says `identical`, `reproducible` (the two runs at the same cut agree, which a
  difference between cuts needs to mean anything), `distinct_cuts` (the two cuts really made a different number of units) and that the stage stored something at all, or the
  comparison proved nothing;
* `flagged_commands()` finds the commands that carry the flag, with no list kept by hand, and `gate_failures()` is what the merge gate asserts: every flagged command has a
  recipe and its verdict passes. A flagged command that cannot be run at two cuts fails until its flag is removed.
"""

import contextlib
import importlib
from dataclasses import dataclass, field

from planetgen.queue import commands

SPLIT_FIELD = "split"
"""str: The `CommandSpec` field that holds a command's split flag (OPS.161); a command without it keeps its suggested split."""


@dataclass
class Context:
    """What a recipe works with: the throwaway database, its `monkeypatch` and any extra fixtures the test passes in (`extras`)."""

    mysql_config: object
    monkeypatch: object
    extras: dict = field(default_factory=dict)


@dataclass(frozen=True)
class Recipe:
    """
    How to run one command kind at two cuts.

    What the stage stored is everything it changed anywhere in the database: the whole database is dumped after `prepare` and again after each run (`dump`), the database is put back
    to the first dump before each run (`restore`), and the dumps are compared table by table, so the recipe names no table and a stage that writes somewhere unexpected is still seen.

    Attributes:
        kind (str): The command kind in `commands.REGISTRY`.
        knobs (tuple): `(dotted module, constant name)` for each constant that cuts the stage; `cuts` gives each one's values, in the same order.
        cuts (tuple): `((first cut values), (second cut values))`, the first small (many units) and the second large (few).
        prepare (callable): `prepare(context)` puts the database in the state the stage works on, once.
        run_params (callable): `run_params(context)` is the params the stage's `expand` is given.
        ignore_columns (tuple): `(table, column)` pairs left out of the comparison, each one a value that is not an output of the stage (an id from a counter the run does not reset);
            every pair needs a reason in `note`.
        views (tuple): `(name, sql)` pairs, a query run after each run whose rows are compared as well; a table whose rows refer to a surrogate key left out of `ignore_columns` is compared
            through the query that joins them, so what a row says about its parent is still checked.
        capture (callable, optional): `capture(context)` gives what the stage stored outside the database, `{name: [sorted strings]}` (the cache a warm-up fills), compared as well.
        reset (callable, optional): `reset(context)` empties what `capture` reads, before each run.
        note (str): Why the recipe is shaped as it is (shown in the report).
    """

    kind: str
    knobs: tuple
    cuts: tuple
    prepare: object
    run_params: object
    ignore_columns: tuple = ()
    views: tuple = ()
    capture: object = None
    reset: object = None
    note: str = ""


@dataclass
class Verdict:
    """What `check` found for one recipe."""

    kind: str
    cuts: tuple
    units: tuple
    identical: bool
    reproducible: bool
    distinct_cuts: bool
    differences: list
    changed_tables: list = field(default_factory=list)

    @property
    def passed(self):
        return self.identical and self.reproducible and self.distinct_cuts and bool(self.changed_tables)

    def line(self):
        words = "IDENTICAL" if self.passed else "FAILS"
        why = [] if self.passed else [text for text, bad in (("the cuts differ in the stored bytes", not self.identical),
                                                             ("the same cut does not reproduce itself", not self.reproducible),
                                                             ("the two cuts make the same number of units", not self.distinct_cuts),
                                                             ("the stage stored nothing, so there is nothing to compare", not self.changed_tables)) if bad]
        return f"{words} {self.kind}: cuts {self.cuts[0]} and {self.cuts[1]} make {self.units[0]} and {self.units[1]} units" + (f" ({'; '.join(why)})" if why else "")


def _module(dotted):
    return importlib.import_module(dotted)


@contextlib.contextmanager
def _knobs(context, recipe, values):
    """Sets the recipe's constants to `values` for the run (the test's `monkeypatch` puts them back)."""
    for (dotted, name), value in zip(recipe.knobs, values):
        context.monkeypatch.setattr(_module(dotted), name, value)
    yield


def run_stage(recipe, context, values):
    """
    Runs the stage of `recipe.kind` at the cut `values`, every unit one after another in unit order, and returns how many units it was cut into.

    Raises:
        AssertionError: The unit's report never advanced or a unit failed (the exception propagates).
    """
    with _knobs(context, recipe, values):
        expansion = commands.expand_function(recipe.kind)(recipe.run_params(context))
        unit = commands.unit_function(recipe.kind)
        for params, _cost in expansion.units:
            unit(params, lambda *args, **kwargs: None)
    return len(expansion.units)


CLOCK_COLUMN_SUFFIXES = ("_at",)
"""tuple: A column that ends like this holds the time a row was written, which is no output of a cut."""


def _tables(conn):
    return [row["table_name"] for row in conn.execute(
        "SELECT table_name AS table_name FROM information_schema.tables WHERE table_schema = DATABASE() AND table_type = 'BASE TABLE' ORDER BY table_name").fetchall()]


def _stored_columns(conn, table):
    """The columns of `table` a row can be written with (not the generated ones), in table order."""
    return [row["column_name"] for row in conn.execute(
        "SELECT column_name AS column_name FROM information_schema.columns WHERE table_schema = DATABASE() AND table_name = ? AND extra NOT LIKE ? "
        "ORDER BY ordinal_position", (table, "%GENERATED%")).fetchall()]


def dump(config):
    """Every row of every table of the database `config` names: `{table: (columns, [row tuples])}`."""
    from planetgen.db import store

    conn = store.get_connection(config)
    try:
        result = {}
        for table in _tables(conn):
            columns = _stored_columns(conn, table)
            rows = conn.execute(f"SELECT {', '.join(f'`{name}`' for name in columns)} FROM `{table}`").fetchall()
            result[table] = (columns, [tuple(row[name] for name in columns) for row in rows])
        return result
    finally:
        conn.close()


def restore(config, dumped):
    """Puts the database back to `dumped` (`dump`): every table emptied (the id counters with it) and its rows written again."""
    from planetgen.db import store

    conn = store.get_connection(config)
    try:
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        for table, (columns, rows) in dumped.items():
            conn.execute(f"TRUNCATE TABLE `{table}`")
            if rows:
                marks = ", ".join("?" for _name in columns)
                conn.executemany(f"INSERT INTO `{table}` ({', '.join(f'`{name}`' for name in columns)}) VALUES ({marks})", rows)
        conn.execute("SET FOREIGN_KEY_CHECKS = 1")
        conn.commit()
    finally:
        conn.close()
    store.forget_id_blocks(config._key())      # the ids this process holds were taken from counters that are back where they were


def _view_rows(config, recipe):
    """`{view name: sorted repr of each row}` of the recipe's `views`, run on the database now."""
    from planetgen.db import store

    if not recipe.views:
        return {}
    conn = store.get_connection(config)
    try:
        return {f"view {name}": sorted(repr(tuple(row.values())) for row in conn.execute(sql).fetchall()) for name, sql in recipe.views}
    finally:
        conn.close()


def _comparable(dumped, ignore):
    """`{table: sorted repr of each row}` with the clock columns and the recipe's ignored columns left out (a row's position is no output)."""
    found = {}
    for table, (columns, rows) in dumped.items():
        keep = [index for index, name in enumerate(columns)
                if not name.endswith(CLOCK_COLUMN_SUFFIXES) and (table, name) not in ignore]
        found[table] = sorted(repr(tuple(row[index] for index in keep)) for row in rows)
    return found


def check(recipe, context):
    """The verdict on `recipe`: the stage at the first cut, at the second and at the first again, each from the same stored state, compared (module note)."""
    recipe.prepare(context)
    base = dump(context.mysql_config)
    ignore = set(recipe.ignore_columns)
    runs = []
    for values in (recipe.cuts[0], recipe.cuts[1], recipe.cuts[0]):
        restore(context.mysql_config, base)
        if recipe.reset is not None:
            recipe.reset(context)
        units = run_stage(recipe, context, values)
        found = {**_comparable(dump(context.mysql_config), ignore), **_view_rows(context.mysql_config, recipe)}
        if recipe.capture is not None:
            found.update({f"stored {name}": rows for name, rows in recipe.capture(context).items()})
        runs.append((units, found))
    (units_a, first), (units_b, second), (_units, again) = runs
    changed = [table for table, rows in _comparable(base, ignore).items() if first[table] != rows] + [name for name, rows in first.items() if name.startswith("stored ") and rows]
    differences = _differences(first, second)
    verdict = Verdict(recipe.kind, recipe.cuts, (units_a, units_b), not differences, _differences(first, again) == [], units_a != units_b, differences)
    verdict.changed_tables = changed
    return verdict


def _differences(left, right):
    """What differs between two comparable dumps, as lines (the first few rows of each table that differs)."""
    found = []
    for name in sorted(set(left) | set(right)):
        a, b = left.get(name, []), right.get(name, [])
        if a == b:
            continue
        only_a, only_b = sorted(set(a) - set(b)), sorted(set(b) - set(a))
        found.append(f"{name}: {len(a)} rows against {len(b)}; only in the first {only_a[:3]}; only in the second {only_b[:3]}")
    return found


def flagged_commands(registry=None):
    """The kinds of the commands that carry the byte-identical flag, from the registry (`CommandSpec.split`), in registry order. None before OPS.161 adds the field."""
    registry = commands.REGISTRY if registry is None else registry
    flag = getattr(commands, "SPLIT_BYTE_IDENTICAL", "byte-identical")
    return [kind for kind, spec in registry.items() if getattr(spec, SPLIT_FIELD, None) == flag]


def gate_failures(recipes, context_for, registry=None):
    """
    Every reason the flags of `registry` are not true (the merge gate asserts this is empty): each flagged command must have a recipe in `recipes` (`{kind: Recipe}`), and the recipe's
    verdict (`context_for(recipe)` gives its context) must pass. A flagged command with no recipe cannot be run at two cuts and fails until the flag is removed.
    """
    failures = []
    for kind in flagged_commands(registry):
        recipe = recipes.get(kind)
        if recipe is None:
            failures.append(f"{kind} is flagged byte identical but has no recipe to run it at two cuts (add one to tests/cut_neutral_recipes.py, or remove the flag)")
            continue
        verdict = check(recipe, context_for(recipe))
        if not verdict.passed:
            failures.append(verdict.line() + "".join(f"\n    {line}" for line in verdict.differences))
    return failures
