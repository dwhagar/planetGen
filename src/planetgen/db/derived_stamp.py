# planetgen/db/derived_stamp.py

"""
Which stored numbers a row holds were worked out by the code that is running (GEN.1479, audit C6 and F6)
=========================================================================================================

Some stored columns are derived from other stored columns by a formula (a planet's habitability score from its conditions). When the
formula changes the old rows are stale, and the jobs that refresh them (`rescore-habitability`) used to recompute every row, whether or
not anything about it had changed. A **stamp** on the row says which formula wrote it:

* a family of derived columns has a stamp column on `star_systems` (`STAMP_COLUMNS`), an unsigned 32-bit number;
* `current(family)` is the stamp of the code that is running, a hash of the source of the modules that decide the numbers, so a changed
  formula is a different stamp and nobody has to remember to bump anything (`RULE_VERSION` is for the code outside those modules that
  reads or writes the family: a change to `store.refresh_habitability` itself; for the ledger, `store.refresh_effect_ledger`);
* whatever writes the family's columns for a system writes the stamp with them (the system save, `store.refresh_habitability`);
* a refresh job takes the rows where the column is NULL or not `current(family)` and leaves the others alone, so a rerun after no change
  to the formula touches no row, a run that was stopped carries on where it stopped, and a freshly generated system is never done twice.

A stamp stands for the formula, not for the data it read. Whatever changes a score's inputs after the save (a cloud's dose, a life
update) re-scores the body itself, as before; a job that must redo everything (`--all`) ignores the stamp.

The stamp of a family is a hash of the *physics* modules its root module imports, however indirectly (parsed from the source, so the
lazy imports inside functions count and nothing needs importing). That is wider than the code that really runs for a body, which can only
make a stamp change when it need not (a row is redone that did not have to be), never miss a change. The modules outside `planetgen.physics`
(the generators, the galaxy, the datasets) are left out: they decide the conditions a score is worked out from, not the score, and a
changed condition is a regenerated body, which is stamped by its own save.
"""

import functools
import hashlib
import importlib.util
import os

RULE_VERSION = {"habitability": 1, "ledger": 1}
"""dict: Per family, a number to bump when the code that reads or writes the family outside its roots' closure changes (for habitability,
`store.refresh_habitability` and the way a system's bodies are scored)."""

ROOTS = {"habitability": ("planetgen.physics.habitability_world", "planetgen.physics.insolation", "planetgen.physics.orbit_elements"),
         "ledger": ("planetgen.physics.effect_ledger", "planetgen.physics.light_ledger", "planetgen.physics.radiation_field")}
"""dict: Per family, the modules whose imports (and theirs) make up the source its stamp hashes."""

STAMP_COLUMNS = {"habitability": "habitability_stamp", "ledger": "ledger_stamp"}
"""dict: Per family, its stamp column on `system_summary`."""

EXTRA_MODULES = ("planetgen.tuning",)
"""tuple: Modules hashed with every family: the tuning knobs the physics reads."""

_PREFIXES = ("planetgen.physics.",)
"""tuple: A module is followed when its name starts with one of these (or is a root or an extra)."""


@functools.lru_cache(maxsize=None)
def _file_of(name):
    """The source file of the module `name`, `None` when it is not a module (a name imported from one, say) or has no source file."""
    try:
        spec = importlib.util.find_spec(name)
    except (ImportError, ValueError):
        return None
    origin = None if spec is None else spec.origin
    return origin if origin and os.path.isfile(origin) else None


def _imports_of(source):
    """The names of the `planetgen` modules a source text imports (every `import` and `from ... import`, nested ones too)."""
    import ast
    names = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
            names.add(node.module)
            names.update(f"{node.module}.{alias.name}" for alias in node.names)
    return {name for name in names if name.startswith("planetgen.")}


@functools.lru_cache(maxsize=None)
def sources(family):
    """The module names the stamp of `family` hashes, sorted: its roots, the extras and every physics module they import."""
    found = set(ROOTS[family]) | set(EXTRA_MODULES)
    todo = list(found)
    while todo:
        path = _file_of(todo.pop())
        if path is None:
            continue
        with open(path, encoding="utf-8") as handle:
            for name in _imports_of(handle.read()):
                if name not in found and name.startswith(_PREFIXES) and _file_of(name):
                    found.add(name)
                    todo.append(name)
    return tuple(sorted(found))


@functools.lru_cache(maxsize=None)
def current(family):
    """The stamp of the running code for `family`: an unsigned 32-bit number (the first four bytes of a SHA-256 over the rule version and
    the source of `sources(family)`)."""
    digest = hashlib.sha256(f"{family} rule {RULE_VERSION[family]}".encode())
    for name in sources(family):
        with open(_file_of(name), "rb") as handle:
            digest.update(name.encode() + b"\0" + handle.read())
    return int.from_bytes(digest.digest()[:4], "big")


def stale_condition(family, alias=""):
    """`(sql, parameters)`: the WHERE condition for the rows of `star_systems` (under table `alias`) whose `family` is not current."""
    column = f"{alias}.{STAMP_COLUMNS[family]}" if alias else STAMP_COLUMNS[family]
    return f"({column} IS NULL OR {column} <> ?)", (current(family),)
