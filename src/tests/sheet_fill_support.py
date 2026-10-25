# tests/sheet_fill_support.py

"""
The fill audit of the property sheet (UX.208): generate sample systems into a scratch database, read each body's sheet
(`describe.sheet.sheet`) and count, for every kind, class and sheet column, how often the column is filled and how often
it is empty. `scripts/sheet_fill_report.py` prints and writes the table; `tests/test_sheet_fill.py` checks the table
against the registry `describe.sheet.APPLIES`.

Systems are made in the process with a fixed seed per system (`draw.set_run_seed`), so a run is repeatable.
"""

import collections
import os
import re

from planetgen.db import store
from planetgen.db.store import MySQLConfig
from planetgen.describe import sheet
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

CLASS_COLUMNS = {"planet": ("planets", "planet_class"), "moon": ("moons", "planet_class"), "star": ("stars", "star_type"),
                 "belt": ("asteroid_belts", None), "comet": ("comets", None)}
"""dict: For each kind of a system's bodies, its table and the column that holds its class (None: one class)."""

BODY_KINDS = tuple(CLASS_COLUMNS)


def server_kwargs():
    """The connection arguments of the test MySQL server (the same variables the tests read)."""
    return dict(
        host=os.environ.get("PLANETGEN_TEST_MYSQL_HOST", "127.0.0.1"),
        port=int(os.environ.get("PLANETGEN_TEST_MYSQL_PORT", 3306)),
        user=os.environ.get("PLANETGEN_TEST_MYSQL_USER", "planetgen"),
        password=os.environ.get("PLANETGEN_TEST_MYSQL_PASSWORD", ""),
    )


def generate(config, systems, first_seed=1000):
    """Generates and saves `systems` systems into the database `config` points at, system `i` from the run seed
    `first_seed + i`. Returns the ids of the saved systems."""
    saved = []
    for index in range(systems):
        draw.set_run_seed(first_seed + index)
        system_config = SystemConfig()
        system_config.PLANETS = True
        system_config.MAX_PLANETS = True
        saved.append(store.save_system(StarSystem(system_config=system_config), system_config, config=config))
    return saved


BARE_NUMBER = re.compile(r"^[-+]?[0-9][0-9,]*\.?[0-9]*([eE][-+]?[0-9]+)?$")
"""re.Pattern: A row text that is only a number: no unit, noun or label says what it is."""


def class_of(kind, raw):
    """The class a body is counted under: the planet class letter, the spectral letter of a star, or the kind itself."""
    if kind == "star":
        return (raw or "?")[:1]
    return raw if raw else kind


def tally(config):
    """
    Counts every sheet column of every body in the database.

    Returns:
        dict: `{(kind, class): {"bodies": n, "columns": {key: [filled, empty where it applies, empty where it does not]}}}`
            where a column is empty when its sheet row was never stored (`stored` False) and applies per
            `describe.sheet.applies` (which reads the body's own stored columns for a conditional rule); the key
            `"bare"` holds a Counter of `(kind, key, label)` for stored values whose text is only a number.
    """
    result = {}
    bare = collections.Counter()
    conn = store.get_connection(config)
    try:
        for kind, (table, class_column) in CLASS_COLUMNS.items():
            for row in conn.execute(f"SELECT * FROM {table}").fetchall():
                row = dict(row)
                klass = class_of(kind, row.get(class_column) if class_column else None)
                entry = result.setdefault((kind, klass), {"bodies": 0, "columns": collections.OrderedDict()})
                entry["bodies"] += 1
                for group in sheet.sheet(conn, kind, row["id"])["groups"]:
                    for item in group["rows"]:
                        if item["key"].startswith("p_"):    # the gases are one row per gas the air holds
                            continue
                        counts = entry["columns"].setdefault(item["key"], [0, 0, 0])
                        if item["stored"]:
                            counts[0] += 1
                            if BARE_NUMBER.match(item["text"].strip()):
                                bare[(kind, item["key"], item["label"])] += 1
                        else:
                            counts[1 if sheet.applies(kind, klass, item["key"], row) else 2] += 1
    finally:
        conn.close()
    result["bare"] = bare
    return result


def lines(result, target=30):
    """The report as text: per kind and class the number of samples (marked when under `target`), then each column that
    is empty in some sample, with how many were empty where it applies and where it does not, then the values printed as
    a bare number."""
    out = []
    for (kind, klass), entry in sorted((k, v) for k, v in result.items() if k != "bare"):
        short = "" if entry["bodies"] >= target else f"   (fewer than {target} samples)"
        out.append(f"{kind} {klass}: {entry['bodies']} bodies{short}")
        for key, (filled, gap, not_applicable) in entry["columns"].items():
            if gap or not_applicable:
                out.append(f"    {key}: {filled} filled, {gap} empty where it applies, {not_applicable} empty where it does not")
    out.append("")
    out.append("Values printed as a bare number (no unit, noun or label), by kind, column and label:")
    for (kind, key, label), count in sorted(result.get("bare", {}).items()):
        out.append(f"    {kind} {key} ({label}): {count} bodies")
    return out
