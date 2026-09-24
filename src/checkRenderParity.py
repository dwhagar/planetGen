#!/usr/bin/env python3
# src/checkRenderParity.py

"""
Pre-upgrade check for schema v29, which drops the stored page text
(`star_systems.wikitext_content`/`markdown_content`) and renders both
formats from the database rows instead (see `stellarObjects/schema.sql`'s
"v29" header note).

Run this against a database that is still at v28 -- before `update.sh`/
`migrateDb.py` applies v29, since that step deletes the stored copies it
compares against. For every system it renders the page fresh
(`stellarObjects/systemRender.py`) and compares it with the stored copy,
line by line, sorting each difference into one of three buckets:

- names: the only differences are names -- a rename, or a planet/moon
  name made unique after the stored text was written. The fresh page is
  the correct one (it matches the rest of the site).
- live values: lines that describe where things are right now (planetary
  and moon wobble, comet positions and activity), which `updateOrbits.py`
  keeps moving and the stored text froze at generation time.
- other: anything else. These are the ones worth a look before upgrading;
  the first few diffs are printed.

`--export-dir DIR` also writes every system's stored text to
`DIR/<database>/<id>.wiki` and `.md` before anything is dropped -- a
readable backup alongside (not instead of) a `mysqldump`.

Usage:
    python3 src/checkRenderParity.py [--export-dir DIR] [--show N]
                                     [--mysql-host HOST] [--mysql-port PORT]
                                     [--mysql-user USER] [--mysql-password PASSWORD]
                                     [--mysql-database DATABASE]

This tool only reads the database (plus writing export files).
"""

import argparse
import difflib
import os
import re
import sys

from stellarObjects._db import (
    add_mysql_connection_args, get_connection, load_star_system, mysql_config_from_args,
)
from stellarObjects.names import COMPANION_SUFFIXES, DIMINUTIVE_PREFIXES
from stellarObjects.systemRender import render_star_system

_LIVE_LINE = re.compile(r"wobble|comet|parabolic|perihelion|apoapsis|aphelion", re.IGNORECASE)
"""Lines describing current positions (`updateOrbits.py` moves them)."""


def _name_variants(name):
    """`name` plus the forms a stored copy may still carry: without a
    companion suffix (planets/moons) or a diminutive prefix (systems)."""
    variants = {name}
    words = name.split(" ")
    if len(words) > 1 and words[-1] in COMPANION_SUFFIXES:
        variants.add(" ".join(words[:-1]))
    if len(words) > 1 and words[0] in DIMINUTIVE_PREFIXES:
        variants.add(" ".join(words[1:]))
    return variants


def _names_in(system, stored_text):
    """Every name that may legitimately differ between the two copies."""
    names = set()
    bodies = list(system.stars) + [system.star]
    for body in system.planets + system.secondary_planets:
        bodies.append(body)
        bodies.extend(getattr(body, "moons", []))
    bodies.extend(system.comets + system.secondary_comets)
    for body in bodies:
        names |= _name_variants(getattr(body, "name", "") or "")
    # The stored copy's own title line holds the system's old name.
    first_line = stored_text.split("\n", 1)[0].strip("=# ").strip()
    if first_line:
        names.add(first_line)
        names.add(f"{first_line} B")
    return sorted((n for n in names if n), key=len, reverse=True)


def _mask_names(line, names):
    for name in names:
        line = line.replace(name, "\x00")
    return line


def classify(stored, rendered, names):
    """Returns `None` when identical, else `"names"`, `"live"` or `"other"`."""
    if stored == rendered:
        return None
    stored_lines, rendered_lines = stored.splitlines(), rendered.splitlines()
    worst = "names"
    matcher = difflib.SequenceMatcher(None, stored_lines, rendered_lines, autojunk=False)
    for tag, a0, a1, b0, b1 in matcher.get_opcodes():
        if tag == "equal":
            continue
        old, new = stored_lines[a0:a1], rendered_lines[b0:b1]
        pairs = zip(old, new) if len(old) == len(new) else [(line, line) for line in old + new]
        for o, n in pairs:
            if o != n and _mask_names(o, names) == _mask_names(n, names):
                continue
            if _LIVE_LINE.search(o) and _LIVE_LINE.search(n):
                worst = "live"
                continue
            return "other"
    return worst


def main():
    parser = argparse.ArgumentParser(
        description="Compare each system's stored wikitext/Markdown with a fresh render, before schema v29 drops it.",
    )
    parser.add_argument("--export-dir", help="Also write every system's stored text to this directory.")
    parser.add_argument("--show", type=int, default=5, help="How many 'other' diffs to print (default 5).")
    add_mysql_connection_args(parser)
    args = parser.parse_args()
    config = mysql_config_from_args(args)

    conn = get_connection(config, ensure_schema=False)
    try:
        columns = {r["Field"] for r in conn.execute("SHOW COLUMNS FROM star_systems").fetchall()}
        if "markdown_content" not in columns:
            print(f"{config.database}: the stored page text is already gone (schema v29 or later) -- nothing to compare.")
            return

        export_dir = os.path.join(args.export_dir, config.database) if args.export_dir else None
        if export_dir:
            os.makedirs(export_dir, exist_ok=True)

        counts = {"identical": 0, "names": 0, "live": 0, "other": 0, "no stored text": 0, "failed": 0}
        shown = 0
        ids = [r["id"] for r in conn.execute("SELECT id FROM star_systems ORDER BY id").fetchall()]
        for system_id in ids:
            row = conn.execute(
                "SELECT name, wikitext_content, markdown_content FROM star_systems WHERE id = ?", (system_id,)
            ).fetchone()
            if export_dir:
                for ext, text in (("wiki", row["wikitext_content"]), ("md", row["markdown_content"])):
                    if text is not None:
                        with open(os.path.join(export_dir, f"{system_id}.{ext}"), "w", encoding="utf-8") as f:
                            f.write(text)
            if row["wikitext_content"] is None and row["markdown_content"] is None:
                counts["no stored text"] += 1
                continue
            try:
                system = load_star_system(conn, system_id)
            except Exception as exc:  # report and keep going -- this is a survey
                counts["failed"] += 1
                print(f"system {system_id} ({row['name']}): could not be loaded: {exc}")
                continue

            worst = None
            for fmt, stored in (("wikitext", row["wikitext_content"]), ("markdown", row["markdown_content"])):
                if stored is None:
                    continue
                rendered = render_star_system(system, fmt)
                verdict = classify(stored, rendered, _names_in(system, stored))
                if verdict == "other" and shown < args.show:
                    shown += 1
                    print(f"--- system {system_id} ({row['name']}), {fmt}: differs beyond names and live values")
                    diff = difflib.unified_diff(stored.splitlines(), rendered.splitlines(),
                                                "stored", "rendered", lineterm="", n=0)
                    for line in list(diff)[:20]:
                        print(f"    {line[:200]}")
                order = (None, "names", "live", "other")
                if order.index(verdict) > order.index(worst):
                    worst = verdict
            counts[worst or "identical"] += 1
    finally:
        conn.close()

    print(f"\n{config.database}: {len(ids)} systems checked")
    for label, count in counts.items():
        if count or label in ("identical", "other"):
            print(f"  {label:>15}: {count}")
    if export_dir:
        print(f"Stored text written to {export_dir}")
    if counts["other"] or counts["failed"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
