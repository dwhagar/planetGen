#!/usr/bin/env python3
# html/index.py

"""
Landing page: lists every `.db` file found in the database directory
(`db/` alongside `html/` by default -- see `lib/dbutil.py`) and links to
`browse.py` for each one.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from dbutil import esc, list_databases, open_readonly, resolve_db_path
from page import run


def _counts(path):
    """Returns (sector_count, system_count) for a quick-glance summary."""
    conn = open_readonly(path)
    try:
        sectors = conn.execute("SELECT COUNT(*) FROM sectors").fetchone()[0]
        systems = conn.execute("SELECT COUNT(*) FROM star_systems").fetchone()[0]
        return sectors, systems
    finally:
        conn.close()


def handler():
    databases = list_databases()

    if not databases:
        body = """
<section class="panel">
<p>No databases found. Generate one first with
<code>sectorGen.py</code> or <code>systemGen.py</code>
(see the project README).</p>
</section>
"""
        return "planetGen Databases", body

    rows = []
    for entry in databases:
        sectors, systems = _counts(resolve_db_path(entry["name"]))
        rows.append(
            "<tr>"
            f'<td><a href="browse.py?db={esc(entry["name"])}">{esc(entry["name"])}</a></td>'
            f'<td>{sectors}</td>'
            f'<td>{systems}</td>'
            f'<td>{entry["size_bytes"]:,} bytes</td>'
            f'<td>{esc(entry["modified_at"])}</td>'
            "</tr>"
        )

    count_badge = f"{len(databases)} database{'s' if len(databases) != 1 else ''} found"
    body = f"""
<p class="badges"><span class="badge">{count_badge}</span></p>
<section class="panel">
<h2>Databases</h2>
<div class="table-scroll"><table>
  <thead>
    <tr><th>Database</th><th>Sectors</th><th>Systems</th><th>Size</th><th>Last modified</th></tr>
  </thead>
  <tbody>
    {''.join(rows)}
  </tbody>
</table></div>
</section>
"""
    return "planetGen Databases", body


run(handler)
