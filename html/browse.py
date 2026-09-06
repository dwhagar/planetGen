#!/usr/bin/env python3
# html/browse.py

"""
Per-database overview: every sector (with its system count) and every
standalone system (one generated with no sector, `sector_id IS NULL`).
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from dbutil import esc, fetch_all, open_readonly, resolve_db_path
from page import query_params, run


def _star_summary(conn, star_system_id, is_binary):
    """One-line star description for a system row: the single/primary
    star's descriptive type, or a binary's combined type."""
    if is_binary:
        row = conn.execute(
            "SELECT binary_type FROM star_systems WHERE id = ?", (star_system_id,)
        ).fetchone()
        return row["binary_type"] if row else ""
    row = conn.execute(
        "SELECT star_type FROM stars WHERE star_system_id = ? AND role = 'single'",
        (star_system_id,),
    ).fetchone()
    return row["star_type"] if row else ""


def handler():
    params = query_params()
    db_name = params.get("db", "")
    path = resolve_db_path(db_name)
    conn = open_readonly(path)
    try:
        sectors = fetch_all(
            conn,
            """
            SELECT s.id, s.name, s.edge_mpc, COUNT(ss.id) AS system_count
            FROM sectors s
            LEFT JOIN star_systems ss ON ss.sector_id = s.id
            GROUP BY s.id
            ORDER BY s.name
            """,
        )
        standalone = fetch_all(
            conn,
            """
            SELECT id, name, is_binary
            FROM star_systems
            WHERE sector_id IS NULL
            ORDER BY name
            """,
        )

        sector_rows = "".join(
            "<tr>"
            f'<td><a href="sector.py?db={esc(db_name)}&id={row["id"]}">{esc(row["name"])}</a></td>'
            f'<td>{row["system_count"]}</td>'
            "</tr>"
            for row in sectors
        ) or '<tr><td colspan="2"><em>None</em></td></tr>'

        standalone_rows = "".join(
            "<tr>"
            f'<td><a href="system.py?db={esc(db_name)}&id={row["id"]}">{esc(row["name"])}</a></td>'
            f'<td>{_esc_bool(row["is_binary"])}</td>'
            f'<td>{esc(_star_summary(conn, row["id"], row["is_binary"]))}</td>'
            "</tr>"
            for row in standalone
        ) or '<tr><td colspan="3"><em>None</em></td></tr>'
    finally:
        conn.close()

    body = f"""
<p><a href="index.py">&larr; Databases</a></p>
<h2>Sectors</h2>
<table>
  <thead><tr><th>Name</th><th>Systems</th></tr></thead>
  <tbody>{sector_rows}</tbody>
</table>

<h2>Standalone Systems</h2>
<table>
  <thead><tr><th>Name</th><th>Binary</th><th>Star type</th></tr></thead>
  <tbody>{standalone_rows}</tbody>
</table>
"""
    return f"Browse: {db_name}", body


def _esc_bool(value):
    return "Yes" if value else "No"


run(handler)
