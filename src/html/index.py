#!/usr/bin/env python3
# html/index.py

"""
Landing page: lists every MySQL schema on the configured server whose
name matches this deployment's prefix (`GET /api/databases`) and links to
`browse.py` for each one.
"""

import os
import sys
from urllib.parse import quote

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from apiclient import ApiError, list_databases
from fmt import esc
from page import query_params, redirect, run


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
        # sector_count/system_count come back None (not 0) for a schema
        # matching the configured prefix but missing this project's own
        # tables -- see GET /api/databases's own comment for why that's
        # reported rather than 500ing the whole listing.
        sector_count = entry["sector_count"] if entry["sector_count"] is not None else "?"
        system_count = entry["system_count"] if entry["system_count"] is not None else "?"
        rows.append(
            "<tr>"
            f'<td><a href="browse.py?db={esc(entry["name"])}">{esc(entry["name"])}</a></td>'
            f'<td>{sector_count}</td>'
            f'<td>{system_count}</td>'
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


# A choice of exactly one database isn't a choice -- skip the picker table
# entirely and go straight to it, same as picking its only row would.
# Zero (nothing to redirect to) and 2+ (an actual choice) both fall through
# to the normal `run(handler)` picker below, unchanged. `?all=1` (used by
# the sidenav's Databases link -- see `lib/page.py`'s `_sidenav_html`)
# forces the picker table even for a single database, since that's the
# only way back to the database-info page (size, sector/system counts,
# last-modified) on the single-database deployments this redirect would
# otherwise strand every other page behind.
try:
    _databases = list_databases()
except ApiError:
    # The planetGen API is unreachable -- fall through to run(handler),
    # whose own ApiError handling renders a clear 502 page instead of a
    # raw traceback (this probe would otherwise crash before run() ever
    # gets a chance to catch anything).
    _databases = None

if _databases is not None and len(_databases) == 1 and not query_params().get("all"):
    redirect(f"browse.py?db={quote(_databases[0]['name'])}")
else:
    run(handler)
