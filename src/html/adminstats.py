#!/usr/bin/env python3
# html/adminstats.py

"""
Admin stats page: server health and statistics about the current
database, from `GET /api/admin/stats` and `GET /api/admin/duplicate-names`
(see `html/api/admin.py`), plus the galaxy tile cache, which lives on
this web server's disk rather than behind the API (`lib/tilecache.py`).

Sections: health (API process, MySQL server, tile cache), the database
itself (size, schema version, exact sector/system counts), when rows were
last created or changed (v27's `created_at`/`modified_at`), per-table
sizes, and every name the uniqueness rules had to decorate (Alpha/Beta...,
Little..., ...Kin) with links to each sector and system carrying one.
Planets and moons have no page of their own, so they link to their system.

The database is the request's `db` when one was posted (the sidenav and
every link here carry it), otherwise the first one `GET /api/databases`
lists -- the same default `index.py` shows.

Same login gate as `admin.py`: not logged in goes to `login.py`, still on
the seeded credentials goes to `changecreds.py`.
"""

import os
import shutil
import sys
import time

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

import tilecache  # noqa: E402
from apiclient import (  # noqa: E402
    ApiError, NotFoundError, admin_duplicate_names, admin_stats, auth_me, list_databases,
)
from fmt import esc, post_link  # noqa: E402
from page import incoming_cookie_header, nav_params, redirect, render, render_error  # noqa: E402

_NAMES_PAGE_SIZE = 100
"""int: Base names per page of the duplicate-names list."""


def _format_bytes(value):
    if value is None:
        return "unknown"
    size = float(value)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024


def _format_duration(seconds):
    if seconds is None:
        return "unknown"
    seconds = int(seconds)
    days, seconds = divmod(seconds, 86400)
    hours, seconds = divmod(seconds, 3600)
    minutes = seconds // 60
    if days:
        return f"{days}d {hours}h"
    if hours:
        return f"{hours}h {minutes}m"
    return f"{minutes}m"


def _format_count(value):
    return "unknown" if value is None else f"{value:,}"


def _parse_offset(raw):
    try:
        return max(0, int(raw))
    except (TypeError, ValueError):
        return 0


def _tile_cache_info():
    """
    Where the galaxy tile cache lives, how much it holds against its
    budget, and how much room is left on that disk. Reads the directory
    without creating it (unlike `tilecache.cache_dir`).
    """
    configured = tilecache.configured_cache_dir()
    info = {
        "enabled": configured is not None,
        "dir": configured,
        "exists": False,
        "size_bytes": 0,
        "files": 0,
        "max_bytes": tilecache.max_cache_bytes(),
        "disk_free_bytes": None,
    }
    if configured is None or not os.path.isdir(configured):
        return info
    info["exists"] = True
    for dirpath, _dirnames, filenames in os.walk(configured):
        for name in filenames:
            try:
                info["size_bytes"] += os.stat(os.path.join(dirpath, name)).st_size
                info["files"] += 1
            except OSError:
                continue
    try:
        info["disk_free_bytes"] = shutil.disk_usage(configured).free
    except OSError:
        pass
    return info


def _tiles(pairs):
    return '<div class="stat-tiles">' + "".join(
        f'<div class="stat-tile"><span class="stat-value">{value}</span>'
        f'<span class="stat-label">{esc(label)}</span></div>'
        for label, value in pairs
    ) + "</div>"


def _kv_table(rows):
    body = "".join(f"<tr><th scope=\"row\">{esc(key)}</th><td>{value}</td></tr>" for key, value in rows)
    return f'<div class="table-scroll"><table><tbody>{body}</tbody></table></div>'


def _health_html(stats, api_ms, cache):
    api = stats["api"]
    mysql = stats.get("mysql") or {}
    database = stats["database"]

    if database["reachable"] and database.get("schema_current"):
        status = '<span class="badge">healthy</span>'
    elif database["reachable"]:
        status = '<span class="badge">schema out of date</span>'
    else:
        status = '<span class="badge">database unreachable</span>'

    memory = api.get("memory")
    memory_text = (
        f"{_format_bytes(memory['available_bytes'])} free of {_format_bytes(memory['total_bytes'])}"
        if memory else "unknown"
    )
    load = api.get("load_average")
    load_text = " / ".join(f"{value:.2f}" for value in load) if load else "unknown"

    if not cache["enabled"]:
        cache_text = "disabled (max size is 0)"
    elif not cache["exists"]:
        cache_text = f"<code>{esc(cache['dir'])}</code> doesn't exist yet"
    else:
        free = _format_bytes(cache["disk_free_bytes"]) if cache["disk_free_bytes"] is not None else "unknown"
        cache_text = (
            f"{_format_bytes(cache['size_bytes'])} of {_format_bytes(cache['max_bytes'])} "
            f"({_format_count(cache['files'])} files) in <code>{esc(cache['dir'])}</code>, {free} free on disk"
        )

    rows = [
        ("Status", status),
        ("API response", f"{api_ms:.0f} ms (queries took {stats.get('query_ms', 0):.0f} ms)"),
        ("API version", esc(api["version"])),
        ("API uptime", esc(_format_duration(api["uptime_seconds"]))),
        ("Python", esc(api["python_version"])),
        ("Load average (1/5/15 min)", esc(load_text)),
        ("Memory", esc(memory_text)),
        ("MySQL version", esc(mysql.get("version") or "unknown")),
        ("MySQL uptime", esc(_format_duration(mysql.get("uptime_seconds")))),
        ("MySQL connections", esc(_format_count(mysql.get("threads_connected")))),
        ("Galaxy tile cache", cache_text),
    ]
    if not database["reachable"]:
        rows.append(("Database error", f'<span class="error">{esc(database.get("detail", ""))}</span>'))
    return f'<section class="panel"><h2>Server Health</h2>{_kv_table(rows)}</section>'


def _database_html(database, db_names):
    counts = database["counts"]
    collisions = database["name_collisions"]
    stamps = {row["table"]: row for row in database["timestamps"]}
    version = database["schema_version"]
    if database["schema_current"]:
        schema_text = f"v{version} (current)"
    elif version is None:
        schema_text = f"not initialized (code expects v{database['schema_expected']}); run migrateDb.py"
    else:
        schema_text = f"v{version}, code expects v{database['schema_expected']}; run migrateDb.py"

    picker = ""
    if len(db_names) > 1:
        options = "".join(
            f'<option value="{esc(name)}"{" selected" if name == database["name"] else ""}>{esc(name)}</option>'
            for name in db_names
        )
        picker = f"""
<form method="post" action="adminstats.py" class="search-form">
  <div class="search-fields">
    <label class="search-field">Database <select name="db">{options}</select></label>
  </div>
  <div class="search-actions"><button type="submit" class="btn">Show</button></div>
</form>
"""

    tiles = _tiles([
        ("Sectors", _format_count(counts.get("sectors"))),
        ("Star systems", _format_count(counts.get("star_systems"))),
        ("Size on disk", _format_bytes(database["size_bytes"])),
        ("Names made unique", _format_count(collisions.get("distinct_base_names"))),
    ])
    rows = [
        ("Database", f"<code>{esc(database['name'])}</code>"),
        ("Schema", esc(schema_text)),
        ("Last system change", esc(stamps.get("star_systems", {}).get("last_modified_at") or "never")),
        ("Newest system", esc(stamps.get("star_systems", {}).get("newest_created_at") or "none")),
        ("Last sector change", esc(stamps.get("sectors", {}).get("last_modified_at") or "never")),
    ]
    return f"""
<section class="panel">
<h2>Database</h2>
{picker}
{tiles}
{_kv_table(rows)}
</section>
"""


def _timestamps_html(database):
    rows = "".join(
        f"<tr><td><code>{esc(row['table'])}</code></td>"
        f"<td>{esc(row['newest_created_at'] or 'none')}</td>"
        f"<td>{esc(row['last_modified_at'] or 'never')}</td></tr>"
        for row in database["timestamps"]
    )
    return f"""
<section class="panel">
<h2>Recent Activity</h2>
<p class="hint">When each kind of top-level row was last created or changed.
A change to a star, planet or moon counts as a change to its system.</p>
<div class="table-scroll"><table>
  <thead><tr><th>Table</th><th>Newest row created</th><th>Last modified</th></tr></thead>
  <tbody>{rows}</tbody>
</table></div>
</section>
"""


def _tables_html(database):
    rows = "".join(
        f"<tr><td><code>{esc(t['name'])}</code></td>"
        f"<td class=\"num\">{_format_count(t['approx_rows'])}</td>"
        f"<td class=\"num\">{_format_bytes(t['data_bytes'])}</td>"
        f"<td class=\"num\">{_format_bytes(t['index_bytes'])}</td></tr>"
        for t in database["tables"]
    )
    return f"""
<section class="panel">
<h2>Tables</h2>
<p class="hint">Row counts here are MySQL's own estimates and can lag behind;
the sector and system totals above are exact.</p>
<div class="table-scroll"><table>
  <thead><tr><th>Table</th><th class="num">Rows (approx.)</th><th class="num">Data</th><th class="num">Indexes</th></tr></thead>
  <tbody>{rows}</tbody>
</table></div>
</section>
"""


_LEVEL_LABELS = {
    "sector": "sector vs sector",
    "system": "system vs system or sector",
    "body": "planet or moon",
}


def _name_link(db_name, row):
    if row["kind"] == "sector":
        return post_link("sector.py", {"db": db_name, "id": row["id"]}, esc(row["name"])) + ' <span class="hint">sector</span>'
    if row["kind"] == "system":
        return post_link("system.py", {"db": db_name, "id": row["id"]}, esc(row["name"])) + ' <span class="hint">system</span>'
    system_link = post_link("system.py", {"db": db_name, "id": row["star_system_id"]}, esc(row["system_name"]))
    return f'{esc(row["name"])} <span class="hint">{row["kind"]} in</span> {system_link}'


def _duplicate_names_html(db_name, database, names, error):
    collisions = database["name_collisions"]
    summary = _tiles([
        ("Names made unique", _format_count(collisions.get("distinct_base_names"))),
        ("Sector collisions", _format_count(collisions.get("sector"))),
        ("System collisions", _format_count(collisions.get("system"))),
        ("Planet/moon collisions", _format_count(collisions.get("body"))),
    ])
    if error:
        listing = f'<p class="error">{esc(error)}</p>'
    elif not names["items"]:
        listing = '<p class="hint">No names have had to be made unique.</p>'
    else:
        rows = []
        for item in names["items"]:
            links = "<br>".join(_name_link(db_name, row) for row in item["rows"]) or '<span class="hint">no rows left</span>'
            levels = ", ".join(_LEVEL_LABELS.get(level, level) for level in item["levels"])
            rows.append(
                f"<tr><td>{esc(item['base_name'])}</td><td>{esc(levels)}</td><td>{links}</td></tr>"
            )
        listing = f"""
<div class="table-scroll"><table>
  <thead><tr><th>Base name</th><th>Collided as</th><th>Now named</th></tr></thead>
  <tbody>{''.join(rows)}</tbody>
</table></div>
{_names_pagination(db_name, names)}
"""
    return f"""
<section class="panel" id="duplicate-names">
<h2>Names Made Unique</h2>
<p class="hint">Names that came up more than once and were given a Greek
letter (Alpha, Beta, ...), a "little" prefix (a system named after a
sector) or a companion suffix (a planet or moon), with every sector,
system, planet and moon that now carries one.</p>
{summary}
{listing}
</section>
"""


def _names_pagination(db_name, names):
    total, offset, limit = names["total"], names["offset"], names["limit"]
    if total <= limit:
        return ""
    shown_to = min(offset + limit, total)
    prev_link = (
        post_link("adminstats.py#duplicate-names", {"db": db_name, "names_offset": max(0, offset - limit)}, "&larr; Prev")
        if offset > 0 else '<span class="hint">&larr; Prev</span>'
    )
    next_link = (
        post_link("adminstats.py#duplicate-names", {"db": db_name, "names_offset": offset + limit}, "Next &rarr;")
        if shown_to < total else '<span class="hint">Next &rarr;</span>'
    )
    return (
        f'<div class="pagination">{prev_link}'
        f"<span>{offset + 1:,}&ndash;{shown_to:,} of {total:,}</span>{next_link}</div>"
    )


def _page_html(db_name, db_names, cookie_header, names_offset):
    started = time.perf_counter()
    stats = admin_stats(cookie_header, db_name)
    api_ms = (time.perf_counter() - started) * 1000
    database = stats["database"]

    sections = [_health_html(stats, api_ms, _tile_cache_info())]
    if database["reachable"]:
        names = None
        names_error = None
        try:
            names = admin_duplicate_names(cookie_header, db_name, limit=_NAMES_PAGE_SIZE, offset=names_offset)
        except ApiError as exc:
            names_error = str(exc)
        sections.extend([
            _database_html(database, db_names),
            _duplicate_names_html(db_name, database, names, names_error),
            _timestamps_html(database),
            _tables_html(database),
        ])
    return "".join(sections)


cookie_header = incoming_cookie_header()
try:
    identity = auth_me(cookie_header)
except ApiError as exc:
    render_error(f"Could not reach the planetGen API ({exc}).", status="502 Bad Gateway")
if identity is None:
    redirect("login.py")
    sys.exit(0)
if identity["must_change_credentials"]:
    redirect("changecreds.py")
    sys.exit(0)

params = nav_params()
try:
    db_names = [entry["name"] for entry in list_databases()]
except (ApiError, NotFoundError) as exc:
    render_error(f"Could not reach the planetGen API ({exc}).", status="502 Bad Gateway")

db_name = params.get("db") or (db_names[0] if db_names else None)
if not db_name:
    render_error("No databases found.", status="404 Not Found")
if db_name not in db_names:
    render_error(f"Unknown database {db_name!r}.", status="404 Not Found")

try:
    body = _page_html(db_name, db_names, cookie_header, _parse_offset(params.get("names_offset")))
except (ApiError, NotFoundError) as exc:
    render_error(f"Could not load stats from the planetGen API ({exc}).", status="502 Bad Gateway")

render("Server & Database Stats", body)
