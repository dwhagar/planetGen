# html/api/admin.py

"""
Admin-only read endpoints behind the admin stats page (`html/adminstats.py`):
server health and statistics about one content database. Both take the
same `?db=` every read route in `routes.py` does, and both require a
logged-in admin past the forced credential change (`require_admin(fresh=True)`),
since table sizes and server details aren't for public visitors.

The queries themselves live in `src/adminStats.py`.
"""

import os
import platform
import time

import pymysql
from flask import Blueprint, jsonify, request

import adminStats
from stellarObjects import _db
from stellarObjects._version import __version__

from .authz import require_admin
from .common import ApiError
from .routes import _paginate, _resolve_requested_db_config, get_db

bp = Blueprint("admin", __name__, url_prefix="/api/admin")

_STARTED_AT = time.time()
"""float: When this API process imported this module -- close enough to
process start for an uptime figure."""


def _memory_info():
    """Total/available memory from `/proc/meminfo`, in bytes, or `None` on
    a system without it."""
    try:
        with open("/proc/meminfo", encoding="ascii") as handle:
            fields = dict(line.split(":", 1) for line in handle if ":" in line)
        total_kb = int(fields["MemTotal"].split()[0])
        available_kb = int(fields["MemAvailable"].split()[0])
    except (OSError, KeyError, ValueError):
        return None
    return {"total_bytes": total_kb * 1024, "available_bytes": available_kb * 1024}


def _api_process_info():
    try:
        load = os.getloadavg()
    except (AttributeError, OSError):
        load = None
    return {
        "version": __version__,
        "python_version": platform.python_version(),
        "pid": os.getpid(),
        "uptime_seconds": int(time.time() - _STARTED_AT),
        "load_average": list(load) if load is not None else None,
        "memory": _memory_info(),
    }


@bp.route("/stats")
@require_admin(fresh=True)
def stats():
    """
    `GET /api/admin/stats[?db=]` -- health of the API process and the
    MySQL server, plus statistics about the selected database: size and
    estimated rows per table, exact sector/system counts, schema version,
    newest/last-modified row times, and how many names the uniqueness
    rules had to decorate (the list itself is `/api/admin/duplicate-names`).

    A database that can't be reached still returns 200 with
    `database.reachable: false`, so the admin page can show what's wrong
    instead of an error page.
    """
    started = time.perf_counter()
    db_name = _resolve_requested_db_config().database
    body = {"api": _api_process_info(), "mysql": None, "database": {"name": db_name, "reachable": False}}

    try:
        conn = get_db()
        conn.execute("SELECT 1")
    except Exception as exc:
        body["database"]["detail"] = str(exc)
        return jsonify(body)

    version = adminStats.schema_version(conn)
    tables = adminStats.table_stats(conn)
    body["mysql"] = adminStats.server_info(conn)
    body["database"] = {
        "name": db_name,
        "reachable": True,
        "schema_version": version,
        "schema_expected": _db.SCHEMA_VERSION,
        "schema_current": version == _db.SCHEMA_VERSION,
        "size_bytes": sum(t["data_bytes"] + t["index_bytes"] for t in tables),
        "counts": adminStats.exact_counts(conn),
        "tables": tables,
        "timestamps": adminStats.timestamp_stats(conn),
        "name_collisions": adminStats.name_collision_summary(conn),
    }
    body["query_ms"] = round((time.perf_counter() - started) * 1000, 1)
    return jsonify(body)


@bp.route("/duplicate-names")
@require_admin(fresh=True)
def duplicate_names():
    """
    `GET /api/admin/duplicate-names[?db=&limit=&offset=]` -- one page of
    the base names that collided and were made unique (Alpha/Beta...,
    Little..., ...Kin), each with every sector, system, planet and moon
    now carrying a form of it. Paginated by base name with the same
    `limit`/`offset` rules as every other listing. See
    `adminStats.duplicate_names` for the response shape.
    """
    limit, offset = _paginate(request.args)
    try:
        result = adminStats.duplicate_names(get_db(), limit=limit, offset=offset)
    except pymysql.err.ProgrammingError as exc:
        # A database from before the v24 name registries has nothing to list.
        raise ApiError(f"could not read the name registries ({exc})", status_code=409)
    return jsonify(result)
