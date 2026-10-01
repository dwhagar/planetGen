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
import sys
import time

import pymysql
from flask import Blueprint, current_app, jsonify, request

import adminStats
from stellarObjects import _db, adminAuth, generationStats, loginThrottle
from stellarObjects._version import __version__

from .authz import audit, require_admin
from .common import ApiError, get_control_db, require_json_body
from .loginguard import with_store
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


def _libraries_dir():
    """Where this process imports Flask from: apt's
    /usr/lib/python3/dist-packages, or /usr/local/.../dist-packages when
    pip installed it. Shows which libraries mod_wsgi really uses."""
    flask = sys.modules.get("flask")
    path = getattr(flask, "__file__", None)
    return os.path.dirname(os.path.dirname(path)) if path else None


def _api_process_info():
    try:
        load = os.getloadavg()
    except (AttributeError, OSError):
        load = None
    return {
        "version": __version__,
        "python_version": platform.python_version(),
        "python_prefix": sys.prefix,
        "libraries_dir": _libraries_dir(),
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
    estimated rows per table, exact sector/system counts, bright-star
    placed/filled counts (`adminStats.bright_star_counts`), schema version,
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
        "bright_stars": adminStats.bright_star_counts(conn),
        "tables": tables,
        "timestamps": adminStats.timestamp_stats(conn),
        "name_collisions": adminStats.name_collision_summary(conn),
    }
    body["query_ms"] = round((time.perf_counter() - started) * 1000, 1)
    return jsonify(body)


LOGIN_FAILURES_DEFAULT_LIMIT = 20
LOGIN_FAILURES_MAX_LIMIT = 200


@bp.route("/login-failures")
@require_admin(fresh=True)
def login_failures():
    """
    `GET /api/admin/login-failures[?limit=]` -- the newest refused sign-ins
    (wrong username or password, a locked login, a wrong current password
    on change-credentials), newest first: `{"items": [{"action",
    "username", "ip", "created_at"}]}`. `limit` defaults to 20, at most
    200. Kept for `adminAuth.LOGIN_FAILURE_RETENTION_DAYS` days.
    """
    raw = request.args.get("limit")
    limit = LOGIN_FAILURES_DEFAULT_LIMIT
    if raw not in (None, ""):
        try:
            limit = int(raw)
        except ValueError:
            raise ApiError("'limit' must be an integer")
        if not 1 <= limit <= LOGIN_FAILURES_MAX_LIMIT:
            raise ApiError(f"'limit' must be between 1 and {LOGIN_FAILURES_MAX_LIMIT}")
    rows = adminAuth.recent_login_failures(get_control_db(), limit=limit)
    return jsonify({"items": [{
        **row,
        # UTC (the connection's zone), with an explicit offset.
        "created_at": row["created_at"].isoformat() + "Z" if row["created_at"] else None,
    } for row in rows]})


def _proxy_warning(locked, failures):
    """
    Whether the site looks like it sits behind a reverse proxy whose
    address isn't unwrapped (`proxy_fix.x_for` is 0): a private address is
    locked, or most recent failed sign-ins share one private or loopback
    address. Then every visitor shares that address, and a lockout of it
    would shut everyone out (SEC.1).
    """
    if (current_app.config.get("PROXY_FIX") or {}).get("x_for"):
        return False
    if any(row["scope"] == loginThrottle.SCOPE_IP and loginThrottle.is_private_address(row["subject"])
           for row in locked):
        return True
    addresses = [row["ip"] for row in failures if row["ip"]]
    if len(addresses) < 10:
        return False
    top = max(set(addresses), key=addresses.count)
    shared = loginThrottle.is_private_address(top) or top in ("127.0.0.1", "::1")
    return shared and addresses.count(top) >= 0.9 * len(addresses)


@bp.route("/generation-stats")
@require_admin(fresh=True)
def generation_stats():
    """
    `GET /api/admin/generation-stats` -- how fast this server generates
    and how big a galaxy gets (PERF.10, `stellarObjects/generationStats.py`):
    `{"buckets": [{"kind", "bucket", "density_low", "density_high",
    "samples", "seconds_per_task", "seconds_per_system",
    "systems_per_task", "stars_per_system", "max_density"}], "sizes":
    {database: {"bytes_per_system", "systems", "total_bytes"}},
    "available": bool}`. `available` is false (and both empty) when the
    control schema is older than v6.
    """
    stats = generationStats.GenerationStats()
    try:
        stats.read(get_control_db())
        available = True
    except Exception:  # noqa: BLE001 -- update.sh not run yet
        available = False
    return jsonify({"buckets": stats.rows(), "sizes": stats.sizes, "available": available})


@bp.route("/lockouts")
@require_admin(fresh=True)
def lockouts():
    """
    `GET /api/admin/lockouts` -- every address and username locked right
    now (SEC.1, SEC.21): `{"items": [{"scope", "subject", "retry_after",
    "locked_until", "level"}], "proxy_warning": bool}`. `scope` is `ip`
    (an IPv6 subject is a /64) or `user` (case-folded). `proxy_warning`
    says the site seems to be behind a reverse proxy without `proxy_fix`.
    """
    rows = with_store(loginThrottle.locked_subjects)
    failures = adminAuth.recent_login_failures(get_control_db(), limit=50)
    return jsonify({
        "items": [{
            "scope": row["scope"],
            "subject": row["subject"],
            "retry_after": row["retry_after"],
            "locked_until": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(row["locked_until"])),
            "level": row["level"],
        } for row in rows],
        "proxy_warning": _proxy_warning(rows, failures),
    })


@bp.route("/lockouts/lift", methods=["POST"])
@require_admin(fresh=True)
def lift_lockout():
    """
    `POST /api/admin/lockouts/lift` `{"scope": "ip"|"user", "subject": str}`
    lifts one lockout and forgets its count (an address's doubling level
    too); `{"all": true}` lifts every one. Returns `{"lifted": n}`.
    Written to the audit and activity logs as `lockout.lift`.
    """
    body = require_json_body()
    if body.get("all") is True:
        lifted = with_store(lambda store: store.lift())
        audit("lockout.lift", target="all", detail=f"lifted={lifted}")
        return jsonify({"lifted": lifted})
    scope = body.get("scope")
    subject = body.get("subject")
    if scope not in loginThrottle.SCOPES:
        raise ApiError("'scope' must be 'ip' or 'user'")
    if not isinstance(subject, str) or not subject.strip() or len(subject) > loginThrottle.MAX_SUBJECT_LENGTH:
        raise ApiError("'subject' is required")
    subject = subject.strip()
    if scope == loginThrottle.SCOPE_USER:
        subject = loginThrottle.normalize_username(subject)
    lifted = with_store(lambda store: store.lift(scope, subject))
    audit("lockout.lift", target=f"{scope}:{subject}", detail=f"lifted={lifted}")
    return jsonify({"lifted": lifted})


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
