# planetgen/api/admin.py

"""
Admin-only read endpoints behind the admin stats page (`html/adminstats.py`):
server health and statistics about one content database. Both take the
same `?db=` every read route in `routes.py` does, and both require a
logged-in admin past the forced credential change (`require_admin(fresh=True)`),
since table sizes and server details aren't for public visitors.

The queries themselves live in `planetgen.db.stats`.
"""

import os
import platform
import secrets
import sys
import time

import pymysql
from flask import Blueprint, current_app, g, jsonify, request

from planetgen.db import stats as adminStats
from planetgen.db import store as _db
from planetgen.admin import auth, throttle
from planetgen.generation import stats as generationStats
from planetgen.api import naming
from planetgen.galaxy import settings_file
from planetgen.names import naming_key
from planetgen._version import __version__

from .authz import audit, require_admin
from .common import ApiError, get_control_db, require_json_body
from .loginguard import with_store
from .schemas import LockoutLift, NamingKeyChange, parse_body
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


def _database_disk(conn, config):
    """The drive holding the database's data directory (never the boot drive, unless the data is there),
    for the stats page: `{"measured", "path", "mount", "total_bytes", "free_bytes", "source", "note"}`."""
    from planetgen.generation import stats as generationStats

    disk = generationStats.database_disk(conn, config.host, config.database)
    if isinstance(disk, generationStats.DiskSpace):
        return {"measured": True, "path": disk.path, "mount": disk.mount, "where": disk.where(),
                "total_bytes": disk.total_bytes, "free_bytes": disk.free_bytes, "source": disk.source, "note": None}
    return {"measured": False, "path": None, "mount": None, "where": None, "total_bytes": None,
            "free_bytes": None, "source": None, "note": disk.reason}


@bp.route("/stats")
@require_admin(fresh=True)
def stats():
    """
    `GET /api/admin/stats[?db=]` -- health of the API process and the
    MySQL server, plus statistics about the selected database: size and
    estimated rows per table, exact sector/system counts, bright-star
    placed/filled counts (`adminStats.bright_star_counts`), the per-sector
    density stats (`adminStats.density_stats`, PERF.11), schema version,
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
        "disk": _database_disk(conn, _resolve_requested_db_config()),
        "counts": adminStats.exact_counts(conn),
        "bright_stars": adminStats.bright_star_counts(conn),
        "sector_stats": adminStats.density_stats(conn),
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
    200. Kept for `auth.LOGIN_FAILURE_RETENTION_DAYS` days.
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
    rows = auth.recent_login_failures(get_control_db(), limit=limit)
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
    if any(row["scope"] == throttle.SCOPE_IP and throttle.is_private_address(row["subject"])
           for row in locked):
        return True
    addresses = [row["ip"] for row in failures if row["ip"]]
    if len(addresses) < 10:
        return False
    top = max(set(addresses), key=addresses.count)
    shared = throttle.is_private_address(top) or top in ("127.0.0.1", "::1")
    return shared and addresses.count(top) >= 0.9 * len(addresses)


@bp.route("/generation-stats")
@require_admin(fresh=True)
def generation_stats():
    """
    `GET /api/admin/generation-stats` -- how fast this server generates
    and how big a galaxy gets (PERF.10, `planetgen/generation/stats.py`):
    `{"buckets": [{"kind", "workers", "bucket", "density_low", "density_high",
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


@bp.route("/generation-stats/reset", methods=["POST"])
@require_admin(fresh=True)
def reset_generation_stats():
    """
    `POST /api/admin/generation-stats/reset` -- deletes every recorded
    generation rate (PERF.32); the next run records afresh and estimates
    from defaults until it has. Returns `{"deleted": n}`. Written to the
    audit and activity logs as `generation_stats.reset`.
    """
    deleted = generationStats.GenerationStats().reset(get_control_db())
    audit("generation_stats.reset", target="generation_stats", detail=f"deleted={deleted}")
    return jsonify({"deleted": deleted})


def _naming_view(database, stored):
    return {
        "database": database, "key": stored["key"] if stored else None,
        "codec_version": stored["codec_version"] if stored else None,
        "current_codec_version": naming_key.CODEC_VERSION,
        "drawn_at": stored["drawn_at"].isoformat() if stored and stored["drawn_at"] else None,
        "changed_at": stored["changed_at"].isoformat() if stored and stored["changed_at"] else None,
        "changed_by": stored["changed_by"] if stored else None,
    }


@bp.route("/naming-key", methods=["GET"])
@require_admin(fresh=True)
def naming_key_get():
    """
    `GET /api/admin/naming-key` -- the galaxy's naming key (GEN.70,
    `planetgen/names/naming_key.py`): `{"database", "key", "codec_version",
    "current_codec_version", "drawn_at", "changed_at", "changed_by"}`.
    `key` is `null` before the galaxy is planned.
    """
    database = _resolve_requested_db_config().database
    return jsonify(_naming_view(database, naming_key.get(get_control_db(), database)))


@bp.route("/naming-key", methods=["POST"])
@require_admin(fresh=True)
def naming_key_set():
    """
    `POST /api/admin/naming-key` `{"key": "0123ABCD"}` sets the naming key
    (8 hex digits), or `{"draw": true}` draws a fresh random one. Every
    object the phoneme codec names takes a new name at once; no row is
    rewritten. 409 before the galaxy is planned. Written to the audit log as
    `naming-key.change`.
    """
    body = parse_body(NamingKeyChange, require_json_body())
    database = _resolve_requested_db_config().database
    key = secrets.token_hex(naming_key.KEY_DIGITS // 2).upper() if body.draw else body.key
    conn = get_control_db()
    try:
        naming_key.change(conn, database, key, g.admin_user["username"])
    except LookupError as exc:
        raise ApiError(str(exc), status_code=409) from exc
    naming.forget(database)
    g.naming_key = key
    audit("naming-key.change", target=database, detail=f"key={key}")
    return jsonify(_naming_view(database, naming_key.get(conn, database)))


@bp.route("/galaxy-settings", methods=["GET"])
@require_admin(fresh=True)
def galaxy_settings():
    """
    `GET /api/admin/galaxy-settings` -- the galaxy's creation-settings files
    (ADM.18, `galaxy/settings_file.py`), newest first: `{"items": [{"name",
    "seed", "version_key", "created_at", "size", "current"}]}`. `current` is
    true for the newest file of each seed; the others are dated backups.
    """
    items, seen = [], set()
    for entry in settings_file.list_files():
        items.append({
            "name": entry["name"], "seed": entry["seed"], "version_key": entry["key"],
            "created_at": entry["when"].isoformat() + "Z", "size": entry["size"],
            "current": entry["seed"] not in seen,
        })
        seen.add(entry["seed"])
    return jsonify({"items": items})


@bp.route("/galaxy-settings/<name>", methods=["GET"])
@require_admin(fresh=True)
def galaxy_settings_file(name):
    """`GET /api/admin/galaxy-settings/<name>` -- one settings file's JSON
    (404 for a name that is not a settings file in the folder)."""
    for entry in settings_file.list_files():
        if entry["name"] == name:
            with open(entry["path"], "r", encoding="utf-8") as handle:
                return current_app.response_class(handle.read(), mimetype="application/json")
    raise ApiError("No such settings file.", status_code=404)


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
    rows = with_store(throttle.locked_subjects)
    failures = auth.recent_login_failures(get_control_db(), limit=50)
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
    body = parse_body(LockoutLift, require_json_body())
    if body.all is True:
        lifted = with_store(lambda store: store.lift())
        audit("lockout.lift", target="all", detail=f"lifted={lifted}")
        return jsonify({"lifted": lifted})
    scope, subject = body.scope, body.subject.strip()
    if scope == throttle.SCOPE_USER:
        subject = throttle.normalize_username(subject)
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
