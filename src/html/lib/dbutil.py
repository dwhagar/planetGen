# html/lib/dbutil.py

"""
Shared helpers for the planetGen CGI web interface.

These scripts run as plain Apache2 CGI scripts on a VPS with a system
Python 3, the `planetGen` package, and its `pymysql`/`DBUtils` runtime
dependencies installed (see `setup.py`'s `install_requires` -- unlike the
pre-MySQL-port version of this module, `stellarObjects` is no longer
optional here: there's no standard-library fallback for talking to
MySQL). Every query here is read-only; the web interface never writes to
a database -- see `queryDb.py`'s module docstring for the same
read-only-by-account-grant convention this module's `open_readonly`
follows.

Not part of the `stellarObjects` package's public API -- this module is
web-plumbing specific to `html/`, kept out of the CGI-mapped document root
(see `examples/apache/planetgen.conf.example`, which denies direct web access to
this `lib/` directory) purely so it can't accidentally be requested and
executed as a script in its own right.

"Multiple databases" (`index.py`'s picker, `?db=` on every other page)
means multiple MySQL schemas on one configured server -- e.g. one
schema per campaign/galaxy (`planetgen`, `planetgen_alpha`, ...) -- rather
than the pre-MySQL-port meaning of multiple `.db` files in a directory.
`DB_PREFIX_ENV_VAR` restricts which schemas on the server are offered,
both so a shared server's unrelated schemas never show up in the picker
and so `?db=` can't be used to probe/open a schema this deployment was
never meant to expose (`resolve_db_name` validates against that same
prefix-filtered list before ever connecting with it).
"""

import html
import os
import re
import sys

DB_PREFIX_ENV_VAR = "PLANETGEN_MYSQL_DATABASE_PREFIX"
"""str: Apache `SetEnv` variable name overriding the default schema-name
prefix -- see the example vhost config."""

DEFAULT_DB_PREFIX = "planetgen"
"""str: Matches `stellarObjects._db.MySQLConfig`'s own default database
name -- a deployment with just one schema names it `planetgen` and never
needs to set `DB_PREFIX_ENV_VAR` at all; one with several names them
`planetgen_<something>` to share the prefix."""

_LIB_DIR = os.path.dirname(os.path.abspath(__file__))
_HTML_DIR = os.path.dirname(_LIB_DIR)
_SRC_DIR = os.path.dirname(_HTML_DIR)
_PROJECT_ROOT = os.path.dirname(_SRC_DIR)

# Falls back to src/ (stellarObjects lives at src/stellarObjects/, src
# layout) so `stellarObjects` is importable even when it hasn't been
# `pip install`-ed system-wide -- true for the default deployment layout
# (`html/` and `src/` as siblings under /var/lib/planetGen). Every caller
# of this module reaches it via its own sys.path setup too, but that isn't
# guaranteed to include this (e.g. browse.py only adds `lib/`), so this
# module makes sure of it independently rather than relying on import
# order.
sys.path.append(os.path.join(_PROJECT_ROOT, "src"))

from stellarObjects.utils import milliparsecs_to_ly
from stellarObjects.physical_constants import LOCAL_STELLAR_DENSITY_LY3
from stellarObjects._db import MySQLConfig, get_connection


class NotFoundError(Exception):
    """Raised for a missing/invalid database name or row id -- callers
    turn this into a 404 response."""


def _server_connection():
    """
    Opens a connection to the configured MySQL server with no specific
    schema selected -- used only to list what schemas exist
    (`information_schema`), never to query planetGen data itself.
    `ensure_schema=False`: there's no default schema selected for
    `_ensure_schema`'s `CREATE TABLE` statements to even apply against,
    and this is a read-only account regardless (see the module docstring).
    """
    return get_connection(MySQLConfig(database=""), ensure_schema=False)


def list_databases():
    """
    Lists every MySQL schema on the configured server whose name starts
    with the configured prefix (`DB_PREFIX_ENV_VAR`, or `DEFAULT_DB_PREFIX`).

    Returns:
        list[dict]: One entry per schema, sorted by name, each with
                    `name`, `size_bytes` (sum of `data_length`/
                    `index_length` across its tables), and `modified_at`
                    (the latest `information_schema.tables.update_time`
                    across its tables, formatted, or `"unknown"` when
                    the storage engine doesn't track it).
    """
    prefix = os.environ.get(DB_PREFIX_ENV_VAR) or DEFAULT_DB_PREFIX
    conn = _server_connection()
    try:
        schema_rows = conn.execute(
            "SELECT schema_name AS name FROM information_schema.schemata "
            "WHERE schema_name LIKE ? ORDER BY schema_name",
            (f"{prefix}%",),
        ).fetchall()

        entries = []
        for schema_row in schema_rows:
            name = schema_row["name"]
            stats = conn.execute(
                "SELECT COALESCE(SUM(data_length + index_length), 0) AS size_bytes, "
                "MAX(update_time) AS modified_at "
                "FROM information_schema.tables WHERE table_schema = ?",
                (name,),
            ).fetchone()
            modified_at = stats["modified_at"]
            entries.append({
                "name": name,
                "size_bytes": int(stats["size_bytes"] or 0),
                "modified_at": modified_at.strftime("%Y-%m-%d %H:%M") if modified_at else "unknown",
            })
        return entries
    finally:
        conn.close()


def resolve_db_name(name):
    """
    Validates a database name supplied via a query string against the
    same prefix-filtered list `list_databases` offers.

    This is what keeps `?db=` from selecting a schema this deployment
    never meant to expose (every other schema on a shared MySQL server,
    `information_schema` itself, etc.) -- only an exact, case-sensitive
    match against a currently-listed schema is accepted.

    Args:
        name (str): The `db` query parameter, e.g. `"planetgen"`.

    Returns:
        MySQLConfig: Ready to pass to `open_readonly`.

    Raises:
        NotFoundError: If `name` is empty or doesn't match a listed schema.
    """
    if not name:
        raise NotFoundError("No database specified.")
    if name not in {entry["name"] for entry in list_databases()}:
        raise NotFoundError(f"No such database: {name!r}")
    return MySQLConfig(database=name)


def open_readonly(config):
    """
    Opens a connection for this read-only web interface -- see the
    module docstring for why "read-only" is enforced by the configured
    account's grants rather than anything this function does itself.

    Args:
        config (MySQLConfig): From `resolve_db_name`.

    Returns:
        stellarObjects._db.Connection: An open connection.
    """
    return get_connection(config, ensure_schema=False)


def fetch_one(conn, query, params=()):
    """Runs `query` and returns the first row, or `None`."""
    return conn.execute(query, params).fetchone()


def fetch_all(conn, query, params=()):
    """Runs `query` and returns every row as a list."""
    return conn.execute(query, params).fetchall()


def esc(value):
    """
    HTML-escapes any value for safe interpolation into a page -- database
    content (system/star/planet names, flavor text, generated wikitext) is
    user-influenced-adjacent (from `--name`, `--system-file`, etc.) and
    must never be trusted as pre-sanitized HTML.

    Args:
        value: Any value; `None` becomes `""`.

    Returns:
        str: The escaped string.
    """
    if value is None:
        return ""
    return html.escape(str(value), quote=True)


_LOCATION_NEIGHBOR_MARKER = " -- nearest: "
_LOCATION_NEIGHBOR_RE = re.compile(r'^(.*) (\([\d.]+ ly\))$')


def linkify_location(db_name, location, name_to_id):
    """
    HTML-escapes a `star_systems.location` string and turns each nearest-
    neighbor name it lists into a link to that system's page.

    `location` is plain text baked in at generation time by
    `stellarObjects._db._format_location_string`, e.g.
    `"Voranthis Kelmoor -- nearest: Alpha Prime (4.2 ly), Beta (5.1 ly)"` --
    the sector name, then up to 3 "Name (distance ly)" entries
    comma-joined after a fixed `" -- nearest: "` marker (empty when the
    sector has no other systems, in which case this is just the sector
    name with nothing to link). This matches that exact format to pull the
    names back out; anything that doesn't fit it (older data predating the
    "-- nearest:" suffix, or a name not found in `name_to_id`) is left as
    plain escaped text rather than guessed at.

    Args:
        db_name (str): The current `?db=` value, for building system.py URLs.
        location (str): The raw `star_systems.location` value.
        name_to_id (dict[str, int]): Every `star_systems.name` -> `id` in
                                     the same sector, for resolving each
                                     neighbor name to a link target.

    Returns:
        str: HTML-safe markup, neighbor names linked where resolvable.
    """
    if not location:
        return ""
    if _LOCATION_NEIGHBOR_MARKER not in location:
        return esc(location)

    prefix, neighbors_part = location.split(_LOCATION_NEIGHBOR_MARKER, 1)
    linked_entries = []
    for entry in neighbors_part.split(", "):
        match = _LOCATION_NEIGHBOR_RE.match(entry)
        name = match.group(1) if match else None
        system_id = name_to_id.get(name) if name is not None else None
        if match and system_id is not None:
            distance = match.group(2)
            linked_entries.append(
                f'<a href="system.py?db={esc(db_name)}&id={system_id}">{esc(name)}</a> {esc(distance)}'
            )
        else:
            linked_entries.append(esc(entry))

    return f"{esc(prefix)}{_LOCATION_NEIGHBOR_MARKER}" + ", ".join(linked_entries)


def format_density(edge_mpc, system_count):
    """
    Formats a sector's star density as systems per cubic light-year, with a
    percentage relative to `LOCAL_STELLAR_DENSITY_LY3` (the real local
    stellar density sector generation targets -- see
    `stellarObjects.spaceSector`'s module docstring) when that comparison
    can be computed.

    Args:
        edge_mpc (float): The sector's cube edge, in milliparsecs
                          (`sectors.edge_mpc`).
        system_count (int): How many systems are placed in the sector.

    Returns:
        str: e.g. `"0.00329 systems/ly&sup3; (116% of local average)"`.
    """
    edge_ly = milliparsecs_to_ly(edge_mpc)
    if not edge_ly:
        return "n/a"

    density_ly3 = system_count / (edge_ly ** 3)
    text = f"{density_ly3:.5f} systems/ly&sup3;"
    if LOCAL_STELLAR_DENSITY_LY3:
        relative_pct = (density_ly3 / LOCAL_STELLAR_DENSITY_LY3) * 100
        text += f" ({relative_pct:,.0f}% of local average)"
    return text
