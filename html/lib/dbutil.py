# html/lib/dbutil.py

"""
Shared helpers for the planetGen CGI web interface.

Deliberately dependency-free (standard library only) -- these scripts are
meant to run as plain Apache2 CGI scripts on a VPS with nothing beyond a
system Python 3 and the `planetGen` package itself installed. Every query
here is read-only; the web interface never writes to a database.

Not part of the `stellarObjects` package's public API -- this module is
web-plumbing specific to `html/`, kept out of the CGI-mapped document root
(see `apache/planetgen.conf.example`, which denies direct web access to
this `lib/` directory) purely so it can't accidentally be requested and
executed as a script in its own right.
"""

import glob
import html
import os
import sqlite3
import sys

DB_DIR_ENV_VAR = "PLANETGEN_DB_DIR"
"""str: Apache `SetEnv` variable name that overrides the default db
directory -- see the example vhost config."""

_LIB_DIR = os.path.dirname(os.path.abspath(__file__))
_HTML_DIR = os.path.dirname(_LIB_DIR)
_PROJECT_ROOT = os.path.dirname(_HTML_DIR)

# Falls back to the project root (html/'s parent) so `stellarObjects` is
# importable even when it hasn't been `pip install`-ed system-wide -- true
# for the default deployment layout (`html/` and `stellarObjects/` as
# siblings under /var/lib/planetGen). Every caller of this module reaches
# it via its own sys.path setup too, but that isn't guaranteed to include
# the project root (e.g. browse.py only adds `lib/`), so this module makes
# sure of it independently rather than relying on import order.
sys.path.append(_PROJECT_ROOT)
try:
    from stellarObjects.utils import milliparsecs_to_ly
    from stellarObjects.physical_constants import LOCAL_STELLAR_DENSITY_LY3
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # callers fall back to showing raw stored units rather than failing.
    milliparsecs_to_ly = None
    LOCAL_STELLAR_DENSITY_LY3 = None

DEFAULT_DB_DIR = os.path.join(_PROJECT_ROOT, "db")
"""str: `db/` alongside `html/`, matching the repo layout and the
recommended deployment layout (`/var/lib/planetGen/db` next to
`/var/lib/planetGen/html`)."""


class NotFoundError(Exception):
    """Raised for a missing/invalid database name or row id -- callers
    turn this into a 404 response."""


def get_db_dir():
    """
    Returns the directory to look for `.db` files in: `PLANETGEN_DB_DIR`
    if the Apache config sets it, otherwise `DEFAULT_DB_DIR`.

    Returns:
        str: Absolute path to the database directory.
    """
    return os.environ.get(DB_DIR_ENV_VAR) or DEFAULT_DB_DIR


def list_databases():
    """
    Lists every `*.db` file directly inside the database directory.

    Returns:
        list[dict]: One entry per file, sorted by name, each with `name`,
                    `size_bytes`, and `modified_at` (ISO 8601 local time).
    """
    db_dir = get_db_dir()
    entries = []
    for path in sorted(glob.glob(os.path.join(db_dir, "*.db"))):
        stat = os.stat(path)
        entries.append({
            "name": os.path.basename(path),
            "size_bytes": stat.st_size,
            "modified_at": _format_mtime(stat.st_mtime),
        })
    return entries


def _format_mtime(epoch_seconds):
    import datetime
    return datetime.datetime.fromtimestamp(epoch_seconds).strftime("%Y-%m-%d %H:%M")


def resolve_db_path(name):
    """
    Validates a database name supplied via a query string and resolves it
    to a real file inside the database directory.

    Only an exact basename match against a file that's actually present is
    accepted -- this is what keeps `?db=` from being used for path
    traversal (`../../etc/passwd` et al.) or for opening arbitrary files
    outside the database directory.

    Args:
        name (str): The `db` query parameter, e.g. `"planetgen.db"`.

    Returns:
        str: Absolute path to the validated `.db` file.

    Raises:
        NotFoundError: If `name` is empty or doesn't match a listed file.
    """
    if not name:
        raise NotFoundError("No database specified.")
    candidate = os.path.basename(name)
    db_dir = get_db_dir()
    path = os.path.join(db_dir, candidate)
    if candidate != name or not os.path.isfile(path) or not candidate.endswith(".db"):
        raise NotFoundError(f"No such database: {name!r}")
    return path


def open_readonly(path):
    """
    Opens a SQLite database strictly read-only, via a `file:` URI with
    `mode=ro` -- this fails outright if the database doesn't already
    exist, rather than silently creating one, and guarantees the web
    interface can never write to it even if a query is buggy.

    Args:
        path (str): Absolute path to the `.db` file (from `resolve_db_path`).

    Returns:
        sqlite3.Connection: A read-only connection with `Row` row access.
    """
    uri = f"file:{path}?mode=ro"
    conn = sqlite3.connect(uri, uri=True)
    conn.row_factory = sqlite3.Row
    return conn


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
        str: e.g. `"0.00329 systems/ly&sup3; (116% of local average)"`, or
             `"unknown (cube edge unit unavailable)"` if `stellarObjects`
             isn't importable in this deployment.
    """
    if milliparsecs_to_ly is None:
        return "unknown (cube edge unit unavailable)"

    edge_ly = milliparsecs_to_ly(edge_mpc)
    if not edge_ly:
        return "n/a"

    density_ly3 = system_count / (edge_ly ** 3)
    text = f"{density_ly3:.5f} systems/ly&sup3;"
    if LOCAL_STELLAR_DENSITY_LY3:
        relative_pct = (density_ly3 / LOCAL_STELLAR_DENSITY_LY3) * 100
        text += f" ({relative_pct:,.0f}% of local average)"
    return text
