# planetgen/galaxy/version_history.py

"""
The version keys an update has run under (OPS.13;
docs/design/reproducible-galaxies.md).

Every update records, for each galaxy database, one row in the control
database's `version_key_history`: the galaxy's seed, the version key
(`version_key.version_key`), the release, the SHA-256 of `requirements.lock`
and the time. Only the last `KEEP` rows per database are kept. The galaxy
seed itself never changes on an update: the key is recorded next to it, since
a changed seed would make a different galaxy.
"""

import datetime
import hashlib
import os

from planetgen._version import __version__
from planetgen.galaxy import version_key

KEEP = 10
"""int: How many history rows are kept per galaxy database."""

LOCK_FILE = "requirements.lock"

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


def requirements_sha256(path=None):
    """SHA-256 hex of the requirements lock file (the checkout's by default),
    or `None` when there is none beside the code."""
    path = path or os.path.join(_REPO_ROOT, LOCK_FILE)
    try:
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    except OSError:
        return None


def record(conn, database, galaxy_seed, lock_sha256=None, key=None):
    """
    Adds a row for `database` (a galaxy database name) from an open control
    connection and drops all but the newest `KEEP` of its rows.

    Args:
        galaxy_seed (bytes): The galaxy's 16-byte seed.
        lock_sha256 (str, optional): Defaults to the checkout's lock file's.
        key (str, optional): Defaults to the running code's version key.

    Returns:
        str: The version key recorded.
    """
    key = key or version_key.version_key()
    lock_sha256 = lock_sha256 if lock_sha256 is not None else requirements_sha256()
    now = datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None)
    conn.execute(
        "INSERT INTO version_key_history (database_name, galaxy_seed, version_key, planetgen_version,"
        " requirements_sha256, recorded_at) VALUES (?, ?, ?, ?, ?, ?)",
        (database, bytes(galaxy_seed).hex().upper(), key, __version__, lock_sha256, now))
    conn.execute(
        "DELETE FROM version_key_history WHERE database_name = ? AND id NOT IN"
        " (SELECT id FROM (SELECT id FROM version_key_history WHERE database_name = ?"
        " ORDER BY id DESC LIMIT ?) newest)",
        (database, database, KEEP))
    conn.commit()
    return key


def history(conn, database):
    """`database`'s rows, newest first, as dicts."""
    rows = conn.execute(
        "SELECT galaxy_seed, version_key, planetgen_version, requirements_sha256, recorded_at"
        " FROM version_key_history WHERE database_name = ? ORDER BY id DESC", (database,)).fetchall()
    return [dict(row) for row in rows]
