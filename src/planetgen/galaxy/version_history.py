# planetgen/galaxy/version_history.py

"""
The version keys an update has run under (OPS.13, OPS.58;
docs/design/generation/reproducible-galaxies.md).

Every update looks, for each galaxy database, at the newest row of the control database's `version_key_history`. When the
galaxy's seed, the version key (`version_key.version_key`), the generator epoch, the generator battery's digest
(`generation.generator_battery`), the SHA-256 of `requirements.lock` and the environment (`environment`) are all as that
row has them, the update only sets the row's `last_seen_at`. When any differs it adds a row, so the history is the list of
changes, not of updates, and the last `KEEP` rows per database hold the last ten changes. The galaxy seed itself never
changes on an update: the key is recorded next to it, since a changed seed would make a different galaxy.
"""

import datetime
import hashlib
import importlib.metadata
import json
import os
import platform

from planetgen._version import __version__
from planetgen.galaxy import version_key
from planetgen.generation import generator_battery
from planetgen.util import generator_epoch

KEEP = 10
"""int: How many history rows are kept per galaxy database."""

LOCK_FILE = "requirements.lock"

LIBRARIES = ("numpy", "astropy", "scipy", "scikit-image")
"""tuple: The distributions whose versions go in the environment record: the ones whose arithmetic the generator's
floats come from."""

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

_FIELDS = ("galaxy_seed", "version_key", "generator_epoch", "battery_digest", "requirements_sha256", "environment")
"""tuple: What an update compares with the newest row: a row is added when any of these differs."""


def requirements_sha256(path=None):
    """SHA-256 hex of the requirements lock file (the checkout's by default),
    or `None` when there is none beside the code."""
    path = path or os.path.join(_REPO_ROOT, LOCK_FILE)
    try:
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    except OSError:
        return None


def environment():
    """
    The environment the generator's floats come from, as a dict: `python`, `system`, `machine`, `libc` (the C library
    and its version, `platform.libc_ver()`; empty where there is none) and the versions of `LIBRARIES` (`None` for one
    that is not installed).
    """
    libc, libc_version = platform.libc_ver()
    found = {"python": platform.python_version(), "system": platform.system(), "machine": platform.machine(),
             "libc": f"{libc} {libc_version}".strip()}
    for name in LIBRARIES:
        try:
            found[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            found[name] = None
    return found


def _decoded(value):
    """A JSON column's value as Python data (the driver gives text)."""
    return json.loads(value) if isinstance(value, (str, bytes)) else value


def _row(conn, database):
    """`database`'s newest row as a dict with its JSON columns decoded, or `None`."""
    row = conn.execute(
        "SELECT id, galaxy_seed, version_key, planetgen_version, requirements_sha256, recorded_at, generator_epoch,"
        " battery_digest, battery_cases, environment, last_seen_at FROM version_key_history"
        " WHERE database_name = ? ORDER BY id DESC LIMIT 1", (database,)).fetchone()
    return None if row is None else _with_decoded(dict(row))


def _with_decoded(row):
    row["battery_cases"] = _decoded(row["battery_cases"])
    row["environment"] = _decoded(row["environment"])
    return row


def record(conn, database, galaxy_seed, lock_sha256=None, key=None, battery=None, epoch=None, env=None):
    """
    Records that an update ran for `database` (a galaxy database name) from an open control connection: sets the newest
    row's `last_seen_at` when nothing it records differs, otherwise adds a row (and drops all but the newest `KEEP` of
    the database's rows).

    Args:
        galaxy_seed (bytes): The galaxy's 16-byte seed.
        lock_sha256 (str, optional): Defaults to the checkout's lock file's.
        key (str, optional): Defaults to the running code's version key.
        battery (generator_battery.BatteryDigest, optional): Defaults to the running code's, worked out once a process.
        epoch (int, optional): Defaults to the running code's generator epoch.
        env (dict, optional): Defaults to `environment()`.

    Returns:
        str: The version key recorded.
    """
    key = key or version_key.version_key()
    battery = battery or generator_battery.running_digest()
    values = {"galaxy_seed": bytes(galaxy_seed).hex().upper(), "version_key": key,
              "generator_epoch": generator_epoch.generator_epoch() if epoch is None else epoch,
              "battery_digest": battery.overall,
              "requirements_sha256": lock_sha256 if lock_sha256 is not None else requirements_sha256(),
              "environment": environment() if env is None else env}
    now = datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None)
    newest = _row(conn, database)
    if newest is not None and all(newest[field] == values[field] for field in _FIELDS):
        conn.execute("UPDATE version_key_history SET last_seen_at = ? WHERE id = ?", (now, newest["id"]))
        conn.commit()
        return key
    conn.execute(
        "INSERT INTO version_key_history (database_name, galaxy_seed, version_key, planetgen_version,"
        " requirements_sha256, recorded_at, generator_epoch, battery_digest, battery_cases, environment, last_seen_at)"
        " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (database, values["galaxy_seed"], key, __version__, values["requirements_sha256"], now,
         values["generator_epoch"], values["battery_digest"],
         json.dumps(generator_battery.cases_of(battery), sort_keys=True),
         json.dumps(values["environment"], sort_keys=True), now))
    conn.execute(
        "DELETE FROM version_key_history WHERE database_name = ? AND id NOT IN"
        " (SELECT id FROM (SELECT id FROM version_key_history WHERE database_name = ?"
        " ORDER BY id DESC LIMIT ?) newest)",
        (database, database, KEEP))
    conn.commit()
    return key


def history(conn, database):
    """`database`'s rows, newest first, as dicts (`last_seen_at` is `recorded_at` for a row nothing has seen again;
    `battery_cases` and `environment` are decoded, `None` for a row recorded before OPS.58)."""
    rows = conn.execute(
        "SELECT galaxy_seed, version_key, planetgen_version, requirements_sha256, recorded_at, generator_epoch,"
        " battery_digest, battery_cases, environment, COALESCE(last_seen_at, recorded_at) AS last_seen_at"
        " FROM version_key_history WHERE database_name = ? ORDER BY id DESC", (database,)).fetchall()
    return [_with_decoded(dict(row)) for row in rows]


def changes(newer, older):
    """
    What differs between a history row and the one before it, as short phrases (`epoch 54 to 55`, `lock file`,
    `numpy 2.1.0 to 2.1.1`, `battery (first case sector-a)`); empty when only the release does, or when `older` is
    `None`.
    """
    if older is None:
        return []
    found = []
    if newer["galaxy_seed"] != older["galaxy_seed"]:
        found.append("galaxy seed")
    if older["generator_epoch"] is not None and newer["generator_epoch"] != older["generator_epoch"]:
        found.append(f"epoch {older['generator_epoch']} to {newer['generator_epoch']}")
    if newer["requirements_sha256"] != older["requirements_sha256"]:
        found.append("lock file")
    old_env, new_env = older["environment"] or {}, newer["environment"] or {}
    for name in sorted(old_env.keys() & new_env.keys()):
        if old_env[name] != new_env[name]:
            found.append(f"{name} {old_env[name]} to {new_env[name]}")
    if older["battery_digest"] is not None and newer["battery_digest"] != older["battery_digest"]:
        first = generator_battery.first_difference(older["battery_cases"] or {}, newer["battery_cases"] or {})
        found.append("battery digest" + (f" (first case {first})" if first else ""))
    return found
