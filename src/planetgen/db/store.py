# planetgen/db/store.py

"""
Database Persistence (private)
===============================

Writes already-generated `StarSystem`/`SpaceSector` objects into the
MySQL database described by `planetgen/db/schema.sql` and
`docs/reference/database-schema.md`. Leading underscore -- this module is an internal
implementation detail of the persistence boundary, not part of the
package's public generation API (`StarSystem`, `SpaceSector`, `Planet`,
etc. remain the public surface).

This module owns every unit conversion at the point of writing a value
into the database -- generation/physics code elsewhere in the package
keeps its own native units (km/AU/ly) throughout and never imports this
module. See `docs/reference/database-schema.md` for the two-tier distance-unit convention
(milliparsecs for sector-scale placement, kilometers for everything else)
and the full column reference.

The read path (`load_star_system`/`load_sector`) reconstructs live objects
back from rows, inverting every unit conversion the write path applies.
Rather than assembling a nested dict shaped like
`StarSystem.to_dict()` (natural for JSON, awkward here since the data is
normalized across many flat tables) and routing through
`StarSystem.from_dict`, each row is mapped directly to the flat dict shape
each *leaf* class's own `from_dict` already accepts (`Star.from_dict`,
`BinaryStarProxy.from_dict`, `Planet.from_dict`, `AsteroidBelt.from_dict`),
and `StarSystem` itself is assembled directly here the same way
`StarSystem.from_dict` assembles it internally -- reusing the leaf
allowlists/reconstruction logic from `planetgen/util/serialization.py`
(Phase 1) without a redundant object-to-dict-to-object round trip for data
that was never nested to begin with.

MySQL port (TODO.md Phase 5): this module used to talk directly to
`sqlite3`. It now goes through `pymysql` (pure-Python driver, no system
libraries to build against on the deployment host) via a small
`Connection` wrapper (below) that keeps every existing call site's
`conn.execute(sql, params)` shape working unchanged -- including this
module's and `planetgen.db.query`'s `?` positional placeholders, which the
wrapper rewrites to `pymysql`'s `%s` at the point of execution, and
`sqlite3.Row`-style `row["column"]` access, which `pymysql`'s
`DictCursor` already provides natively. Real concurrent access uses a
connection pool (SQLAlchemy's `QueuePool`) rather than opening a fresh TCP
connection per call, per TODO.md's "add real connection pooling" note.
"""

import argparse
import collections
import contextlib
import dataclasses
import functools
import hashlib
import json
import math
import os
import random
import re
import threading
import time
import unicodedata
from types import SimpleNamespace
from collections import namedtuple

import pymysql
import pymysql.cursors
import sqlalchemy.event
import sqlalchemy.exc
import sqlalchemy.pool

from planetgen.admin import activity_log
from planetgen.db import alembic_runner, kind_tables, login, object_rows, orbit_trees, parity, point_columns, scatter_systems, sector_mass, storage_checksum
from planetgen.db import fingerprint_tree as fingerprintTree
from planetgen.names import object_id as objectId
from planetgen import tuning
from planetgen.galaxy import component_field as componentField, components as galaxyComponents, containment, mass_chain as massChain, nebula_field as nebulaField, object_uid as objectUid, scales, scatter_limits as scatterLimits, seed as galaxySeed, shape_draw, uid as galaxyUid, version_key as versionKey
from planetgen.physics import activity, atmosphere, bulk_mass, constants as physical_constants, climate_draws, habitability, habitability_world, hydrosphere, insolation, interior_heat, kepler, magnetism, pre_main_sequence, radiation, radiogenic, spin
from planetgen.physics import aerosols, body_colour, bulk_composition, daylight, giant_interior, life_classes, moon_composition, ocean_tides, orbit_elements, rings, surface_composition
from planetgen.generation import generator_battery as generatorBattery
from planetgen.util import fingerprint_encoding as fingerprintEncoding, generator_epoch as generatorEpoch, log
from planetgen.util.settings import env_name, get_settings
from planetgen.generation import cloud as cloudContext, multistar, steps, stored_trees
from planetgen.generation.belt import AsteroidBelt
from planetgen.generation.phenomena.asteroid_field import AsteroidField, asteroid_field_designation
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.population import facilities as facility_rules, facility_attrs
from planetgen.generation.comet import Comet, comet_designation, rename_comet_designation
from planetgen.generation.config import SystemConfig
from planetgen.generation.prevalence import FEATURES as PREVALENCE_FEATURES
from planetgen.generation.binary import BinaryStarProxy
from planetgen.galaxy.density import GalaxyShape, with_family_columns, with_system_richness
from planetgen.galaxy.frame import galaxy_to_local_pc, local_to_galaxy_pc
from planetgen.galaxy.spherical import (
    CELL_REACH_EDGES, MAX_ADDRESSABLE_SHELL, SectorCell, provisional_sector_designation, sector_address_at,
)
from planetgen.names.bodies import rename_prefix, wide_pair_first_word
from planetgen.names.wordlists import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from planetgen.names.uniqueness import (
    MAX_SYSTEM_NAME_WORDS, fits_word_limit, has_diminutive, resolve_diminutive, resolve_greek_roman_collision,
    word_limit_for,
)
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.planet import Planet
from planetgen.physics.position import MAX_UPDATE_INTERVAL_S, THRESHOLDS_M, update_interval_s
from planetgen.physics import rings, rogue_habitability
from planetgen.physics.rogue_surface import ROGUE_SURFACE_FIELDS
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet, interstellar_comet_designation
from planetgen.galaxy.sector import SectorSystemEntry, SpaceSector, classify_octant, distance_between
from planetgen.generation.star import Star, compressed_heliosphere_radius
from planetgen.generation.phenomena.quasar import Quasar
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.generation.system import StarSystem
from planetgen.galaxy.galactic_orbit import calculate_galactic_orbit
from planetgen.galaxy.straight_line import straight_line_position_pc
from planetgen.galaxy.system_position import galactic_velocity_ms, system_velocity_ms
from planetgen.names.wordsalad import generate_phoneme_salad_name, generate_sector_name
from planetgen.physics.units import ly_to_milliparsecs, ly_to_pc, milliparsecs_to_ly, mpc_to_pc, pc_to_ly
from planetgen.generation.wide_binary import WideBinaryPair
from planetgen.util import draw
from planetgen.util import detmath
from planetgen.physics import unit_registry as units

_au_to_km = units.converter("au", "km")
_km_to_au = units.converter("km", "au")
_yr_to_s = units.converter("yr", "s")
_s_to_yr = units.converter("s", "yr")
_au_yr_to_kms = units.converter("au_yr", "km_s")
_SECONDS_PER_YEAR = _yr_to_s(1.0)

_SYSTEM_POINT = point_columns.LAYOUTS["star_systems"]
"""point_columns.PointLayout: The point of a star system (GEN.335), whose SQL below takes its column lists from it."""

SCHEMA_VERSION = alembic_runner.head_version()
"""int: Matches `star_systems.schema_version` and the highest row in the
`schema_migrations` table (see `planetgen/db/schema.sql`'s header
comment). Also the target version `migrate_database` brings a database's
`schema_migrations` bookkeeping up to. It is the number of the newest
Alembic revision (`planetgen/db/migrations/`)."""


class SchemaTooNewError(RuntimeError):
    """
    Raised when a database's schema version is newer than this code's
    (`SCHEMA_VERSION`, or `CONTROL_SCHEMA_VERSION` for the control
    schema): a newer planetGen has already migrated it, and this older
    code must not write to it or try to migrate it (TEST.10). The fix is
    to update the code (`update.sh`), not the database.
    """

    def __init__(self, database, version, expected, what="schema"):
        self.database = database
        self.version = version
        self.expected = expected
        super().__init__(
            f"database '{database}' is at {what} v{version}, newer than this code's v{expected}: "
            f"a newer planetGen has already migrated it. Update planetGen (git pull, or update.sh) "
            f"before using this database; it is refused rather than migrated or written to."
        )


class SchemaTooOldError(RuntimeError):
    """
    Raised when a database's schema version is older than the oldest this
    code upgrades from (`alembic_runner.BASELINE_VERSION`, DB.11): the
    migration steps below it were removed once Alembic took over. The fix
    is to bring the database to the baseline with an earlier planetGen
    (the last release before DB.11's cleanup), then update.
    """

    def __init__(self, database, version):
        self.database = database
        self.version = version
        super().__init__(
            f"database '{database}' is at schema v{version}, older than v{alembic_runner.BASELINE_VERSION}, "
            f"the oldest this planetGen upgrades from. Run update.sh from an earlier planetGen checkout "
            f"(before the Alembic cleanup) to bring it to v{alembic_runner.BASELINE_VERSION} first."
        )


class SchemaVersionMissingError(RuntimeError):
    """
    Raised when a database has the galaxy tables but no schema version
    recorded (`schema_migrations` empty or missing, DB.97): nothing says
    which migrations it needs, and this code no longer guesses from the
    shape of its tables. The fix is to restore it from a backup, or to
    drop the database and let `update.sh` create a new one.
    """

    def __init__(self, database):
        self.database = database
        super().__init__(
            f"database '{database}' has galaxy tables but no schema version recorded (schema_migrations is empty "
            f"or missing), so planetGen cannot tell which migrations it needs. Restore it from a backup, or drop "
            f"the database and run update.sh to create a new one."
        )


def _has_galaxy_tables(conn):
    """Whether the connected database already holds the galaxy tables (it is not a new, empty one)."""
    row = conn.execute(
        "SELECT COUNT(*) AS n FROM information_schema.tables"
        " WHERE table_schema = DATABASE() AND table_name = 'star_systems'").fetchone()
    return bool(row["n"])


def _refuse_unversioned(conn, stored):
    """Raises `SchemaVersionMissingError` for an existing database with no recorded version."""
    if stored is None and _has_galaxy_tables(conn):
        database = conn._config.database if getattr(conn, "_config", None) is not None else "?"
        raise SchemaVersionMissingError(database)


def _refuse_older(conn, version):
    """Raises `SchemaTooOldError` when `version` is below the baseline."""
    if version is not None and version < alembic_runner.BASELINE_VERSION:
        database = conn._config.database if getattr(conn, "_config", None) is not None else "?"
        raise SchemaTooOldError(database, version)


def _refuse_newer(conn, version, expected, what="schema"):
    """Raises `SchemaTooNewError` when `version` is above `expected`."""
    if version is not None and version > expected:
        database = conn._config.database if getattr(conn, "_config", None) is not None else "?"
        raise SchemaTooNewError(database, version, expected, what)

_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

SCHEMA_PATH = os.path.join(_PACKAGE_DIR, "schema.sql")
"""str: Path to the DDL file applied by `_ensure_schema`."""

NAMED_LOCK_TIMEOUT_S = 50
"""int: How long `Connection.lock_until_commit` waits for a named lock,
matching InnoDB's default `innodb_lock_wait_timeout`."""

LOCK_HOLDER_ROW_WAIT_S = 3
"""int: How long a named lock's holder waits for a row lock
(`innodb_lock_wait_timeout`) before giving up, rolling back and letting
`save_sector` start it over (PERF.21). A session waiting for the named
lock still holds the rows it wrote first, and InnoDB can't see a wait on
a named lock, so when the holder needs one of those rows (a name
collision renames an existing system, which the holder's nearest-neighbor
rows then reference) neither would move until InnoDB's own 50 s timeout.
Giving up early breaks that wait at once."""

CONTROL_SCHEMA_VERSION = 21
"""int: Version counter for `control_schema.sql`, independent of
`SCHEMA_VERSION` above -- see that file's header comment for why the
control plane (admin identities/sessions/API keys/audit log) is a
separate schema with its own versioning. v2 added `login_throttle`
(SEC.1, SEC.21; dropped in v8, SEC.30: the counters live in Redis), v5 the work queue's `work_jobs`/`work_tasks`/
`work_lease` (PERF.8), v6 `generation_stats`/`generation_size` (PERF.3,
PERF.10), v7 the job tree's columns on `work_jobs` and the queue pause
on `work_lease` (ADM.12, ADM.10), v8 the drop of `login_throttle`
(SEC.30), v9 `galaxy_naming` (GEN.70), v10 `work_job_args` and the drop
of the JSON columns `work_jobs.argv`, `work_tasks.result` (DB.13), v11 `version_key_history` (OPS.13), v12 `generation_stats` keyed by worker count and
version key (PERF.32), v17 `generate_recipes` and `generate_recipe_fields` (ADM.125), v18 the epoch, battery digest, environment and last-seen columns of `version_key_history` (OPS.58), v19 `admin_sessions.password_confirmed_at` (USR.31), v20 `admin_users.avatar` (SEC.36), v21 the drop of the work queue's four tables, `work_jobs`, `work_job_args`, `work_tasks` and `work_lease` (PERF.100). New tables need nothing more than
`CREATE TABLE IF NOT EXISTS`; new columns on an existing table are
added by `_add_control_columns`."""

_CONTROL_COLUMNS = {
    "admin_api_keys": [
        ("key_prefix", "CHAR(8) NULL"),
        ("expires_at", "TIMESTAMP NULL"),
    ],
    "version_key_history": [
        ("generator_epoch", "SMALLINT UNSIGNED NULL"),
        ("battery_digest", "CHAR(64) NULL"),
        ("battery_cases", "JSON NULL"),
        ("environment", "JSON NULL"),
        ("last_seen_at", "DATETIME(6) NULL"),
    ],
    "admin_sessions": [
        ("password_confirmed_at", "DATETIME NULL"),
    ],
    "admin_users": [
        ("email", "VARCHAR(254) NULL"),
        ("role", "VARCHAR(16) NOT NULL DEFAULT 'admin'"),
        ("status", "VARCHAR(16) NOT NULL DEFAULT 'active'"),
        ("owner_flag", "TINYINT AS (CASE WHEN role = 'owner' THEN 1 END) VIRTUAL"),
        ("avatar", "VARCHAR(40) NULL"),
    ],
}
"""dict: Columns added to existing control tables after they were first
created (v7), as `table -> [(column, definition)]`, in order. Must match
`control_schema.sql`'s `CREATE TABLE` statements."""

_CONTROL_CONSTRAINTS = (
    ("admin_users", "uq_admin_users_email", "ADD UNIQUE KEY uq_admin_users_email (email)"),
    ("admin_users", "uq_admin_users_owner", "ADD UNIQUE KEY uq_admin_users_owner (owner_flag)"),
    ("admin_users", "chk_admin_users_role",
     "ADD CONSTRAINT chk_admin_users_role CHECK (role IN ('user', 'moderator', 'admin', 'owner'))"),
    ("admin_users", "chk_admin_users_status",
     "ADD CONSTRAINT chk_admin_users_status CHECK (status IN ('active', 'disabled'))"),
)
"""tuple: `(table, name, ALTER clause)`: the keys and CHECKs added to an existing control table (v15), each only when
missing. Must match `control_schema.sql`'s `CREATE TABLE` statements."""

CONTROL_SCHEMA_PATH = os.path.join(_PACKAGE_DIR, "control_schema.sql")
"""str: Path to the DDL file applied by `_ensure_control_schema`."""

CONTROL_DB_ENV_VAR = env_name("control_database")
"""str: Env var naming the one MySQL schema the control plane lives in
(admin identities are global to a deployment, not per-galaxy -- see
`control_schema.sql`'s header comment); the settings model's
`control_database` (`planetgen.util.settings`)."""


def _stored_version(conn, table):
    """`MAX(version)` of a migrations table (`schema_migrations` or
    `control_schema_migrations`), or `None` while it doesn't exist yet or
    is empty -- read before any DDL runs, so a newer database is refused
    before this code's older `schema.sql` (its `CREATE OR REPLACE VIEW`s
    among them) touches it."""
    row = conn.execute(
        "SELECT COUNT(*) AS n FROM information_schema.tables WHERE table_schema = DATABASE() AND table_name = ?",
        (table,),
    ).fetchone()
    if not row["n"]:
        return None
    return conn.execute(f"SELECT MAX(version) AS v FROM {table}").fetchone()["v"]


class MySQLConfig:
    """
    MySQL connection parameters, read from environment variables --
    mirrors every other entry point in this project (`generation/run_sector.py`,
    `generation/run_system.py`, `planetgen.db.query`, `planetgen/api/config.py`) reading its own
    `PLANETGEN_*` variable rather than hardcoding a value, so a deployment
    points every tool at the same server via its process environment
    (e.g. the Apache vhost's `SetEnv`, or a `systemd`/`gunicorn` unit's
    environment file) without editing code. No password default -- unlike
    host/port/user/database, a blank password is a real (if unusual)
    credential, not an obviously-safe placeholder, so getting it wrong by
    omission fails loudly at connect time instead of silently connecting
    as some other account.

    A caller needing a different database than the process-wide default
    (every test in this project's own suite included -- see
    `src/tests/conftest.py`) builds its own `MySQLConfig` instance
    directly rather than going through environment variables at all.

    Precedence for each field left unset here is: the matching
    `PLANETGEN_MYSQL_*` environment variable, then `config.json`'s
    `mysql` section (see `planetgen.util.settings`), then the
    hardcoded default below.
    """

    def __init__(self, host=None, port=None, user=None, password=None, database=None):
        defaults = get_settings().mysql
        self.host = host if host is not None else defaults.host
        self.port = int(port if port is not None else defaults.port)
        self.user = user if user is not None else defaults.user
        self.password = password if password is not None else defaults.password
        self.database = database if database is not None else defaults.database

    def _key(self):
        """A hashable identity for this config, used to key the pool
        cache below -- two `MySQLConfig` instances with the same
        connection parameters should share one pool rather than each
        opening their own."""
        return (self.host, self.port, self.user, self.password, self.database)


DEFAULT_MYSQL_CONFIG = MySQLConfig()
"""MySQLConfig: The process-wide default, built from `PLANETGEN_MYSQL_*`
env vars (or their defaults) at import time. Every entry point that
doesn't need a different database (i.e. everything except this project's
own test suite) uses this implicitly by passing `config=None` through to
`get_connection`."""


def _mysql_port(text):
    """argparse `type` for `--mysql-port`: a whole number from 1 to 65535,
    so an impossible port is a usage error up front instead of a
    connection error later (OPS.6)."""
    try:
        port = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"invalid int value: {text!r}")
    if not 1 <= port <= 65535:
        raise argparse.ArgumentTypeError("must be between 1 and 65535")
    return port


def add_mysql_connection_args(parser):
    """
    Adds the `--mysql-host`/`--mysql-port`/`--mysql-user`/
    `--mysql-password`/`--mysql-database` optional overrides to `parser`,
    shared by every CLI entry point in this project (`generation/run_sector.py`,
    `generation/run_system.py`, `generation/run_galaxy.py`, `planetgen.db.query`, `planetgen.cli.migrate`) instead
    of each one re-declaring the same five arguments -- pair with
    `mysql_config_from_args` to turn the parsed result into a
    `MySQLConfig`.

    Every flag defaults to `None` (falls through to `MySQLConfig`'s own
    `PLANETGEN_MYSQL_*` env var default) rather than duplicating that
    default in the `--help` text here, which would drift out of sync with
    `MySQLConfig.__init__`'s actual defaults over time.

    Args:
        parser (argparse.ArgumentParser or argparse._ActionsContainer):
            The parser (or subparser) to add the arguments to.
    """
    parser.add_argument('--mysql-host', type=str,
                         help="MySQL host. Defaults to $PLANETGEN_MYSQL_HOST, or 127.0.0.1.")
    parser.add_argument('--mysql-port', type=_mysql_port,
                         help="MySQL port. Defaults to $PLANETGEN_MYSQL_PORT, or 3306.")
    parser.add_argument('--mysql-user', type=str,
                         help="MySQL user. Defaults to $PLANETGEN_MYSQL_USER, or 'planetgen'.")
    parser.add_argument('--mysql-password', type=str,
                         help="MySQL password. Defaults to $PLANETGEN_MYSQL_PASSWORD, or empty.")
    parser.add_argument('--mysql-database', type=str,
                         help="MySQL database name. Defaults to $PLANETGEN_MYSQL_DATABASE, or 'planetgen'.")


def mysql_config_from_args(args) -> MySQLConfig:
    """
    Builds a `MySQLConfig` from the `--mysql-*` arguments
    `add_mysql_connection_args` adds -- any flag left unset (`None`) falls
    through to `MySQLConfig`'s own env-var/hardcoded default, exactly like
    omitting that argument to `MySQLConfig()` directly.

    Args:
        args (argparse.Namespace): Parsed arguments, from a parser that
            called `add_mysql_connection_args`.

    Returns:
        MySQLConfig: Ready to pass as `get_connection`/`save_sector`/
            `save_system`/`migrate_database`'s `config` argument.
    """
    return MySQLConfig(
        host=args.mysql_host, port=args.mysql_port, user=args.mysql_user,
        password=args.mysql_password, database=args.mysql_database,
    )

_pools = {}
"""dict: `MySQLConfig._key() -> sqlalchemy `QueuePool`, one pool per distinct set of
connection parameters seen so far in this process. A WSGI worker or a
single CLI invocation only ever needs one (the process-wide default), but
keying by config rather than keeping a single module-level pool lets
tests (and any future multi-database use) point at a second database
within the same process without the two stomping on each other's pooled
connections."""
_pools_lock = threading.Lock()


STATEMENT_TIMEOUT_ERRORS = (3024,)
"""tuple: The errors a statement stopped by `statement_timeout_s` raises: MySQL's `MAX_EXECUTION_TIME`."""

LOCK_WAIT_ERRORS = (1205, 1213)
"""tuple: The errors of a statement that waited too long for a row lock (1205) or lost a deadlock (1213): the
transaction holding the lock is somebody else's, and a short web write gives up on it rather than wait (ADM.53)."""

def _statement_timeout_sql(seconds):
    return f"SESSION MAX_EXECUTION_TIME = {int(float(seconds) * 1000)}"


SQL_MODE_ENV_VAR = "PLANETGEN_MYSQL_SQL_MODE"
"""str: Env var that, when set, pins every pooled connection's session
`sql_mode` to its value (TEST.7). The test suite sets it to MySQL 8's
default (`STRICT_SQL_MODE`), so a GROUP BY or a truncation a looser
`sql_mode` forgives fails in the tests too. Unset (every deployment): the server's own `sql_mode` applies."""

STRICT_SQL_MODE = ("ONLY_FULL_GROUP_BY,STRICT_TRANS_TABLES,NO_ZERO_IN_DATE,NO_ZERO_DATE,"
                   "ERROR_FOR_DIVISION_BY_ZERO,NO_ENGINE_SUBSTITUTION")
"""str: MySQL 8.0's default `sql_mode`."""


def _session_sql_mode():
    mode = os.environ.get(SQL_MODE_ENV_VAR, "").strip()
    if mode and not re.fullmatch(r"[A-Za-z_,]+", mode):
        raise ValueError(f"{SQL_MODE_ENV_VAR} must be a comma-separated list of sql_mode names, not {mode!r}")
    return mode


def _get_pool(config, statement_timeout_s=None, lock_wait_s=None):
    key = config._key() if not (statement_timeout_s or lock_wait_s) else (
        config._key(), float(statement_timeout_s or 0), float(lock_wait_s or 0))
    pool = _pools.get(key)
    if pool is not None:
        return pool
    # One pool per key even when several threads are the first to ask (TEST.218): each would otherwise build its own,
    # and all but the last stored would be dropped with a connection open on the server.
    with _pools_lock:
        if key not in _pools:
            init = "SET time_zone = '+00:00'"
            if statement_timeout_s:
                init += ", " + _statement_timeout_sql(statement_timeout_s)
            if lock_wait_s:
                init += f", SESSION innodb_lock_wait_timeout = {max(int(lock_wait_s), 1)}"
            sql_mode = _session_sql_mode()
            if sql_mode:
                init += f", SESSION sql_mode = '{sql_mode}'"
            connect_args = dict(
                host=config.host,
                port=config.port,
                user=config.user,
                password=config.password,
                database=config.database,
                charset="utf8mb4",
                cursorclass=pymysql.cursors.DictCursor,
                autocommit=False,
                # TIMESTAMP columns are stored in UTC but read back in the
                # session's zone; pinning it to UTC makes every time read the
                # same whatever the server's own zone is (the pages convert to
                # each viewer's zone in the browser).
                # PERF.17: a web pool also stops any statement that runs
                # past `statement_timeout_s`.
                init_command=init,
            )
            # Never fail over to a fresh connection mid-use: SQLAlchemy's pool
            # hands out one DBAPI connection and never re-runs a statement on
            # another, so a deadlock leaves the caller to see it (the error
            # 1452s of PERF.14's parallel test came from a silent re-run in a
            # transaction MySQL had already rolled back). A dead pooled
            # connection is still replaced when it's checked out
            # (`_ping_on_checkout`); one lost mid-use raises.
            # A bare pool, not an Engine: an Engine's dialect would probe the
            # server through these DictCursor connections and read rows by
            # position.
            pool = sqlalchemy.pool.QueuePool(
                lambda: login.connect(**connect_args),
                pool_size=10,
                max_overflow=10,
            )
            sqlalchemy.event.listen(pool, "checkout", _ping_on_checkout)
            # Open one connection now, so an unreachable server or a wrong
            # password fails here rather than at the first query.
            pool.connect().close()
            _pools[key] = pool
        return _pools[key]


def _ping_on_checkout(dbapi_connection, _record, _proxy):
    """Replaces a pooled connection the server has dropped, before the
    caller gets it."""
    try:
        dbapi_connection.ping(reconnect=False)
    except pymysql.err.Error:
        raise sqlalchemy.exc.DisconnectionError() from None


def close_pool(config):
    """
    Closes and discards the cached engine for `config`, if one exists.

    `_get_pool` never evicts an entry from `_pools` on its own -- fine for
    the handful of long-lived databases a real deployment or WSGI worker
    ever points at, but a caller that opens many short-lived, uniquely
    named databases in one process (this project's own `mysql_config` test
    fixture, one throwaway schema per test) would otherwise leave every
    prior pool's `mincached` connections open for the rest of the
    process's life, eventually exhausting the server's `max_connections`.
    Call this once such a database is being dropped for good.
    """
    key = config._key()
    for pool_key in [k for k in _pools if k == key or (isinstance(k, tuple) and k and k[0] == key)]:
        _pools.pop(pool_key).dispose()
    _schema_ensured.discard(key)
    forget_id_blocks(key)


class _Cursor:
    """
    Thin proxy around a real `pymysql` cursor that normalizes
    `fetchall()` to always return a `list` -- confirmed by testing,
    `pymysql`'s own cursor returns `()` (a tuple) for zero matching rows
    but a `list` when rows exist, an inconsistency `sqlite3`'s cursor
    (always a `list`, regardless of row count) never had, and that a
    caller comparing a query result against `[]` (rather than checking
    `len(...)` or truthiness) would otherwise trip over. Every other
    attribute (`.lastrowid`, `.fetchone()`, ...) is forwarded to the real
    cursor unchanged.
    """

    def __init__(self, cursor):
        self._cursor = cursor

    def fetchall(self):
        return list(self._cursor.fetchall())

    def __getattr__(self, name):
        return getattr(self._cursor, name)


_SECRET_SQL = re.compile(r"password|token|key_hash|secret", re.IGNORECASE)
_SQL_LOG_LIMIT = 1000


def _sql_for_log(sql, params):
    """One-line SQL plus its parameters for the debug log, with the
    parameters withheld on statements touching credentials."""
    text = " ".join(sql.split())
    shown = "<withheld: credentials>" if _SECRET_SQL.search(text) else repr(params)
    line = f"{text} | params={shown}"
    return line if len(line) <= _SQL_LOG_LIMIT else line[:_SQL_LOG_LIMIT] + f"... ({len(line)} chars)"


ID_BLOCK_TABLES = frozenset(kind_tables.tables("id_block"))
"""frozenset: Tables whose new rows take their `id` from `id_blocks`
(schema v45, PERF.13) rather than from AUTO_INCREMENT. Every plain
`INSERT INTO <one of these> (...) VALUES (...)` that names no `id`
column gets one added by `Connection.execute`, so a sector's children
know their parents' ids before anything is written, and `cur.lastrowid`
still answers as before."""

BATCH_CHILD_TABLES = frozenset((
    "system_config_slots", "planet_reflection_spectrum", "moon_reflection_spectrum", "planet_interior_layers", "moon_interior_layers",
    "planet_gases", "moon_gases", "planet_life_class_scores", "moon_life_class_scores", "planet_life", "moon_life",
    "asteroid_belt_composition",
    "comet_composition", "interstellar_comet_composition", "asteroid_field_composition", "orbit_nodes", "dark_members",
))
"""frozenset: Child tables nothing refers to by id, so `Connection.batched`
may hold their INSERTs back and write them many rows at a time with no
`id_blocks` involvement (their own AUTO_INCREMENT ids are never read)."""

_ID_BLOCK_MIN = 64
_ID_BLOCK_MAX = 4096
_BATCH_ROWS = 500
"""int: Most rows per multi-row INSERT when `Connection.flush` writes a
batch; fewer when their bytes would pass `max_allowed_packet`."""

_PACKET_MARGIN = 1024
"""int: Bytes of `max_allowed_packet` `Connection.flush` leaves for the
protocol around a statement."""

_INSERT_RE = re.compile(
    r"^\s*INSERT\s+INTO\s+`?(\w+)`?\s*\(([^()]*)\)\s*VALUES\s*(\(.*\))\s*;?\s*$",
    re.IGNORECASE | re.DOTALL,
)
_InsertShape = namedtuple("_InsertShape", "table columns values has_id")
_insert_shapes = {}

_schema_ensured = set()
"""set: `MySQLConfig._key()`s this process already ran `_ensure_schema`
against (PERF.12): replaying `schema.sql` on every checkout cost about
122 statements per connection, twice per filled sector."""

_id_blocks = {}
_FLAGGING_TABLES = frozenset(storage_checksum.CONTENT_TABLES) - {"sectors"}
"""frozenset: The tables an INSERT into marks a sector's `parity_dirty` (`Connection._flag_insert`, DB.45)."""

_id_blocks_off = set()
_id_lock = threading.Lock()
_fk_ranks = {}
_fk_same_rank_parents = {}
_max_packets = {}


def _insert_shape(sql):
    """The parsed shape of a plain single-row `INSERT ... VALUES (...)`,
    or `None` for anything else (an upsert, `INSERT ... SELECT`, a
    statement that isn't an INSERT). Cached per SQL string -- every
    query in this module is a literal."""
    shape = _insert_shapes.get(sql)
    if shape is None and sql not in _insert_shapes:
        match = _INSERT_RE.match(sql)
        if match and not re.search(r"\bON\s+DUPLICATE\b|\bSELECT\b", sql, re.IGNORECASE):
            columns = [c.strip().strip("`") for c in match.group(2).split(",")]
            shape = _InsertShape(match.group(1).lower(), columns, match.group(3),
                                 any(c.lower() == "id" for c in columns))
        _insert_shapes[sql] = shape
    return shape


def forget_id_blocks(key):
    """Drops this process's cached id blocks for one database
    (`MySQLConfig._key()`), so the next id is reserved afresh -- after
    `planetgen.cli.reset` empties the tables (it keeps `id_blocks`, DB.3), or a
    migration adds `id_blocks`."""
    with _id_lock:
        for block_key in [k for k in _id_blocks if k[0] == key]:
            del _id_blocks[block_key]
        _id_blocks_off.discard(key)


def _allocate_id(config, table):
    """
    The next `id` for a new `table` row (PERF.13), from this process's
    current block of ids, fetching a new block from `id_blocks` when it
    runs out. A block is reserved on its own short autocommitted
    connection (`UPDATE ... SET next_id = LAST_INSERT_ID(next_id + n)`),
    so writers never wait on each other's sector transactions, and a
    rolled-back sector only leaves a gap, as AUTO_INCREMENT would. Each
    reservation starts no lower than the table's current `MAX(id) + 1`,
    so rows written before v45 (or with an explicit id) are never
    reused. Blocks start at 64 ids and double up to 4,096 per table.

    Returns:
        int or None: The id, or `None` when this database has no
            `id_blocks` table yet (not migrated to v45); the INSERT then
            falls back to AUTO_INCREMENT.
    """
    key = config._key()
    with _id_lock:
        if key in _id_blocks_off:
            return None
        block = _id_blocks.get((key, table))
        if block is None or block[0] >= block[1]:
            size = _ID_BLOCK_MIN if block is None else min(_ID_BLOCK_MAX, block[2] * 2)
            try:
                start = _reserve_id_block(config, table, size, 1 if block is None else block[1])
            except pymysql.err.ProgrammingError as exc:
                if exc.args and exc.args[0] == 1146:  # no id_blocks table yet
                    _id_blocks_off.add(key)
                    return None
                raise
            block = _id_blocks[(key, table)] = [start, start + size, size]
        block[0] += 1
        return block[0] - 1


def _reserve_id_block(config, table, size, at_least):
    """Reserves `size` ids for `table` and returns the first. `at_least`
    is past every id this process already handed out for it, which the
    table's `MAX(id)` can't show while those rows are uncommitted."""
    raw = _get_pool(config).connect()
    try:
        cur = raw.cursor()
        cur.execute(f"SELECT COALESCE(MAX(id), 0) + 1 AS floor_id FROM {table}")
        floor_id = max(cur.fetchone()["floor_id"], at_least)
        # One upsert, so the row lock is taken exclusive from the start:
        # INSERT IGNORE then UPDATE had two processes each hold a shared
        # lock on a new table's row and deadlock (1213) upgrading it
        # (TEST.81).
        cur.execute(
            "INSERT INTO id_blocks (table_name, next_id) VALUES (%s, LAST_INSERT_ID(%s + %s))"
            " ON DUPLICATE KEY UPDATE next_id = LAST_INSERT_ID(GREATEST(next_id, %s) + %s)",
            (table, floor_id, size, floor_id, size),
        )
        cur.execute("SELECT LAST_INSERT_ID() AS end_id")
        end = cur.fetchone()["end_id"]
        raw.commit()
    except Exception:
        raw.rollback()
        raise
    finally:
        raw.close()
    return end - size


def _table_ranks(conn, key):
    """
    `{table: rank}` with every table ranked after the tables its foreign
    keys point at, so `Connection.flush` writes parents before children
    (cached per database for the life of the process). Tables that point
    at each other in a circle (star systems, stars, black holes and
    supernova remnants, through the v39 containment columns) share a
    rank; `_same_rank_parents` keeps their rows in insertion order (see
    `Connection._batch_level`).
    """
    ranks = _fk_ranks.get(key)
    if ranks is None or key not in _fk_same_rank_parents:
        cur = conn.cursor()
        cur.execute(
            "SELECT TABLE_NAME AS child, REFERENCED_TABLE_NAME AS parent FROM information_schema.KEY_COLUMN_USAGE"
            " WHERE TABLE_SCHEMA = DATABASE() AND REFERENCED_TABLE_NAME IS NOT NULL"
        )
        parents = {}
        for row in cur.fetchall():
            child, parent = row["child"].lower(), row["parent"].lower()
            parents.setdefault(child, set())
            parents.setdefault(parent, set())
            if child != parent:
                parents[child].add(parent)
        ranks = _condensed_ranks(parents)
        _fk_same_rank_parents[key] = {
            table: {parent for parent in above if ranks[parent] == ranks[table]} for table, above in parents.items()
        }
        _fk_ranks[key] = ranks
    return ranks


def _same_rank_parents(conn, key):
    """`{table: the tables it refers to that share its rank}` -- non-empty
    only for the FK cycle's tables (`_table_ranks`)."""
    _table_ranks(conn, key)
    return _fk_same_rank_parents[key]


def _max_packet(connection):
    """The server's `max_allowed_packet` for a `Connection` (cached per
    database), the most bytes one statement may take."""
    key = connection._config._key() if connection._config is not None else None
    packet = _max_packets.get(key)
    if packet is None:
        cur = connection._conn.cursor()
        cur.execute("SELECT @@max_allowed_packet AS n")
        packet = int(cur.fetchone()["n"])
        if key is not None:
            _max_packets[key] = packet
    return packet


def _condensed_ranks(parents):
    """Ranks for `_table_ranks`: Tarjan's strongly connected components of
    the child -> parent graph, each ranked one past its highest parent
    component."""
    index, low, stack, on_stack, component = {}, {}, [], set(), {}
    counter = [0]

    def visit(node):
        index[node] = low[node] = counter[0]
        counter[0] += 1
        stack.append(node)
        on_stack.add(node)
        for parent in parents[node]:
            if parent not in index:
                visit(parent)
                low[node] = min(low[node], low[parent])
            elif parent in on_stack:
                low[node] = min(low[node], index[parent])
        if low[node] == index[node]:
            members = []
            while True:
                member = stack.pop()
                on_stack.discard(member)
                members.append(member)
                if member == node:
                    break
            for member in members:
                component[member] = node

    for node in sorted(parents):
        if node not in index:
            visit(node)
    component_ranks = {}

    def rank(root):
        if root not in component_ranks:
            component_ranks[root] = 0
            members = [table for table, owner in component.items() if owner == root]
            above = {component[parent] for table in members for parent in parents[table]} - {root}
            component_ranks[root] = 1 + max((rank(other) for other in above), default=-1)
        return component_ranks[root]

    return {table: rank(component[table]) for table in parents}


_NEEDS_ESCAPING = re.compile("[\\\\'\"\0\n\r\x1a]").search
"""The characters pymysql escapes in a string (`converters._escape_table`); a string with none of them is quoted as it is."""

_ALL_MARKS = re.compile(r"\((\?, )*\?\)")
"""A `VALUES (...)` shape made only of `?` marks, which `_row_literal` can write without splitting it."""


def _row_literal(escape, pieces, row):
    """
    One row of a batched INSERT as SQL text: the statement's `VALUES (?, ?, 0)` shape split at its `?` marks
    (`pieces`; `None` for a shape that is only marks, `(?, ?, ?)`) with each value written out the way pymysql writes
    it (`Connection.escape`), so the text is the one the driver would have made. Plain ints, floats and `None`, and
    strings with nothing to escape, are formatted here without the driver's per-value dispatch, which was most of the
    time of an insert (PERF.82); everything else (strings with quotes, bools, numpy numbers, dates) goes through
    `escape` itself.
    """
    parts = []
    append = parts.append
    for value in row:
        kind = type(value)
        if kind is float:
            text = repr(value)
            if "e" not in text:
                if text in ("inf", "-inf", "nan"):
                    raise pymysql.err.ProgrammingError(f"{text} can not be used with MySQL")
                text += "e0"
            append(text)
        elif kind is int:
            append(str(value))
        elif kind is str:
            append(escape(value) if _NEEDS_ESCAPING(value) else "'" + value + "'")
        elif value is None:
            append("NULL")
        else:
            append(escape(value))
    if pieces is None:
        return "(" + ", ".join(parts) + ")"
    out = [pieces[0]]
    for piece, part in zip(pieces[1:], parts):
        out.append(part)
        out.append(piece)
    return "".join(out)


class Connection:
    """
    Thin wrapper around a pooled `pymysql` connection that keeps this
    module's (and `planetgen.db.query`'s) existing `conn.execute(sql, params)` /
    `cur.lastrowid` / `row["column"]` call sites working unchanged after
    the SQLite -> MySQL port, instead of touching every one of their SQL
    strings and call sites individually:

      - `execute` rewrites `?` positional placeholders to `pymysql`'s
        `%s` before handing the query to a real DB-API cursor (every
        query in this codebase is a literal, module-level string -- never
        built from request input -- so this rewrite is safe; no query
        here contains a literal `?` character it wasn't meant as a
        placeholder, or a literal `%` in the query template itself).
      - Rows come back from `pymysql.cursors.DictCursor` already
        supporting `row["column"]`, matching `sqlite3.Row` -- no
        additional translation needed there.
      - The returned cursor is wrapped in `_Cursor` so `.fetchall()`
        always returns a `list` -- see that class's own docstring.
      - `with conn:` commits on a clean exit and rolls back on an
        exception, same as `sqlite3.Connection`'s context-manager
        behavior -- and, same as `sqlite3.Connection`, does NOT close the
        connection either way.
      - An INSERT into one of `ID_BLOCK_TABLES` gets its id from
        `_allocate_id` (PERF.13), and inside `with conn.batched():` such
        INSERTs (and those into `BATCH_CHILD_TABLES`) are held back and
        written as multi-row INSERTs, one per table and statement shape,
        the next time anything else runs or the batch ends.
    """

    def __init__(self, pooled_conn, config=None):
        self._conn = pooled_conn
        self._config = config
        self._batch = None
        self._batch_depth = 0
        self._batch_levels = {}
        self.prereserved_names = None
        self.deferred_name_confirmations = None
        self.claim_connection = None
        """`Connection` or `None`: a second connection of the same database on which `reserve_system_names` claims
        names in their own short transactions (PERF.49), so the registry's row locks are not held to this
        transaction's commit. Set by `_save_with_retries`."""
        self.name_claims = []
        """list: `(base_name, count, count_after)` per claim made on `claim_connection`, for
        `_release_name_claims` if this transaction fails."""
        self.later_name_keys = set()
        """set: Registry keys this transaction claimed after another holder had; it is not their first holder."""
        self.burned_name_keys = {}
        """dict: Registry key -> claims of this transaction's own that a redrawn name left counted."""
        self.sector_name_keys = set()
        """set: Registry keys it claimed that an existing sector also uses as its base name."""
        self._txn_locks = []
        self._short_row_waits = False
        self.uid_issuer = None
        """`_UidIssuer` or `None`: while a sector is being saved, gives each new row its unique ID in the INSERT
        (PERF.44)."""
        self.sector_centers_pc = {}
        """dict: `sector_id -> (x, y, z)` or `None`, the sector centers `_system_velocity_kms` has read inside the
        current `batched()` scope (a read flushes the held-back INSERTs, so each is read once)."""
        self.orbit_epoch_unix = None
        """int or None: the orbit epoch a binary's tree is stored at, read once inside a `batched()` scope (GEN.311)."""
        self._after_commit = []
        """list: Callables run once the transaction has committed (`after_commit`); dropped by a rollback."""
        self._parity_touched = set()
        """set: The sector ids this transaction marked `parity_dirty`, for the parity file to catch up after the commit
        (`parity_touch`); empty when `parity.enabled` is off."""
        self._flagged = set()
        """set: The `(kind, id)` of sectors and systems this transaction has already marked `parity_dirty` for an INSERT
        (DB.45), and the systems it inserted itself (their sector was marked by the INSERT); cleared at commit and rollback."""

    def _run(self, sql, params):
        cur = self._conn.cursor()
        traced = log.debug_log_active()
        start = time.perf_counter()
        try:
            cur.execute(sql.replace("?", "%s"), params)
        except Exception as exc:
            if isinstance(exc, pymysql.err.OperationalError) and exc.args and exc.args[0] in STATEMENT_TIMEOUT_ERRORS:
                # PERF.64: the page only says "took too long"; the log says which statement.
                activity_log.event("DB", "statement_timeout", db=self._config.database if self._config else "?",
                                   seconds=f"{time.perf_counter() - start:.1f}", sql=_sql_for_log(sql, params)[:800])
            if traced:
                log.trace(f"SQL failed after {(time.perf_counter() - start) * 1000:.2f}ms: "
                          f"{_sql_for_log(sql, params)} -> {type(exc).__name__}: {exc}", area="sql", stacklevel=4)
            raise
        if traced:
            log.trace(f"SQL {(time.perf_counter() - start) * 1000:.2f}ms, {cur.rowcount} row(s): "
                      f"{_sql_for_log(sql, params)}", area="sql", stacklevel=4)
        return _Cursor(cur)

    @contextlib.contextmanager
    def statement_limit(self, seconds):
        """Lowers this connection's statement time limit to `seconds` for the block (PERF.76), restoring the
        limit it had; never raises a limit already lower, and does nothing for a connection with no server."""
        if self._config is None:
            yield
            return
        cur = self._conn.cursor()
        cur.execute("SELECT @@SESSION.max_execution_time")
        row = cur.fetchone()
        current = float((list(row.values())[0] if isinstance(row, dict) else row[0]) or 0) / 1000
        if current and current <= seconds:
            yield
            return
        cur.execute("SET " + _statement_timeout_sql(seconds))
        try:
            yield
        finally:
            cur.execute("SET " + (_statement_timeout_sql(current) if current else "SESSION max_execution_time = 0"))

    def execute(self, sql, params=()):
        shape = _insert_shape(sql)
        if shape is not None:
            new_id = None
            if shape.table in ID_BLOCK_TABLES and not shape.has_id and self._config is not None:
                new_id = _allocate_id(self._config, shape.table)
            if new_id is not None:
                shape = shape._replace(columns=["id", *shape.columns], values="(?, " + shape.values[1:], has_id=True)
                params = (new_id, *params)
            if self.uid_issuer is not None and shape.table in _UID_ISSUED_TABLES:
                uid = self.uid_issuer.issue(shape.table, shape.columns, params)
                if uid is not None:
                    shape = shape._replace(columns=[*shape.columns, "uid"], values=shape.values[:-1] + ", ?)")
                    params = (*params, uid)
            self._flag_insert(shape, params, new_id)
            if self._batch is not None and (new_id is not None or shape.table in BATCH_CHILD_TABLES):
                level = self._batch_level(shape.table)
                self._batch.setdefault((shape.table, tuple(shape.columns), shape.values, level), []).append(
                    tuple(params))
                return _InsertedCursor(new_id)
            self.flush()
            if new_id is not None:
                cur = self._run(f"INSERT INTO {shape.table} ({', '.join(shape.columns)}) VALUES {shape.values}", params)
                return _InsertedCursor(new_id, cur)
            return self._run(sql, params)
        self.flush()
        return self._run(sql, params)

    def _flag_insert(self, shape, params, new_id):
        """DB.45: an INSERT of a row a sector's storage checksum covers marks that sector's `parity_dirty` (once per sector
        or system per transaction). A row filed by `sector_id` marks that sector; one filed by `star_system_id` marks the
        system's sector, unless this transaction inserted the system (its own INSERT marked the sector)."""
        table = shape.table
        if table not in _FLAGGING_TABLES:
            return
        columns = shape.columns
        if table == "star_systems" and new_id is not None:
            self._flagged.add(("fresh", new_id))
        for kind, column in (("sector", "sector_id"), ("system", "star_system_id")):
            if column in columns and params[columns.index(column)] is not None:
                key = (kind, params[columns.index(column)])
                if key in self._flagged or (kind == "system" and ("fresh", key[1]) in self._flagged):
                    continue
                self._flagged.add(key)
                self._run(storage_checksum.MARK_SECTOR_SQL if kind == "sector" else storage_checksum.MARK_SYSTEM_SQL, (key[1],))
                self.parity_touch("id = ?" if kind == "sector" else "id = (SELECT sector_id FROM star_systems WHERE id = ?)",
                                  (key[1],))

    def batched(self):
        """
        `with conn.batched():` -- holds back plain INSERTs into
        `ID_BLOCK_TABLES` and `BATCH_CHILD_TABLES` and writes them as
        multi-row INSERTs (PERF.13). Any other statement, a commit, or the
        end of the block writes what's held first, so reads inside the
        block still see every row added before them. Errors from a held
        row surface at that write, not at its `execute` call. Nests.
        """
        return _BatchScope(self)

    def _batch_level(self, table):
        """
        The round of `flush`'s writes a newly held `table` row goes in:
        after every held row of a table it refers to that shares its rank
        (the FK cycle, `_table_ranks`), so those rows keep their insertion
        order -- a held star never goes ahead of the held system it points
        at -- while each table's rows still share as few statements as
        that order allows (every system in round 0, their stars in 1).
        Always 0 for any other table.
        """
        if not self._batch:
            self._batch_levels = {}
        levels = self._batch_levels
        level = levels.get(table, 0)
        key = self._config._key() if self._config is not None else None
        for parent in _same_rank_parents(self._conn, key).get(table, ()):
            if parent in levels:
                level = max(level, levels[parent] + 1)
        levels[table] = level
        return level

    def flush(self):
        """
        Writes every held-back INSERT, parents before children (by
        `_table_ranks`, then `_batch_level`), each statement as many rows
        as fit both `_BATCH_ROWS` and the server's `max_allowed_packet`.

        Raises:
            pymysql.err.OperationalError: 1153 for a single row too big to
                send, before anything of it reaches the server (which would
                drop the connection instead).
        """
        batch = self._batch
        if not batch:
            return
        self._batch = {}
        key = self._config._key() if self._config is not None else None
        ranks = _table_ranks(self._conn, key)
        order = sorted(enumerate(batch.items()),
                       key=lambda item: (ranks.get(item[1][0][0], 0), item[1][0][3], item[0]))
        budget = _max_packet(self) - _PACKET_MARGIN
        cur = self._conn.cursor()
        escape = self._conn.escape
        start = time.perf_counter()
        statements = 0
        for _position, ((table, columns, values, _level), rows) in order:
            head = f"INSERT INTO {table} ({', '.join(columns)}) VALUES "
            pieces = None if _ALL_MARKS.fullmatch(values) else values.split("?")
            chunk, size = [], len(head)
            for row in rows:
                text = _row_literal(escape, pieces, row)
                # The statement's bytes: a character is 1 byte unless the text is not plain ASCII (the names).
                row_size = (len(text) if text.isascii() else len(text.encode())) + 1
                if len(head) + row_size > budget:
                    raise pymysql.err.OperationalError(
                        1153, f"A {table} row of {row_size} bytes doesn't fit the server's "
                              f"max_allowed_packet ({budget + _PACKET_MARGIN} bytes); nothing was sent for it")
                if chunk and (len(chunk) >= _BATCH_ROWS or size + row_size > budget):
                    cur.execute(head + ",".join(chunk))
                    statements += 1
                    chunk, size = [], len(head)
                chunk.append(text)
                size += row_size
            if chunk:
                cur.execute(head + ",".join(chunk))
                statements += 1
        if log.debug_log_active():
            log.trace(f"SQL batch {(time.perf_counter() - start) * 1000:.2f}ms: "
                      f"{sum(len(rows) for rows in batch.values())} row(s) in {statements} INSERT(s)", area="sql", stacklevel=3)

    def executemany(self, sql, seq_of_params):
        """
        Runs one parameterized `INSERT`/`UPDATE`/`DELETE` against every
        tuple in `seq_of_params` -- `pymysql.cursors.Cursor.executemany`
        batches these into as few round trips as the driver can manage,
        rather than this module looping `execute` once per row itself.
        `seq_of_params` may be a generator (`insert_sector`'s per-vertex
        rows, in particular, are built as one) -- consumed exactly once,
        same as `sqlite3.Connection.executemany`.
        """
        self.flush()
        cur = self._conn.cursor()
        rows = list(seq_of_params)
        start = time.perf_counter()
        cur.executemany(sql.replace("?", "%s"), rows)
        if log.debug_log_active():
            log.trace(f"SQL executemany {(time.perf_counter() - start) * 1000:.2f}ms, {len(rows)} parameter "
                      f"row(s): {_sql_for_log(sql, rows[:3])}{' ...' if len(rows) > 3 else ''}", area="sql",
                      stacklevel=3)
        return _Cursor(cur)

    def execute_literal(self, sql):
        """
        Runs one statement that already holds its values (no `?` marks, no parameter formatting, no insert
        rewriting), for the bulk writes that build their own multi-row `INSERT` text (PERF.73).
        """
        self.flush()
        cur = self._conn.cursor()
        start = time.perf_counter()
        cur.execute(sql)
        if log.debug_log_active():
            log.trace(f"SQL literal {(time.perf_counter() - start) * 1000:.2f}ms, {cur.rowcount} row(s): {sql[:120]}",
                      area="sql", stacklevel=3)
        return _Cursor(cur)

    def executescript(self, script):
        """
        Runs a `;`-separated sequence of DDL statements -- `pymysql` has
        no `sqlite3.Connection.executescript` equivalent (no
        multi-statement single `execute` call). Comment lines (`--`,
        matching this project's own SQL style throughout `schema.sql`)
        are stripped before splitting; safe specifically for this
        project's own DDL, which never embeds a `;` inside a string
        literal.
        """
        self.flush()
        statements = []
        for line in script.splitlines():
            stripped = line.strip()
            if stripped.startswith("--"):
                continue
            statements.append(line)
        for statement in "\n".join(statements).split(";"):
            statement = statement.strip()
            if statement:
                if log.debug_log_active():
                    log.trace(f"SQL script statement: {statement[:200]}", area="sql")
                self._conn.cursor().execute(statement)
        self._conn.commit()

    def lock_until_commit(self, name, timeout_s=NAMED_LOCK_TIMEOUT_S):
        """
        Takes the MySQL named lock `name` (`GET_LOCK`) and holds it until
        this connection's next commit or rollback, so whatever the
        transaction reads after this is read with every other holder's
        work committed. A no-op when already held.

        Raises:
            pymysql.err.OperationalError: 1205 when `timeout_s` passes
                first (retried by `save_sector`, like a row lock wait).
        """
        if name in self._txn_locks:
            return
        row = self.execute("SELECT GET_LOCK(?, ?) AS ok", (name, timeout_s)).fetchone()
        if row["ok"] != 1:
            raise pymysql.err.OperationalError(1205, f"Lock wait timeout exceeded waiting for lock {name!r}")
        self._txn_locks.append(name)
        if not self._short_row_waits:
            # The holder must not wait long on a row: see
            # `LOCK_HOLDER_ROW_WAIT_S`.
            self._run("SET @planetgen_row_wait = @@SESSION.innodb_lock_wait_timeout,"
                      " SESSION innodb_lock_wait_timeout = ?", (LOCK_HOLDER_ROW_WAIT_S,))
            self._short_row_waits = True

    def _release_txn_locks(self):
        while self._txn_locks:
            name = self._txn_locks.pop()
            try:
                self._run("SELECT RELEASE_LOCK(?)", (name,))
            except pymysql.err.MySQLError:  # a lost session has already dropped it
                pass
        if self._short_row_waits:
            self._short_row_waits = False
            try:
                self._run("SET SESSION innodb_lock_wait_timeout = @planetgen_row_wait", ())
            except pymysql.err.MySQLError:  # a lost session has no setting left to restore
                pass

    def after_commit(self, callback):
        """Runs `callback()` once this transaction has committed (never when it rolls back). A callback that raises is
        logged and the others still run: the transaction is already durable."""
        self._after_commit.append(callback)

    def parity_touch(self, where, params=()):
        """
        Notes the sectors `WHERE where` (with `params`) selects as changed by this transaction, so that once it commits the
        parity file is brought up to date for them (`parity.protect_sectors`). Does nothing unless `parity.enabled` is on.
        """
        if not parity.enabled():
            return
        found = {row["id"] for row in self._run(f"SELECT id FROM sectors WHERE {where}", tuple(params)).fetchall()}
        if found and not self._parity_touched:
            self._after_commit.append(self._protect_touched)
        self._parity_touched |= found

    def _protect_touched(self):
        touched, self._parity_touched = self._parity_touched, set()
        parity.protect_sectors(self._config, touched)

    def commit(self):
        self.flush()
        self._conn.commit()
        self._flagged.clear()
        self._release_txn_locks()
        if log.debug_log_active():
            log.trace("SQL commit", area="sql")
        callbacks, self._after_commit = self._after_commit, []
        for callback in callbacks:
            try:
                callback()
            except Exception as exc:  # noqa: BLE001 -- the transaction is durable; what follows it must not undo that
                log.error(f"After-commit step failed: {type(exc).__name__}: {exc}")

    def rollback(self):
        if log.debug_log_active():
            log.trace(f"SQL rollback ({sum(len(rows) for rows in self._batch.values()) if self._batch else 0} "
                      f"held-back row(s) dropped)", area="sql")
        if self._batch:
            self._batch = {}
        self._flagged.clear()
        self._after_commit = []
        self._parity_touched = set()
        try:
            self._conn.rollback()
        finally:
            self._release_txn_locks()

    def close(self):
        self._batch = None
        self._batch_depth = 0
        self._release_txn_locks()
        self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is None:
            self.commit()
        else:
            self.rollback()
        return False


class _BatchScope:
    def __init__(self, conn):
        self._conn = conn

    def __enter__(self):
        conn = self._conn
        if conn._batch_depth == 0:
            conn._batch = {}
            conn.sector_centers_pc = {}
            conn.orbit_epoch_unix = None
        conn._batch_depth += 1
        return conn

    def __exit__(self, exc_type, exc_value, traceback):
        conn = self._conn
        conn._batch_depth -= 1
        if conn._batch_depth == 0:
            try:
                if exc_type is None:
                    conn.flush()
            finally:
                conn._batch = None
                conn.sector_centers_pc = {}
                conn.orbit_epoch_unix = None
        return False


class _InsertedCursor:
    """The cursor `Connection.execute` returns for an INSERT whose id came
    from `_allocate_id` (or that `batched` held back): `.lastrowid` is
    that id."""

    def __init__(self, lastrowid, cursor=None):
        self.lastrowid = lastrowid
        self._cursor = cursor
        self.rowcount = 1

    def fetchone(self):
        return None

    def fetchall(self):
        return []

    def __getattr__(self, name):
        if self._cursor is None:
            raise AttributeError(name)
        return getattr(self._cursor, name)


def get_connection(config=None, ensure_schema=True, statement_timeout_s=None):
    """
    Opens a pooled MySQL connection (creating the pool for `config` on
    first use), applying the schema if needed.

    Safe to call repeatedly against the same database -- `_ensure_schema`
    uses `CREATE TABLE IF NOT EXISTS`/`CREATE OR REPLACE VIEW` throughout
    (every index and foreign key is declared inline within its table, per
    `schema.sql`'s own header note on why -- MySQL's `CREATE INDEX` has no
    `IF NOT EXISTS` form), so an existing database is left untouched
    beyond having any missing tables/views added.

    Args:
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.
        statement_timeout_s (float, optional): Stop any statement on this
            connection that runs longer than this many seconds (PERF.17;
            the web's read-only connections), raising an `OperationalError`
            whose code is in `STATEMENT_TIMEOUT_ERRORS`. Such connections
            come from their own pool. `None` (the default): no limit.
        ensure_schema (bool): Whether to run `_ensure_schema` (DDL) on
            this connection. `True` (the default) suits every read-write
            caller (`save_sector`/`save_system`/the generation CLIs) --
            the whole point of `CREATE ... IF NOT EXISTS` is that a fresh
            database gets its schema the first time anything connects.
            Every read-only caller (`planetgen.db.query`'s `open_readonly`, the
            Flask API's `get_db`, `html/lib/dbutil.py`) should instead
            pass `False`: this project's own docs recommend pointing
            those tools at a database account with `SELECT`-only grants
            (see `planetgen.db.query`'s module docstring), and DDL requires
            `CREATE`, which such an account deliberately doesn't have --
            attempting `_ensure_schema` there would fail every single
            connection with a permissions error instead of just skipping
            a step that a read-write caller has already done once.

    Returns:
        Connection: An open connection (schema-initialized, if
                    `ensure_schema`). Supports `.execute(sql, params)`
                    (returning a cursor with `.lastrowid`/`.fetchone()`/
                    `.fetchall()`, and `row["column"]` access on each
                    row), `with conn:` for a commit-on-success/
                    rollback-on-exception block, and `.close()` (returns
                    the underlying connection to its pool rather than
                    truly closing a socket).
    """
    config = config or DEFAULT_MYSQL_CONFIG
    conn = Connection(_get_pool(config, statement_timeout_s).connect(), config)
    if log.debug_log_active():
        log.trace(f"SQL connection from the pool: {config.database}@{config.host}:{config.port} as {config.user}"
                  f"{'' if statement_timeout_s is None else f', statement limit {statement_timeout_s}s'}"
                  f"{', schema checked' if ensure_schema and config._key() not in _schema_ensured else ''}",
                  area="sql")
    if ensure_schema and config._key() not in _schema_ensured:
        try:
            if not _schema_is_current(conn):
                _ensure_schema(conn)
        except BaseException:
            conn.close()
            raise
        _schema_ensured.add(config._key())
    return conn


def escape_like(value):
    """
    Escapes `\\`, `%` and `_` in `value` so it matches only itself inside
    a SQL `LIKE` pattern -- paired with `ESCAPE '\\\\'` in the query (the
    same convention as `db.query._search_like_pattern`). Without it a
    prefix like `planetgen_` would also match `planetgenX...`, and a `%`
    would match anything.
    """
    return value.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


def configured_control_database():
    """The control schema's name: `CONTROL_DB_ENV_VAR`, else
    `config.json`'s `control_database`
    (what `control_mysql_config` connects to)."""
    return get_settings().control_database


def list_databases(base_config=None, prefix=None, with_counts=False):
    """
    Lists every MySQL schema on `base_config`'s server whose name starts
    with `prefix` (default: `mysql.database_prefix`) --
    "multiple databases" here means multiple MySQL schemas on one
    configured server (e.g. one schema per campaign/galaxy: `planetgen`,
    `planetgen_alpha`, ...), the way both the `html/` CGI browser's
    database picker and the Flask API's own `?db=`/`/api/databases`
    offer a choice among them.

    Args:
        base_config (MySQLConfig, optional): Connection parameters
            (host/port/user/password) to list schemas from -- its own
            `database` field is ignored (this opens a connection with no
            specific schema selected). Defaults to `DEFAULT_MYSQL_CONFIG`.
        prefix (str, optional): Overrides the env-var-derived default.
        with_counts (bool): Also count each schema's sectors and star
            systems (`sector_count`/`system_count`, `None` when it lacks
            this project's tables), over the same one connection: a
            connection per schema would leave a pool open for each
            (TEST.93: in the test suite, one per other worker's database,
            until the server ran out of connections).

    The control schema (`configured_control_database`: admin logins,
    sessions, API keys) is never listed, even when its name shares the
    prefix (`planetgen_control` does by default), so no `?db=` can select
    it either (`resolve_database`).

    Returns:
        list[dict]: One entry per matching schema, sorted by name, each
                    with `name`, `size_bytes` (sum of `data_length`/
                    `index_length` across its tables), and `modified_at`
                    (the latest `information_schema.tables.update_time`
                    across its tables, formatted `"%Y-%m-%d %H:%M"`, or
                    `"unknown"` when the storage engine doesn't track it).
    """
    base_config = base_config or DEFAULT_MYSQL_CONFIG
    prefix = prefix or get_settings().mysql.database_prefix
    conn = get_connection(
        MySQLConfig(
            host=base_config.host, port=base_config.port,
            user=base_config.user, password=base_config.password, database="",
        ),
        ensure_schema=False,
    )
    try:
        schema_rows = conn.execute(
            "SELECT schema_name AS name FROM information_schema.schemata "
            "WHERE schema_name LIKE ? ESCAPE '\\\\' AND schema_name <> ? ORDER BY schema_name",
            (f"{escape_like(prefix)}%", configured_control_database()),
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
            entry = {
                "name": name,
                "size_bytes": int(stats["size_bytes"] or 0),
                "modified_at": modified_at.strftime("%Y-%m-%d %H:%M") if modified_at else "unknown",
            }
            if with_counts:
                entry.update(_schema_counts(conn, name))
            entries.append(entry)
        return entries
    finally:
        conn.close()


def _holds_galaxy(base_config, name):
    """Whether schema `name` holds this project's galaxy tables (`sectors`). A schema that shares the prefix and holds
    something else (the real `planetgen_control` on a server whose configured control schema has another name) is no
    database to select: every route would answer 500 from its first missing table (API.32)."""
    conn = get_connection(
        MySQLConfig(host=base_config.host, port=base_config.port,
                    user=base_config.user, password=base_config.password, database=""),
        ensure_schema=False,
    )
    try:
        return conn.execute("SELECT 1 FROM information_schema.tables WHERE table_schema = ? AND table_name = 'sectors'"
                            " LIMIT 1", (name,)).fetchone() is not None
    finally:
        conn.close()


def _schema_counts(conn, name):
    """`list_databases`'s `sector_count`/`system_count` for schema `name`,
    both `None` when it can't be read (no such tables: mid-migration, or
    an unrelated schema sharing the prefix; dropped since the listing)."""
    schema = "`" + name.replace("`", "``") + "`"
    try:
        sectors = conn.execute(f"SELECT COUNT(*) AS n FROM {schema}.sectors").fetchone()["n"]
        systems = conn.execute(f"SELECT COUNT(*) AS n FROM {schema}.star_systems").fetchone()["n"]
    except pymysql.MySQLError:
        return {"sector_count": None, "system_count": None}
    return {"sector_count": sectors, "system_count": systems}


def resolve_database(base_config, name, prefix=None):
    """
    Validates a database name supplied by a caller (e.g. a web request's
    `?db=`/`db` parameter) against the same prefix-filtered list
    `list_databases` offers, and returns a ready-to-use `MySQLConfig`.

    This is what keeps an externally-supplied database name from
    selecting a schema this deployment never meant to expose (every other
    schema on a shared MySQL server, `information_schema` itself, etc.) --
    only an exact, case-sensitive match against a currently-listed schema
    is accepted.

    Args:
        base_config (MySQLConfig): Connection parameters (host/port/user/
            password) to validate/connect against.
        name (str): The requested database name.
        prefix (str, optional): Passed through to `list_databases`.

    Returns:
        MySQLConfig: `base_config` with `database` set to `name`.

    Raises:
        ValueError: If `name` is empty, doesn't match a listed schema, or the schema holds no galaxy (API.32).
    """
    if not name:
        raise ValueError("No database specified.")
    if name not in {entry["name"] for entry in list_databases(base_config, prefix=prefix)}:
        raise ValueError(f"No such database: {name!r}")
    if not _holds_galaxy(base_config, name):
        raise ValueError(f"No such database: {name!r}")
    return MySQLConfig(
        host=base_config.host, port=base_config.port,
        user=base_config.user, password=base_config.password, database=name,
    )


SCHEMA_LOCK_WAIT_S = 600
"""int: How long a connection waits for another one to finish creating or
migrating the same database's schema (`_schema_lock`) before it gives up.
Creating a new database's schema takes a few seconds; a migration that
rewrites big tables can take minutes."""


@contextlib.contextmanager
def _schema_lock(conn, what="schema"):
    """
    Holds the connected database's schema lock (a MySQL named lock, one per
    database and per `what`) while schema DDL and its version row are
    written (DB.5): when several first connections reach an empty database
    at once, one creates the schema and the others wait, then find it
    current. Without it each ran `schema.sql` and inserted the baseline
    `schema_migrations` row, and all but one failed on a duplicate key.

    Held on the session, not the transaction: the DDL commits implicitly.

    Raises:
        pymysql.err.OperationalError: 1205 when `SCHEMA_LOCK_WAIT_S` passes
            before the lock is free.
    """
    database = conn.execute("SELECT DATABASE() AS db").fetchone()["db"] or ""
    name = f"planetgen-{what}:" + hashlib.sha256(database.encode("utf-8")).hexdigest()[:32]
    got = conn.execute("SELECT GET_LOCK(?, ?) AS ok", (name, SCHEMA_LOCK_WAIT_S)).fetchone()["ok"]
    if got != 1:
        raise pymysql.err.OperationalError(
            1205, f"Lock wait timeout exceeded waiting for {database!r}'s {what} lock")
    # A fresh transaction, so what follows reads what the last holder committed.
    conn.commit()
    try:
        yield
    finally:
        try:
            conn.execute("SELECT RELEASE_LOCK(?)", (name,))
        except pymysql.err.MySQLError:  # a lost session has already dropped it
            pass


def _schema_is_current(conn):
    """
    Whether the connected database already holds the schema `SCHEMA_VERSION` names, so `_ensure_schema` has nothing
    to add (PERF.174): its version row is the code's version. Every change to `schema.sql` comes with a migration and a
    new version, and a database behind it is refused by `_ensure_schema`, so a database at the current version already
    has every table and view. A work unit is a new process, and without this each one ran the schema's ten-odd
    `CREATE TABLE IF NOT EXISTS` and its views again under the schema lock (about 32,000 statements in one four-minute
    scatter) to find nothing to do.
    """
    return _stored_version(conn, "schema_migrations") == SCHEMA_VERSION


def _ensure_schema(conn):
    """
    Applies `schema.sql` to `conn`, then bootstraps `schema_migrations`
    (inserting the baseline row) if it's empty -- idempotent, safe to call
    on a database that already has some, all, or none of the schema, and
    from several connections at once (`_schema_lock`).

    The baseline is `SCHEMA_VERSION` for a new database. An existing one
    whose `schema_migrations` was emptied or lost is refused
    (`SchemaVersionMissingError`, DB.97) before any DDL runs.

    Args:
        conn (Connection): The connection to apply the schema to.

    Raises:
        SchemaTooNewError: The database is already past `SCHEMA_VERSION`
            (checked before any DDL runs).
        SchemaVersionMissingError: The database has tables but no recorded
            version.
    """
    with _schema_lock(conn):
        _apply_schema(conn)


def _apply_schema(conn):
    """`_ensure_schema`'s work, for a caller already holding the schema
    lock."""
    stored = _stored_version(conn, "schema_migrations")
    _refuse_newer(conn, stored, SCHEMA_VERSION)
    _refuse_unversioned(conn, stored)
    _refuse_older(conn, stored)
    with open(SCHEMA_PATH, "r", encoding="utf-8") as f:
        conn.executescript(f.read())
    if conn._config is not None:
        forget_id_blocks(conn._config._key())

    row = conn.execute("SELECT COUNT(*) AS n FROM schema_migrations").fetchone()
    if row["n"] == 0:
        conn.execute("INSERT INTO schema_migrations (version) VALUES (?)", (SCHEMA_VERSION,))
    conn.commit()


def open_write(config=None):
    """
    Opens a connection for a write-capable caller against an already-
    existing content database (`ensure_schema=False` -- same reasoning as
    `db.query.open_readonly`: if the account passed in (see
    `docs/deployment/README.md`'s "MySQL accounts") lacks `CREATE`/`ALTER`
    grants, attempting `_ensure_schema`'s DDL here would fail every
    connection instead of just skipping a step a full-access account has
    already done once, via `planetgen.cli.migrate`).

    Args:
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.

    Returns:
        Connection: An open connection.
    """
    return get_connection(config, ensure_schema=False)


def control_mysql_config(base_config=None):
    """
    Builds a `MySQLConfig` pointed at the control schema (see
    `control_schema.sql`'s header comment), reusing `base_config`'s
    host/port/user/password -- typically the same account `open_write`
    uses (`docs/deployment/README.md`'s "MySQL accounts"), since the
    control schema needs the same `SELECT`/`INSERT`/`UPDATE`/`DELETE`
    grants, just on a different schema name.

    Args:
        base_config (MySQLConfig, optional): Connection parameters to
            reuse host/port/user/password from. Defaults to
            `DEFAULT_MYSQL_CONFIG`.

    Returns:
        MySQLConfig: `base_config` with `database` replaced by
            `CONTROL_DB_ENV_VAR`.
    """
    base_config = base_config or DEFAULT_MYSQL_CONFIG
    database = configured_control_database()
    return MySQLConfig(
        host=base_config.host, port=base_config.port,
        user=base_config.user, password=base_config.password, database=database,
    )


def get_control_connection(config=None, ensure_schema=False, lock_wait_s=None):
    """
    Opens a pooled connection to the control schema.

    Args:
        config (MySQLConfig, optional): Connection parameters. Defaults
            to `control_mysql_config()`.
        ensure_schema (bool): Whether to run `_ensure_control_schema`
            (DDL) on this connection -- `False` by default (the normal
            runtime case: the Flask API's account may have no `CREATE`
            grant, same reasoning as `open_write` above).
            `admin.auth.bootstrap_control_schema` passes `True`, using a
            full-access account (the same one `planetgen.cli.migrate` already
            uses), to create/update the control schema once per deploy.
        lock_wait_s (float, optional): How many seconds a statement on this
            connection waits for a row lock before it fails with
            `LOCK_WAIT_ERRORS` (the server's own limit, 50 s, when `None`):
            the web gives up on a lock sooner than that so a request does
            not hold one of the app's few threads for the whole wait (ADM.53).

    Returns:
        Connection: An open connection.
    """
    config = config or control_mysql_config()
    conn = Connection(_get_pool(config, lock_wait_s=lock_wait_s).connect(), config)
    if ensure_schema:
        try:
            _ensure_control_schema(conn)
        except BaseException:
            conn.close()
            raise
    return conn


def _ensure_control_schema(conn):
    """
    Applies `control_schema.sql` to `conn`, then bootstraps
    `control_schema_migrations` (inserting `CONTROL_SCHEMA_VERSION` when
    the newest row is older, or the table is empty) -- idempotent, mirrors `_ensure_schema`
    above exactly, just against the control schema's own DDL file/version
    counter.

    Args:
        conn (Connection): The connection to apply the schema to.
    """
    with _schema_lock(conn, "control"):
        _apply_control_schema(conn)


def _apply_control_schema(conn):
    """`_ensure_control_schema`'s work, under the control schema's lock."""
    _refuse_newer(conn, _stored_version(conn, "control_schema_migrations"), CONTROL_SCHEMA_VERSION,
                  "control schema")
    _drop_old_generation_stats(conn)
    with open(CONTROL_SCHEMA_PATH, "r", encoding="utf-8") as f:
        conn.executescript(f.read())
    _add_control_columns(conn)

    row = conn.execute("SELECT MAX(version) AS v FROM control_schema_migrations").fetchone()
    if row["v"] is not None and row["v"] < 13:
        _grandfather_api_key_scopes(conn)
    if row["v"] is not None and row["v"] < 15:
        _promote_control_owner(conn)
    if row["v"] is not None and row["v"] < 21:
        _drop_work_queue_tables(conn)
    if row["v"] is None or row["v"] < CONTROL_SCHEMA_VERSION:
        conn.execute("INSERT INTO control_schema_migrations (version) VALUES (?)", (CONTROL_SCHEMA_VERSION,))
        conn.commit()
        if row["v"] is not None:
            activity_log.event("DB", "migrate", db=configured_control_database(), from_version=row["v"],
                              to_version=CONTROL_SCHEMA_VERSION)


def _grandfather_api_key_scopes(conn):
    """
    v13 (API.9): every API key made before scopes existed was an admin key,
    so each gets the `admin` scope (which implies all the others).
    """
    conn.execute("INSERT IGNORE INTO admin_api_key_scopes (key_id, scope) SELECT id, 'admin' FROM admin_api_keys")
    conn.commit()


def _promote_control_owner(conn):
    """
    v15 (DB.58): the accounts that existed before roles are admins (the column's default), and the lowest id, the
    account the install seeded, becomes the one owner. Does nothing when an owner is already there or there are no
    accounts (the install seeds the owner).
    """
    if conn.execute("SELECT 1 FROM admin_users WHERE role = 'owner' LIMIT 1").fetchone() is not None:
        return
    first = conn.execute("SELECT id, username FROM admin_users ORDER BY id LIMIT 1").fetchone()
    if first is None:
        return
    conn.execute("UPDATE admin_users SET role = 'owner' WHERE id = ?", (first["id"],))
    conn.commit()
    log.normal(f"Control schema v15: account '{first['username']}' (id {first['id']}), the lowest, is now the owner.")


def _drop_old_generation_stats(conn):
    """
    v12 (PERF.32): `generation_stats` gained `workers` and `version_key`
    and a new primary key. Its rows are deleted by the first run of any
    other version anyway, so an older table is dropped (the schema file
    then creates the new one) rather than rebuilt.
    """
    old = conn.execute(
        "SELECT COUNT(*) AS n FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'generation_stats'").fetchone()["n"] and not conn.execute(
        "SELECT COUNT(*) AS n FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'generation_stats' AND COLUMN_NAME = 'workers'").fetchone()["n"]
    if old:
        conn.execute("DROP TABLE generation_stats")
        conn.commit()


def _drop_work_queue_tables(conn):
    """
    v21 (PERF.100): drops the work queue's tables, `work_job_args`, `work_tasks`, `work_jobs` and `work_lease` (the
    children first, they reference `work_jobs`). The job tree has had no writer since PERF.98, and a run and the pause of
    the queues are on Redis.
    """
    for table in ("work_job_args", "work_tasks", "work_jobs", "work_lease"):
        conn.execute(f"DROP TABLE IF EXISTS {table}")
    conn.commit()


def _add_control_columns(conn):
    """
    Adds `_CONTROL_COLUMNS` and `_CONTROL_CONSTRAINTS` to
    control tables created before them -- each only when missing, so it
    is safe to run on every deploy (MySQL has no `ADD COLUMN IF NOT EXISTS`).
    """
    def existing(sql, table):
        return {row["name"] for row in conn.execute(sql, (table,)).fetchall()}

    for table, columns in _CONTROL_COLUMNS.items():
        have = existing("SELECT COLUMN_NAME AS name FROM information_schema.COLUMNS"
                        " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", table)
        missing = [f"ADD COLUMN {name} {definition}" for name, definition in columns if name not in have]
        if missing:
            conn.execute(f"ALTER TABLE {table} {', '.join(missing)}")
    for table, name, clause in _CONTROL_CONSTRAINTS:
        if name not in existing("SELECT CONSTRAINT_NAME AS name FROM information_schema.TABLE_CONSTRAINTS"
                                " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", table):
            conn.execute(f"ALTER TABLE {table} {clause}")
    conn.commit()


def _tristate(value):
    """
    Converts a `SystemConfig` tri-state value (`True`/`False`/`None`) to
    the schema's nullable `INTEGER` `0`/`1`/`NULL` representation.

    Args:
        value (bool or None): The tri-state value.

    Returns:
        int or None: `1`, `0`, or `None`.
    """
    return None if value is None else int(bool(value))


def _lifespan_gy(value):
    """
    Converts a star's `lifespan` (a float, or `float('inf')` for white
    dwarfs) to the schema's convention: `NULL` means infinite. Never the
    non-standard JSON `Infinity` token -- pymysql itself rejects a raw
    `float('inf')` before it ever reaches the server ("inf can not be
    used with MySQL"), so this must run before any `INSERT`.

    Args:
        value (float): The lifespan in billions of years, or `float('inf')`.

    Returns:
        float or None: The lifespan, or `None` if infinite.
    """
    return None if value == float("inf") else value


def insert_system_config(conn, config: SystemConfig) -> int:
    """
    Inserts a `SystemConfig` "recipe" row (plus its `SLOTS` child rows, if
    any) and returns the new `system_configs.id`.

    Always inserts a new row -- configs aren't deduplicated, since each
    generated `StarSystem` has its own config instance regardless of
    whether its values happen to match another system's.

    Args:
        conn (Connection): An open, schema-initialized connection.
        config (SystemConfig): The recipe to persist.

    Returns:
        int: The new `system_configs.id`.
    """
    cur = conn.execute(
        """
        INSERT INTO system_configs (
            habitable_world, asteroid_belt, comets, large_star, moons,
            max_planets, planets, star_type, name, age, intelligent_life,
            binary_system, wide_binary, num_orbits, stars, max_stars, {prevalence_columns}
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, {prevalence_marks})
        """.format(prevalence_columns=", ".join(f"prevalence_{feature}" for feature in PREVALENCE_FEATURES),
                   prevalence_marks=", ".join("?" for _ in PREVALENCE_FEATURES)),
        (
            _tristate(config.HABITABLE_WORLD),
            _tristate(config.ASTEROID_BELT),
            _tristate(config.COMETS),
            _tristate(config.LARGE_STAR),
            _tristate(config.MOONS),
            _tristate(config.MAX_PLANETS),
            _tristate(config.PLANETS),
            config.STAR_TYPE,
            config.NAME,
            config.AGE,
            _tristate(config.INTELLIGENT_LIFE),
            _tristate(config.BINARY_SYSTEM),
            _tristate(config.WIDE_BINARY),
            config.NUM_ORBITS,
            config.STARS,
            config.MAX_STARS,
            *(float((config.PREVALENCE or {}).get(feature) or 0.0) for feature in PREVALENCE_FEATURES),
        ),
    )
    config_id = cur.lastrowid

    for orbit_index, slot in enumerate(config.SLOTS or []):
        if slot is None:
            continue
        conn.execute(
            """
            INSERT INTO system_config_slots (config_id, orbit_index, type, planet_class, moons)
            VALUES (?, ?, ?, ?, ?)
            """,
            (config_id, orbit_index, slot.get("type"), slot.get("planet_class"), slot.get("moons")),
        )

    return config_id


# ---------------------------------------------------------------------------
# Name-uniqueness reservation (planetgen/names/uniqueness.py, v24) --
# one reserve/confirm pair per naming level (sector, system), called from
# insert_sector/insert_star_system below. "Reserve" runs *before* that
# level's own INSERT (it doesn't know the new row's id yet, but may need
# to rename an existing row via a plain UPDATE); "confirm" runs *after*,
# now that the id is known, to upsert that level's registry row. Stars,
# planets and moons are named from their system (names/bodies.py, v34), so
# they need no registry of their own. See names/uniqueness.py's own module
# docstring for the full sector > system hierarchy.
# ---------------------------------------------------------------------------

def _rename_bodies_with_prefix(conn, star_system_id, old_prefix, new_prefix):
    """
    Renames every star, planet and moon in one system whose name carries
    `old_prefix` (`names.bodies.rename_prefix`), so derived names such as
    `Voranthis II` follow their system or star. Names set by hand that
    don't carry the prefix are left alone.
    """
    for table in ("stars", "planets", "moons"):
        rows = conn.execute(f"SELECT id, name FROM {table} WHERE star_system_id = ?", (star_system_id,)).fetchall()
        for row in rows:
            renamed = rename_prefix(row["name"], old_prefix, new_prefix)
            if renamed is not None and renamed != row["name"]:
                conn.execute(f"UPDATE {table} SET name = ? WHERE id = ?", (renamed, row["id"]))
    _rename_comets(conn, "star_system_id = ?", (star_system_id,), old_prefix, new_prefix)


def _rename_comets(conn, where, params, old_host, new_host):
    """Moves the comets matching `where` whose designation
    (`cometData.comet_designation`) names `old_host` over to `new_host`."""
    for row in conn.execute(f"SELECT id, name FROM comets WHERE {where}", params).fetchall():
        renamed = rename_comet_designation(row["name"], old_host, new_host)
        if renamed is not None and renamed != row["name"]:
            conn.execute("UPDATE comets SET name = ? WHERE id = ?", (renamed, row["id"]))


SYSTEM_NAME_MAX_LENGTH = 200
"""int: The longest name a star system or star may be given (TEST.12).
Every name column is VARCHAR(255), but these names grow after they're
chosen -- a collision prefixes them ("Omicron ", "Malutki "), and planets,
moons and comets are named after them (`<name> XII ab`) -- so a 255-character
system name failed to save, or to rename, with MySQL error 1406."""


def check_system_name_length(name):
    """Raises ValueError when `name` is longer than `SYSTEM_NAME_MAX_LENGTH`.
    Checked where a name is chosen (`insert_star_system`, the API, the
    CLI), not in the renames, which also add collision prefixes."""
    if name is not None and len(name) > SYSTEM_NAME_MAX_LENGTH:
        raise ValueError(f"a star system or star name may be at most {SYSTEM_NAME_MAX_LENGTH} characters "
                         f"(this one has {len(name)})")


def rename_star_system(conn, star_system_id, new_name):
    """
    Renames one system and every star, planet and moon still named after
    it (`_rename_bodies_with_prefix`). The one place a system's name
    changes after it's saved: the name-uniqueness decorations below, and
    `PATCH /api/systems/<id>`.

    Returns:
        bool: `False` if no such system exists.
    """
    row = conn.execute("SELECT name, binary_configuration FROM star_systems WHERE id = ?",
                       (star_system_id,)).fetchone()
    if row is None:
        return False
    conn.execute(
        "UPDATE star_systems SET name = ?, modified_at = CURRENT_TIMESTAMP(3) WHERE id = ?",
        (new_name, star_system_id),
    )
    storage_checksum.mark_system(conn, star_system_id)
    if row["binary_configuration"] == "wide":
        # A wide pair's stars and its primary's planets carry only the
        # system name's first word (GEN.62).
        _rename_bodies_with_prefix(conn, star_system_id, wide_pair_first_word(row["name"]),
                                   wide_pair_first_word(new_name))
    else:
        _rename_bodies_with_prefix(conn, star_system_id, row["name"], new_name)
    return True


def rename_star(conn, star_id, new_name):
    """
    Renames one star. A single star shares its system's name, so renaming
    it renames the system (`rename_star_system`). A binary's star is
    renamed on its own, along with the planets and moons named after it
    (a wide pair's).

    Returns:
        bool: `False` if no such star exists.
    """
    row = conn.execute("SELECT star_system_id, role, name FROM stars WHERE id = ?", (star_id,)).fetchone()
    if row is None:
        return False
    if row["role"] == "single":
        return rename_star_system(conn, row["star_system_id"], new_name)
    conn.execute("UPDATE stars SET name = ? WHERE id = ?", (new_name, star_id))
    storage_checksum.mark_system(conn, row["star_system_id"])
    # Since GEN.62 a wide pair's planets carry one word of their star's
    # name: the primary's the shared first word, the secondary's its own
    # last. Planets named before carry the whole star name.
    prefixes = [(row["name"], new_name)]
    if row["name"] and new_name and len(row["name"].split()) == 2:
        pick = 0 if row["role"] == "primary" else -1
        prefixes.append((row["name"].split()[pick], new_name.split()[pick]))
    for table in ("planets", "moons"):
        rows = conn.execute(
            f"SELECT id, name FROM {table} WHERE star_system_id = ? AND star_id = ?",
            (row["star_system_id"], star_id),
        ).fetchall()
        for body in rows:
            renamed = next((r for r in (rename_prefix(body["name"], old, new) for old, new in prefixes)
                            if r is not None), None)
            if renamed is not None:
                conn.execute(f"UPDATE {table} SET name = ? WHERE id = ?", (renamed, body["id"]))
    _rename_comets(conn, "star_system_id = ? AND star_id = ?", (row["star_system_id"], star_id),
                   row["name"], new_name)
    touch_star_system(conn, row["star_system_id"])
    return True


def rename_body(conn, table, body_id, new_name):
    """
    Renames one planet or moon (`table` is `'planets'` or `'moons'`).
    A planet's moons keep their names.

    Returns:
        bool: `False` if no such row exists.
    """
    if table not in ("planets", "moons"):
        raise ValueError(f"not a body table: {table!r}")
    row = conn.execute(f"SELECT star_system_id FROM {table} WHERE id = ?", (body_id,)).fetchone()
    if row is None:
        return False
    conn.execute(f"UPDATE {table} SET name = ? WHERE id = ?", (new_name, body_id))
    touch_star_system(conn, row["star_system_id"])
    return True


UNIQUE_NAME_TABLES = (
    "sectors", "star_systems", "stars",
    "black_holes", "neutron_stars", "nebulae", "supernova_remnants", "rogue_planets", "quasars",
)
"""tuple: The tables whose names must be unique across the galaxy --
what `name_in_use` searches (the uniquely named phenomena since v40).
Planets, moons, comets and asteroid fields never are."""


def name_in_use(conn, name, exclude=None):
    """
    Whether any uniquely named object -- a sector, system, star or
    uniquely named phenomenon (`NAMED_PHENOMENON_TABLES`) -- is
    already called exactly `name`, the check a rename runs before it goes
    ahead. Planets and moons are left out (Boss, 2026-09-30): their names
    come from their star's (`names/bodies.py`), so they're unique whenever
    it is, and a body renamed by hand may share another body's name.
    Each table's `name` column is indexed, so this is one index lookup
    per table.

    Args:
        conn (Connection): An open connection.
        name (str): The proposed name.
        exclude (iterable, optional): `(table, id)` pairs for the rows
            being renamed, which may already carry `name`.

    Returns:
        str or None: The table holding the clash (one of
            `UNIQUE_NAME_TABLES`), or `None` when the name is free.
    """
    exclude = {tuple(pair) for pair in (exclude or ())}
    for table in UNIQUE_NAME_TABLES:
        rows = conn.execute(f"SELECT id FROM {table} WHERE name = ? LIMIT 3", (name,)).fetchall()
        if any((table, row["id"]) not in exclude for row in rows):
            return table
    return None


def _rename_existing_system_for_diminutive(conn, base_name):
    """
    Called while reserving a *sector* name that collides with an existing
    system's base name -- decorates that system (never the sector) with
    the next diminutive prefix, applied on top of whatever its current
    name already is (it may already carry its own Greek/Roman decoration
    from an unrelated system-vs-system collision).

    Args:
        conn (Connection): Part of the same transaction as the caller's
            own sector INSERT.
        base_name (str): The sector's own (undecorated) candidate name.

    Returns:
        bool: `True` if resolved (including "no colliding system exists
            at all", a no-op). `False` if `names.DIMINUTIVE_PREFIXES` is
            exhausted for this base name, or the prefixed system name
            would go past `MAX_SYSTEM_NAME_WORDS` (GEN.46) -- the caller
            must draw an entirely fresh sector name instead.
    """
    row = conn.execute(
        "SELECT first_star_system_id, first_object_table, first_object_id, diminutive_index "
        "FROM system_name_registry WHERE base_name = ? FOR UPDATE",
        (base_name,),
    ).fetchone()
    if row is None:
        return True

    prefix, next_index = resolve_diminutive(row["diminutive_index"])
    if prefix is None:
        return False

    current = _registry_holder_name(conn, row)
    if current is not None and not fits_word_limit(f"{prefix} {current}", word_limit_for(base_name)):
        # A third word isn't allowed (GEN.46). A holder that already
        # carries a diminutive can never meet a sector name (those only
        # get Greek letters), so it can stay as it is; any other holder
        # keeps its name and the sector draws a fresh one.
        return has_diminutive(current)
    if current is not None:
        _rename_registry_holder(conn, row, f"{prefix} {current}")
    # else: that holder was since deleted -- nothing left to rename,
    # but diminutive_index still advances below so a later collision on
    # this base name doesn't reuse the same prefix.
    conn.execute(
        "UPDATE system_name_registry SET diminutive_index = ? WHERE base_name = ?",
        (next_index, base_name),
    )
    return True


def reserve_sector_name(conn, candidate_name, sector_id):
    """
    Sector name-uniqueness reservation, run just after the sector's own
    INSERT (with its candidate name) in the same transaction -- resolves a
    sector-vs-sector collision (`names.uniqueness.resolve_greek_roman_collision`),
    then a cross-level collision against an existing system's base name
    (renaming *that* system instead of this sector's own name -- see
    `_rename_existing_system_for_diminutive`). Draws an entirely fresh
    candidate (`planetgen.names.wordsalad.generate_sector_name`) and starts
    over whenever any of those mechanisms is exhausted.

    The registry row is written first (`INSERT ... ON DUPLICATE KEY
    UPDATE occurrence_count = occurrence_count + 1`) and read back
    after (PERF.14): the upsert's row lock makes a second writer with the
    same base name wait for this transaction, then count past it, where
    the old `SELECT ... FOR UPDATE` of a missing name took a gap lock that
    two writers could deadlock on. A drawn-again name leaves its base's
    count one high, which only means the next holder gets a decoration
    it could have skipped.

    Args:
        conn (Connection): Part of the same transaction as the caller's
            own sector INSERT.
        candidate_name (str): The freshly generated name to reserve.
        sector_id (int): The new `sectors.id` -- recorded as the base
            name's first holder when it is one.

    Returns:
        tuple: `(final_name, base_name)` -- `final_name` is what the
            `sectors` row's `name` column should hold.
    """
    while True:
        base = candidate_name
        conn.execute(
            "INSERT INTO sector_name_registry (base_name, occurrence_count, first_sector_id) VALUES (?, 1, ?) "
            "ON DUPLICATE KEY UPDATE occurrence_count = occurrence_count + 1",
            (base, sector_id),
        )
        row = conn.execute(
            "SELECT occurrence_count, base_name FROM sector_name_registry WHERE base_name = ?", (base,),
        ).fetchone()
        existing_count = row["occurrence_count"] - 1

        new_name, rename = resolve_greek_roman_collision(base, existing_count)
        if new_name is None:
            candidate_name = generate_sector_name()
            continue
        # The system side first: if it can't make room, this sector draws
        # a fresh name and the sector holder below keeps its own.
        if not _rename_existing_system_for_diminutive(conn, base):
            candidate_name = generate_sector_name()
            continue
        if rename is not None:
            # The holder keeps the spelling it was registered under
            # (`WHERE name = ?` still finds it: the collation ignores case
            # and accents).
            old_name, renamed_to = resolve_greek_roman_collision(row["base_name"], existing_count)[1]
            conn.execute("UPDATE sectors SET name = ?, parity_dirty = 1 WHERE name = ? AND id <> ?", (renamed_to, old_name, sector_id))

        return new_name, base


def _regenerate_star_name():
    """A fresh star-name candidate, for `reserve_system_names`'
    exhaustion fallback -- same generator `star.Star` itself uses."""
    return generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)


NAMED_PHENOMENON_TABLES = (
    "black_holes", "neutron_stars", "nebulae", "supernova_remnants", "rogue_planets", "quasars",
)
"""tuple: The phenomenon tables whose names are drawn through
`system_name_registry` alongside star systems (v40, Boss 2026-09-30):
the uniquely named objects. Comets and asteroid fields get designations
(`comet_designation`, `asteroid_field_designation`) instead, and planets
and moons are named from their star."""


def _registry_holder_name(conn, row):
    """The current name of a `system_name_registry` row's first holder --
    a star system, or a phenomenon (v40) -- or `None` once it's deleted."""
    if row["first_star_system_id"] is not None:
        table, holder_id = "star_systems", row["first_star_system_id"]
    elif row["first_object_table"] in NAMED_PHENOMENON_TABLES:
        table, holder_id = row["first_object_table"], row["first_object_id"]
    else:
        return None
    current = conn.execute(f"SELECT name FROM {table} WHERE id = ?", (holder_id,)).fetchone()
    return None if current is None else current["name"]


def _rename_registry_holder(conn, row, new_name):
    """Renames a `system_name_registry` row's first holder: a star system
    (`rename_star_system`) or a phenomenon (`rename_phenomenon`)."""
    if row["first_star_system_id"] is not None:
        rename_star_system(conn, row["first_star_system_id"], new_name)
    elif row["first_object_table"] in NAMED_PHENOMENON_TABLES:
        rename_phenomenon(conn, row["first_object_table"], row["first_object_id"], new_name)


def rename_phenomenon(conn, table, object_id, new_name):
    """
    Renames one uniquely named phenomenon (`NAMED_PHENOMENON_TABLES`). A
    supernova remnant's collapsed core, named `"<remnant> Core"`, follows
    it.

    Returns:
        bool: `False` if no such row exists.
    """
    if table not in NAMED_PHENOMENON_TABLES:
        raise ValueError(f"not a named phenomenon table: {table!r}")
    row = conn.execute(f"SELECT name FROM {table} WHERE id = ?", (object_id,)).fetchone()
    if row is None:
        return False
    conn.execute(f"UPDATE {table} SET name = ? WHERE id = ?", (new_name, object_id))
    if table == "supernova_remnants":
        core = conn.execute(
            "SELECT compact_remnant_black_hole_id, compact_remnant_neutron_star_id "
            "FROM supernova_remnants WHERE id = ?", (object_id,),
        ).fetchone()
        for core_table, core_id in (("black_holes", core["compact_remnant_black_hole_id"]),
                                    ("neutron_stars", core["compact_remnant_neutron_star_id"])):
            if core_id is not None:
                conn.execute(
                    f"UPDATE {core_table} SET name = ? WHERE id = ? AND name = ?",
                    (f"{new_name} Core", core_id, f"{row['name']} Core"),
                )
    return True


_NAME_BATCH = 500
"""int: Base names per registry statement in `reserve_system_names`."""


def _name_key(name):
    """`name` folded the way the name columns' `utf8mb4_unicode_ci`
    compares it -- case, accents and trailing spaces ignored ("Vega",
    "VEGA" and "Véga" are one name to the registries), so the Python side
    groups names exactly as the unique keys do."""
    decomposed = unicodedata.normalize("NFKD", name)
    return "".join(c for c in decomposed if not unicodedata.combining(c)).casefold().rstrip(" ")


def _registry_rows(conn, bases):
    """`system_name_registry` rows for `bases`, keyed by `_name_key` of
    their stored base name."""
    rows = {}
    for first in range(0, len(bases), _NAME_BATCH):
        chunk = bases[first:first + _NAME_BATCH]
        for row in conn.execute(
            "SELECT id, base_name, occurrence_count, diminutive_index, first_star_system_id, first_object_table, "
            f"first_object_id FROM system_name_registry WHERE base_name IN ({', '.join('?' * len(chunk))})",
            tuple(chunk),
        ).fetchall():
            rows[_name_key(row["base_name"])] = row
    return rows


def _claim_names(c, keys, counts, spelled):
    """
    The registry statements of one `reserve_system_names` pass on connection
    `c`: which of `keys` an existing sector already uses as its base name,
    then one multi-row upsert in sorted order adding each key's `counts`,
    then the rows read back. Returns `(sector_hits, rows)`, `rows` keyed by
    `_name_key` (a key stored under a spelling `_name_key` doesn't match is
    looked up by its spelling).
    """
    sector_hits = set()
    for first in range(0, len(keys), _NAME_BATCH):
        chunk = [spelled[key] for key in keys[first:first + _NAME_BATCH]]
        sector_hits.update(_name_key(row["base_name"]) for row in c.execute(
            f"SELECT base_name FROM sector_name_registry WHERE base_name IN ({', '.join('?' * len(chunk))})",
            tuple(chunk),
        ).fetchall())
    for first in range(0, len(keys), _NAME_BATCH):
        chunk = keys[first:first + _NAME_BATCH]
        c.execute(
            "INSERT INTO system_name_registry (base_name, occurrence_count) VALUES "
            + ", ".join(["(?, ?)"] * len(chunk))
            + " ON DUPLICATE KEY UPDATE occurrence_count = occurrence_count + VALUES(occurrence_count)",
            tuple(value for key in chunk for value in (spelled[key], counts[key])),
        )
    rows = _registry_rows(c, [spelled[key] for key in keys])
    for key in keys:
        if key not in rows:
            rows[key] = c.execute(
                "SELECT id, base_name, occurrence_count, diminutive_index, first_star_system_id, "
                "first_object_table, first_object_id FROM system_name_registry WHERE base_name = ?",
                (spelled[key],),
            ).fetchone()
    return sector_hits, rows


def _release_name_claims(claim, claims):
    """
    Gives back the counts `reserve_system_names` claimed on `claim` when the
    save they were for rolled back (PERF.49), so a retry or a resumed run
    draws the same decorations it would have. A count is only taken back
    while it is still what the claim left (nobody has claimed after it):
    otherwise the count stays one high, which just skips a decoration.
    """
    try:
        with claim:
            for base, count, after in sorted(reversed(claims), key=lambda item: _name_key(item[0])):
                claim.execute(
                    "UPDATE system_name_registry SET occurrence_count = occurrence_count - ? "
                    "WHERE base_name = ? AND occurrence_count = ?", (count, base, after))
                claim.execute(
                    "DELETE FROM system_name_registry WHERE base_name = ? AND occurrence_count <= 0 "
                    "AND first_star_system_id IS NULL AND first_object_id IS NULL", (base,))
    except pymysql.err.Error as exc:
        log.debug(f"Could not give back {len(claims)} name claims: {exc}")


def reserve_system_names(conn, candidate_names):
    """
    System name-uniqueness reservation for many systems and uniquely
    named phenomena at once (PERF.14) -- a whole sector's worth in a few
    statements instead of three per name. Each name resolves a
    same-registry collision (`resolve_greek_roman_collision`, renaming
    the first holder where that calls for it), then a cross-level
    collision against an existing sector's base name (a diminutive prefix
    on *this* name, per `resolve_diminutive` -- never the sector). Planet
    and moon names are never searched: they're derived from the system
    name (`names/bodies.py`), so they're unique whenever it is. A name whose
    decorations are exhausted is drawn again (`_regenerate_star_name`)
    and goes round once more.

    Locking: the sector registry is read without locks (a locking read
    there waited on, and deadlocked with, other writers' uncommitted new
    sector names, a few times per 30 sectors with three workers), then one
    multi-row `INSERT ... ON DUPLICATE KEY UPDATE occurrence_count =
    occurrence_count + k` in sorted base-name order claims every base
    name (row locks, always taken in the same order; no gap locks, which
    is what deadlocked four parallel writers before), then the rows are
    read back to learn how many holders came before. A sector and a
    system given the same new base name by two writers in the same
    instant can both miss the diminutive; anything else that still
    deadlocks is retried by `save_sector`. The registry's first
    holder and diminutive index are written by `confirm_system_names`
    once the new rows' ids are known. A drawn-again name leaves its
    base's count one high (harmless: the next holder just gets a
    decoration it could have skipped).

    Args:
        conn (Connection): Part of the same transaction as the callers'
            own INSERTs.
        candidate_names (list): Freshly generated names, in insertion
            order. Two equal names here are told apart as two holders.

    Returns:
        list: `(final_name, base_name, diminutive_index)` per candidate,
            in order -- what the row's name should be, and what
            `confirm_system_names` should write.
    """
    names = list(candidate_names)
    results = [None] * len(names)
    todo = list(range(len(names)))
    own_claims = {}
    while todo:
        # Each name's key as this pass inserts it, taken once: a name
        # redrawn below (`_regenerate_star_name`) belongs to the next pass,
        # and must not count as a holder of whatever row its new spelling
        # happens to match later in this one (TEST.85: that row was then
        # short of holders, an existing count of -1).
        key_of = {i: _name_key(names[i]) for i in todo}
        counts = {}
        for i in todo:
            counts[key_of[i]] = counts.get(key_of[i], 0) + 1
        spelled = {}
        for i in todo:
            spelled.setdefault(key_of[i], names[i])
        keys = sorted(counts)
        claim = getattr(conn, "claim_connection", None)
        if claim is not None:
            # PERF.49: the registry rows are locked only for this short
            # transaction, not until the whole sector is saved.
            with claim:
                sector_hits, rows = _claim_names(claim, keys, counts, spelled)
        else:
            sector_hits, rows = _claim_names(conn, keys, counts, spelled)

        # Keys the collation still counts as one name (a spelling
        # `_name_key` folds differently) land on one registry row: they
        # are one name's holders, counted together.
        by_row = {}
        for key in keys:
            row = rows[key]
            by_row.setdefault(row["id"], (row, set()))[1].add(key)
        if claim is not None:
            conn.name_claims.extend((row["base_name"], sum(counts[key] for key in row_keys), row["occurrence_count"])
                                    for row, row_keys in by_row.values())
            conn.sector_name_keys.update(sector_hits)
        todo_by_key = {}
        for i in todo:
            todo_by_key.setdefault(key_of[i], []).append(i)

        retry = []
        for row, row_keys in by_row.values():
            uses = sorted(i for key in row_keys for i in todo_by_key.get(key, ()))
            existing_before = row["occurrence_count"] - len(uses)
            # (A name this call drew again counts the earlier pass's claim too: that is no other holder.)
            if claim is not None and existing_before - own_claims.get(row["id"], 0) > 0:
                conn.later_name_keys.update(row_keys)
            own_claims[row["id"]] = own_claims.get(row["id"], 0) + len(uses)
            diminutive_index = row["diminutive_index"]
            holder = "db" if existing_before > 0 else None
            for offset, i in enumerate(uses):
                base = names[i]
                new_name, rename = resolve_greek_roman_collision(
                    base, existing_before + offset, max_words=MAX_SYSTEM_NAME_WORDS)
                prefix, next_index = (None, diminutive_index)
                if new_name is not None and row_keys & sector_hits:
                    prefix, next_index = resolve_diminutive(diminutive_index)
                    if prefix is None or not fits_word_limit(f"{prefix} {new_name}", word_limit_for(base)):
                        new_name = None
                if new_name is None:
                    # Every decoration for this base is used, or would
                    # make the name longer than two words (GEN.46).
                    conn.burned_name_keys[key_of[i]] = conn.burned_name_keys.get(key_of[i], 0) + 1
                    names[i] = _regenerate_star_name()
                    retry.append(i)
                    continue
                diminutive_index = next_index
                if rename is not None:
                    # Matched by id, not by the rename tuple's assumed
                    # old name: the holder may carry a diminutive (an
                    # earlier system-vs-sector collision), which the
                    # Greek/Roman decoration replaces; uniqueness holds
                    # either way. The holder keeps its own spelling
                    # ("Vega" -> "Alpha Vega" when "vega" arrives).
                    if holder == "db":
                        _rename_registry_holder(
                            conn, row, resolve_greek_roman_collision(row["base_name"], existing_before + offset)[1][1])
                    elif holder is not None:
                        results[holder][0] = resolve_greek_roman_collision(names[holder], existing_before + offset)[1][1]
                if prefix is not None:
                    new_name = f"{prefix} {new_name}"
                results[i] = [new_name, base, diminutive_index]
                if holder is None:
                    holder = i
        todo = sorted(retry)
    return [tuple(result) for result in results]


def reserve_system_name(conn, candidate_name):
    """
    `reserve_system_names` for one name -- see there.

    Args:
        conn (Connection): Part of the same transaction as the caller's
            own system INSERT.
        candidate_name (str): The freshly generated name to reserve.

    Returns:
        tuple: `(final_name, base_name, diminutive_index)`.
    """
    return reserve_system_names(conn, [candidate_name])[0]


def confirm_system_names(conn, confirmations):
    """
    Second half of system name reservation (`reserve_system_names`):
    records each base name's first holder, now that the new rows' ids
    are known, and its diminutive index, in one multi-row upsert. A
    holder is only written to a row that has none yet (the row this
    transaction just made), so passing every new row's id is safe.

    Args:
        conn (Connection): The same transaction.
        confirmations (list): `(base_name, star_system_id, object_table,
            object_id, diminutive_index)` -- `star_system_id` for a
            system, `object_table`/`object_id` for a phenomenon (v40),
            the others `None`. The first entry per base name is its
            holder; the last sets its diminutive index.
    """
    claimed = getattr(conn, "claim_connection", None) is not None
    by_base = {}
    held = {}
    for base, system_id, table, object_id, diminutive_index in confirmations:
        key = _name_key(base)
        held[key] = held.get(key, 0) + 1
        if key in by_base:
            by_base[key][4] = diminutive_index
        else:
            if claimed and key in conn.later_name_keys:
                # Claimed after another holder (`reserve_system_names`): the first claimer is the holder.
                system_id = table = object_id = None
            by_base[key] = [base, system_id, table, object_id, diminutive_index]
    ordered = sorted(by_base.items())
    vacant = "first_star_system_id IS NULL AND first_object_id IS NULL"
    for first in range(0, len(ordered), _NAME_BATCH):
        chunk = [row for _key, row in ordered[first:first + _NAME_BATCH]]
        conn.execute(
            "INSERT INTO system_name_registry "
            "(base_name, occurrence_count, first_star_system_id, first_object_table, first_object_id, diminutive_index) "
            "VALUES " + ", ".join(["(?, 1, ?, ?, ?, ?)"] * len(chunk))
            + f" ON DUPLICATE KEY UPDATE first_object_table = IF({vacant}, VALUES(first_object_table), first_object_table),"
            f" first_star_system_id = IF({vacant}, VALUES(first_star_system_id), first_star_system_id),"
            f" first_object_id = IF({vacant}, VALUES(first_object_id), first_object_id),"
            " diminutive_index = VALUES(diminutive_index)",
            tuple(value for row in chunk for value in row),
        )
    return _rename_unconfirmed_holders(conn, ordered, held) if claimed else {}


_HOLDER_RENAME_SPAN = 40
"""int: Claims past which a first holder is never renamed again (the Greek tier hands off to the roman at 24)."""


def _rename_unconfirmed_holders(conn, ordered, held):
    """
    PERF.49: a name claimed on its own short transaction can be claimed
    again before its first holder saves, and the later claimer then has no
    holder row to rename ("Vega" to "Alpha Vega"). The holder renames itself
    here instead, once it is confirmed: whatever renames the claims since its
    own would have made. Returns `{key: new name}` for the holders renamed.
    """
    mine = [(key, row) for key, row in ordered if row[1] is not None or row[3] is not None]
    if not mine:
        return {}
    registry = _registry_rows(conn, [row[0] for _key, row in mine])
    healed = {}
    for key, (base, system_id, table, object_id, _diminutive_index) in mine:
        registered = registry.get(key)
        spent = held[key] + conn.burned_name_keys.get(key, 0)  # this transaction's own claims
        if registered is None or registered["occurrence_count"] <= spent:
            continue
        if (registered["first_star_system_id"], registered["first_object_id"]) != (system_id, object_id):
            continue
        target = None
        for later in range(spent, min(registered["occurrence_count"], _HOLDER_RENAME_SPAN)):
            new_name, rename = resolve_greek_roman_collision(base, later, max_words=MAX_SYSTEM_NAME_WORDS)
            if new_name is not None and key in conn.sector_name_keys:
                # A claimer that shares the base with a sector gets a diminutive as well, and
                # redraws its name instead when that is too long (`reserve_system_names`).
                prefix, _next_index = resolve_diminutive(None)
                if prefix is None or not fits_word_limit(f"{prefix} {new_name}", word_limit_for(base)):
                    new_name = None
            if new_name is not None and rename is not None:
                target = rename[1]
        if target is None:
            continue
        if system_id is not None:
            rename_star_system(conn, system_id, target)
        else:
            rename_phenomenon(conn, table, object_id, target)
        healed[key] = target
    return healed


def confirm_system_name(conn, base_name, star_system_id, diminutive_index):
    """`confirm_system_names` for one star system (or, when the caller
    holds a deferral list -- `insert_sector` -- adds it to that)."""
    _confirm_name(conn, (base_name, star_system_id, None, None, diminutive_index))


def confirm_object_name(conn, base_name, table, object_id, diminutive_index):
    """`confirm_system_name`'s counterpart for a uniquely named phenomenon
    (`NAMED_PHENOMENON_TABLES`, v40): the registry row, when this is the
    base name's first holder, points at `table`/`object_id`."""
    _confirm_name(conn, (base_name, None, table, object_id, diminutive_index))


def _confirm_name(conn, confirmation):
    if confirmation[0] is None:
        return  # named by its object ID (GEN.64): nothing to register
    deferred = getattr(conn, "deferred_name_confirmations", None)
    if deferred is not None:
        deferred.append(confirmation)
    else:
        confirm_system_names(conn, [confirmation])


def _take_name(conn, obj):
    """Reserves `obj.name` (`reserve_system_name`) and sets it to the
    final name -- or, when `insert_sector` already reserved it with the
    rest of its sector (`prereserved_names`), takes that. Returns
    `(base_name, diminutive_index)` for the confirm call."""
    prereserved = getattr(conn, "prereserved_names", None)
    if prereserved is not None and id(obj) in prereserved:
        return prereserved.pop(id(obj))
    final_name, base, diminutive_index = reserve_system_name(conn, obj.name)
    obj.name = final_name
    return base, diminutive_index


def _reserve_phenomenon_name(conn, obj, kind=None, placement=None):
    """
    Names a phenomenon and returns `(base_name, diminutive_index)` for
    `confirm_object_name`. A placed one (`placement` with a center) is
    named by its object ID (GEN.64, `_claim_object_ids`), or takes the ID
    `insert_sector` already claimed for it, and returns `(None, None)`:
    nothing to register. An unplaced one, or one whose name was given by
    hand (`name_given`), reserves its name (`_take_name`) as before.
    """
    prereserved = getattr(conn, "prereserved_names", None)
    if prereserved is not None and id(obj) in prereserved:
        return prereserved.pop(id(obj))
    if kind is not None and _placed(placement) and not getattr(obj, "name_given", False):
        _claim_object_ids(conn, [(obj, kind, _placement_center(placement))])
        return None, None
    return _take_name(conn, obj)


OBJECT_ID_TABLES = {
    **kind_tables.generator_type_tables(order="object_ids"),
    "black-hole-core": "black_holes",
    "neutron-star-core": "neutron_stars",
}
"""dict: `objectId.KIND_CODES` kind -> the table its objects (and their
ID names) live in."""


def _named_by_object_id(conn, obj, kind, placement):
    """Names a placed comet or asteroid field by its object ID (GEN.64),
    or takes the one `insert_sector` claimed for it. `False` for an
    unplaced one, which keeps its designation."""
    prereserved = getattr(conn, "prereserved_names", None)
    if prereserved is not None and id(obj) in prereserved:
        prereserved.pop(id(obj))
        return True
    if _placed(placement):
        _claim_object_ids(conn, [(obj, kind, _placement_center(placement))])
        return True
    return False


def _take_core_id(conn, core, kind, placement):
    """Names a placed supernova remnant core by its own object ID
    (GEN.64), or takes the one `insert_sector` claimed for it. An
    unplaced core keeps `"<remnant> Core"`."""
    prereserved = getattr(conn, "prereserved_names", None)
    if prereserved is not None and id(core) in prereserved:
        prereserved.pop(id(core))
    elif _placed(placement):
        _claim_object_ids(conn, [(core, kind, _placement_center(placement))])


def _placed(placement):
    """Whether `placement` gives a galaxy-frame center."""
    return bool(placement) and placement.get("center_x_pc") is not None


def _placement_center(placement):
    """`placement`'s `(x, y, z)` center in parsecs."""
    return tuple(placement[f"center_{axis}_pc"] for axis in "xyz")


def _claim_object_ids(conn, items, address=None):
    """
    Names each object by its 64-bit object ID (GEN.64, `objectId`):
    `items` is `(obj, kind, center_pc)` triples, in generation order, from
    the sector at `address` (`(shell, band, wedge)`, `None` for an object
    saved on its own). Objects of one kind that share an ID take the
    collision numbers 0, 1, 2 ... in the order of their sectors' addresses
    and then their generation order (GEN.307), whichever sector saved
    first: when a stored object of a higher address holds a number this
    one outranks, the stored one moves up to the next (`_rank_chain`). Sets
    each `obj.name`, and when `insert_sector` holds a reservation map
    (`prereserved_names`), records each object there as needing no
    registry entry.
    """
    by_table = {}
    for obj, kind, center in items:
        by_table.setdefault(OBJECT_ID_TABLES[kind], []).append((obj, objectId.pack(kind, center)))
    prereserved = getattr(conn, "prereserved_names", None)
    for table, entries in by_table.items():
        groups = {}
        for obj, object_id in entries:
            groups.setdefault(object_id, []).append(obj)
        taken = set()
        for base, objs in groups.items():
            for obj, name in zip(objs, _rank_chain(conn, table, base, len(objs), address, taken)):
                obj.name = name
                if prereserved is not None:
                    prereserved[id(obj)] = (None, None)


def _address_of(galaxy_position):
    """`(shell, band, wedge)` of a sector's galaxy position, `None` when it has none."""
    if not galaxy_position or galaxy_position.get("shell_index") is None:
        return None
    return (galaxy_position["shell_index"], galaxy_position["band_index"], galaxy_position["wedge_index"])


def _address_key(address):
    """The order of sector addresses: `(shell, band, wedge)`, an object saved on its own after every sector."""
    return (0,) + tuple(address) if address is not None else (1, 0, 0, 0)


def _chain(base, length):
    """The first `length` IDs of `base`'s collision chain (`objectId.bump`) as names."""
    names, object_id = [], base
    for _ in range(length):
        names.append(objectId.format_id(object_id))
        object_id = objectId.bump(object_id)
    return names


def _rank_chain(conn, table, base, count, address, taken):
    """
    The names of `count` new objects whose ID is `base`, from one sector at `address`: the first free numbers of the
    chain, except where a stored object of a higher address holds one: the numbers are then given out in the order of
    the addresses (and, in one sector, of the order the objects were stored in), so the stored objects of the higher
    addresses move up past the new ones. `taken` holds the names this call has already handed out; this adds its own.
    """
    length = count + 16
    while True:
        names = [name for name in _chain(base, length) if name not in taken]
        holders = conn.execute(
            f"SELECT t.id, t.name, s.shell_index, s.band_index, s.wedge_index FROM {table} t"
            f" LEFT JOIN sectors s ON s.id = t.sector_id WHERE t.name IN ({', '.join('?' * len(names))})",
            tuple(names)).fetchall()
        if len(holders) + count < len(names):
            break
        length *= 2

    def rank(row):
        stored = None if row["shell_index"] is None else (row["shell_index"], row["band_index"], row["wedge_index"])
        return _address_key(stored), row["id"]

    mine = _address_key(address)
    lower = sorted((row for row in holders if rank(row)[0] <= mine), key=rank)
    higher = sorted((row for row in holders if rank(row)[0] > mine), key=rank)
    if not higher:
        held = {row["name"].upper() for row in holders}
        given = [name for name in names if name not in held][:count]
    else:
        given = names[len(lower):len(lower) + count]
        for row, name in zip(lower + higher, names[:len(lower)] + names[len(lower) + count:]):
            if row["name"].upper() != name:
                conn.execute(f"UPDATE {table} SET name = ? WHERE id = ?", (name, row["id"]))
    taken.update(given)
    return given


_AU_KM = _au_to_km(1.0)


_SPIN_FIELD = object_rows.computed(spin.SPIN_FIELDS, lambda obj, context: spin.spin_values(obj),
                                   lambda row: {name: row[name] for name in spin.SPIN_FIELDS})
"""Field: the spin axis and axial tilt."""

_GALACTIC_ORBIT_FIELD = object_rows.computed(
    ("galactic_orbital_period_gy", "galactic_orbital_phase_deg", "galactic_min_update_interval_years"),
    lambda obj, context: (obj.galactic_orbital_period_gy, obj.galactic_orbital_phase_deg,
                          obj.galactic_min_update_interval_years),
    lambda row: {name: row[name] for name in ("galactic_orbital_period_gy", "galactic_orbital_phase_deg",
                                              "galactic_min_update_interval_years")})
"""Field: the galactic orbit of a star or free object."""

_PLACEMENT_FIELD = object_rows.computed(
    ("center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc"),
    lambda obj, context: _placement_values(context.get("placement")))
"""Field: where a free object sits in the galaxy (all `None` unplaced)."""

_FREE_VECTOR_FIELD = object_rows.computed(
    ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
    lambda obj, context: _free_vector(obj, context.get("placement")))
"""Field: the movement vector a free object is stored with."""


def _au_pair(inner_column, outer_column, name):
    """A field storing the AU pair `obj.<name>` (a habitable zone) as two kilometre columns."""
    return object_rows.computed(
        (inner_column, outer_column),
        lambda obj, context: (getattr(obj, name)[0] * _AU_KM, getattr(obj, name)[1] * _AU_KM),
        lambda row: {name: [row[inner_column] / _AU_KM, row[outer_column] / _AU_KM]})


def _flag(column):
    """A column holding the attribute of its name as 0/1; it loads back as a `bool`."""
    return object_rows.computed(column, lambda obj, context: int(getattr(obj, column)),
                                lambda row: {column: bool(row[column])})


def _whole_number(value):
    """A stored `REAL` that is a whole number as the `int` the generator makes it (`5800`, not `5800.0`).

    `generate_star` sets the temperature with `int(round(...))`, but a `DOUBLE` column returns every number as a
    float and `f"{temperature} K"` renders the two differently. A fractional value is kept (an accretion disk's
    temperature is drawn as a float and truncating it would print 3,676,064.7 K as 3,676,064 K).
    """
    return int(value) if value is not None and float(value).is_integer() else value


FREE_PARENT = ("free", None)
"""tuple: The `(parent_type, parent_id)` of an object that orbits nothing (DB.73)."""


def star_parent(parent_node_id, star_system_id):
    """`(parent_type, parent_id)` of a star: the orbit node it hangs on, or the system for the head of one."""
    return ("node", parent_node_id) if parent_node_id is not None else ("system", star_system_id)


def anchored_parent(star_id):
    """`(parent_type, parent_id)` of a remnant: the star it is anchored to, else free."""
    return ("star", star_id) if star_id is not None else FREE_PARENT


STAR_FIELDS = (
    object_rows.context("star_system_id"), object_rows.context("role"),
    object_rows.attr("name"), object_rows.attr("star_type", "type"), object_rows.attr("yerkes_class"),
    object_rows.attr("mass_kg", "mass"), object_rows.attr("radius_km", "radius"),
    object_rows.computed("temperature_k", lambda star, context: star.temperature,
                         lambda row: {"temperature": _whole_number(row["temperature_k"])}),
    object_rows.attr("luminosity_w", "luminosity"),
    object_rows.attr("age_gy", "age"),
    object_rows.computed("lifespan_gy", lambda star, context: _lifespan_gy(star.lifespan),
                         lambda row: {"lifespan": row["lifespan_gy"]}),
    _au_pair("habitable_zone_inner_km", "habitable_zone_outer_km", "habitable_zone"),
    object_rows.scaled("system_perimeter_km", "system_perimeter", _AU_KM),
    object_rows.scaled("heliosphere_radius_km", "heliosphere_radius", _AU_KM),
    _GALACTIC_ORBIT_FIELD,
    object_rows.scaled("wide_binary_a_crit_km", "a_crit_au", _AU_KM, nullable=True),
    object_rows.scaled("reflex_offset_x_km", "reflex_offset_x", _AU_KM, null=0.0, nullable=True),
    object_rows.scaled("reflex_offset_y_km", "reflex_offset_y", _AU_KM, null=0.0, nullable=True),
    object_rows.scaled("reflex_offset_z_km", "reflex_offset_z", _AU_KM, null=0.0, nullable=True),
    object_rows.optional("rotation_period_hours"), _SPIN_FIELD,
    object_rows.each(activity.STAR_ACTIVITY_FIELDS), object_rows.optional("metallicity_fe_h"),
    object_rows.optional("state"), object_rows.optional("decay_age_yr"),
    object_rows.scaled("snow_line_km", "snow_line_au", _AU_KM),
    object_rows.each(pre_main_sequence.DISK_FIELDS),
    object_rows.context("parent_node_id"), object_rows.context("slot"),
    object_rows.computed(("parent_type", "parent_id"),
                         lambda star, context: star_parent(context["parent_node_id"], context["star_system_id"])),
    object_rows.computed(("offset_x_au", "offset_y_au", "offset_z_au"), lambda star, context: tuple(context["offset_au"])),
    object_rows.computed("flags", lambda star, context: context["flags"]),
)
"""tuple: The `Field`s of a `stars` row (`object_rows`)."""


def insert_star(conn, star, star_system_id, role, parent_node_id=None, offset_au=(0.0, 0.0, 0.0), slot=None, flags=0) -> int:
    """
    Inserts a `stars` row for one individual `Star` (never a
    `BinaryStarProxy` -- see `insert_star_system`).

    Args:
        conn (Connection): An open, schema-initialized connection.
        star: The `Star` instance.
        star_system_id (int): The owning `star_systems.id`.
        role (str): `'primary'`, `'secondary'`, `'member'` (a star beyond the pair of a system of more than two, GEN.630)
            or `'single'`.
        parent_node_id (int, optional): The `orbit_nodes.node_no` the star is a child of; `None` for a single star.
        offset_au (tuple, optional): The star's `(x, y, z)` offset from its system's centre of mass, AU (GEN.311).
        slot (int, optional): Which child of its parent node the star is, 0 or 1; `None` for a single star (GEN.630).
        flags (int, optional): The star's member flag bits (`generation.multistar.MEMBER_FLAGS`); a runaway star's primary
            carries the bit of its origin (GEN.329); 0 for none.

    Returns:
        int: The new `stars.id`.
    """
    return object_rows.insert_row(conn, "stars", STAR_FIELDS, star, {
        "star_system_id": star_system_id, "role": role, "parent_node_id": parent_node_id, "offset_au": offset_au,
        "slot": slot, "flags": flags})


def _insert_reflection_spectrum(conn, table, id_column, owner_id, visible, non_visible):
    """Shared body for `planet_reflection_spectrum`/
    `moon_reflection_spectrum` -- identical shape, different owning
    table/column (always this module's own literal constants, never request-derived)."""
    for spectrum_type, values in (("visible", visible), ("non_visible", non_visible)):
        for position, value in enumerate(values or []):
            conn.execute(
                f"INSERT INTO {table} ({id_column}, spectrum_type, position, value) VALUES (?, ?, ?, ?)",
                (owner_id, spectrum_type, position, value),
            )


def _radiation_values(body):
    return tuple(getattr(body, name, None) for name in radiation.BODY_FIELDS)


def _radiation_load(row):
    data = {name: row[name] for name in radiation.BODY_FIELDS}
    data["ozone_loss_flag"] = _tristate_from_db(row["ozone_loss_flag"])
    return data


_POSITION_COLUMNS = point_columns.LAYOUTS["planets"].position
_VELOCITY_COLUMNS = point_columns.LAYOUTS["planets"].velocity

BODY_FIELDS = (
    object_rows.attr("body_type"), object_rows.attr("name"),
    object_rows.attr("planet_class"), object_rows.scaled("distance_km", "distance", _AU_KM),
    object_rows.attr("radius_km", "radius"), object_rows.attr("mass_kg", "mass"),
    object_rows.attr("volume_km3", "volume"), object_rows.attr("period_years", "period"), object_rows.attr("zone"),
    object_rows.attr("description"), object_rows.attr("gravity_g", "gravity"),
    object_rows.attr("surface_temperature_k", "surface_temperature"),
    object_rows.attr("density_g_cm3", "density"), object_rows.attr("atmosphere"),
    object_rows.attr("atm_density"), object_rows.attr("atm_molar_density"),
    object_rows.attr("atmospheric_pressure_pa", "atmospheric_pressure"), object_rows.attr("composition"),
    object_rows.attr("scale_height_km", "scale_height"), object_rows.attr("hill_radius_km", "hill_radius"),
    object_rows.scaled("min_orbit_distance_km", "min_orbit_distance", _AU_KM, nullable=True),
    _au_pair("habitable_zone_inner_km", "habitable_zone_outer_km", "habitable_zone"),
    object_rows.each(life_classes.BODY_FIELDS),
    object_rows.attr("life_stage"), object_rows.attr("evolutionary_pace"),
    object_rows.attr("life_milestone_age_gyr"), object_rows.attr("life_system_age_gyr"),
    *(object_rows.attr(field) for field in tuning.LIFE_SCHEDULE_FIELDS),
    object_rows.attr("flavor_key"), object_rows.attr("flavor_text_count"),
    object_rows.attr("orbital_inclination_deg"), object_rows.attr("orbital_ascending_node_deg"), object_rows.attr("eccentricity"),
    *(object_rows.scaled(column, f"position_{axis}", _AU_KM) for axis, column in zip("xyz", _POSITION_COLUMNS)),
    object_rows.attr("orbital_speed_kms"),
    *(object_rows.attr(column, f"velocity_{axis}_kms") for axis, column in zip("xyz", _VELOCITY_COLUMNS)),
    object_rows.attr("min_update_interval_years"),
    object_rows.attr("rotation_period_hours"), _SPIN_FIELD,
    object_rows.each(atmosphere.ATMOSPHERE_FIELDS),
    object_rows.each(climate_draws.CLIMATE_FIELDS),
    object_rows.each(magnetism.BODY_FIELDS),
    object_rows.each(hydrosphere.HYDROSPHERE_FIELDS),
    object_rows.each(radiogenic.BODY_FIELDS),
    object_rows.each((*interior_heat.BODY_FIELDS, *interior_heat.OUTPUT_FIELDS, *interior_heat.REGIME_FIELDS)),
    object_rows.each(bulk_composition.COLUMN_FIELDS),
    object_rows.each(surface_composition.COLUMN_FIELDS),
    object_rows.each(daylight.FIELDS),
    object_rows.each(ocean_tides.FIELDS),
    object_rows.each(moon_composition.BODY_FIELDS),
    object_rows.computed("carbon_rock", lambda body, context: _tristate(body.carbon_rock),
                         lambda row: {"carbon_rock": _tristate_from_db(row["carbon_rock"])}),
    object_rows.each(tuple(name for name in giant_interior.BODY_FIELDS if name != "has_helium_shell")),
    object_rows.each(rings.BODY_FIELDS),
    object_rows.each(body_colour.BODY_FIELDS),
    object_rows.computed("has_helium_shell", lambda body, context: _tristate(body.has_helium_shell),
                         lambda row: {"has_helium_shell": _tristate_from_db(row["has_helium_shell"])}),
    object_rows.computed(radiation.BODY_FIELDS, lambda body, context: _radiation_values(body), _radiation_load),
    object_rows.each(habitability_world.BODY_FIELDS),
    object_rows.each(aerosols.BODY_FIELDS),
    object_rows.computed("born_beyond_snow_line", lambda body, context: _tristate(body.born_beyond_snow_line),
                         lambda row: {"born_beyond_snow_line": _tristate_from_db(row["born_beyond_snow_line"])}),
)
"""tuple: The `Field`s `planets` and `moons` share (`insert_planet`/`insert_moon` write them, `db.edits` updates
them in place, `_planet_or_moon_row_to_dict` reads them back)."""

PLANET_ONLY_FIELDS = tuple(
    object_rows.scaled(f"reflex_offset_{axis}_km", f"reflex_offset_{axis}", _AU_KM, null=0.0, nullable=True)
    for axis in "xyz")
"""tuple: The `planets` fields `moons` lacks (a moon hosts no moons)."""

PLANET_BODY_FIELDS = BODY_FIELDS + PLANET_ONLY_FIELDS


INTERIOR_LAYER_COLUMNS = (*giant_interior.LAYER_FIELDS, *giant_interior.GAS_MATERIALS, *giant_interior.HEAVY_MATERIALS)
"""tuple: The value columns of `planet_interior_layers`/`moon_interior_layers` (DB.89), after the owner id and `layer_index`."""


def _insert_interior_layers(conn, table, id_column, owner_id, layers):
    """Shared body for `planet_interior_layers`/`moon_interior_layers` (this module's own literal table names)."""
    columns = ", ".join(INTERIOR_LAYER_COLUMNS)
    marks = ", ".join("?" for _ in range(len(INTERIOR_LAYER_COLUMNS) + 2))
    for index, layer in enumerate(layers or []):
        conn.execute(
            f"INSERT INTO {table} ({id_column}, layer_index, {columns}) VALUES ({marks})",
            (owner_id, index, *(layer[name] for name in INTERIOR_LAYER_COLUMNS)),
        )


def body_child_rows(conn, body, body_id):
    """Writes a planet's or moon's reflection spectrum, interior layer rows, gas rows and life class score rows (after any old ones were deleted)."""
    prefix, id_column = ("moon", "moon_id") if body.is_moon else ("planet", "planet_id")
    _insert_reflection_spectrum(
        conn, f"{prefix}_reflection_spectrum", id_column, body_id,
        body.reflection_spectrum_visible, body.reflection_spectrum_non_visible,
    )
    _insert_interior_layers(conn, f"{prefix}_interior_layers", id_column, body_id, body.interior_layers)
    for formula, partial_pa in sorted(atmosphere.gas_partials_pa(body).items()):
        if partial_pa > 0.0:
            conn.execute(f"INSERT INTO {prefix}_gases ({id_column}, gas, partial_pa) VALUES (?, ?, ?)",
                         (body_id, formula, partial_pa))
    for class_id, score in sorted((getattr(body, life_classes.SCORES_FIELD, None) or {}).items()):
        if score > 0:
            conn.execute(f"INSERT INTO {prefix}_life_class_scores ({id_column}, life_class, score) VALUES (?, ?, ?)",
                         (body_id, class_id, score))
    for entry in getattr(body, life_classes.LIVING_FIELD, None) or ():
        ages = list(entry.get("ages") or [])
        ages += [None] * (len(tuning.LIFE_SCHEDULE_FIELDS) - len(ages))
        kinds = list(entry.get("kinds") or [None] * 4)
        conn.execute(f"INSERT INTO {prefix}_life ({id_column}, {', '.join(LIVING_COLUMNS)}) VALUES (?, {', '.join('?' * len(LIVING_COLUMNS))})",
                     (body_id, entry["life_class"], entry["share"], entry.get("stage"), entry.get("milestone_age_gyr"), entry.get("system_age_gyr"),
                      entry.get("animal_grade"), entry.get("complexity"), *ages, *kinds))


def body_gases(conn, body_ids, is_moon):
    """`{body id: {gas formula: partial pressure Pa}}` of the planets or moons `body_ids`: the air their rows keep in
    `planet_gases` or `moon_gases` (a body with none has no key)."""
    table, id_column = ("moon_gases", "moon_id") if is_moon else ("planet_gases", "planet_id")
    air = {}
    ids = list(body_ids)
    for start in range(0, len(ids), 1000):
        chunk = ids[start:start + 1000]
        marks = ", ".join("?" * len(chunk))
        for row in conn.execute(
                f"SELECT {id_column}, gas, partial_pa FROM {table} WHERE {id_column} IN ({marks})", tuple(chunk)).fetchall():
            air.setdefault(row[id_column], {})[row["gas"]] = row["partial_pa"]
    return air


def body_class_scores(conn, body_ids, is_moon):
    """`{body id: {class id: score 1 to 9}}` of the planets or moons `body_ids`: the life classes their rows keep in
    `planet_life_class_scores` or `moon_life_class_scores` (a body that supports none has no key)."""
    table, id_column = ("moon_life_class_scores", "moon_id") if is_moon else ("planet_life_class_scores", "planet_id")
    found = {}
    ids = list(body_ids)
    for start in range(0, len(ids), 1000):
        chunk = ids[start:start + 1000]
        marks = ", ".join("?" * len(chunk))
        for row in conn.execute(
                f"SELECT {id_column}, life_class, score FROM {table} WHERE {id_column} IN ({marks})", tuple(chunk)).fetchall():
            found.setdefault(row[id_column], {})[row["life_class"]] = row["score"]
    return found


LIVING_COLUMNS = ("life_class", "share", "life_stage", "life_milestone_age_gyr", "life_system_age_gyr", "life_animal_grade", "life_complexity",
                  *tuning.LIFE_SCHEDULE_FIELDS, "life_microbial", "life_flora", "life_fauna", "life_sentient")
"""tuple: The columns of `planet_life` and `moon_life` after the body's id, in the order `body_child_rows` writes them (GEN.1209)."""


def living_entry(row):
    """The `life_living` entry (`physics.life_classes.LIVING_FIELD`) of a `planet_life` or `moon_life` row."""
    return {"life_class": row["life_class"], "share": row["share"], "stage": row["life_stage"], "milestone_age_gyr": row["life_milestone_age_gyr"],
            "system_age_gyr": row["life_system_age_gyr"], "animal_grade": row["life_animal_grade"], "complexity": row["life_complexity"],
            "ages": [row[field] for field in tuning.LIFE_SCHEDULE_FIELDS],
            "kinds": [row["life_microbial"], row["life_flora"], row["life_fauna"], row["life_sentient"]]}


def body_living(conn, body_ids, is_moon):
    """`{body id: [life_living entries, the dominant class first]}` of the planets or moons `body_ids`: the life classes their rows keep in `planet_life`
    or `moon_life` (a body with none has no key). The dominant class is the one with the largest share; the rest follow by share."""
    table, id_column = ("moon_life", "moon_id") if is_moon else ("planet_life", "planet_id")
    found = {}
    ids = list(body_ids)
    for start in range(0, len(ids), 1000):
        chunk = ids[start:start + 1000]
        marks = ", ".join("?" * len(chunk))
        for row in conn.execute(
                f"SELECT {id_column}, {', '.join(LIVING_COLUMNS)} FROM {table} WHERE {id_column} IN ({marks}) ORDER BY {id_column}, share DESC, life_class",
                tuple(chunk)).fetchall():
            found.setdefault(row[id_column], []).append(living_entry(row))
    return found


def dominant_first(entries, number, variant):
    """The `life_living` entries `entries` (or `None`) with the dominant class, from the body's `life_class` number and `life_class_variant` letter, first."""
    if not entries:
        return None
    dominant = life_classes.class_id(number, variant)
    return sorted(entries, key=lambda entry: entry["life_class"] != dominant)


def with_life_inputs(conn, rows, is_moon):
    """`with_gases` rows that also hold the temperature of the star each orbits as `star_teff_k` (`None` for a circumbinary planet), the
    one input of a life class re-score (`habitability_world.compute`) a `planets` or `moons` row lacks."""
    rows = with_gases(conn, rows, is_moon)
    star_ids = sorted({row["star_id"] for row in rows if row.get("star_id") is not None})
    teff = {}
    for start in range(0, len(star_ids), 1000):
        chunk = star_ids[start:start + 1000]
        for star in conn.execute(f"SELECT id, temperature_k FROM stars WHERE id IN ({', '.join('?' * len(chunk))})", tuple(chunk)).fetchall():
            teff[star["id"]] = star["temperature_k"]
    return [{**row, "star_teff_k": teff.get(row.get("star_id"))} for row in rows]


def with_gases(conn, rows, is_moon):
    """`rows` (dicts from `SELECT *` on `planets` or `moons`) each as a copy that also holds its air under `atmosphere.GAS_FIELD`, the
    key the habitability, class and colour code read (the same as a `Planet` attribute)."""
    air = body_gases(conn, [row["id"] for row in rows], is_moon)
    return [{**dict(row), atmosphere.GAS_FIELD: air.get(row["id"], {})} for row in rows]


def planet_parent(star_id, star_system_id):
    """`(parent_type, parent_id)` of a planet (GEN.289): the star it orbits, or, for a close pair's circumbinary
    planet (`star_id` None), the system, whose barycenter it orbits."""
    return ("star", star_id) if star_id is not None else ("system", star_system_id)


PLANET_FIELDS = (
    object_rows.context("star_system_id"), object_rows.context("star_id"), object_rows.context("orbital_index"),
    object_rows.computed(("parent_type", "parent_id"),
                         lambda planet, context: planet_parent(context["star_id"], context["star_system_id"])),
    *PLANET_BODY_FIELDS,
)
"""tuple: The `Field`s of a `planets` row."""

def moon_parent(context):
    """`(parent_type, parent_id)` of a moon: its planet, or the rogue planet it orbits when it has no planet."""
    if context.get("rogue_planet_id") is not None:
        return ("rogue_planet", context["rogue_planet_id"])
    return ("planet", context["planet_id"])


MOON_FIELDS = (
    object_rows.context("planet_id"), object_rows.context("star_system_id"), object_rows.context("star_id"),
    object_rows.context("orbital_index"),
    object_rows.Field(("rogue_planet_id",), lambda moon, context: context.get("rogue_planet_id")),
    object_rows.computed(("parent_type", "parent_id"), lambda moon, context: moon_parent(context)),
    *BODY_FIELDS,
)
"""tuple: The `Field`s of a `moons` row."""


def insert_planet(conn, planet, star_system_id, star_id, orbital_index) -> int:
    """
    Inserts a `planets` row for one top-level `Planet` (never a moon --
    schema v2 gives moons their own table, see `insert_moon`), then
    inserts each of its moons.

    Args:
        conn (Connection): An open, schema-initialized connection.
        planet (Planet): The top-level planet to persist.
        star_system_id (int): The owning `star_systems.id`.
        star_id (int or None): The specific `stars.id` this planet orbits,
                               or `None` for a binary system (see the
                               `planets.star_id` column comment in
                               `schema.sql`). Threaded unchanged into every
                               `insert_moon` call for this planet's moons.
        orbital_index (int): This planet's position in the star's ordered
                             `planets` list.

    Returns:
        int: The new `planets.id`.
    """
    planet_id = object_rows.insert_row(conn, "planets", PLANET_FIELDS, planet, {
        "star_system_id": star_system_id, "star_id": star_id, "orbital_index": orbital_index})

    body_child_rows(conn, planet, planet_id)

    for position, moon in enumerate(planet.moons or []):
        insert_moon(conn, moon, star_system_id, star_id, planet_id, position)

    return planet_id


def insert_moon(conn, moon, star_system_id, star_id, planet_id, orbital_index) -> int:
    """
    Inserts a `moons` row for one moon (a `Planet` instance with
    `is_moon=True`). Unlike `insert_planet`, this never recurses --
    `Planet.__init__` only calls `generate_moons` `if not self.is_moon`,
    so a moon never has moons of its own.

    Args:
        conn (Connection): An open, schema-initialized connection.
        moon (Planet): The moon to persist.
        star_system_id (int): The owning `star_systems.id` (same value the
                              parent planet was inserted with).
        star_id (int or None): Same value the parent planet was inserted
                               with -- see `insert_planet`'s docstring.
        planet_id (int): The owning `planets.id` -- the planet this moon
                         orbits.
        orbital_index (int): This moon's position in its parent planet's
                             `moons` list.

    Returns:
        int: The new `moons.id`.
    """
    moon_id = object_rows.insert_row(conn, "moons", MOON_FIELDS, moon, {
        "planet_id": planet_id, "star_system_id": star_system_id, "star_id": star_id, "orbital_index": orbital_index,
        "rogue_planet_id": None})

    body_child_rows(conn, moon, moon_id)

    return moon_id


PARENT_TYPES = ("star", "planet", "system", "node", "rogue_planet", "free")
"""tuple: What a `parent_type` can name (DB.73): the table of its `parent_id` is `star` -> `stars`, `planet` -> `planets`,
`system` -> `star_systems`, `rogue_planet` -> `rogue_planets`; a `node` parent is an `orbit_nodes.node_no` of the row's own
system, and a `free` one has no `parent_id`."""

PARENT_LINK_TABLES = ("stars", "planets", "moons", "asteroid_belts", "comets", "black_holes", "neutron_stars",
                      "rogue_planets", "interstellar_comets")
"""tuple: The kind tables that say what they orbit in `parent_type` and `parent_id`."""

_PARENT_ROW_TABLES = {"star": "stars", "planet": "planets", "system": "star_systems", "rogue_planet": "rogue_planets"}


def set_parent(conn, table, row_id, parent_type, parent_id=None):
    """
    Changes what one object orbits (DB.73): its row keeps its `id` and its `uid`, and so does everything else about it,
    so a planet ejected from its star is the same planet, now `free`.

    Args:
        conn (Connection): Part of the caller's transaction.
        table (str): One of `PARENT_LINK_TABLES`.
        row_id (int): The row.
        parent_type (str): One of `PARENT_TYPES`.
        parent_id (int, optional): The parent's row id (a `node`'s `node_no`); `None` exactly for `free`.

    Raises:
        ValueError: An unknown table, row or type, a parent that is not there, a `free` parent with an id or any other
            without one, or an object made its own parent.
    """
    if table not in PARENT_LINK_TABLES:
        raise ValueError(f"{table} has no parent link")
    if parent_type not in PARENT_TYPES:
        raise ValueError(f"no parent type {parent_type!r}; one of {', '.join(PARENT_TYPES)}")
    if (parent_type == "free") != (parent_id is None):
        raise ValueError("a free object has no parent id and every other parent has one")
    row = conn.execute(f"SELECT * FROM {table} WHERE id = ?", (row_id,)).fetchone()
    if row is None:
        raise ValueError(f"no {table} row with id {row_id}")
    if parent_type in _PARENT_ROW_TABLES:
        parent_table = _PARENT_ROW_TABLES[parent_type]
        if parent_table == table and parent_id == row_id:
            raise ValueError("an object cannot orbit itself")
        if conn.execute(f"SELECT 1 FROM {parent_table} WHERE id = ?", (parent_id,)).fetchone() is None:
            raise ValueError(f"no {parent_table} row with id {parent_id}")
    elif parent_type == "node":
        system_id = row["star_system_id"] if "star_system_id" in row.keys() else None
        if system_id is None or conn.execute("SELECT 1 FROM orbit_nodes WHERE star_system_id = ? AND node_no = ?",
                                             (system_id, parent_id)).fetchone() is None:
            raise ValueError(f"no orbit node {parent_id} in this object's system")
    extra, params = "", [parent_type, parent_id]
    if table == "moons":
        extra = ", rogue_planet_id = ?"
        params.append(parent_id if parent_type == "rogue_planet" else None)
    conn.execute(f"UPDATE {table} SET parent_type = ?, parent_id = ?{extra} WHERE id = ?", (*params, row_id))


def insert_moon_of_rogue_planet(conn, moon, rogue_planet_id, orbital_index) -> int:
    """
    Inserts a `moons` row for a moon of a free planet (DB.73): one with no star, system or planet, whose parent is the
    `rogue_planets` row `rogue_planet_id`. Deleting that rogue planet, or its sector, deletes the moon.

    Args:
        conn (Connection): An open, schema-initialized connection.
        moon (Planet): The moon to persist (a `Planet` with `is_moon=True`).
        rogue_planet_id (int): The `rogue_planets.id` it orbits.
        orbital_index (int): Its position in the rogue planet's moons.

    Returns:
        int: The new `moons.id`.
    """
    moon_id = object_rows.insert_row(conn, "moons", MOON_FIELDS, moon, {
        "planet_id": None, "star_system_id": None, "star_id": None, "orbital_index": orbital_index,
        "rogue_planet_id": rogue_planet_id})
    body_child_rows(conn, moon, moon_id)
    return moon_id


def load_moons_of_rogue_planet(conn, rogue_planet_id, system_config=None):
    """
    The moons of the rogue planet `rogue_planet_id`, in orbital order, rebuilt as `insert_moon_of_rogue_planet` stored
    them, each with its row id as `db_id`. A moon has no star to draw light from, so none is threaded in.
    """
    rows = conn.execute("SELECT * FROM moons WHERE rogue_planet_id = ? ORDER BY orbital_index, id",
                        (rogue_planet_id,)).fetchall()
    config = system_config if system_config is not None else SystemConfig()
    moons = []
    for row in rows:
        moon = Planet.from_dict(_planet_or_moon_row_to_dict(conn, row, is_moon=True), None, config)
        moon.db_id = row["id"]
        moons.append(moon)
    return moons


BELT_FIELDS = (
    object_rows.context("star_system_id"), object_rows.context("star_id"), object_rows.context("orbital_index"),
    object_rows.computed(("parent_type", "parent_id"),
                         lambda belt, context: planet_parent(context["star_id"], context["star_system_id"])),
    object_rows.scaled("distance_km", "distance", _AU_KM), object_rows.scaled("lower_limit_km", "lower_limit", _AU_KM),
    object_rows.scaled("upper_limit_km", "upper_limit", _AU_KM),
    object_rows.attr("density"), object_rows.attr("mass_kg"), object_rows.optional("chondrite_group"),
    object_rows.computed("composition_summary", lambda belt, context: belt.get_composition_summary()),
    # GEN.337: a belt is a ring about its star, so its point is the star's (the system centre) and it carries no vector
    # of its own.
    object_rows.computed(("position_x_km", "position_y_km", "position_z_km", "velocity_x_kms", "velocity_y_kms",
                          "velocity_z_kms"), lambda belt, context: (0.0,) * 6),
)
"""tuple: The `Field`s of an `asteroid_belts` row."""


def insert_asteroid_belt(conn, belt: AsteroidBelt, star_system_id, orbital_index, star_id=None) -> int:
    """
    Inserts an `asteroid_belts` row (plus its `asteroid_belt_composition`
    child rows).

    Args:
        conn (Connection): An open, schema-initialized connection.
        belt (AsteroidBelt): The belt to persist.
        star_system_id (int): The owning `star_systems.id`.
        orbital_index (int): This belt's position in its owning star's
                             `planets`/`secondary_planets` list (shared
                             index space with `Planet` entries within that
                             same list, so orbital order across both types
                             is preserved -- see `insert_star_system` for
                             how a 'wide' binary's two lists each restart
                             this index at 0, disambiguated by `star_id`).
        star_id (int, optional): The specific owning `stars.id`, same
            semantics as `planets.star_id` (see that column's own schema
            comment) -- `None` for a single star's or a 'close' binary's
            belt, set for a 'wide' binary's.

    Returns:
        int: The new `asteroid_belts.id`.
    """
    belt_id = object_rows.insert_row(conn, "asteroid_belts", BELT_FIELDS, belt, {
        "star_system_id": star_system_id, "star_id": star_id, "orbital_index": orbital_index})

    for position, (component, concentration) in enumerate(belt.composition):
        conn.execute(
            """
            INSERT INTO asteroid_belt_composition (belt_id, position, component, concentration)
            VALUES (?, ?, ?, ?)
            """,
            (belt_id, position, component, concentration),
        )

    return belt_id


COMET_FIELDS = (
    object_rows.context("star_system_id"), object_rows.context("star_id"),
    object_rows.computed(("parent_type", "parent_id"),
                         lambda comet, context: planet_parent(context["star_id"], context["star_system_id"])),
    object_rows.attr("name"), object_rows.attr("orbit_type"), object_rows.attr("period_class"),
    object_rows.attr("nucleus_diameter_km"),
    object_rows.computed("composition_summary", lambda comet, context: comet.get_composition_summary()),
    object_rows.scaled("perihelion_distance_km", "perihelion_distance_au", _AU_KM),
    object_rows.attr("eccentricity"), object_rows.attr("inclination_deg"), object_rows.attr("arg_periapsis_deg"),
    object_rows.attr("ascending_node_deg"), object_rows.attr("orbital_period_years"),
    object_rows.attr("mean_anomaly_deg"), object_rows.attr("parabolic_mean_anomaly"),
    object_rows.attr("min_update_interval_years"), object_rows.attr("primary_mass_solar"),
    _flag("is_active"),
    object_rows.scaled("distance_km", "distance_au", _AU_KM), object_rows.scaled("position_x_km", "position_x_au", _AU_KM),
    object_rows.scaled("position_y_km", "position_y_au", _AU_KM), object_rows.scaled("position_z_km", "position_z_au", _AU_KM),
    object_rows.attr("orbital_speed_kms"),
    object_rows.attr("velocity_x_kms"), object_rows.attr("velocity_y_kms"), object_rows.attr("velocity_z_kms"),
    object_rows.optional("rotation_period_hours"), _SPIN_FIELD,
    object_rows.attr("mass_kg"),
)
"""tuple: The `Field`s of a `comets` row."""


def insert_comet(conn, comet: Comet, star_system_id, star_id=None) -> int:
    """
    Inserts a `comets` row (plus its `comet_composition` child rows).

    Modeled directly on `insert_asteroid_belt` above -- a `Comet` has no
    `orbital_index` (it isn't part of the planets/belts orbital-slot
    ordering; see `system.StarSystem._generate_comets`'s docstring),
    so this only needs `star_system_id`/`star_id`, not a position within
    another list.

    Args:
        conn (Connection): An open, schema-initialized connection.
        comet (Comet): The comet to persist.
        star_system_id (int): The owning `star_systems.id`.
        star_id (int, optional): The specific owning `stars.id`, same
            `planets.star_id` real semantics (see that column's own schema
            comment) -- a real `stars.id` for a single star or a 'wide'
            binary's comet, `None` only for a 'close' binary's (which
            orbits the merged pair, not one individually-stored star row).

    Returns:
        int: The new `comets.id`.
    """
    comet_id = object_rows.insert_row(conn, "comets", COMET_FIELDS, comet, {
        "star_system_id": star_system_id, "star_id": star_id})

    for position, component in enumerate(comet.composition):
        conn.execute(
            "INSERT INTO comet_composition (comet_id, position, component) VALUES (?, ?, ?)",
            (comet_id, position, component),
        )

    return comet_id


def _remnant_galactic(remnant, context):
    """The movement vector and galactic orbit of a standalone remnant; an anchored one (`star_id` set) leaves them `NULL`
    (its motion already lives on its own `stars` row, so there is one source of truth for its position)."""
    if context["star_id"] is not None:
        return (None,) * 6
    return (*_free_vector(remnant, context["placement"]), remnant.galactic_orbital_period_gy,
            remnant.galactic_orbital_phase_deg, remnant.galactic_min_update_interval_years)


def _remnant_spin(remnant, context):
    """GEN.104: like the galactic fields, an anchored remnant's spin is on its stars row."""
    return (None,) * len(spin.SPIN_FIELDS) if context["star_id"] is not None else spin.spin_values(remnant)


_REMNANT_GALACTIC_FIELD = object_rows.computed(
    ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms", "galactic_orbital_period_gy", "galactic_orbital_phase_deg",
     "galactic_min_update_interval_years"), _remnant_galactic)
_REMNANT_SPIN_FIELD = object_rows.computed(spin.SPIN_FIELDS, _remnant_spin)

BLACK_HOLE_FIELDS = (
    object_rows.context("star_id"), object_rows.context("sector_id"), object_rows.attr("name"),
    object_rows.computed(("parent_type", "parent_id"), lambda remnant, context: anchored_parent(context["star_id"])),
    object_rows.attr("mass_class"), object_rows.attr("mass_solar"), object_rows.attr("event_horizon_radius_km"),
    object_rows.attr("spin"), _flag("has_accretion_disk"),
    object_rows.attr("temperature_k", "temperature"), object_rows.attr("luminosity_w", "luminosity"),
    object_rows.attr("age_gy", "age"), _REMNANT_GALACTIC_FIELD, _PLACEMENT_FIELD, _REMNANT_SPIN_FIELD,
)
"""tuple: The `Field`s of a `black_holes` row."""

NEUTRON_STAR_FIELDS = (
    object_rows.context("star_id"), object_rows.context("sector_id"), object_rows.attr("name"),
    object_rows.computed(("parent_type", "parent_id"), lambda remnant, context: anchored_parent(context["star_id"])),
    object_rows.attr("mass_solar"), object_rows.attr("radius_km", "radius"), object_rows.attr("spin_period_ms"),
    object_rows.attr("magnetic_field_gauss"), object_rows.attr("pulsar_type"),
    object_rows.attr("surface_temperature_k"), object_rows.attr("luminosity_w", "luminosity"),
    object_rows.attr("age_gy", "age"), _REMNANT_GALACTIC_FIELD, _PLACEMENT_FIELD, _REMNANT_SPIN_FIELD,
)
"""tuple: The `Field`s of a `neutron_stars` row."""

_BLACK_HOLE_OWN_COLUMNS = ("mass_class", "mass_solar", "event_horizon_radius_km", "spin", "has_accretion_disk")
_NEUTRON_STAR_OWN_COLUMNS = ("mass_solar", "spin_period_ms", "magnetic_field_gauss", "pulsar_type",
                             "surface_temperature_k")
"""tuple: The columns of a remnant's satellite row that are not already on its `stars` row."""


def insert_black_hole(conn, black_hole: BlackHole, star_id=None, sector_id=None, placement=None,
                      register_name=True) -> int:
    """
    Inserts a `black_holes` row (see `schema.sql`'s "v16"/"v17"/"v21"
    header notes).

    Args:
        conn (Connection): An open, schema-initialized connection.
        black_hole (BlackHole): The black hole to persist.
        star_id (int, optional): The owning `stars.id`, when this black
            hole anchors a `StarSystem` (`phenomenonGen.py --anchor-system`
            -- see `insert_star_system`'s single-star branch, the only
            caller that passes this). `None` for one generated standalone.
            When set, the `galactic_orbital_*` columns are left `NULL` --
            an anchored remnant's motion already lives on its own `stars`
            row (written by `insert_star`), so this avoids two sources of
            truth for the same object's position. `sector_id`/`placement`
            below are only ever given for a standalone black hole -- an
            anchored one belongs to its own `StarSystem`, not directly to
            a sector.
        sector_id (int, optional): The sector this standalone black hole
            was generated as part of (`run_sector.generate_sector_phenomena`)
            or placed near (`phenomenonGen.py --sector-id`). `None` (the
            default) leaves it unplaced -- `placement` must also be `None`
            in that case.
        placement (dict, optional): `center_x_pc`/`center_y_pc`/
            `center_z_pc`/`galactic_radius_pc` -- either the exact
            galaxy-frame position this black hole's own Hill-sphere-aware
            in-sector placement converts to
            (`_galaxy_placement_from_sector_offset`, via `insert_sector`),
            or `compute_phenomenon_placement`'s independent random jitter
            (`phenomenonGen.py`'s own standalone `--sector-id` use, which
            has no specific in-sector position to convert). `None` to
            leave this black hole unplaced.
        register_name (bool): Reserve the name through
            `system_name_registry` (v40). Skipped for an anchored black
            hole (its system holds the name) and a supernova remnant's
            core (`insert_supernova_remnant` passes `False`). A placed
            one is named by its object ID instead (GEN.64), a core by its
            own core-type ID; an unplaced core keeps `"<remnant> Core"`.

    Returns:
        int: The new `black_holes.id`.
    """
    name_base = None
    if register_name and star_id is None:
        name_base, diminutive_index = _reserve_phenomenon_name(conn, black_hole, "black-hole", placement)
    elif star_id is None:
        _take_core_id(conn, black_hole, "black-hole-core", placement)
    row_id = object_rows.insert_row(conn, "black_holes", BLACK_HOLE_FIELDS, black_hole, {
        "star_id": star_id, "sector_id": sector_id, "placement": placement})
    if name_base is not None:
        confirm_object_name(conn, name_base, "black_holes", row_id, diminutive_index)
    return row_id


def insert_neutron_star(conn, neutron_star: NeutronStar, star_id=None, sector_id=None, placement=None,
                        register_name=True) -> int:
    """
    Inserts a `neutron_stars` row (see `schema.sql`'s "v16"/"v17"/"v21"
    header notes).

    Args:
        conn (Connection): An open, schema-initialized connection.
        neutron_star (NeutronStar): The neutron star to persist.
        star_id (int, optional): The owning `stars.id`, when this neutron
            star anchors a `StarSystem` (`phenomenonGen.py --anchor-system`
            -- see `insert_star_system`'s single-star branch, the only
            caller that passes this). `None` for one generated standalone.
            When set, the `galactic_orbital_*` columns are left `NULL` --
            see `insert_black_hole`'s identical reasoning; `sector_id`/
            `placement` below likewise only ever apply to a standalone
            neutron star.
        sector_id (int, optional): As in `insert_black_hole`.
        placement (dict, optional): As in `insert_black_hole`.
        register_name (bool): As in `insert_black_hole`.

    Returns:
        int: The new `neutron_stars.id`.
    """
    name_base = None
    if register_name and star_id is None:
        name_base, diminutive_index = _reserve_phenomenon_name(conn, neutron_star, "neutron-star", placement)
    elif star_id is None:
        _take_core_id(conn, neutron_star, "neutron-star-core", placement)
    row_id = object_rows.insert_row(conn, "neutron_stars", NEUTRON_STAR_FIELDS, neutron_star, {
        "star_id": star_id, "sector_id": sector_id, "placement": placement})
    if name_base is not None:
        confirm_object_name(conn, name_base, "neutron_stars", row_id, diminutive_index)
    return row_id


def insert_nebula(conn, nebula: Nebula, sector_id=None, placement=None) -> int:
    """
    Inserts a `nebulae` row (see `schema.sql`'s "v16"/"v18" header notes).

    Args:
        conn (Connection): An open, schema-initialized connection.
        nebula (Nebula): The nebula to persist.
        sector_id (int, optional): The nearest already-generated sector to
            `placement`'s own center -- a convenience "home" link, not this
            phenomenon's real geometry (see `schema.sql`'s "v18" note).
            `None` (the default) for a phenomenon never placed in the
            galaxy at all -- `placement` must also be `None` in that case.
        placement (dict, optional): `compute_phenomenon_placement`'s
            return shape (`center_x_pc`/`center_y_pc`/`center_z_pc`/
            `galactic_radius_pc`), or `None` to leave this nebula
            unplaced (the v16/v17 default behavior).

    Returns:
        int: The new `nebulae.id`.
    """
    name_base, diminutive_index = _reserve_phenomenon_name(conn, nebula, "nebula", placement)
    placement = placement or {}
    cur = conn.execute(
        """
        INSERT INTO nebulae (
            sector_id, name, nebula_class, nebula_type, radius_ly, composition, formation_cause,
            dominant_species, density_cm3, temperature_k, extinction_av, age_myr, q_total, hosts,
            velocity_x_kms, velocity_y_kms, velocity_z_kms, mass_kg, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years,
            center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, nebula.name, nebula.nebula_class, nebula.nebula_type, nebula.radius_ly,
            nebula.composition, nebula.formation_cause,
            nebula.dominant_species, nebula.density_cm3, nebula.temperature_k, nebula.extinction_av,
            nebula.age_myr, nebula.q_total, nebula.hosts,
            *_free_vector(nebula, placement), bulk_mass.gas_cloud_mass_kg(nebula.density_cm3, nebula.radius_ly), nebula.galactic_orbital_period_gy,
            nebula.galactic_orbital_phase_deg, nebula.galactic_min_update_interval_years,
            placement.get("center_x_pc"), placement.get("center_y_pc"),
            placement.get("center_z_pc"), placement.get("galactic_radius_pc"),
        ),
    )
    confirm_object_name(conn, name_base, "nebulae", cur.lastrowid, diminutive_index)
    _store_nebula_shape(conn, cur.lastrowid, nebula)
    _refresh_containment_around(conn, placement, nebula.radius_ly)
    return cur.lastrowid


def _store_nebula_shape(conn, nebula_id, nebula):
    """Writes `nebula`'s shape (GEN.75) to its row's `shape_*` columns and to
    `nebula_shape_balls`."""
    from planetgen.galaxy import nebula_shape

    scalars, balls = nebula_shape.to_columns(nebula.get_shape())
    sets = ", ".join(f"{column} = ?" for column in scalars)
    conn.execute(f"UPDATE nebulae SET {sets} WHERE id = ?", (*scalars.values(), nebula_id))
    storage_checksum.mark_rows(conn, "nebulae", [nebula_id])
    conn.execute("DELETE FROM nebula_shape_balls WHERE nebula_id = ?", (nebula_id,))
    for ball in balls:
        conn.execute(
            "INSERT INTO nebula_shape_balls (nebula_id, ball_index, center_x, center_y, center_z, radius)"
            " VALUES (?, ?, ?, ?, ?, ?)", (nebula_id, *ball))


def _placement_values(placement):
    """
    `placement`'s four `center_x_pc`/`center_y_pc`/`center_z_pc`/
    `galactic_radius_pc` values in column order, all `None` for a `None`
    placement -- the tail of every placeable phenomenon's `INSERT`.
    """
    placement = placement or {}
    return (
        placement.get("center_x_pc"), placement.get("center_y_pc"),
        placement.get("center_z_pc"), placement.get("galactic_radius_pc"),
    )

def _default_place_pc():
    """tuple: Where an object with no place in the galaxy is taken to be for its movement vector: the radius that
    holds half of the galaxy's disc stars on the +x axis (`galaxy.scales`), the place
    `galactic_orbit.generate_galactic_orbit_fields` works its speed out for."""
    return (scales.default_galactocentric_pc(), 0.0, 0.0)


def _speed_kms(row):
    """km/s: the length of the movement vector a row stores in `velocity_x/y/z_kms` -- the galactic speed every
    reader used to take from the scalar `galactic_orbital_speed_kms` (GEN.336). `None` when the row has no vector."""
    parts = (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"])
    if any(part is None for part in parts):
        return None
    return math.sqrt(sum(part * part for part in parts))


def _free_vector(phenomenon, placement):
    """`(vx, vy, vz)`, km/s, galactic axes: the movement a free object is stored with (GEN.336) -- its
    `galactic_orbital_speed_kms` along the tangent of the galaxy's rotation at its place (`galactic_velocity_ms`), or at
    `_default_place_pc()` when it has none (a phenomenon never put in the galaxy). A remnant the galaxy ejected (GEN.238)
    carries the velocity it was given."""
    ejected = getattr(phenomenon, "ejecta_velocity_kms", None)
    if ejected is not None:
        return tuple(ejected)
    if not placement or placement.get("center_x_pc") is None:
        point = _default_place_pc()
    else:
        point = (placement["center_x_pc"], placement["center_y_pc"], placement["center_z_pc"])
    return tuple(c / 1000.0 for c in galactic_velocity_ms(point, phenomenon.galactic_orbital_speed_kms))


def insert_supernova_remnant(conn, remnant: SupernovaRemnant, sector_id=None, placement=None) -> int:
    """
    Inserts a `supernova_remnants` row (see `schema.sql`'s "v16"/"v28"
    header notes), plus (for a core-collapse progenitor whose collapsed
    core is still detectable) its embedded `black_holes`/`neutron_stars`
    row, which shares the remnant's own `sector_id` and sits at its
    `placement` moved by the core's birth kick (`compact_offset_ly`).

    Args:
        conn (Connection): An open, schema-initialized connection.
        remnant (SupernovaRemnant): The remnant to persist.
        sector_id (int, optional): The sector this remnant was generated
            as part of (`run_sector.generate_sector_phenomena`) or placed
            near (`phenomenonGen.py --sector-id`). `None` for one never
            generated as part of any sector.
        placement (dict, optional): As in `insert_black_hole` -- the
            galaxy-frame center, or `None` to leave this remnant unplaced.

    Returns:
        int: The new `supernova_remnants.id`.
    """
    old_name = remnant.name
    name_base, diminutive_index = _reserve_phenomenon_name(conn, remnant, "supernova-remnant", placement)
    core = remnant.compact_remnant
    if core is not None and core.name == f"{old_name} Core":
        core.name = f"{remnant.name} Core"
    compact_remnant_kind = None
    black_hole_id = None
    neutron_star_id = None
    core_placement = placement
    offset_ly = getattr(remnant, "compact_offset_ly", None)
    if placement is not None and offset_ly is not None:
        # The core has drifted off-center by its birth kick.
        x, y, z = (placement[f"center_{axis}_pc"] + ly_to_pc(offset_ly[i]) for i, axis in enumerate("xyz"))
        core_placement = {"center_x_pc": x, "center_y_pc": y, "center_z_pc": z,
                          "galactic_radius_pc": math.sqrt(x * x + y * y + z * z)}
    if isinstance(remnant.compact_remnant, BlackHole):
        compact_remnant_kind = "black_hole"
        black_hole_id = insert_black_hole(
            conn, remnant.compact_remnant, sector_id=sector_id, placement=core_placement, register_name=False,
        )
    elif isinstance(remnant.compact_remnant, NeutronStar):
        compact_remnant_kind = "neutron_star"
        neutron_star_id = insert_neutron_star(
            conn, remnant.compact_remnant, sector_id=sector_id, placement=core_placement, register_name=False,
        )

    cur = conn.execute(
        """
        INSERT INTO supernova_remnants (
            sector_id, name, remnant_class, morphology, age_years, radius_ly, progenitor_type,
            compact_remnant_kind, compact_remnant_black_hole_id, compact_remnant_neutron_star_id,
            dominant_species, density_cm3, temperature_k, extinction_av,
            velocity_x_kms, velocity_y_kms, velocity_z_kms, mass_kg, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years,
            center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, remnant.name, remnant.remnant_class, remnant.morphology, remnant.age_years,
            remnant.radius_ly, remnant.progenitor_type, compact_remnant_kind, black_hole_id, neutron_star_id,
            remnant.dominant_species, remnant.density_cm3, remnant.temperature_k, remnant.extinction_av,
            *_free_vector(remnant, placement), bulk_mass.gas_cloud_mass_kg(remnant.density_cm3, remnant.radius_ly), remnant.galactic_orbital_period_gy,
            remnant.galactic_orbital_phase_deg, remnant.galactic_min_update_interval_years,
            *_placement_values(placement),
        ),
    )
    confirm_object_name(conn, name_base, "supernova_remnants", cur.lastrowid, diminutive_index)
    _refresh_containment_around(conn, placement, remnant.radius_ly)
    return cur.lastrowid


def _rogue_surface_values(planet):
    """`ROGUE_SURFACE_FIELDS`' values from `planet`, as stored."""
    return tuple(
        int(value) if field == "has_liquid_water" and value is not None else value
        for field, value in ((field, getattr(planet, field, None)) for field in ROGUE_SURFACE_FIELDS)
    )


ROGUE_PLANET_FIELDS = (
    object_rows.context("sector_id"), object_rows.attr("name"), object_rows.attr("planet_type"),
    object_rows.optional("planet_class"), object_rows.attr("mass_bin"), object_rows.attr("mass_kg"),
    object_rows.attr("radius_km"), object_rows.attr("composition"),
    object_rows.computed("has_internal_heat", lambda planet, context: int(planet.has_internal_heat)),
    object_rows.computed("has_moons", lambda planet, context: int(planet.has_moons)),
    object_rows.computed(tuple(ROGUE_SURFACE_FIELDS), lambda planet, context: _rogue_surface_values(planet)),
    object_rows.each(rogue_habitability.ROGUE_FIELDS),
    object_rows.each(rings.BODY_FIELDS),  # GEN.836
    _FREE_VECTOR_FIELD, _GALACTIC_ORBIT_FIELD, _PLACEMENT_FIELD,
    object_rows.optional("rotation_period_hours"), _SPIN_FIELD,
)
"""tuple: The `Field`s of a `rogue_planets` row."""


def insert_rogue_planet(conn, planet: RoguePlanet, sector_id=None, placement=None) -> int:
    """
    Inserts a `rogue_planets` row (see `schema.sql`'s "v16"/"v28" header
    notes).

    Args:
        conn (Connection): An open, schema-initialized connection.
        planet (RoguePlanet): The rogue planet to persist.
        sector_id (int, optional): As in `insert_supernova_remnant`.
        placement (dict, optional): As in `insert_supernova_remnant`.

    Returns:
        int: The new `rogue_planets.id`.
    """
    name_base, diminutive_index = _reserve_phenomenon_name(conn, planet, "rogue-planet", placement)
    row_id = object_rows.insert_row(conn, "rogue_planets", ROGUE_PLANET_FIELDS, planet, {
        "sector_id": sector_id, "placement": placement})
    confirm_object_name(conn, name_base, "rogue_planets", row_id, diminutive_index)
    return row_id


def sector_code(conn, sector_id):
    """
    How a designation names a sector (v40): its grid designation
    (`spherical.provisional_sector_designation`) when it sits in the
    spherical grid, else its name. `None` for no sector.
    """
    if sector_id is None:
        return None
    row = conn.execute(
        "SELECT name, shell_index, band_index, wedge_index FROM sectors WHERE id = ?", (sector_id,),
    ).fetchone()
    if row is None:
        return None
    if None not in (row["shell_index"], row["band_index"], row["wedge_index"]):
        try:
            return provisional_sector_designation(row["shell_index"], row["band_index"], row["wedge_index"])
        except ValueError:
            pass
    return row["name"]


def _next_in_sector(conn, table, sector_id):
    """How many `table` rows `sector_id` (or no sector, for `None`)
    already has, plus one -- the next designation number."""
    if sector_id is None:
        row = conn.execute(f"SELECT COUNT(*) AS n FROM {table} WHERE sector_id IS NULL").fetchone()
    else:
        row = conn.execute(f"SELECT COUNT(*) AS n FROM {table} WHERE sector_id = ?", (sector_id,)).fetchone()
    return row["n"] + 1


INTERSTELLAR_COMET_FIELDS = (
    object_rows.context("sector_id"), object_rows.attr("name"), object_rows.attr("nucleus_diameter_km"),
    object_rows.attr("velocity_kms"), object_rows.computed("is_active", lambda comet, context: int(comet.is_active)),
    object_rows.computed("composition_summary", lambda comet, context: comet.get_composition_summary()),
    _FREE_VECTOR_FIELD,
    object_rows.attr("mass_kg"),
    _GALACTIC_ORBIT_FIELD, _PLACEMENT_FIELD, object_rows.optional("rotation_period_hours"), _SPIN_FIELD,
)
"""tuple: The `Field`s of an `interstellar_comets` row."""


def insert_interstellar_comet(conn, comet: InterstellarComet, sector_id=None, placement=None) -> int:
    """
    Inserts an `interstellar_comets` row (plus its
    `interstellar_comet_composition` child rows; see `schema.sql`'s
    "v16"/"v28" header notes).

    Args:
        conn (Connection): An open, schema-initialized connection.
        comet (InterstellarComet): The comet to persist.
        sector_id (int, optional): As in `insert_supernova_remnant`.
        placement (dict, optional): As in `insert_supernova_remnant`.

    Returns:
        int: The new `interstellar_comets.id`.
    """
    if not _named_by_object_id(conn, comet, "comet", placement):
        comet.name = interstellar_comet_designation(
            sector_code(conn, sector_id), _next_in_sector(conn, "interstellar_comets", sector_id),
        )
    comet_id = object_rows.insert_row(conn, "interstellar_comets", INTERSTELLAR_COMET_FIELDS, comet, {
        "sector_id": sector_id, "placement": placement})

    for position, component in enumerate(comet.composition):
        conn.execute(
            "INSERT INTO interstellar_comet_composition (comet_id, position, component) VALUES (?, ?, ?)",
            (comet_id, position, component),
        )

    return comet_id


ASTEROID_FIELD_FIELDS = (
    object_rows.context("sector_id"), object_rows.attr("name"), object_rows.attr("field_class"),
    object_rows.attr("composition_family"), object_rows.attr("density"), object_rows.attr("radius_ly"),
    object_rows.optional("chondrite_group"),
    object_rows.computed("composition_summary", lambda field, context: field.get_composition_summary()),
    _FREE_VECTOR_FIELD, object_rows.computed("mass_kg", lambda field, context: 0.0), _GALACTIC_ORBIT_FIELD,
    _PLACEMENT_FIELD,
)
"""tuple: The `Field`s of an `asteroid_fields` row."""


def insert_asteroid_field(conn, field: AsteroidField, sector_id=None, placement=None) -> int:
    """
    Inserts an `asteroid_fields` row (plus its `asteroid_field_composition`
    child rows; see `schema.sql`'s "v17"/"v18" header notes) -- the
    standalone counterpart to `insert_asteroid_belt`, following the
    identical belt-plus-child-rows shape.

    Args:
        conn (Connection): An open, schema-initialized connection.
        field (AsteroidField): The asteroid field to persist.
        sector_id (int, optional): The nearest already-generated sector to
            `placement`'s own center -- see `insert_nebula`'s identical
            parameter for the full explanation. `None` for a field never
            placed in the galaxy.
        placement (dict, optional): `compute_phenomenon_placement`'s
            return shape, or `None` to leave this field unplaced (the
            v16/v17 default behavior).

    Returns:
        int: The new `asteroid_fields.id`.
    """
    if not _named_by_object_id(conn, field, "asteroid-field", placement):
        field.name = asteroid_field_designation(
            field.field_class, sector_code(conn, sector_id), _next_in_sector(conn, "asteroid_fields", sector_id),
        )
    field_id = object_rows.insert_row(conn, "asteroid_fields", ASTEROID_FIELD_FIELDS, field, {
        "sector_id": sector_id, "placement": placement})

    for position, (component, concentration) in enumerate(field.composition):
        conn.execute(
            """
            INSERT INTO asteroid_field_composition (field_id, position, component, concentration)
            VALUES (?, ?, ?, ?)
            """,
            (field_id, position, component, concentration),
        )

    return field_id


GALACTIC_CENTER_PLACEMENT = {
    "center_x_pc": 0.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 0.0,
}
"""dict: The galactic origin, in `_placement_values`' shape -- where every
placed quasar sits."""


def insert_quasar(conn, quasar: Quasar, sector_id=None, placement=None) -> int:
    """
    Inserts a `quasars` row (see `schema.sql`'s "v31" header note).

    A quasar is a galaxy's nucleus, so any placement at all is snapped to
    the galactic center (`GALACTIC_CENTER_PLACEMENT`): `insert_sector`'s
    converted in-sector offset already lands there up to rounding, and
    `save_phenomenon`'s `--sector-id` jitter would otherwise scatter it.

    Args:
        conn (Connection): An open, schema-initialized connection.
        quasar (Quasar): The quasar to persist.
        sector_id (int, optional): The core sector it belongs to.
        placement (dict, optional): Any non-`None` value places it at the
            galactic center; `None` leaves it unplaced.

    Returns:
        int: The new `quasars.id`.
    """
    name_base, diminutive_index = _reserve_phenomenon_name(
        conn, quasar, "quasar", GALACTIC_CENTER_PLACEMENT if placement is not None else None)
    cur = conn.execute(
        """
        INSERT INTO quasars (
            sector_id, name, black_hole_mass_solar, event_horizon_radius_km, eddington_ratio,
            luminosity_w, accretion_rate_solar_per_year, broad_line_region_light_days,
            is_radio_loud, jet_length_ly, active_age_years,
            center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc,
            velocity_x_kms, velocity_y_kms, velocity_z_kms, epoch_unix
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, quasar.name, quasar.black_hole_mass_solar, quasar.event_horizon_radius_km,
            quasar.eddington_ratio, quasar.luminosity_w, quasar.accretion_rate_solar_per_year,
            quasar.broad_line_region_light_days, int(quasar.is_radio_loud), quasar.jet_length_ly,
            quasar.active_age_years,
            *_placement_values(GALACTIC_CENTER_PLACEMENT if placement is not None else None),
            # GEN.337: the nucleus is the origin of the galaxy frame and does not move (a stated zero vector), and its
            # point holds from the database's orbit epoch. It is never due: the orbit update skips it.
            0.0, 0.0, 0.0, get_orbit_epoch_unix(conn),
        ),
    )
    confirm_object_name(conn, name_base, "quasars", cur.lastrowid, diminutive_index)
    return cur.lastrowid


def _check_quasar_sector(conn, sector_id):
    """
    Refuses to place a quasar "at" a sector that doesn't host the galactic
    nucleus (the core, shell 0), or into a
    galaxy that already has one -- a galaxy has a single nucleus.

    Raises:
        ValueError: If either rule would be broken.
    """
    row = conn.execute("SELECT shell_index, band_index FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    if row is None or row["shell_index"] != 0 or row["band_index"] != 0:
        raise ValueError(
            f"a quasar can only be placed in the core sector (shell 0); "
            f"sector {sector_id} is not one"
        )
    existing = conn.execute("SELECT id FROM quasars WHERE center_x_pc IS NOT NULL LIMIT 1").fetchone()
    if existing is not None:
        raise ValueError(f"this galaxy already has a quasar at its center (quasars.id={existing['id']})")


_PHENOMENON_INSERTERS = {
    "black-hole": insert_black_hole,
    "neutron-star": insert_neutron_star,
    "nebula": insert_nebula,
    "asteroid-field": insert_asteroid_field,
    "supernova-remnant": insert_supernova_remnant,
    "rogue-planet": insert_rogue_planet,
    "comet": insert_interstellar_comet,
    "quasar": insert_quasar,
}
"""dict: `tuning.PHENOMENON_TYPE_CHOICES` value -> the
`insert_*` function for its table. Every one takes `sector_id=` and
`placement=` keywords (v28 gave the last three types placement columns),
so `insert_sector` and `save_phenomenon` share this one lookup."""


def save_phenomenon(phenomenon, system_config: SystemConfig, phenomenon_type: str, config=None,
                     sector_id=None) -> int:
    """
    Opens the database and persists a single generated exotic phenomenon
    (`phenomenonGen.py`) in one transaction -- the exotic-phenomenon
    counterpart to `save_system`.

    Args:
        phenomenon: The generated phenomenon -- a `BlackHole`/`NeutronStar`
            (standalone), a `StarSystem` (an anchored compact remnant, see
            `StarSystem.__init__`'s `compact_remnant` parameter and
            `phenomenonGen.py --anchor-system`), a `Nebula`, a
            `SupernovaRemnant`, a `RoguePlanet`, an `InterstellarComet`, or
            an `AsteroidField`.
        system_config (SystemConfig): The config it was generated from --
            only actually persisted (via `insert_star_system`) when
            `phenomenon` is a `StarSystem`; no other phenomenon type has a
            `SystemConfig` row of its own.
        phenomenon_type (str): One of
            `tuning.PHENOMENON_TYPE_CHOICES`, naming which
            table `phenomenon` belongs in (ignored when `phenomenon` is a
            `StarSystem`, which always goes through `insert_star_system`
            regardless of whether it's anchored by a black hole or a
            neutron star).
        config (MySQLConfig, optional): Connection parameters. Defaults to
            `DEFAULT_MYSQL_CONFIG`.
        sector_id (int, optional): The sector to associate this phenomenon
            with (`phenomenonGen.py --sector-id`, or `insert_sector`'s own
            call for a `run_sector.generate_sector_phenomena`-generated
            one). This also places the phenomenon in the galaxy near this
            already-placed sector via `compute_phenomenon_placement` (every
            type has placement columns since v28). `None` (the default)
            leaves it unplaced/unlinked, exactly like every earlier schema
            version.

    Returns:
        int: The new row's id, in whichever table `phenomenon`/
            `phenomenon_type` maps to.

    Raises:
        ValueError: If `phenomenon_type` isn't one of the recognized
                   choices (and `phenomenon` isn't a `StarSystem`), or if
                   `sector_id` is given but that sector has no galaxy
                   placement of its own (see `compute_phenomenon_placement`).
    """
    conn = get_connection(config)
    try:
        with conn:
            if isinstance(phenomenon, StarSystem):
                # An anchored compact remnant: insert_star_system's own
                # single-star branch already inserts the satellite
                # black_holes/neutron_stars row -- see that function.
                system_id = insert_star_system(conn, phenomenon, system_config)
                assign_uids(conn, system_ids=[system_id])
                return system_id
            inserter = _PHENOMENON_INSERTERS.get(phenomenon_type)
            if inserter is None:
                raise ValueError(f"Unknown phenomenon type: {phenomenon_type!r}")
            if phenomenon_type == "quasar" and sector_id is not None:
                _check_quasar_sector(conn, sector_id)
            placement = compute_phenomenon_placement(conn, sector_id) if sector_id is not None else None
            row_id = inserter(conn, phenomenon, sector_id=sector_id, placement=placement)
            assign_uids(conn, sector_id=sector_id,
                        phenomenon=None if sector_id is not None else (OBJECT_ID_TABLES[phenomenon_type], row_id))
            if sector_id is not None:
                refresh_containment(conn, [sector_id])
                refresh_nearest_systems(conn, [sector_id])
                # The galaxy's content stamp (`db.query.galaxy_content_state`)
                # only sees sectors and systems; without this, a phenomenon
                # added from the command line never reaches a cached node
                # or page (TEST.42).
                touch_sector(conn, sector_id)
            return row_id
    finally:
        conn.close()


_NULL_PROXY_ONLY_BINARY_FIELDS = (None,) * 14
"""Placeholder for the 14 `star_systems.binary_*` columns that only ever
describe a merged `BinaryStarProxy` (a 'close'/P-type pair) -- always NULL
for a 'wide'/S-type pair or a single star, since no merged effective star
exists to describe in either of those cases. See `schema.sql`'s "v15"
header note and `_proxy_only_binary_fields`."""

_NULL_MUTUAL_ORBIT_FIELDS = (None,) * 16
"""Placeholder for the 16 `star_systems.binary_mutual_*`/
`binary_primary_position_*`/`binary_secondary_position_*`/
`binary_secondary_mass_fraction` columns (excluding `binary_separation_km`,
handled separately in `insert_star_system` since it sits earlier in column
order, alongside the eccentricity/periapsis/apoapsis columns it's grouped
with) -- NULL for a single (non-binary) star. See
`_mutual_orbit_fields_from_proxy`/`_mutual_orbit_fields_from_wide_binary`."""


def _proxy_only_binary_fields(proxy: BinaryStarProxy):
    """
    Extracts the 14 `star_systems.binary_*` column values that only ever
    apply to a merged `BinaryStarProxy` -- `binary_type` (the pair's
    combined spectral-summary string) through
    `binary_galactic_min_update_interval_years` -- in the exact order
    `insert_star_system`'s `INSERT` lists them. NOT used for a 'wide'
    (S-type) pair, which has no merged effective star (see `schema.sql`'s
    "v15" header note); its two stars' own equivalent data lives on their
    own `stars` rows instead.

    Args:
        proxy (BinaryStarProxy): The system's combined-pair proxy.

    Returns:
        tuple: 15 values, ready to splice into the `INSERT` parameters.
    """
    return (
        proxy.type,
        proxy.temperature,
        proxy.radius,
        proxy.mass,
        proxy.luminosity,
        proxy.age,
        _lifespan_gy(proxy.lifespan),
        _au_to_km(proxy.habitable_zone[0]),
        _au_to_km(proxy.habitable_zone[1]),
        _au_to_km(proxy.system_perimeter),
        _au_to_km(proxy.heliosphere_radius),
        proxy.galactic_orbital_period_gy,
        proxy.galactic_orbital_phase_deg,
        proxy.galactic_min_update_interval_years,
    )


def _mutual_orbit_fields_from_proxy(proxy: BinaryStarProxy):
    """
    Extracts the 16 `star_systems.binary_mutual_*`/`binary_primary_position_*`/
    `binary_secondary_position_*`/`binary_secondary_mass_fraction` column
    values (period, speed, inclination, ascending node, phase,
    update-guard interval, x/y/z position, each star's own barycenter
    offset, and the mass fraction -- everything except
    `binary_separation_km`, handled separately) from a 'close' pair's
    `BinaryStarProxy`.

    Returns:
        tuple: 16 values, ready to splice into the `INSERT` parameters.
    """
    return (
        proxy.binary_mutual_orbital_period_years,
        proxy.binary_mutual_orbital_speed_kms,
        proxy.binary_mutual_orbital_inclination_deg,
        proxy.binary_mutual_orbital_ascending_node_deg,
        proxy.binary_mutual_orbital_phase_deg,
        proxy.binary_mutual_min_update_interval_years,
        _au_to_km(proxy.binary_mutual_position_x),
        _au_to_km(proxy.binary_mutual_position_y),
        _au_to_km(proxy.binary_mutual_position_z),
        _au_to_km(proxy.binary_primary_position_x),
        _au_to_km(proxy.binary_primary_position_y),
        _au_to_km(proxy.binary_primary_position_z),
        _au_to_km(proxy.binary_secondary_position_x),
        _au_to_km(proxy.binary_secondary_position_y),
        _au_to_km(proxy.binary_secondary_position_z),
        proxy.binary_secondary_mass_fraction,
    )


def _mutual_orbit_fields_from_wide_binary(pair):
    """
    The same 16 `star_systems.binary_mutual_*`/`binary_primary_position_*`/
    `binary_secondary_position_*`/`binary_secondary_mass_fraction` columns
    as `_mutual_orbit_fields_from_proxy`, from a 'wide' pair's
    `wide_binary.WideBinaryPair` instead -- these columns are reused
    unchanged across both binary configurations (see `schema.sql`'s "v15"
    header note): a wide pair's own (circular-approximation) mutual orbit
    fits the exact same shape a close pair's already occupies.

    Args:
        pair (WideBinaryPair): The system's wide-binary orbital pair.

    Returns:
        tuple: 16 values, ready to splice into the `INSERT` parameters.
    """
    return (
        pair.period_years,
        pair.speed_kms,
        pair.inclination_deg,
        pair.ascending_node_deg,
        pair.phase_deg,
        pair.min_update_interval_years,
        _au_to_km(pair.position_x_au),
        _au_to_km(pair.position_y_au),
        _au_to_km(pair.position_z_au),
        _au_to_km(pair.primary_position_x_au),
        _au_to_km(pair.primary_position_y_au),
        _au_to_km(pair.primary_position_z_au),
        _au_to_km(pair.secondary_position_x_au),
        _au_to_km(pair.secondary_position_y_au),
        _au_to_km(pair.secondary_position_z_au),
        pair.secondary_mass_fraction,
    )


def _format_location_string(sector_name, neighbors):
    """
    Builds the `star_systems.location` display string: the owning sector's
    name, followed by up to 3 nearest neighbors and their distances in
    light-years, nearest first -- see `schema.sql`'s "v3" header note.

    Used by the live-object write path (`_location_for_entry`, used from
    `insert_sector`).

    Args:
        sector_name (str): The owning `sectors.name`.
        neighbors (list): Up to 3 `(name, distance_ly)` tuples, nearest
                          first. Empty when the sector has no other systems.

    Returns:
        str: e.g. `"Voranthis Kelmoor -- nearest: Alpha Vesta (4.2 ly), ..."`,
             or just `sector_name` when `neighbors` is empty.
    """
    if not neighbors:
        return sector_name
    parts = [f"{name} ({distance_ly:.1f} ly)" for name, distance_ly in neighbors]
    return f"{sector_name} -- nearest: " + ", ".join(parts)


def _location_for_entry(sector: SpaceSector, entry: SectorSystemEntry, neighbors=None) -> str:
    """
    Computes `entry`'s `star_systems.location` string from the live
    `sector` it belongs to -- up to 3 nearest neighbors
    (`SpaceSector.nearest_neighbors`), nearest first, with their distances
    in light-years (`distance_between`; `SectorSystemEntry.position` is
    already light-years, so no `milliparsecs_to_ly` conversion is needed
    here -- that only applies to the `position_x/y/z_mpc` storage columns,
    not this in-memory computation).

    Args:
        sector (SpaceSector): The sector `entry` belongs to.
        entry (SectorSystemEntry): The system to compute a location for.
        neighbors (list, optional): `entry`'s 3 nearest neighbors when the caller already has them for every
            entry (`SpaceSector.all_nearest_neighbors`, PERF.81); worked out here otherwise.

    Returns:
        str: See `_format_location_string`.
    """
    if neighbors is None:
        neighbors = sector.nearest_neighbors(entry, count=3)
    neighbor_info = [
        (neighbor.star_system.name, distance_between(entry, neighbor))
        for neighbor in neighbors
    ]
    return _format_location_string(sector.name, neighbor_info)


def _designate_comets(star_system):
    """
    Gives every comet in `star_system` its designation
    (`cometData.comet_designation`, v40), numbered per host star: the
    system's name for a single star or close pair, each star's own name
    for a wide pair.
    """
    if getattr(star_system, "binary_type", None) == "wide":
        hosts = [(star_system.primary_star.name, star_system.comets),
                 (star_system.secondary_star.name, star_system.secondary_comets)]
    else:
        hosts = [(star_system.name, star_system.comets)]
    for host_name, comets in hosts:
        for index, comet in enumerate(comets, start=1):
            comet.name = comet_designation(host_name, index, comet)


def _galaxy_point_pc(conn, sector_id, offset_mpc):
    """A system's point in the galaxy frame, parsecs (`star_systems.galaxy_*_pc`), from its sector-local offset in
    milliparsecs; `(None, None, None)` while it has no position or its sector has no center."""
    center = _sector_center_pc(conn, sector_id)
    if center is None or offset_mpc is None:
        return (None, None, None)
    return local_to_galaxy_pc(center, tuple(c / 1000.0 for c in offset_mpc))


def _sector_center_pc(conn, sector_id):
    """A sector's center `(x, y, z)`, parsecs, or `None` when it has none (a standalone sector, or no sector at all).
    Inside a `batched()` scope each is read once (a read flushes the held-back INSERTs)."""
    if sector_id is None:
        return None
    if conn._batch is not None and sector_id in conn.sector_centers_pc:
        return conn.sector_centers_pc[sector_id]
    row = conn.execute("SELECT center_x_pc, center_y_pc, center_z_pc FROM sectors WHERE id = ?",
                       (sector_id,)).fetchone()
    center = None if row is None or row["center_x_pc"] is None else tuple(row[f"center_{axis}_pc"] for axis in "xyz")
    if conn._batch is not None:
        conn.sector_centers_pc[sector_id] = center
    return center


def _system_velocity_kms(conn, star_system, sector_id, position_pc):
    """
    The galactic velocity, km/s, `star_system` has at `position_pc`
    (`(x, y, z)` parsecs from the center of sector `sector_id`): the rotation
    curve's tangent at its place plus its runaway motion
    (`system_position.system_velocity_ms`). While it has no place it is taken at `_default_place_pc()`, where the
    generators draw its speed, so its vector keeps that speed.
    """
    center = None if sector_id is None or position_pc is None else _sector_center_pc(conn, sector_id)
    if center is None:
        place = _default_place_pc()
    else:
        place = tuple(c + offset for c, offset in zip(center, position_pc))
    return tuple(v / 1000.0 for v in system_velocity_ms(star_system, place))


def insert_star_system(conn, star_system: StarSystem, system_config: SystemConfig,
                        sector_id=None, position=None, location=None) -> int:
    """
    Inserts a full `StarSystem` -- the `star_systems` row, its `stars` row(s),
    and every planet/moon/asteroid belt/comet it contains -- into the database.

    No page text is stored (v29): Markdown and plain text are rendered on demand
    from these rows by `planetgen/db/render.py` -- see
    `schema.sql`'s "v29" header note.

    `star_system.name` is reserved via `reserve_system_name` (v24,
    `names/uniqueness.py`), so every other system and sector anywhere in
    the database stays distinct from this one, then every star, planet and
    moon is named from the final name (`StarSystem.assign_names`, v34).
    See that function's own docstring for the full mechanism.

    Args:
        conn (Connection): An open, schema-initialized connection.
        star_system (StarSystem): The generated system to persist.
        system_config (SystemConfig): The config it was generated from.
        sector_id (int, optional): The owning `sectors.id`, if this system
                                   belongs to a sector. `None` for a
                                   standalone system.
        position (tuple, optional): `(x, y, z)` in light-years, relative to
                                    the sector's center (as stored on
                                    `SectorSystemEntry.position`). `None`
                                    if the system isn't placed in a sector
                                    -- `position_x/y/z_mpc`, `quadrant`, and
                                    `location` are then left `NULL`
                                    regardless of the `location` argument.
        location (str, optional): The precomputed `star_systems.location`
                                  string (see `schema.sql`'s "v3" header
                                  note and `_location_for_entry`) -- callers
                                  with a full `SpaceSector` (`insert_sector`)
                                  compute this once per entry and pass it in,
                                  since deriving it needs sibling systems
                                  this function doesn't otherwise see.
                                  Ignored (forced `None`) when `position` is
                                  `None`.

    Returns:
        int: The new `star_systems.id`.
    """
    config_id = insert_system_config(conn, system_config)

    # binary_type (None | "close" | "wide") is the authoritative
    # discriminator (see system.StarSystem.__init__); is_binary is kept
    # for the schema's own pre-existing column and now means "this system
    # has two stars", true for either configuration. proxy_like gates
    # exactly the columns that only ever describe a merged BinaryStarProxy
    # -- see schema.sql's "v15" header note.
    binary_configuration = getattr(star_system, "binary_type", None)
    is_binary = binary_configuration is not None
    proxy_like = isinstance(star_system.star, BinaryStarProxy)

    if proxy_like:
        proxy_only_fields = _proxy_only_binary_fields(star_system.star)
        mutual_orbit_fields = _mutual_orbit_fields_from_proxy(star_system.star)
        separation_km = _au_to_km(star_system.star.binary_separation_au)
        # A close pair's mutual orbit is treated as circular (tidal
        # circularization is a legitimate simplification at its 0.05-0.25
        # AU separations) -- periapsis == apoapsis == separation, eccentricity 0.
        eccentricity = 0.0
        periapsis_km = apoapsis_km = separation_km
    elif binary_configuration == "wide":
        proxy_only_fields = _NULL_PROXY_ONLY_BINARY_FIELDS
        mutual_orbit_fields = _mutual_orbit_fields_from_wide_binary(star_system.wide_binary)
        separation_km = _au_to_km(star_system.wide_binary.separation_au)
        eccentricity = star_system.wide_binary.eccentricity
        periapsis_km = _au_to_km(star_system.wide_binary.periapsis_au)
        apoapsis_km = _au_to_km(star_system.wide_binary.apoapsis_au)
    else:
        proxy_only_fields = _NULL_PROXY_ONLY_BINARY_FIELDS
        mutual_orbit_fields = _NULL_MUTUAL_ORBIT_FIELDS
        separation_km = eccentricity = periapsis_km = apoapsis_km = None

    # v20: circumbinary (P-type) planets' combined pull on the whole pair
    # -- only meaningful for a 'close' pair (see schema.sql's "v20" header
    # note); NULL for 'wide'/single, same "meaningful only when applicable"
    # convention `black_holes`' galactic columns already use.
    if binary_configuration == "close":
        planetary_wobble_fields = (
            _au_to_km(star_system.binary_planetary_wobble_x),
            _au_to_km(star_system.binary_planetary_wobble_y),
            _au_to_km(star_system.binary_planetary_wobble_z),
        )
    else:
        planetary_wobble_fields = (None, None, None)

    if position is not None:
        position_x_mpc = ly_to_milliparsecs(position[0])
        position_y_mpc = ly_to_milliparsecs(position[1])
        position_z_mpc = ly_to_milliparsecs(position[2])
        quadrant, _magnitudes = classify_octant(position)
    else:
        position_x_mpc = position_y_mpc = position_z_mpc = None
        quadrant = None
        location = None
    velocity_kms = _system_velocity_kms(conn, star_system, sector_id,
                                        None if position is None else tuple(ly_to_pc(c) for c in position))
    galaxy_pc = _galaxy_point_pc(conn, sector_id, None if position is None else
                                 (position_x_mpc, position_y_mpc, position_z_mpc))

    # Name-uniqueness (v24, names/uniqueness.py) -- the stars, planets and
    # moons are renamed from the final name in place, so the caller's
    # object matches what's stored.
    check_system_name_length(star_system.name)
    name_base, diminutive_index = _take_name(conn, star_system)
    star_system.assign_names(star_system.name)
    _designate_comets(star_system)

    cur = conn.execute(
        """
        INSERT INTO star_systems (
            sector_id, system_config_id, name,
            position_x_mpc, position_y_mpc, position_z_mpc, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc, quadrant, location,
            is_binary, binary_configuration,
            binary_separation_km, binary_eccentricity, binary_periapsis_km, binary_apoapsis_km,
            binary_type, binary_temperature_k, binary_radius_km,
            binary_effective_mass_kg, binary_effective_luminosity_w, binary_age_gy, binary_lifespan_gy,
            binary_habitable_zone_inner_km, binary_habitable_zone_outer_km,
            binary_system_perimeter_km, binary_heliosphere_radius_km,
            binary_galactic_orbital_period_gy,
            binary_galactic_orbital_phase_deg, binary_galactic_min_update_interval_years,
            binary_mutual_orbital_period_years, binary_mutual_orbital_speed_kms,
            binary_mutual_orbital_inclination_deg, binary_mutual_orbital_ascending_node_deg,
            binary_mutual_orbital_phase_deg, binary_mutual_min_update_interval_years,
            binary_mutual_position_x_km, binary_mutual_position_y_km, binary_mutual_position_z_km,
            binary_primary_position_x_km, binary_primary_position_y_km, binary_primary_position_z_km,
            binary_secondary_position_x_km, binary_secondary_position_y_km, binary_secondary_position_z_km,
            binary_secondary_mass_fraction,
            binary_planetary_wobble_x_km, binary_planetary_wobble_y_km, binary_planetary_wobble_z_km,
            system_flavor_text, runaway_class, runaway_speed_kms,
            velocity_x_kms, velocity_y_kms, velocity_z_kms, schema_version
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, config_id, star_system.name,
            position_x_mpc, position_y_mpc, position_z_mpc, *galaxy_pc, quadrant, location,
            int(is_binary), binary_configuration,
            separation_km, eccentricity, periapsis_km, apoapsis_km,
            *proxy_only_fields,
            *mutual_orbit_fields,
            *planetary_wobble_fields,
            star_system.system_flavor_text,
            getattr(star_system, "runaway_class", None), getattr(star_system, "runaway_speed_kms", None),
            *velocity_kms,
            SCHEMA_VERSION,
        ),
    )
    star_system_id = cur.lastrowid
    confirm_system_name(conn, name_base, star_system_id, diminutive_index)
    system_epoch = getattr(star_system, "epoch_unix", None)
    if system_epoch is not None:  # a scattered hypervelocity star holds at the time it was drawn (GEN.137)
        _update_by_id(conn, "star_systems", ("epoch_unix",), [(star_system_id, system_epoch)])

    # GEN.311: the system's tree of orbits. A generated binary is one node, made from its pair; a tree the system was
    # given (`orbit_nodes`) is stored as it is. The phases hold at the database's orbit epoch (now, before any update).
    nodes = tuple(star_system.orbit_nodes)
    generated = bool(nodes) and nodes[0].epoch_unix == stored_trees.GENERATED_EPOCH_UNIX    # GEN.631: a drawn tree
    if (not nodes and star_system.binary_type is not None) or generated:
        orbit_epoch = conn.orbit_epoch_unix if conn._batch is not None else None
        if orbit_epoch is None:
            orbit_epoch = get_orbit_epoch_unix(conn)
            if orbit_epoch is None:
                orbit_epoch = int(time.time())
            if conn._batch is not None:
                conn.orbit_epoch_unix = orbit_epoch
        if generated:
            nodes = tuple(dataclasses.replace(node, epoch_unix=float(orbit_epoch)) for node in nodes)
        else:
            nodes = orbit_trees.binary_nodes(star_system, orbit_epoch)
    orbit_trees.write_nodes(conn, star_system_id, nodes, replace=False)
    orbit_trees.write_dark_members(conn, star_system_id, star_system.dark_members, replace=False)    # GEN.1013
    star_system.orbit_nodes = nodes
    places = orbit_trees.star_places(star_system)
    # GEN.329: how a runaway star came to move so fast, in the flags of its (first) star row
    origin_flags = multistar.RUNAWAY_ORIGIN_FLAGS.get(getattr(star_system, "runaway_origin", None), 0)

    if proxy_like:
        insert_star(conn, star_system.primary_star, star_system_id, "primary", *places[0], flags=origin_flags)
        insert_star(conn, star_system.secondary_star, star_system_id, "secondary", *places[1])
        primary_star_id = secondary_star_id = None  # planets orbit the proxy, not a stored star row -- see planets.star_id
    elif binary_configuration == "wide":
        primary_star_id = insert_star(conn, star_system.primary_star, star_system_id, "primary", *places[0],
                                      flags=origin_flags)
        secondary_star_id = insert_star(conn, star_system.secondary_star, star_system_id, "secondary", *places[1])
    else:
        primary_star_id = insert_star(conn, star_system.star, star_system_id, "single", *places[0], flags=origin_flags)
        secondary_star_id = None
        # v16: a StarSystem anchored by a compact remnant (generation/phenomena/compact_remnant.py,
        # phenomenonGen.py --anchor-system) -- insert_star above already
        # wrote the stars row (BlackHole/NeutronStar's Star-compatible
        # attribute surface, with yerkes_class set to the literal 'BH'/'NS'
        # marker), so only the satellite black_holes/neutron_stars row
        # (this remnant's own extra fields) remains.
        if isinstance(star_system.star, BlackHole):
            insert_black_hole(conn, star_system.star, star_id=primary_star_id)
        elif isinstance(star_system.star, NeutronStar):
            insert_neutron_star(conn, star_system.star, star_id=primary_star_id)

    member_ids = [insert_star(conn, member, star_system_id, "member", *place)
                  for member, place in zip(star_system.member_stars, places[2:])]

    for orbital_index, obj in enumerate(star_system.planets):
        if obj.body_type == "a":
            insert_asteroid_belt(conn, obj, star_system_id, orbital_index, star_id=primary_star_id)
        else:
            insert_planet(conn, obj, star_system_id, primary_star_id, orbital_index)

    # Only ever non-empty for a 'wide' binary (see StarSystem.__init__) --
    # its own, independent index space, disambiguated from star_system.planets'
    # by secondary_star_id (see load_star_system's own per-star_id grouping).
    for orbital_index, obj in enumerate(star_system.secondary_planets):
        if obj.body_type == "a":
            insert_asteroid_belt(conn, obj, star_system_id, orbital_index, star_id=secondary_star_id)
        else:
            insert_planet(conn, obj, star_system_id, secondary_star_id, orbital_index)

    # A member star's planets (GEN.332): each member has its own orbital-index space, told apart by its star_id.
    for member_id, bodies in zip(member_ids, star_system.member_planets):
        for orbital_index, obj in enumerate(bodies):
            if obj.body_type == "a":
                insert_asteroid_belt(conn, obj, star_system_id, orbital_index, star_id=member_id)
            else:
                insert_planet(conn, obj, star_system_id, member_id, orbital_index)

    # Comets have their own list, own table, and no orbital_index -- see
    # insert_comet's own docstring for why they don't share the
    # planets/asteroid_belts orbital-slot loop above.
    for comet in star_system.comets:
        insert_comet(conn, comet, star_system_id, star_id=primary_star_id)
    for comet in star_system.secondary_comets:
        insert_comet(conn, comet, star_system_id, star_id=secondary_star_id)

    return star_system_id


# --- Containment (v39): what sits inside a nebula or supernova remnant ---

CONTAINABLE_TABLES = kind_tables.tables("containable", order="containable")
"""tuple: Every table with `inside_nebula_id`/`inside_remnant_id` (schema
v39). Asteroid fields can sit inside a cloud but never contain anything;
nebulae nest (a small one inside a larger nebula or remnant)."""

CONTAINER_TABLES = (("nebulae", "inside_nebula_id"), ("supernova_remnants", "inside_remnant_id"))
"""tuple: `(container table, the column that points at it)`."""


def surrounding_cloud(conn, row):
    """
    The nebula or supernova remnant a row sits inside (schema v39), from
    its `inside_nebula_id`/`inside_remnant_id`.

    Args:
        conn (Connection): An open connection.
        row (Mapping): Any row carrying those two columns.

    Returns:
        dict or None: The cloud dict of `generation.cloud` -- `type`
            (`"nebula"` or `"supernova_remnant"`), `id`, `name`, `class`,
            `descriptor` (a nebula's type, a remnant's morphology),
            `density_cm3`, `temperature_k`, `extinction_av`, `radius_ly`,
            `age_myr` (a remnant's age in Myr), `center` (galaxy frame,
            parsecs, `None` when the cloud is not placed) and `host_stars`
            -- or `None` in open space (GEN.357 added the last five).
    """
    if row["inside_nebula_id"] is not None:
        found = conn.execute("SELECT id, name, nebula_class AS class, nebula_type AS descriptor, density_cm3,"
                             " temperature_k, extinction_av, radius_ly, age_myr, hosts, center_x_pc, center_y_pc,"
                             " center_z_pc FROM nebulae WHERE id = ?", (row["inside_nebula_id"],)).fetchone()
        kind, age_myr = "nebula", None if found is None else found["age_myr"]
        hosts = None if found is None else found["hosts"]
    elif row["inside_remnant_id"] is not None:
        found = conn.execute("SELECT id, name, remnant_class AS class, morphology AS descriptor, density_cm3,"
                             " temperature_k, extinction_av, radius_ly, age_years, center_x_pc, center_y_pc,"
                             " center_z_pc FROM supernova_remnants WHERE id = ?",
                             (row["inside_remnant_id"],)).fetchone()
        kind, hosts = "supernova_remnant", None
        age_myr = None if found is None or found["age_years"] is None else found["age_years"] / 1.0e6
    else:
        return None
    if found is None:
        return None
    center = None if found["center_x_pc"] is None else (found["center_x_pc"], found["center_y_pc"], found["center_z_pc"])
    return cloudContext.cloud_dict(kind, found["id"], found["name"], found["class"], found["descriptor"],
                                   found["density_cm3"], found["temperature_k"], found["extinction_av"],
                                   found["radius_ly"], age_myr, center, hosts)


def _sector_half_diagonal_pc(edge_mpc):
    """Half a sector cube's space diagonal, parsecs -- how far any point in
    the sector can be from its center."""
    return (edge_mpc / 1000.0) * math.sqrt(3) / 2


def _placed_containers(conn, low, high, reach_pc=0.0):
    """
    Every placed nebula and supernova remnant whose sphere can reach the box
    `low`-`high` (galaxy-frame parsecs, each an `(x, y, z)`), padded by
    `reach_pc`, as dicts with `column`, `id`, `center` and `radius_pc`.
    """
    containers = []
    for table, column in CONTAINER_TABLES:
        max_row = conn.execute(f"SELECT MAX(radius_ly) AS r FROM {table} WHERE center_x_pc IS NOT NULL").fetchone()
        if max_row["r"] is None:
            continue
        pad = ly_to_pc(max_row["r"]) + reach_pc
        rows = conn.execute(
            f"SELECT id, center_x_pc, center_y_pc, center_z_pc, radius_ly FROM {table}"
            " WHERE center_x_pc BETWEEN ? AND ? AND center_y_pc BETWEEN ? AND ? AND center_z_pc BETWEEN ? AND ?",
            (low[0] - pad, high[0] + pad, low[1] - pad, high[1] + pad, low[2] - pad, high[2] + pad),
        ).fetchall()
        for row in rows:
            container = {
                "column": column, "id": row["id"],
                "center": (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]),
                "radius_pc": ly_to_pc(row["radius_ly"]),
            }
            if table == "nebulae":
                container["shape"] = _lazy_nebula_shape(conn, row["id"])
            containers.append(container)
    return containers


def _lazy_nebula_shape(conn, nebula_id):
    """A function returning nebula `nebula_id`'s `NebulaShape` (GEN.75),
    loaded on the first call: most containers never see a point inside
    their sphere."""
    loaded = []

    def shape():
        if not loaded:
            from planetgen.db import query

            loaded.append(query.nebula_shape(conn, nebula_id)[1])
        return loaded[0]

    return shape


def refresh_containment(conn, sector_ids, on_progress=None):
    """
    Sets `inside_nebula_id`/`inside_remnant_id` (schema v39) on every star
    system and phenomenon filed under `sector_ids`: a 3D distance test
    against every nebula and supernova remnant that reaches those sectors,
    keeping the innermost (smallest) container. Rows whose container
    didn't change aren't written. Called for a newly generated sector, for
    every sector a newly placed nebula or remnant reaches, and by the v39
    migration.

    Args:
        conn (Connection): Part of the caller's transaction.
        sector_ids (iterable): `sectors.id` values; unplaced ones are skipped.
        on_progress (callable, optional): `on_progress("Sectors", done, total)`
            after each batch of sectors.

    Returns:
        dict: How many objects entered and left a container (GEN.107):
            `CONTAINMENT_COUNTS` -> count. Moving from one nebula straight
            into another counts as leaving one and entering the other.
    """
    counts = dict.fromkeys(CONTAINMENT_COUNTS, 0)
    sector_ids = sorted(set(sector_ids))
    with steps.hooked(on_progress, "Containment: sectors", "containment", len(sector_ids)) as report:
        for start in range(0, len(sector_ids), 500):
            _refresh_containment_batch(conn, sector_ids[start:start + 500], counts)
            _report(report, "Sectors", min(start + 500, len(sector_ids)), len(sector_ids))
    return counts


CONTAINMENT_COUNTS = ("entered_nebula", "left_nebula", "entered_remnant", "left_remnant")
"""tuple: The keys of `refresh_containment`'s counts."""


def _count_containment_change(counts, old, new):
    """Adds one row's change of `(nebula, remnant)` container to `counts`."""
    for (old_id, new_id), kind in zip(zip(old, new), ("nebula", "remnant")):
        if old_id == new_id:
            continue
        if old_id is not None:
            counts[f"left_{kind}"] += 1
        if new_id is not None:
            counts[f"entered_{kind}"] += 1


def _refresh_containment_batch(conn, sector_ids, counts):
    marks = ", ".join("?" * len(sector_ids))
    sectors = conn.execute(
        f"SELECT id, center_x_pc, center_y_pc, center_z_pc, edge_mpc FROM sectors"
        f" WHERE id IN ({marks}) AND center_x_pc IS NOT NULL",
        tuple(sector_ids),
    ).fetchall()
    if not sectors:
        return
    centers = {row["id"]: (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]) for row in sectors}
    reach = max(_sector_half_diagonal_pc(row["edge_mpc"]) for row in sectors)
    low = tuple(min(c[i] for c in centers.values()) for i in range(3))
    high = tuple(max(c[i] for c in centers.values()) for i in range(3))
    containers = _placed_containers(conn, low, high, reach)
    marks = ", ".join("?" * len(centers))
    placed_ids = tuple(centers)

    changes = {}

    def apply(table, row, point, own_radius_pc=0.0, own=None):
        best = containment.innermost_container(point, containers, own_radius_pc, own) if point is not None else None
        new = (best["id"] if best and best["column"] == "inside_nebula_id" else None,
               best["id"] if best and best["column"] == "inside_remnant_id" else None)
        old = (row["inside_nebula_id"], row["inside_remnant_id"])
        if new != old:
            changes.setdefault(table, []).append((row["id"], *new))
            _count_containment_change(counts, old, new)

    for row in conn.execute(
        f"SELECT id, sector_id, position_x_mpc, position_y_mpc, position_z_mpc, inside_nebula_id, inside_remnant_id"
        f" FROM star_systems WHERE sector_id IN ({marks})",
        placed_ids,
    ).fetchall():
        point = None
        if row["position_x_mpc"] is not None:
            offset = (row["position_x_mpc"] / 1000.0, row["position_y_mpc"] / 1000.0, row["position_z_mpc"] / 1000.0)
            point = local_to_galaxy_pc(centers[row["sector_id"]], offset)
        apply("star_systems", row, point)

    for table in CONTAINABLE_TABLES[1:]:
        radius = ", radius_ly" if table == "nebulae" else ""
        for row in conn.execute(
            f"SELECT id, center_x_pc, center_y_pc, center_z_pc, inside_nebula_id, inside_remnant_id{radius}"
            f" FROM {table} WHERE sector_id IN ({marks})",
            placed_ids,
        ).fetchall():
            point = None if row["center_x_pc"] is None else (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
            if table == "nebulae":
                apply(table, row, point, ly_to_pc(row["radius_ly"]), ("inside_nebula_id", row["id"]))
            else:
                apply(table, row, point)

    for table, rows in changes.items():
        _update_by_id(conn, table, ("inside_nebula_id", "inside_remnant_id"), rows)
    _refresh_cloud_doses(conn, changes.get("star_systems", ()))


def refresh_system_dose(conn, system_id):
    """Re-derives one system's cloud-dependent doses (`_refresh_cloud_doses`)
    from the container it is stored inside: after an edit rewrote its
    bodies with the open-space multiplier."""
    row = conn.execute("SELECT inside_nebula_id, inside_remnant_id FROM star_systems WHERE id = ?",
                       (system_id,)).fetchone()
    if row is not None:
        _refresh_cloud_doses(conn, [(system_id, row["inside_nebula_id"], row["inside_remnant_id"])])


def compact_host_of_star(conn, star_id):
    """The `habitability_world.COMPACT_HOSTS` name (`pulsar` for a pulsing neutron star) of the star `star_id`
    when it is a neutron star or black hole, else `None`: the host a planet of it is rated against."""
    if star_id is None:
        return None
    row = conn.execute(
        "SELECT st.yerkes_class, (SELECT ns.pulsar_type FROM neutron_stars ns WHERE ns.star_id = st.id LIMIT 1)"
        " AS pulsar_type FROM stars st WHERE st.id = ?", (star_id,)).fetchone()
    if row is None:
        return None
    return habitability_world.compact_host_name(
        SimpleNamespace(yerkes_class=row["yerkes_class"], pulsar_type=row["pulsar_type"] or "non-pulsing"))


def _stored_native_life(row):
    """The `habitability.LIFE_STAGES` stage of a stored planet or moon `row` (`planets` or `moons`), read
    from its life class and the milestone its timeline reached (`habitability_world.native_life`), for a re-score."""
    return habitability_world.native_life(row["life_class"], row["life_stage"], row["life_system_age_gyr"] is not None)


def stored_world(conn, table, row, with_life=True):
    """
    The `habitability.World` of a stored rocky planet or moon `row` (`planets` or `moons`), as its stored score
    was worked out (a pulsar's, neutron star's or black hole's planet at its lethal dose), or `None` for a gas
    giant or a body saved without a score (POP.12). `with_life` false leaves out the body's native life
    (two queries a body) when only the physical conditions matter.
    """
    if row["body_type"] == "g" or row["surface_temperature_k"] is None or row["energy_flux_w_m2"] is None:
        return None
    values = {**with_gases(conn, [row], table == "moons")[0], "is_moon": table == "moons",
              "native_life": _stored_native_life(row) if with_life else "none"}
    return habitability_world.world_from(values, lethal_host=compact_host_of_star(conn, row["star_id"]))


def _refresh_cloud_doses(conn, system_changes):
    """
    GEN.87: re-derives the cosmic-ray multiplier (`dose_helio_mult`) and the
    surface dose of every rocky world in each system whose container
    changed. A nebula or remnant presses the star's heliosphere in and the
    star's wind stops shielding the system (`radiation.heliosphere_compression`);
    leaving it restores 1.

    Args:
        conn (Connection): Part of the caller's transaction.
        system_changes (iterable): `(star_system_id, inside_nebula_id,
            inside_remnant_id)` for each system that changed container.
    """
    for system_id, nebula_id, remnant_id in system_changes:
        compression = 0.0
        cloud = surrounding_cloud(conn, {"inside_nebula_id": nebula_id, "inside_remnant_id": remnant_id})
        if cloud is not None:
            open_km = conn.execute(
                "SELECT COALESCE(s.binary_heliosphere_radius_km, (SELECT MAX(heliosphere_radius_km) FROM stars"
                " WHERE star_system_id = s.id)) AS radius_km FROM star_systems s WHERE s.id = ?",
                (system_id,)).fetchone()["radius_km"]
            if open_km:
                open_au = _km_to_au(open_km)
                squeezed = compressed_heliosphere_radius(open_au, cloud["density_cm3"], cloud["temperature_k"])
                compression = radiation.heliosphere_compression(open_au, squeezed)
        compact = conn.execute(
            "SELECT 1 FROM stars WHERE star_system_id = ? AND yerkes_class IN ('NS', 'BH') LIMIT 1",
            (system_id,)).fetchone() is not None
        for table in ("planets", "moons"):
            rows = []
            for row in with_life_inputs(conn, conn.execute(
                f"SELECT * FROM {table} WHERE star_system_id = ? AND dose_gcr_msv_yr IS NOT NULL", (system_id,)
            ).fetchall(), table == "moons"):
                column = radiation.column_g_cm2(row["atmospheric_pressure_pa"], row["gravity_g"])
                multiplier = radiation.helio_multiplier(column, compression)
                total = radiation.total_dose_msv_yr(
                    row["dose_gcr_msv_yr"], multiplier, row["dose_sep_msv_yr"], row["dose_ground_msv_yr"], row["dose_belt_msv_yr"])
                scored = ()
                if not compact and row["energy_flux_w_m2"] is not None:
                    # GEN.89: the dose is one of the score's inputs. (A compact host's planet is
                    # rated at a lethal dose whatever the cosmic rays do.)
                    scores = habitability_world.compute({**row, "surface_dose_msv_yr": total,
                                                         "dose_helio_mult": multiplier, "is_moon": table == "moons",
                                                         "native_life": _stored_native_life(row)})
                    scored = tuple(scores[name] for name in habitability_world.SCORE_FIELDS)
                rows.append((row["id"], multiplier, total, *scored))
            if rows:
                columns = ("dose_helio_mult", "surface_dose_msv_yr")
                updates = [r for r in rows if len(r) == 3]
                scored_updates = [r for r in rows if len(r) > 3]
                if updates:
                    _update_by_id(conn, table, columns, updates, touch=True)
                if scored_updates:
                    _update_by_id(conn, table, columns + habitability_world.SCORE_FIELDS, scored_updates, touch=True)
    from planetgen.db import block_stats    # MAP.143: a new score can move a sector's habitable count

    block_stats.refresh_sectors(conn, recount_habitable_worlds(conn, [change[0] for change in system_changes]))


def refresh_habitability(conn, system_id):
    """
    GEN.217: re-derives the stored habitability numbers of every scored rocky planet and moon of one system from
    its stored conditions: `flux_variation` and the `habitability_world.SCORE_FIELDS` (the index, the six domains
    with their tiers and the equipment tier). It is how a body saved before a change to the rules
    (GEN.148's flux swing, GEN.207's six domains, GEN.215's airless Air) gets the numbers a new one would, without
    regenerating it. A system of one star has the swing of the planet's own orbit (`orbit_elements.flux_swing`, 0 for
    a circle); a system of two is loaded and its swing worked out as
    generation does (`insolation.refresh_system`). A planet of a pulsar, neutron star or black hole is rated against
    its host (`compact_host_of_star`). A gas giant, and a body saved before v76 (no `energy_flux_w_m2`), are left alone.

    Args:
        conn (Connection): Part of the caller's transaction.
        system_id (int): The `star_systems` row.

    Returns:
        int: Bodies updated.
    """
    swings = {}
    if conn.execute("SELECT COUNT(*) AS n FROM stars WHERE star_system_id = ?", (system_id,)).fetchone()["n"] > 1:
        system = load_star_system(conn, system_id)
        insolation.refresh_system(system)
        for planet in system.all_planets:
            if planet.body_type == "a":
                continue
            swings[("planets", planet.db_id)] = planet.flux_variation
            for moon in planet.moons:
                swings[("moons", moon.db_id)] = moon.flux_variation
    else:
        # a lone star: the swing of the planet's own orbit (GEN.727), which a moon shares with its planet
        planet_swings = {}
        for row in conn.execute("SELECT id, eccentricity FROM planets WHERE star_system_id = ?", (system_id,)).fetchall():
            planet_swings[row["id"]] = swings[("planets", row["id"])] = orbit_elements.flux_swing(row["eccentricity"])
        for row in conn.execute("SELECT id, planet_id FROM moons WHERE star_system_id = ?", (system_id,)).fetchall():
            swings[("moons", row["id"])] = planet_swings.get(row["planet_id"], 0.0)
    columns = habitability_world.SCORE_FIELDS + ("flux_variation",)
    updated = 0
    for table in ("planets", "moons"):
        rows = []
        for row in with_life_inputs(conn, conn.execute(
                f"SELECT * FROM {table} WHERE star_system_id = ? AND energy_flux_w_m2 IS NOT NULL AND body_type <> 'g'",
                (system_id,)).fetchall(), table == "moons"):
            swing = swings.get((table, row["id"]), 0.0)
            scores = habitability_world.compute(
                {**row, "flux_variation": swing, "is_moon": table == "moons",
                 "native_life": _stored_native_life(row)},
                compact_host_of_star(conn, row["star_id"]))
            rows.append((row["id"], *(scores[name] for name in habitability_world.SCORE_FIELDS), swing))
        if rows:
            _update_by_id(conn, table, columns, rows, touch=True)
            updated += len(rows)
    return updated


PLACED_PHENOMENON_TABLES = kind_tables.tables("phenomenon", order="placed")
"""tuple: The phenomenon tables with a galaxy-frame center, each with a
`quadrant` column since v41."""

NEAREST_SYSTEMS_COUNT = 3
"""int: How many nearest star systems `nearest_systems` keeps per object."""

NEAREST_SYSTEMS_SEARCH_PC = 4.0
"""float: How far (parsecs, one sector edge) `refresh_nearest_systems`
looks for an object's nearest systems. An object with fewer systems than
`NEAREST_SYSTEMS_COUNT` inside that distance keeps fewer rows."""

_NEAREST_GRID_CELL_PC = 1.0


class _SystemGrid:
    """Star systems bucketed into `_NEAREST_GRID_CELL_PC` cubes, for
    nearest-neighbor searches without numpy."""

    def __init__(self, systems):
        self.cells = {}
        self.systems = list(systems)
        for system_id, point in self.systems:
            self.cells.setdefault(self._key(point), []).append((system_id, point))

    @staticmethod
    def _key(point):
        return tuple(math.floor(c / _NEAREST_GRID_CELL_PC) for c in point)

    def nearest(self, point, count=NEAREST_SYSTEMS_COUNT, limit_pc=NEAREST_SYSTEMS_SEARCH_PC, exclude=None):
        """The `count` nearest `(distance_pc, system_id)` within
        `limit_pc` of `point`, nearest first, leaving out `exclude`."""
        cx, cy, cz = self._key(point)
        best = []
        max_shell = int(math.ceil(limit_pc / _NEAREST_GRID_CELL_PC)) + 1
        if len(self.systems) < (2 * max_shell + 1) ** 3:
            # Fewer systems than cells to visit (one new sector's systems,
            # merged into its neighbors' lists): checking each is cheaper
            # than walking mostly empty shells.
            best = sorted(
                (distance, system_id)
                for system_id, other in self.systems
                if system_id != exclude and (distance := detmath.dist(point, other)) <= limit_pc
            )
            return best[:count]
        for shell in range(max_shell + 1):
            if len(best) >= count and best[count - 1][0] <= (shell - 1) * _NEAREST_GRID_CELL_PC:
                break
            for dx in range(-shell, shell + 1):
                for dy in range(-shell, shell + 1):
                    for dz in range(-shell, shell + 1):
                        if max(abs(dx), abs(dy), abs(dz)) != shell:
                            continue
                        for system_id, other in self.cells.get((cx + dx, cy + dy, cz + dz), ()):
                            if system_id == exclude:
                                continue
                            distance = detmath.dist(point, other)
                            if distance <= limit_pc:
                                best.append((distance, system_id))
            best.sort()
            del best[count:]
        return best


def _sector_centers(conn, sector_ids):
    """`{sector_id: (x, y, z)}` for the placed sectors among `sector_ids`."""
    centers = {}
    sector_ids = sorted(set(sector_ids))
    for start in range(0, len(sector_ids), 500):
        batch = sector_ids[start:start + 500]
        marks = ", ".join("?" * len(batch))
        for row in conn.execute(
            f"SELECT id, center_x_pc, center_y_pc, center_z_pc FROM sectors"
            f" WHERE id IN ({marks}) AND center_x_pc IS NOT NULL",
            tuple(batch),
        ).fetchall():
            centers[row["id"]] = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
    return centers


def _placed_systems(conn, centers):
    """`(system_id, sector_id, galaxy_point)` for every positioned system
    in the sectors `centers` maps."""
    systems = []
    ids = sorted(centers)
    for start in range(0, len(ids), 500):
        batch = ids[start:start + 500]
        marks = ", ".join("?" * len(batch))
        for row in conn.execute(
            f"SELECT id, sector_id, position_x_mpc, position_y_mpc, position_z_mpc FROM star_systems"
            f" WHERE sector_id IN ({marks}) AND position_x_mpc IS NOT NULL",
            tuple(batch),
        ).fetchall():
            offset = (row["position_x_mpc"] / 1000.0, row["position_y_mpc"] / 1000.0, row["position_z_mpc"] / 1000.0)
            systems.append((row["id"], row["sector_id"], local_to_galaxy_pc(centers[row["sector_id"]], offset)))
    return systems


def _placed_objects(conn, centers):
    """`(table, object_id, sector_id, galaxy_point)` for every placed star
    system and phenomenon in the sectors `centers` maps."""
    objects = [("star_systems", system_id, sector_id, point)
               for system_id, sector_id, point in _placed_systems(conn, centers)]
    ids = sorted(centers)
    for start in range(0, len(ids), 500):
        batch = ids[start:start + 500]
        marks = ", ".join("?" * len(batch))
        for table in PLACED_PHENOMENON_TABLES:
            for row in conn.execute(
                f"SELECT id, sector_id, center_x_pc, center_y_pc, center_z_pc FROM {table}"
                f" WHERE sector_id IN ({marks}) AND center_x_pc IS NOT NULL",
                tuple(batch),
            ).fetchall():
                objects.append((table, row["id"], row["sector_id"],
                                (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])))
    return objects


def _stored_nearest(conn, sector_ids):
    """`{(table, object_id): [(distance_pc, system_id), ...]}`, nearest
    first, from `nearest_systems` for objects in `sector_ids`."""
    stored = {}
    sector_ids = sorted(set(sector_ids))
    for start in range(0, len(sector_ids), 500):
        batch = sector_ids[start:start + 500]
        marks = ", ".join("?" * len(batch))
        for row in conn.execute(
            f"SELECT object_table, object_id, neighbor_system_id, distance_pc FROM nearest_systems"
            f" WHERE sector_id IN ({marks}) ORDER BY object_table, object_id, neighbor_rank",
            tuple(batch),
        ).fetchall():
            stored.setdefault((row["object_table"], row["object_id"]), []).append(
                (row["distance_pc"], row["neighbor_system_id"]))
    return stored


def _same_neighbors(old, new):
    return [s for _d, s in old] == [s for _d, s in new] and all(
        abs(a - b) < 1e-9 for (a, _s), (b, _t) in zip(old, new))


def _write_nearest(conn, rows_by_object, sector_by_object):
    """Replaces the stored neighbor list of each object in
    `rows_by_object` (`{(table, id): [(distance_pc, system_id), ...]}`)."""
    if not rows_by_object:
        return
    keys = list(rows_by_object)
    for table in {table for table, _id in keys}:
        ids = [object_id for t, object_id in keys if t == table]
        for start in range(0, len(ids), 500):
            batch = ids[start:start + 500]
            marks = ", ".join("?" * len(batch))
            conn.execute(f"DELETE FROM nearest_systems WHERE object_table = ? AND object_id IN ({marks})",
                         (table, *batch))
    values = [
        (sector_by_object[key], key[0], key[1], key[1] if key[0] == "star_systems" else None,
         rank, system_id, distance)
        for key, neighbors in rows_by_object.items()
        for rank, (distance, system_id) in enumerate(neighbors, start=1)
    ]
    if values:
        conn.executemany(
            "INSERT INTO nearest_systems (sector_id, object_table, object_id, star_system_id,"
            " neighbor_rank, neighbor_system_id, distance_pc) VALUES (?, ?, ?, ?, ?, ?, ?)",
            values,
        )


def _octant_updates(conn, objects, centers):
    """Sets each phenomenon's `quadrant` (v41) from where its center sits
    in its sector, writing only rows that change."""
    by_table = {}
    for table, object_id, sector_id, point in objects:
        if table == "star_systems":
            continue
        label, _magnitudes = classify_octant(galaxy_to_local_pc(centers[sector_id], point))
        by_table.setdefault(table, []).append((label, object_id))
    for table, pairs in by_table.items():
        ids = [object_id for _label, object_id in pairs]
        current = {}
        for start in range(0, len(ids), 500):
            batch = ids[start:start + 500]
            marks = ", ".join("?" * len(batch))
            for row in conn.execute(f"SELECT id, quadrant FROM {table} WHERE id IN ({marks})", tuple(batch)).fetchall():
                current[row["id"]] = row["quadrant"]
        changed = [(object_id, label) for label, object_id in pairs if current.get(object_id) != label]
        _update_by_id(conn, table, ("quadrant",), changed)


def _update_by_id(conn, table, columns, rows, touch=False):
    """
    Sets `columns` on many rows of `table` in as few statements as
    possible (PERF.13): one `UPDATE ... SET col = CASE id WHEN ? THEN ?
    ... END WHERE id IN (...)` per 500 rows, rather than one UPDATE per
    row. `modified_at` is left alone unless `touch`.

    Args:
        conn (Connection): Part of the caller's transaction.
        table (str): One of this module's own table names.
        columns (tuple): Column names to set.
        rows (list): `(id, value, ...)` with one value per column.
        touch (bool): Let `modified_at` update as usual.
    """
    for first in range(0, len(rows), _BATCH_ROWS):
        chunk = rows[first:first + _BATCH_ROWS]
        cases = " ".join(["WHEN ? THEN ?"] * len(chunk))
        sets = [f"{column} = CASE id {cases} END" for column in columns]
        if not touch:
            sets.append("modified_at = modified_at")
        params = [value for index in range(len(columns)) for row in chunk for value in (row[0], row[index + 1])]
        conn.execute(
            f"UPDATE {table} SET {', '.join(sets)} WHERE id IN ({', '.join('?' * len(chunk))})",
            (*params, *(row[0] for row in chunk)),
        )


def _sectors_near(conn, centers, reach_pc):
    """The placed sectors whose centers lie within `reach_pc` of any of
    `centers`' values, as `{sector_id: center}` (bounding box, then exact)."""
    if not centers:
        return {}
    low = tuple(min(c[i] for c in centers.values()) - reach_pc for i in range(3))
    high = tuple(max(c[i] for c in centers.values()) + reach_pc for i in range(3))
    rows = conn.execute(
        "SELECT id, center_x_pc, center_y_pc, center_z_pc FROM sectors"
        " WHERE center_x_pc BETWEEN ? AND ? AND center_y_pc BETWEEN ? AND ? AND center_z_pc BETWEEN ? AND ?",
        (low[0], high[0], low[1], high[1], low[2], high[2]),
    ).fetchall()
    near = {}
    points = list(centers.values())
    for row in rows:
        center = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
        if any(detmath.dist(center, p) <= reach_pc for p in points):
            near[row["id"]] = center
    return near


def _edge_pc(conn):
    row = conn.execute("SELECT MAX(edge_mpc) AS edge FROM sectors").fetchone()
    return (row["edge"] or 0) / 1000.0


def refresh_nearest_systems(conn, sector_ids, on_progress=None):
    """
    Recomputes the stored nearest star systems (`nearest_systems`, v41)
    and each phenomenon's `quadrant` for every placed system and
    phenomenon filed under `sector_ids`, searching every sector within
    `NEAREST_SYSTEMS_SEARCH_PC`. Rows that didn't change aren't written.
    Called by the v41 migration and the correlative update, which move
    things; `insert_sector` uses the cheaper `_add_sector_to_nearest`.

    Args:
        conn (Connection): Part of the caller's transaction.
        sector_ids (iterable): `sectors.id` values; unplaced ones are skipped.
        on_progress (callable, optional): `on_progress("Sectors", done, total)`
            after each batch of sectors.

    Returns:
        set: The `(table, id)` of every object whose list changed.
    """
    sector_ids = sorted(set(sector_ids))
    half_diagonal = _edge_pc(conn) * math.sqrt(3) / 2
    all_changed = set()
    with steps.hooked(on_progress, "Nearest systems: sectors", "nearest", len(sector_ids)) as report:
        for start in range(0, len(sector_ids), 200):
            centers = _sector_centers(conn, sector_ids[start:start + 200])
            _report(report, "Sectors", start, len(sector_ids))
            if not centers:
                continue
            near = _sectors_near(conn, centers, 2 * half_diagonal + NEAREST_SYSTEMS_SEARCH_PC)
            grid = _SystemGrid((system_id, point) for system_id, _sector, point in _placed_systems(conn, near))
            objects = _placed_objects(conn, centers)
            _octant_updates(conn, objects, centers)
            stored = _stored_nearest(conn, centers)
            changed, sectors = {}, {}
            for table, object_id, sector_id, point in objects:
                exclude = object_id if table == "star_systems" else None
                neighbors = grid.nearest(point, exclude=exclude)
                key = (table, object_id)
                if not _same_neighbors(stored.get(key, []), neighbors):
                    changed[key] = neighbors
                    sectors[key] = sector_id
            _write_nearest(conn, changed, sectors)
            all_changed.update(changed)
        _report(report, "Sectors", len(sector_ids), len(sector_ids))
    return all_changed


def _add_sector_to_nearest(conn, sector_id):
    """
    Fills `nearest_systems` and `quadrant` for a newly generated sector
    (`refresh_nearest_systems`), then merges its systems into the lists of
    objects in neighboring sectors that already had theirs -- adding
    systems can only bring neighbors nearer, so the stored lists plus the
    new sector's systems are enough.
    """
    refresh_nearest_systems(conn, [sector_id])
    _merge_into_neighbors(conn, [sector_id])


def _merge_into_neighbors(conn, sector_ids, skip=()):
    """
    Merges the systems of `sector_ids` into the nearest-system lists of the
    objects in the sectors around them (not those in `sector_ids` or
    `skip`, which `refresh_nearest_systems` computes whole).
    """
    own = _sector_centers(conn, sector_ids)
    if not own:
        return
    half_diagonal = _edge_pc(conn) * math.sqrt(3) / 2
    near = _sectors_near(conn, own, 2 * half_diagonal + NEAREST_SYSTEMS_SEARCH_PC)
    for sector_id in (*sector_ids, *skip):
        near.pop(sector_id, None)
    if not near:
        return
    new_systems = [(system_id, point) for system_id, _sector, point in _placed_systems(conn, own)]
    if not new_systems:
        return
    grid = _SystemGrid(new_systems)
    stored = _stored_nearest(conn, near)
    changed, sectors = {}, {}
    for table, object_id, other_sector, point in _placed_objects(conn, near):
        key = (table, object_id)
        old = stored.get(key, [])
        merged = sorted(old + grid.nearest(point))[:NEAREST_SYSTEMS_COUNT]
        if not _same_neighbors(old, merged):
            changed[key] = merged
            sectors[key] = other_sector
    _write_nearest(conn, changed, sectors)


LINK_BATCH_SECTORS = 20
"""int: Sectors `link_sector_neighbors` links in one transaction, which holds the neighbor lock and moves the
progress bar once. 200 was a transaction of about ten minutes on a dense run, a bar that stood still for as
long and a lock nothing else could take (PERF.84); batches of 2 to 12 cost the same in total."""

LINK_CELL_PC = 24.0
"""float: Edge of the cells `link_sector_neighbors` sorts sectors into, so a batch is a neighbourhood and the
systems it loads around itself are few."""

LINK_PHASES = 3
"""int: Steps `link_sector_neighbors` takes per sector (its progress units)."""


def link_order(config, sector_ids):
    """`sector_ids` ordered by the `LINK_CELL_PC` cell of their centers (then id), so consecutive sectors sit near
    each other: a batch of them loads the systems around one place, not around the whole run (PERF.84)."""
    conn = get_connection(config)
    try:
        centers = _sector_centers(conn, sector_ids)
    finally:
        conn.close()

    def cell(sector_id):
        center = centers.get(sector_id)
        if center is None:
            return (0, 0, 0, sector_id)
        return (*(math.floor(c / LINK_CELL_PC) for c in center), sector_id)

    return sorted(sector_ids, key=cell)


def link_sector_neighbors(config, sector_ids, on_progress=None, skip=None):
    """
    PERF.45: the neighbour step of saving sectors, done once for many (a
    `galaxy` run saves its sectors with `link_neighbors=False` and links
    them here at its end): containment in nebulae and remnants, each
    object's nearest systems and quadrant, and the new systems merged into
    the lists of the already-linked sectors around. Done in batches of
    `LINK_BATCH_SECTORS`, each a transaction under the neighbor lock (so an
    on-demand sector saved meanwhile is not missed), retried on a deadlock.
    The result equals linking each sector as it was saved.

    Args:
        config (MySQLConfig): Connection parameters.
        sector_ids (iterable): The `sectors.id` values to link.
        on_progress (callable, optional): `on_progress(done, total, step)`
            in steps (`LINK_PHASES` per sector: containment, nearest
            systems, the neighbours' lists), `step` naming the one being
            done (UX.83).
        skip (iterable, optional): The sectors of the run that are linked whole by their own call (PERF.113: a
            unit of the link stage links a chunk of the run): systems are not merged into their lists. Defaults to
            `sector_ids`.

    Returns:
        int: How many sectors were linked.
    """
    sector_ids = sorted(set(sector_ids))
    skip = set(sector_ids) if skip is None else set(skip) | set(sector_ids)
    total = len(sector_ids) * LINK_PHASES
    sector_ids = link_order(config, sector_ids)

    def report(base, phase, label, size):
        if on_progress is not None:
            on_progress(base * LINK_PHASES + phase * size, total, label)

    def link(batch, base):
        def work(conn):
            conn.lock_until_commit(_neighbor_lock_name(conn))
            report(base, 0, "containment", len(batch))
            refresh_containment(conn, batch)
            report(base, 1, "nearest systems", len(batch))
            refresh_nearest_systems(conn, batch)
            report(base, 2, "neighbours' lists", len(batch))
            _merge_into_neighbors(conn, batch, skip=skip)
        return work

    for start in range(0, len(sector_ids), LINK_BATCH_SECTORS):
        batch = sector_ids[start:start + LINK_BATCH_SECTORS]
        _save_with_retries(config, [], link(batch, start))
        report(start + len(batch), 0, "", 0)
    return len(sector_ids)


# ---------------------------------------------------------------------------
# Facilities (schema v42) -- see planetgen/population/facilities.py
# for the rules and the orbit math.
# ---------------------------------------------------------------------------

class FacilityError(ValueError):
    """A facility that breaks the placement rules, or whose host doesn't
    exist (`not_found`)."""

    def __init__(self, message, not_found=False):
        super().__init__(message)
        self.not_found = not_found


def _facility_host(conn, host_type, host_id):
    """
    Reads one facility host: `(columns, mass_kg, radius_km, body_type,
    sphere_km)`, where `columns` are the `facilities` host columns to set.
    Mass, radius and sphere (the edge of its sphere of influence:
    `facilities.star_host`'s heliosphere for a star, the Hill sphere for a
    planet or moon) describe what an orbit circles: the star for an
    asteroid belt (the pair, for a close pair's belt), `None` for hosts
    nothing orbits.

    Raises:
        FacilityError: If no such host exists (`not_found`).
    """
    missing = FacilityError(f"no such {host_type.replace('_', ' ')}: {host_id}", not_found=True)
    if host_type in ("star", "asteroid_belt"):
        if host_type == "star":
            row = conn.execute("SELECT star_system_id, id AS star_id FROM stars WHERE id = ?", (host_id,)).fetchone()
        else:
            row = conn.execute("SELECT star_system_id, star_id FROM asteroid_belts WHERE id = ?",
                               (host_id,)).fetchone()
        if row is None:
            raise missing
        system = conn.execute("SELECT binary_configuration, binary_separation_km, binary_heliosphere_radius_km"
                              " FROM star_systems WHERE id = ?", (row["star_system_id"],)).fetchone()
        stars = conn.execute("SELECT id, role, mass_kg, radius_km, heliosphere_radius_km FROM stars"
                             " WHERE star_system_id = ? ORDER BY id", (row["star_system_id"],)).fetchall()
        star_id = row["star_id"]
        if star_id is None:
            # A single star's or a close pair's belt: around the primary
            # (the pair as one, for a close pair).
            star_id = next((star["id"] for star in stars if star["role"] in ("single", "primary")), stars[0]["id"])
        mass, radius, sphere = facility_rules.star_host(
            [dict(star) for star in stars], star_id, system["binary_configuration"],
            system["binary_separation_km"], system["binary_heliosphere_radius_km"])
        if host_type == "star":
            return {"star_system_id": row["star_system_id"], "star_id": host_id}, mass, radius, None, sphere
        return ({"star_system_id": row["star_system_id"], "asteroid_belt_id": host_id},
                mass, radius, None, sphere)
    if host_type in ("planet", "moon"):
        table = "planets" if host_type == "planet" else "moons"
        row = conn.execute(f"SELECT star_system_id, mass_kg, radius_km, body_type, hill_radius_km FROM {table}"
                           " WHERE id = ?", (host_id,)).fetchone()
        if row is None:
            raise missing
        return ({"star_system_id": row["star_system_id"], f"{host_type}_id": host_id},
                row["mass_kg"], row["radius_km"], row["body_type"], row["hill_radius_km"])
    if host_type == "asteroid_field":
        if conn.execute("SELECT 1 FROM asteroid_fields WHERE id = ?", (host_id,)).fetchone() is None:
            raise missing
        return {"asteroid_field_id": host_id}, None, None, None, None
    if host_type == "space":
        if conn.execute("SELECT 1 FROM sectors WHERE id = ?", (host_id,)).fetchone() is None:
            raise FacilityError(f"no such sector: {host_id}", not_found=True)
        return {"sector_id": host_id}, None, None, None, None
    raise FacilityError(f"unknown host type {host_type!r}")


def facility_orbit(conn, host_type, host_id, distance_km=None):
    """
    The circular orbit an orbital facility would have around a star,
    planet or moon (`facilities.orbit_for`), without saving anything --
    what the web form shows before saving -- plus the range the form's
    slider covers (`facilities.orbit_limits`): `min_distance_km` and
    `max_distance_km`.

    Raises:
        FacilityError: If the host can't be orbited, doesn't exist, or the
            distance is inside it or outside its sphere of influence.
    """
    if host_type not in ("star", "planet", "moon"):
        raise FacilityError(f"nothing orbits a {host_type.replace('_', ' ')}")
    _columns, mass, radius, _body_type, sphere = _facility_host(conn, host_type, host_id)
    lowest, highest = facility_rules.orbit_limits(radius, sphere)
    try:
        orbit = facility_rules.orbit_for(mass, radius, distance_km, highest)
    except ValueError as exc:
        raise FacilityError(str(exc)) from exc
    return {**orbit, "min_distance_km": lowest, "max_distance_km": highest}


def _belt_orbit(conn, belt_id, mass, radius):
    """A spot in an asteroid belt for an asteroid facility
    (`facilities.belt_position`) and the circular orbit around its star
    from there: `orbit_for`'s dict plus `phase_deg`."""
    belt = conn.execute("SELECT lower_limit_km, upper_limit_km FROM asteroid_belts WHERE id = ?",
                        (belt_id,)).fetchone()
    distance_km, phase_deg = facility_rules.belt_position(belt["lower_limit_km"], belt["upper_limit_km"])
    try:
        orbit = facility_rules.orbit_for(mass, radius, max(distance_km, radius * 1.01))
    except ValueError as exc:
        raise FacilityError(str(exc)) from exc
    orbit["phase_deg"] = phase_deg
    return orbit


def add_facility(conn, name, kind, placement, host_type, host_id, distance_km=None, phase_deg=None,
                 offset_ly=None, description=None, attrs=None):
    """
    Stores one facility after checking it against the placement rules
    (`facilities.check_facility`).

    Args:
        conn (Connection): Part of the caller's transaction.
        name (str): What it's called.
        kind (str): A `tuning.FACILITY_KINDS` key.
        placement (str): `terrestrial`, `orbital`, `asteroid` or `standalone`.
        host_type (str): `star`, `planet`, `moon`, `asteroid_belt`,
            `asteroid_field`, or `space` (then `host_id` is a sector).
        host_id (int): The host row's id.
        distance_km (float, optional): An orbital facility's orbit radius
            (`facilities.orbit_for`'s default otherwise), inside its
            host's sphere of influence (`facilities.orbit_limits`).
        phase_deg (float, optional): Where along its orbit it starts;
            random otherwise.
        offset_ly (tuple, optional): A stand-alone facility's `(x, y, z)`
            from its sector's center, light-years, along the sector's own
            axes (`frame.sector_orientation`). The center if left
            out.
        description (str, optional): Free text.
        attrs (dict, optional): `type_id`, `affiliation_polity_id` or `affiliation_name`, `ratings`,
            `ratings_note` and `fields`, as `facility_attrs.resolve` reads them (POP.10).

    An asteroid facility in a belt takes no distance: it gets a random
    spot in the belt (`facilities.belt_position`) and the circular orbit
    around its star from there, which `advance_facility_orbits` moves
    like an orbital facility's.

    Returns:
        int: The new `facilities.id`.

    Raises:
        FacilityError: If the rules refuse it or the host doesn't exist.
        FacilityTypeError: If the type, affiliation, ratings or field values are refused.
    """
    columns, mass, radius, body_type, sphere = _facility_host(conn, host_type, host_id)
    problem = facility_rules.check_facility(kind, placement, host_type, body_type)
    if problem:
        raise FacilityError(problem)
    plan = facility_attrs.resolve(conn, kind, attrs) if attrs else None

    orbit = {}
    if placement == "orbital":
        try:
            orbit = facility_rules.orbit_for(mass, radius, distance_km, facility_rules.orbit_limits(radius, sphere)[1])
        except ValueError as exc:
            raise FacilityError(str(exc)) from exc
        if phase_deg is None:
            phase_deg = draw.uniform(0.0, 360.0)
        orbit["phase_deg"] = phase_deg % 360.0
    elif distance_km is not None or phase_deg is not None:
        raise FacilityError("only an orbital facility takes an orbit distance or phase")
    elif host_type == "asteroid_belt":
        orbit = _belt_orbit(conn, host_id, mass, radius)

    placement_values = (None, None, None, None)
    velocity_kms = (0.0, 0.0, 0.0)
    if host_type == "space":
        sector = conn.execute("SELECT center_x_pc, center_y_pc, center_z_pc, edge_mpc FROM sectors WHERE id = ?",
                              (host_id,)).fetchone()
        if sector["center_x_pc"] is None:
            raise FacilityError("that sector isn't placed in the galaxy")
        offset_ly = tuple(offset_ly) if offset_ly is not None else (0.0, 0.0, 0.0)
        half_edge_ly = milliparsecs_to_ly(sector["edge_mpc"]) / 2
        if len(offset_ly) != 3 or any(not math.isfinite(c) or abs(c) > half_edge_ly for c in offset_ly):
            raise FacilityError(f"the offset must be three numbers within {half_edge_ly:.2f} ly of the center")
        center = (sector["center_x_pc"], sector["center_y_pc"], sector["center_z_pc"])
        point = local_to_galaxy_pc(center, tuple(ly_to_pc(c) for c in offset_ly))
        placement_values = (*point, math.sqrt(sum(c * c for c in point)))
        # GEN.125: on the galaxy's rotation curve at its place.
        speed_kms, _period_gy = calculate_galactic_orbit(pc_to_ly(detmath.hypot(point[0], point[1])))
        velocity_kms = tuple(v / 1000.0 for v in galactic_velocity_ms(point, speed_kms))
    elif offset_ly is not None:
        raise FacilityError("only a stand-alone facility has a position in space")

    birth_sector_id = columns.get("sector_id")
    if birth_sector_id is None and columns.get("star_system_id") is not None:
        row = conn.execute("SELECT sector_id FROM star_systems WHERE id = ?", (columns["star_system_id"],)).fetchone()
        birth_sector_id = None if row is None else row["sector_id"]
    [facility_uid] = _runtime_uids(conn, _sector_birth_address(conn, birth_sector_id), 1)   # born now (GEN.172)
    cur = conn.execute(
        """
        INSERT INTO facilities (
            name, kind, placement, host_type, star_system_id, star_id, planet_id, moon_id,
            asteroid_belt_id, asteroid_field_id, sector_id,
            center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc,
            velocity_x_kms, velocity_y_kms, velocity_z_kms,
            orbit_distance_km, orbit_period_years, orbital_speed_kms, orbit_phase_deg, description, uid, mass_kg
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            name, kind, placement, host_type, columns.get("star_system_id"), columns.get("star_id"),
            columns.get("planet_id"), columns.get("moon_id"), columns.get("asteroid_belt_id"),
            columns.get("asteroid_field_id"), columns.get("sector_id"),
            *placement_values, *velocity_kms,
            orbit.get("distance_km"), orbit.get("period_years"), orbit.get("orbital_speed_kms"),
            orbit.get("phase_deg"), description, facility_uid,
            # GEN.337: the mass of what it sits on (a star, planet or moon); a belt, a field and open space have none
            # that a facility could be said to hold, a stated zero.
            mass if host_type in ("star", "planet", "moon") and mass else 0.0,
        ),
    )
    if plan is not None:
        facility_attrs.write(conn, cur.lastrowid, plan)
    if columns.get("star_system_id") is not None:
        touch_star_system(conn, columns["star_system_id"])
        record_system_interests(conn, "ss.id = ?", (columns["star_system_id"],))    # a colony makes its world inhabited (MAP.321)
    return cur.lastrowid


def update_facility(conn, facility_id, name=None, description=None, attrs=None):
    """
    Changes a facility's name, description, type, affiliation, ratings or field values (POP.10); its place never
    changes. Only what is given changes: `attrs` is read like `add_facility`'s, and a key sent as `None` clears it.

    Returns:
        bool: `False` when there is no such facility.

    Raises:
        FacilityTypeError: If the type, affiliation, ratings or field values are refused.
    """
    row = conn.execute("SELECT kind, star_system_id FROM facilities WHERE id = ?", (facility_id,)).fetchone()
    if row is None:
        return False
    plan = facility_attrs.resolve(conn, row["kind"], attrs or {}, current=facility_attrs.current_attrs(conn, facility_id))
    if name is not None:
        plan["columns"]["name"] = name
    if description is not None:
        plan["columns"]["description"] = description or None
    facility_attrs.write(conn, facility_id, plan)
    if row["star_system_id"] is not None:
        touch_star_system(conn, row["star_system_id"])
        record_system_interests(conn, "ss.id = ?", (row["star_system_id"],))
    return True


def delete_facility(conn, facility_id):
    """Deletes one facility. Returns `False` if there was none."""
    row = conn.execute("SELECT star_system_id FROM facilities WHERE id = ?", (facility_id,)).fetchone()
    if row is None:
        return False
    storage_checksum.mark_rows(conn, "facilities", [facility_id])
    conn.execute("DELETE FROM facilities WHERE id = ?", (facility_id,))
    if row["star_system_id"] is not None:
        touch_star_system(conn, row["star_system_id"])
        record_system_interests(conn, "ss.id = ?", (row["star_system_id"],))
    return True


# ---------------------------------------------------------------------------
# Galactic motion (GEN.6): the correlative update moves every placed
# star system, phenomenon and stand-alone facility along its galactic orbit,
# refiles it under whichever generated sector it drifted into, then refreshes
# containment, octants, nearest systems and location text.
# ---------------------------------------------------------------------------

def _rotate_about_axis(point, angle_rad):
    """`point` turned counterclockwise (seen from galactic north) about the
    galactic axis -- the direction galactic phase increases."""
    x, y, z = point
    cos_a, sin_a = math.cos(angle_rad), math.sin(angle_rad)
    return (x * cos_a - y * sin_a, x * sin_a + y * cos_a, z)


def _galactic_turn(elapsed_years, period_gy):
    """The angle (radians) an orbit of `period_gy` sweeps in `elapsed_years`,
    whole orbits dropped first (GEN.108: a step longer than an orbit wraps
    exactly instead of losing the fraction to rounding)."""
    if not period_gy or period_gy <= 0 or elapsed_years <= 0:
        return 0.0
    period_years = period_gy * 1e9
    return 2 * math.pi * math.fmod(elapsed_years, period_years) / period_years


class _SectorIndex:
    """Every placed grid sector by address, for refiling moved objects."""

    def __init__(self, conn):
        self.edge_pc = _edge_pc(conn)
        self.by_address = {}
        self.info = {}
        for row in conn.execute(
            "SELECT id, name, shell_index, band_index, wedge_index, center_x_pc, center_y_pc, center_z_pc"
            " FROM sectors WHERE center_x_pc IS NOT NULL"
        ).fetchall():
            center = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
            self.info[row["id"]] = (row["name"], center)
            if row["shell_index"] is not None:
                self.by_address[(row["shell_index"], row["band_index"], row["wedge_index"])] = row["id"]

    def sector_at(self, point, current):
        """The generated sector holding `point`, or `current` when the
        cell it drifted into hasn't been generated (it stays filed where
        it was until that sector exists)."""
        if self.edge_pc <= 0:
            return current
        return self.by_address.get(sector_address_at(point, self.edge_pc), current)


def advance_galactic_positions(conn, clock, on_progress=None, budget=None):
    """
    Moves every due star system, standalone phenomenon, stand-alone
    facility and unbuilt scatter point (`advance_scatter_points`) along its galactic orbit, from its epoch to `clock.end_unix`:
    a system by its close pair's or primary star's period (its stars'
    `galactic_orbital_phase_deg` and a close pair's
    `binary_galactic_orbital_phase_deg` advance with it, each by its own
    period), a phenomenon by its own (its `galactic_orbital_phase_deg`
    with it), a facility by the rotation curve at its radius
    (`galactic_orbit.calculate_galactic_orbit`). Sectors are fixed cells,
    so an object that drifts into another generated sector is refiled
    there (`sector_id`, its sector-relative position, and later its
    octant, location text, containment and nearest systems -- see
    `refresh_after_motion`). A pure move keeps `modified_at`; a change of
    sector bumps it.

    Due (GEN.106): an object's `next_update_due` (indexed) is at or before
    the update's end, so it has moved at least 0.01 mpc on its galactic
    orbit since its epoch. A row without one (new or edited since the last
    update) gets one first, from its epoch and galactic orbital speed. An
    object that isn't due is neither moved nor counted. A system with no
    place in a generated sector still advances its phases on the same
    clock.

    Args:
        conn (Connection): An open, read-write connection.
        clock (OrbitClock): The span to move over (`orbit_clock`).
        on_progress (callable, optional): `on_progress(label, done, total)`
            after the star systems, each phenomenon table and the
            facilities are moved.
        budget (RowBudget, optional): The run's row ceiling (PERF.103); a
            table over its limit moves its oldest-due objects and leaves
            the rest due.

    Returns:
        dict: `moved` (objects whose galaxy position moved) and `refiled`
            counts, `sectors` -- the ids of every sector something moved
            in or out of -- and `counts`, `{name: rows_moved}`:
            `"star_systems"`, `"binary_galactic_orbits"`, each orbiting
            phenomenon table and `"space_facilities"`.
    """
    index = _SectorIndex(conn)
    moved = refiled = 0
    touched = set()
    counts = {}
    seconds_per_year = _yr_to_s(1.0)
    # The star systems, each phenomenon table (the quasars, which don't orbit, report nothing), the facilities, the
    # unbuilt scatter points.
    steps = len(PLACED_PHENOMENON_TABLES) + 1 + len(SCATTER_MOTION_TABLES)

    systems = conn.execute(
        f"""
        SELECT ss.id, ss.sector_id, {_SYSTEM_POINT.select("position", "velocity", "epoch", "due", alias="ss")},
               ss.binary_configuration, ss.runaway_class, ss.runaway_speed_kms,
               COALESCE(ss.binary_galactic_orbital_period_gy, s.galactic_orbital_period_gy) AS period_gy
        FROM star_systems ss
        LEFT JOIN stars s ON s.star_system_id = ss.id AND s.role IN ('single', 'primary')
        WHERE ss.next_update_due IS NULL OR ss.next_update_due <= ?
        ORDER BY ss.next_update_due, ss.id
        """,
        (clock.end_unix,),
    ).fetchall()
    scheduled, clocks, phases, updates, refile_updates = [], [], [], [], []
    close_pairs = 0
    for row in systems:
        ballistic = row["runaway_class"] == "hypervelocity" and row["runaway_speed_kms"]
        epoch, interval, due = _due_interval_s(
            row, clock, row["runaway_speed_kms"] if ballistic else _speed_kms(row), "galactic")
        if due > clock.end_unix or (budget is not None and not budget.take("star_systems")):
            scheduled.append((epoch, due, row["id"]))
            continue
        elapsed_years = (clock.end_unix - epoch) / seconds_per_year
        clocks.append((elapsed_years, clock.end_unix, clock.end_unix + interval, row["id"]))
        phases.append((elapsed_years, row["id"]))
        close_pairs += row["binary_configuration"] == "close"
        angle = _galactic_turn(elapsed_years, row["period_gy"])
        if ballistic:
            angle = None  # not on an orbit: a straight line (GEN.137)
        if angle == 0.0 or row["position_x_mpc"] is None or row["sector_id"] not in index.info:
            continue
        _name, center = index.info[row["sector_id"]]
        offset = (row["position_x_mpc"] / 1000.0, row["position_y_mpc"] / 1000.0, row["position_z_mpc"] / 1000.0)
        if ballistic:
            point = straight_line_position_pc(
                local_to_galaxy_pc(center, offset), (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"]),
                clock.end_unix - epoch)
            velocity = (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"])
        else:
            point = _rotate_about_axis(local_to_galaxy_pc(center, offset), angle)
            velocity = _rotate_about_axis((row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"]), angle)
        sector_id = index.sector_at(point, row["sector_id"])
        local_mpc = tuple(c * 1000.0 for c in galaxy_to_local_pc(index.info[sector_id][1], point))
        quadrant, _magnitudes = classify_octant(local_mpc)
        moved += 1
        if sector_id != row["sector_id"]:
            refiled += 1
            touched.update((sector_id, row["sector_id"]))
            refile_updates.append((sector_id, *local_mpc, *point, quadrant, *velocity, row["id"]))
        else:
            updates.append((*local_mpc, *point, quadrant, *velocity, row["id"]))
    if scheduled:
        conn.executemany("UPDATE star_systems SET epoch_unix = ?, next_update_due = ?, modified_at = modified_at"
                         " WHERE id = ?", scheduled)
    if clocks:
        conn.executemany(
            "UPDATE star_systems SET binary_galactic_orbital_phase_deg = CASE"
            " WHEN binary_configuration = 'close' AND binary_galactic_orbital_period_gy > 0"
            " THEN MOD(binary_galactic_orbital_phase_deg + 360 * MOD(?, binary_galactic_orbital_period_gy * 1e9)"
            " / (binary_galactic_orbital_period_gy * 1e9), 360)"
            " ELSE binary_galactic_orbital_phase_deg END,"
            " epoch_unix = ?, next_update_due = ?, modified_at = modified_at WHERE id = ?", clocks)
        conn.executemany(
            "UPDATE stars SET galactic_orbital_phase_deg ="
            " MOD(galactic_orbital_phase_deg + 360 * MOD(?, galactic_orbital_period_gy * 1e9)"
            " / (galactic_orbital_period_gy * 1e9), 360)"
            " WHERE star_system_id = ? AND galactic_orbital_period_gy > 0", phases)
    if updates:
        conn.executemany(
            f"UPDATE star_systems SET {_SYSTEM_POINT.assignments('position', 'galaxy_position')}, quadrant = ?,"
            f" {_SYSTEM_POINT.assignments('velocity')}, modified_at = modified_at WHERE id = ?",
            updates)
    if refile_updates:
        conn.executemany(
            f"UPDATE star_systems SET sector_id = ?, {_SYSTEM_POINT.assignments('position', 'galaxy_position')},"
            f" quadrant = ?, {_SYSTEM_POINT.assignments('velocity')} WHERE id = ?",
            refile_updates)
        _refile_nearest_rows(conn, "star_systems", [(update[0], update[-1]) for update in refile_updates])
    counts["star_systems"] = len(clocks)
    counts["binary_galactic_orbits"] = close_pairs
    _report(on_progress, "star_systems", 1, steps)

    for done, table in enumerate(PLACED_PHENOMENON_TABLES, start=2):
        if table == "quasars":
            continue  # the galaxy's nucleus sits at the center and doesn't orbit it
        # An anchored black hole or neutron star moves with its system's star row.
        anchored = " AND star_id IS NULL" if table in ("black_holes", "neutron_stars") else ""
        scheduled, clocks, updates, refile_updates = [], [], [], []
        for row in conn.execute(
            f"SELECT id, sector_id, center_x_pc, center_y_pc, center_z_pc, galactic_orbital_period_gy,"
            f" velocity_x_kms, velocity_y_kms, velocity_z_kms, epoch_unix, next_update_due FROM {table}"
            f" WHERE (next_update_due IS NULL OR next_update_due <= ?){anchored} ORDER BY next_update_due, id",
            (clock.end_unix,)
        ).fetchall():
            epoch, interval, due = _due_interval_s(row, clock, _speed_kms(row), "galactic")
            if due > clock.end_unix or (budget is not None and not budget.take(table)):
                scheduled.append((epoch, due, row["id"]))
                continue
            elapsed_years = (clock.end_unix - epoch) / seconds_per_year
            clocks.append((elapsed_years, clock.end_unix, clock.end_unix + interval, row["id"]))
            # A black hole or neutron star with a velocity and no galactic orbit is one the galaxy ejected (GEN.238):
            # a straight line, like a hypervelocity star's (GEN.137). The galaxy's own black hole has no velocity.
            ballistic = (table in ("black_holes", "neutron_stars") and not row["galactic_orbital_period_gy"]
                         and bool(_speed_kms(row)))
            angle = None if ballistic else _galactic_turn(elapsed_years, row["galactic_orbital_period_gy"])
            if angle == 0.0 or row["center_x_pc"] is None or (elapsed_years <= 0 and ballistic):
                continue
            velocity = (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"])
            if ballistic:
                point = straight_line_position_pc((row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]), velocity,
                                                  clock.end_unix - epoch)
            else:
                point = _rotate_about_axis((row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]), angle)
                velocity = _rotate_about_axis(velocity, angle)
            radius = math.sqrt(sum(c * c for c in point))
            moved += 1
            sector_id = index.sector_at(point, row["sector_id"]) if row["sector_id"] is not None else None
            if sector_id != row["sector_id"]:
                refiled += 1
                touched.update((sector_id, row["sector_id"]))
                refile_updates.append((sector_id, *point, radius, *velocity, row["id"]))
            else:
                updates.append((*point, radius, *velocity, row["id"]))
        if scheduled:
            conn.executemany(f"UPDATE {table} SET epoch_unix = ?, next_update_due = ?, modified_at = modified_at"
                             f" WHERE id = ?", scheduled)
        if clocks:
            conn.executemany(
                f"UPDATE {table} SET galactic_orbital_phase_deg = CASE WHEN galactic_orbital_period_gy > 0"
                f" THEN MOD(galactic_orbital_phase_deg + 360 * MOD(?, galactic_orbital_period_gy * 1e9)"
                f" / (galactic_orbital_period_gy * 1e9), 360)"
                f" ELSE galactic_orbital_phase_deg END, epoch_unix = ?, next_update_due = ?,"
                f" modified_at = modified_at WHERE id = ?", clocks)
        if updates:
            conn.executemany(
                f"UPDATE {table} SET center_x_pc = ?, center_y_pc = ?, center_z_pc = ?, galactic_radius_pc = ?,"
                f" velocity_x_kms = ?, velocity_y_kms = ?, velocity_z_kms = ?,"
                f" modified_at = modified_at WHERE id = ?", updates)
        if refile_updates:
            conn.executemany(
                f"UPDATE {table} SET sector_id = ?, center_x_pc = ?, center_y_pc = ?, center_z_pc = ?,"
                f" galactic_radius_pc = ?, velocity_x_kms = ?, velocity_y_kms = ?, velocity_z_kms = ? WHERE id = ?",
                refile_updates)
            _refile_nearest_rows(conn, table, [(update[0], update[-1]) for update in refile_updates])
        counts[table] = len(clocks)
        _report(on_progress, table, done, steps)

    scheduled, updates, refile_updates = [], [], []
    for row in conn.execute(
        "SELECT id, sector_id, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc,"
        " velocity_x_kms, velocity_y_kms, velocity_z_kms, epoch_unix, next_update_due FROM facilities"
        " WHERE host_type = 'space' AND center_x_pc IS NOT NULL"
        " AND (next_update_due IS NULL OR next_update_due <= ?) ORDER BY next_update_due, id", (clock.end_unix,)
    ).fetchall():
        speed_kms, period_gy = calculate_galactic_orbit(pc_to_ly(detmath.hypot(row["center_x_pc"], row["center_y_pc"])))
        epoch, interval, due = _due_interval_s(row, clock, speed_kms, "galactic")
        if due > clock.end_unix or (budget is not None and not budget.take("space_facilities")):
            scheduled.append((epoch, due, row["id"]))
            continue
        angle = _galactic_turn((clock.end_unix - epoch) / seconds_per_year, period_gy)
        next_due = clock.end_unix + interval
        if angle == 0.0:
            scheduled.append((clock.end_unix, next_due, row["id"]))
            continue
        point = _rotate_about_axis((row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]), angle)
        velocity = _rotate_about_axis((row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"]), angle)
        moved += 1
        sector_id = index.sector_at(point, row["sector_id"])
        if sector_id != row["sector_id"]:
            refiled += 1
            touched.update((sector_id, row["sector_id"]))
            refile_updates.append((sector_id, *point, row["galactic_radius_pc"], *velocity,
                                   clock.end_unix, next_due, row["id"]))
        else:
            updates.append((*point, *velocity, clock.end_unix, next_due, row["id"]))
    if scheduled:
        conn.executemany("UPDATE facilities SET epoch_unix = ?, next_update_due = ?, modified_at = modified_at"
                         " WHERE id = ?", scheduled)
    if updates:
        conn.executemany("UPDATE facilities SET center_x_pc = ?, center_y_pc = ?, center_z_pc = ?,"
                         " velocity_x_kms = ?, velocity_y_kms = ?, velocity_z_kms = ?,"
                         " epoch_unix = ?, next_update_due = ?, modified_at = modified_at WHERE id = ?", updates)
    if refile_updates:
        conn.executemany("UPDATE facilities SET sector_id = ?, center_x_pc = ?, center_y_pc = ?, center_z_pc = ?,"
                         " galactic_radius_pc = ?, velocity_x_kms = ?, velocity_y_kms = ?, velocity_z_kms = ?,"
                         " epoch_unix = ?, next_update_due = ? WHERE id = ?", refile_updates)
    counts["space_facilities"] = len(updates) + len(refile_updates)
    _report(on_progress, "facilities", steps - len(SCATTER_MOTION_TABLES), steps)

    points = advance_scatter_points(conn, clock, on_progress, steps=(steps - len(SCATTER_MOTION_TABLES), steps),
                                    budget=budget)
    counts.update(points["counts"])
    moved += points["moved"]
    refiled += points["refiled"]

    touched.discard(None)
    for sector_id in touched:
        touch_sector(conn, sector_id)
    return {"moved": moved, "refiled": refiled, "sectors": touched, "counts": counts}


SCATTER_MOTION_TABLES = ("bright_stars", "phenomenon_scatter")
"""tuple: The tables of objects not yet built (GEN.223): the star and phenomenon scatters. Each row is a point in
space of its own, moved with everything else; once the sector builds it, the built object carries on."""

_SCATTER_ROWS = {
    # table -> the rows that are still only a point: not built, and not the nucleus (which sits at the center)
    "bright_stars": "star_system_id IS NULL",
    "phenomenon_scatter": "built_at IS NULL AND kind <> 'quasar'",
}
_SCATTER_COLUMNS = {
    "bright_stars": "id, NULL AS kind",
    "phenomenon_scatter": "id, kind",
}


def advance_scatter_points(conn, clock, on_progress=None, steps=(1, 1), budget=None):
    """
    Moves every due scatter point (`SCATTER_MOTION_TABLES`) along its galactic
    orbit from its epoch to `clock.end_unix`, the way `advance_galactic_positions`
    moves a built star or phenomenon (GEN.223): a point on the rotation curve
    turns about the galactic axis by the angle its period sweeps, a
    hypervelocity star flies a straight line. A point's own columns hold all of
    it, one per component: position (milliparsecs, galaxy frame), velocity
    (km/s), `epoch_unix`, `next_update_due`, and the mass it was drawn with. A
    row with no velocity yet is on the rotation curve (the update gives it the
    circular velocity at its radius), and one with no due time gets one first.
    A point that crosses into another cell is refiled there (`shell_index`,
    `band_index`, `wedge_index`), so the sector that is generated there later
    builds what is in it now.

    Args:
        conn (Connection): An open, read-write connection.
        clock (OrbitClock): The span to move over (`orbit_clock`).
        on_progress (callable, optional): `on_progress(label, done, total)`
            after each table.
        steps (tuple): `(done, total)` of the steps `on_progress` already
            counts, so the caller's bar continues.
        budget (RowBudget, optional): The run's row ceiling (PERF.103); a
            table over its limit moves its oldest-due points.

    Returns:
        dict: `moved` and `refiled` counts, and `counts`, `{table: rows_moved}`.
    """
    shape = get_galaxy_shape(conn)
    edge_pc = shape.edge_pc if shape is not None else 0.0
    moved = refiled = 0
    counts = {}
    seconds_per_year = _yr_to_s(1.0)
    done, total = steps
    for table in SCATTER_MOTION_TABLES:
        scheduled, updates = [], []
        rows = conn.execute(
            f"SELECT {_SCATTER_COLUMNS[table]}, shell_index, band_index, wedge_index,"
            f" {point_columns.LAYOUTS[table].select('position', 'velocity', 'epoch', 'due')}"
            f" FROM {table} WHERE {_SCATTER_ROWS[table]}"
            f" AND (next_update_due IS NULL OR next_update_due <= ?) ORDER BY next_update_due, id",
            (clock.end_unix,)).fetchall()
        for row in rows:
            point = (row["position_x_mpc"] / 1000.0, row["position_y_mpc"] / 1000.0, row["position_z_mpc"] / 1000.0)
            speed_kms, period_gy = calculate_galactic_orbit(pc_to_ly(detmath.hypot(point[0], point[1])))
            stored = None not in (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"])
            velocity = ((row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"]) if stored
                        else tuple(c / 1000.0 for c in galactic_velocity_ms(point, speed_kms)))
            ballistic = row["kind"] == "hypervelocity-star" and stored
            speed = math.sqrt(sum(c * c for c in velocity)) if ballistic else speed_kms
            epoch, interval, due = _due_interval_s(row, clock, speed, "galactic")
            if due > clock.end_unix or (budget is not None and not budget.take(table)):
                scheduled.append((epoch, due, row["id"]))
                continue
            next_due = clock.end_unix + interval
            elapsed_years = (clock.end_unix - epoch) / seconds_per_year
            if ballistic:
                new_point = straight_line_position_pc(point, velocity, clock.end_unix - epoch)
                new_velocity = velocity
            else:
                angle = _galactic_turn(elapsed_years, period_gy)
                new_point = _rotate_about_axis(point, angle)
                new_velocity = _rotate_about_axis(velocity, angle)
            address = (row["shell_index"], row["band_index"], row["wedge_index"])
            if edge_pc > 0:
                found = sector_address_at(new_point, edge_pc)
                if found[0] <= MAX_ADDRESSABLE_SHELL:
                    address = found
            moved += 1
            refiled += address != (row["shell_index"], row["band_index"], row["wedge_index"])
            updates.append((*address, *(c * 1000.0 for c in new_point), *new_velocity, clock.end_unix, next_due,
                            row["id"]))
        if scheduled:
            conn.executemany(f"UPDATE {table} SET epoch_unix = ?, next_update_due = ? WHERE id = ?", scheduled)
        if updates:
            conn.executemany(
                f"UPDATE {table} SET shell_index = ?, band_index = ?, wedge_index = ?,"
                f" {point_columns.LAYOUTS[table].assignments('position', 'velocity', 'epoch', 'due')} WHERE id = ?",
                updates)
        counts[table] = len(updates)
        done += 1
        _report(on_progress, table, done, total)
    return {"moved": moved, "refiled": refiled, "counts": counts}


def _refile_nearest_rows(conn, table, moves):
    """Files an object's stored `nearest_systems` rows under its new
    sector (`moves` is `(sector_id, object_id)` pairs); an object that
    left every sector loses them, since the table needs a sector."""
    kept = [(sector_id, table, object_id) for sector_id, object_id in moves if sector_id is not None]
    dropped = [(table, object_id) for sector_id, object_id in moves if sector_id is None]
    if kept:
        conn.executemany("UPDATE nearest_systems SET sector_id = ? WHERE object_table = ? AND object_id = ?", kept)
    if dropped:
        conn.executemany("DELETE FROM nearest_systems WHERE object_table = ? AND object_id = ?", dropped)


def advance_facility_orbits(conn, clock, budget=None):
    """Advances every due orbital facility's `orbit_phase_deg` by its period
    from its epoch to `clock.end_unix`, the way `advance_orbital_phases`
    does for moons, and every asteroid facility in a belt around its star
    the same way. Returns the count."""
    where = "placement IN ('orbital', 'asteroid') AND orbit_period_years > 0"
    _schedule_updates(conn, clock, "facilities", _FACILITY_INTERVAL_SQL, where, keep_modified=True)
    limited, tail = _ordered_limit(budget, "facilities")
    moved = conn.execute(
        f"UPDATE facilities SET orbit_phase_deg = MOD(orbit_phase_deg + {_phase_turn_sql('orbit_period_years')},"
        f" 360.0), epoch_unix = ?, next_update_due = ? + {_FACILITY_INTERVAL_SQL},"
        f" modified_at = modified_at WHERE {where} AND next_update_due <= ?{tail}",
        (clock.end_unix, clock.end_unix, clock.end_unix, clock.end_unix),
    ).rowcount
    if limited is not None:
        budget.defer_remaining(conn, clock, "facilities")
    return moved


def refresh_after_motion(conn, refiled_sectors=(), on_progress=None):
    """
    Brings everything that depends on where things are up to date after
    `advance_galactic_positions`: containment (`refresh_containment`),
    each phenomenon's octant and every object's nearest systems
    (`refresh_nearest_systems`) for every placed sector, then the
    `star_systems.location` text of every system whose nearest systems
    changed or that sits in a sector something moved in or out of, built
    from its sector's name and the stored nearest systems.

    Args:
        conn (Connection): An open, read-write connection.
        refiled_sectors (iterable): The sectors something moved in or out of.
        on_progress (callable, optional): `on_progress(label, done, total)`
            as containment (`"Containment: sectors"`), nearest systems
            (`"Nearest systems: sectors"`) and locations (`"Locations:
            systems"`) are refreshed.

    Returns:
        dict: `"locations"`, how many location texts were rewritten, and
            how many objects entered or left a nebula or supernova remnant
            (`refresh_containment`'s counts, GEN.107).
    """
    placed = [row["id"] for row in conn.execute("SELECT id FROM sectors WHERE center_x_pc IS NOT NULL").fetchall()]
    containment = refresh_containment(
        conn, placed, lambda _label, done, total: _report(on_progress, "Containment: sectors", done, total))
    changed = refresh_nearest_systems(
        conn, placed, lambda _label, done, total: _report(on_progress, "Nearest systems: sectors", done, total))
    system_ids = {object_id for table, object_id in changed if table == "star_systems"}
    refiled_sectors = sorted(set(refiled_sectors) - {None})
    for start in range(0, len(refiled_sectors), 500):
        batch = refiled_sectors[start:start + 500]
        marks = ", ".join("?" * len(batch))
        system_ids.update(row["id"] for row in conn.execute(
            f"SELECT id FROM star_systems WHERE sector_id IN ({marks})", tuple(batch)).fetchall())
    system_ids = sorted(system_ids)
    rewritten = 0
    for start in range(0, len(system_ids), 500):
        batch = system_ids[start:start + 500]
        marks = ", ".join("?" * len(batch))
        names = {row["id"]: row["name"] for row in conn.execute(
            f"SELECT ss.id, sec.name FROM star_systems ss JOIN sectors sec ON sec.id = ss.sector_id"
            f" WHERE ss.id IN ({marks})", tuple(batch)).fetchall()}
        neighbors = {}
        for row in conn.execute(
            f"SELECT n.object_id, n.distance_pc, ss.name FROM nearest_systems n"
            f" JOIN star_systems ss ON ss.id = n.neighbor_system_id"
            f" WHERE n.object_table = 'star_systems' AND n.object_id IN ({marks})"
            f" ORDER BY n.object_id, n.neighbor_rank", tuple(batch)).fetchall():
            neighbors.setdefault(row["object_id"], []).append((row["name"], pc_to_ly(row["distance_pc"])))
        updates = [(_format_location_string(names[system_id], neighbors.get(system_id, [])), system_id)
                   for system_id in batch if system_id in names]
        if updates:
            conn.executemany("UPDATE star_systems SET location = ?, modified_at = modified_at WHERE id = ?", updates)
            rewritten += len(updates)
        _report(on_progress, "Locations: systems", min(start + 500, len(system_ids)), len(system_ids))
    return {"locations": rewritten, **containment}


# ---------------------------------------------------------------------------
# Bright-star pre-placement (schema v43) -- storage. The scatter itself
# (planetgen plan) and the fill draw stars from Physics' sampling API.
# ---------------------------------------------------------------------------

BRIGHT_STAR_COLUMNS = (
    "shell_index", "band_index", "wedge_index", *point_columns.LAYOUTS["bright_stars"].position,
    "population", "star_type", "yerkes_class", "mass_kg", "radius_km", "temperature_k", "luminosity_w",
    "age_gy", "lifespan_gy", "initial_mass_sol", "phase_end_age_gy", "seed", "metallicity_fe_h",
    "rotation_period_hours",
)
"""tuple: The `bright_stars` columns a scatter writes, in the order
`insert_bright_stars` expects each row's values."""


def clear_bright_stars(conn):
    """
    Empties `bright_stars`, sets every sector's own bright-star level
    (`sector_stats`, v53, GEN.44) back to untouched (-1; a filled sector
    stays at 0, and the level a delete would put back is forgotten too)
    and forgets the scatter's threshold and seed -- a plan re-run or a new
    galaxy starts over. `TRUNCATE` (an implicit commit), since a real
    scatter leaves tens of millions of rows.
    """
    conn.execute("TRUNCATE TABLE bright_stars")
    scatter_systems.delete_all_member_trees(conn, "bright_stars")
    conn.execute("UPDATE sector_stats SET bright_level_sol = -1 WHERE bright_level_sol > 0")
    conn.execute("UPDATE sector_stats SET level_before_fill_sol = -1 WHERE level_before_fill_sol > 0")
    conn.execute("UPDATE sector_stats SET bright_mass_sol = NULL WHERE bright_mass_sol IS NOT NULL")
    conn.execute("UPDATE galaxy_shape SET bright_star_min_luminosity_sol = NULL, bright_star_seed = NULL,"
                 " bright_star_mass_limit_sol = NULL")
    forget_scatter_epoch(conn, ("mass", "luminosity", "band"))
    conn.commit()


def delete_unbuilt_bright_stars(conn, address):
    """Deletes the stars of one cell that no system was built around, with their members: what a run that failed
    before it recorded a level left behind, drawn again whole. `address` is `(shell, band, wedge)`."""
    where = "shell_index = ? AND band_index = ? AND wedge_index = ? AND star_system_id IS NULL"
    scatter_systems.delete_member_trees_where(
        conn, "bright_stars", "t.shell_index = ? AND t.band_index = ? AND t.wedge_index = ? AND t.star_system_id IS NULL",
        tuple(address))
    conn.execute(f"DELETE FROM bright_stars WHERE {where}", tuple(address))


def _cell_tiers(limits, cells):
    """`{(shell, band, wedge): tier index}` of each of `cells`, one tier lookup per shell, of the bands the cells are in
    (`ScatterLimits.tiers_of`)."""
    by_shell = {}
    for cell in cells:
        by_shell.setdefault(cell[0], []).append(cell)
    found = {}
    for shell_index, group in by_shell.items():
        bands = sorted({cell[1] for cell in group})
        tiers = limits.tiers_of(shell_index, bands)
        place = {band: index for index, band in enumerate(bands)}
        for cell in group:
            found[cell] = int(tiers[place[cell[1]]])
    return found


def _delete_bright_stars_by_id(conn, ids, batch_size=500):
    """Deletes the bright stars with these ids and their member trees."""
    for start in range(0, len(ids), batch_size):
        chunk = ids[start:start + batch_size]
        marks = ", ".join("?" for _ in chunk)
        scatter_systems.delete_member_trees_where(conn, "bright_stars", f"t.id IN ({marks})", tuple(chunk))
        conn.execute(f"DELETE FROM bright_stars WHERE id IN ({marks})", tuple(chunk))


NO_MASS_LIMIT_SOL = 1.0e300
"""float: A mass limit above any star, for a clear that finds none recorded."""


def clear_bright_star_pass(conn, part, limits):
    """
    Deletes one pass of the star scatter and leaves the other (GEN.196):
    `part` "mass" removes the stars born at their band's mass limit or more (the mass
    pass's), "luminosity" the lighter ones (the luminosity pass's). The
    sectors' own bright-star levels go back to untouched, as
    `clear_bright_stars` does, and what the removed pass recorded is
    forgotten: the mass limit for "mass", the floor and seed for "luminosity".

    `limits` (`galaxy.scatter_limits.ScatterLimits`, GEN.333) gives each band's mass limit, which the core's
    differs from the baseline's: a star at or above the highest limit is the mass pass's wherever it is, one under
    the lowest the luminosity pass's, and the few between are told apart by their band's tier.
    """
    if part not in ("mass", "luminosity"):
        raise ValueError(f"part must be 'mass' or 'luminosity', not {part!r}")
    if limits.tiers is None:
        # Only the mass pass ran, or only stars with no recorded split: the one limit stored (none puts every star
        # under it, so the luminosity pass clears them all and the mass pass nothing).
        stored = bright_star_mass_limit(conn)
        masses = [NO_MASS_LIMIT_SOL if stored is None else stored]
    else:
        masses = [found.mass_limit_sol for found in limits.tiers]
    low, high = min(masses), max(masses)
    if part == "mass":
        scatter_systems.delete_member_trees_where(conn, "bright_stars", "t.initial_mass_sol >= ?", (high,))
        conn.execute("DELETE FROM bright_stars WHERE initial_mass_sol >= ?", (high,))
        forget = "bright_star_mass_limit_sol = NULL"
    else:
        scatter_systems.delete_member_trees_where(conn, "bright_stars", "t.initial_mass_sol < ?", (low,))
        conn.execute("DELETE FROM bright_stars WHERE initial_mass_sol < ?", (low,))
        forget = "bright_star_min_luminosity_sol = NULL, bright_star_seed = NULL"
    if low < high:
        between = conn.execute(
            "SELECT id, shell_index, band_index, wedge_index, initial_mass_sol FROM bright_stars"
            " WHERE initial_mass_sol >= ? AND initial_mass_sol < ?", (low, high)).fetchall()
        tiers = _cell_tiers(limits, {(row["shell_index"], row["band_index"], row["wedge_index"]) for row in between})
        _delete_bright_stars_by_id(conn, [
            row["id"] for row in between
            if (row["initial_mass_sol"] >= masses[tiers[(row["shell_index"], row["band_index"], row["wedge_index"])]])
            == (part == "mass")])
    conn.execute("UPDATE sector_stats SET bright_level_sol = -1 WHERE bright_level_sol > 0")
    conn.execute("UPDATE sector_stats SET level_before_fill_sol = -1 WHERE level_before_fill_sol > 0")
    conn.execute("UPDATE sector_stats SET bright_mass_sol = NULL WHERE bright_mass_sol IS NOT NULL")
    conn.execute(f"UPDATE galaxy_shape SET {forget}")
    forget_scatter_epoch(conn, ("mass",) if part == "mass" else ("luminosity", "band"))
    conn.commit()


def scatter_stored(conn, part):
    """
    Whether the scatter of `part` has rows stored: "mass" (the massive-star pass, which records the mass limit it placed
    from), "luminosity" (stars lighter than that limit; with no mass pass recorded, any star at all, since nothing then
    says where to split them), "phenomena" or "nebula-field". `planetgen scatter` draws only where none is stored and
    `planetgen clear-scatter` empties it.
    """
    if part == "mass":
        return bright_star_mass_limit(conn) is not None
    if part == "luminosity":
        limit = bright_star_mass_limit(conn)
        if limit is None:
            return conn.execute("SELECT 1 FROM bright_stars LIMIT 1").fetchone() is not None
        low = min([found.mass_limit_sol for found in scatter_limits_of(conn).tiers or ()
                   if found.mass_limit_sol is not None] or [limit])
        return conn.execute("SELECT 1 FROM bright_stars WHERE initial_mass_sol < ? LIMIT 1", (low,)).fetchone() is not None
    if part == "phenomena":
        return conn.execute("SELECT 1 FROM phenomenon_scatter LIMIT 1").fetchone() is not None
    if part == "nebula-field":
        return conn.execute("SELECT 1 FROM nebula_field LIMIT 1").fetchone() is not None
    raise ValueError(f"part must be 'mass', 'luminosity', 'phenomena' or 'nebula-field', not {part!r}")


SCATTER_KINDS = ("mass", "luminosity", "band", "phenomena", "nebula-field")
"""tuple: The kinds of scatter that have a row in `scatter_epochs` (GEN.764)."""


def record_scatter_epoch(conn, kinds, run_id=None):
    """
    Records that the scatter of each of `kinds` (`SCATTER_KINDS`) was drawn under the running generator epoch, by the
    run `run_id` (None when it is not known), replacing what an earlier draw of the kind recorded. Not committed.
    """
    kinds = (kinds,) if isinstance(kinds, str) else tuple(kinds)
    epoch = generatorEpoch.generator_epoch()
    for kind in kinds:
        if kind not in SCATTER_KINDS:
            raise ValueError(f"kind must be one of {', '.join(SCATTER_KINDS)}, not {kind!r}")
        conn.execute(
            "INSERT INTO scatter_epochs (kind, generator_epoch, run_id) VALUES (?, ?, ?)"
            " ON DUPLICATE KEY UPDATE generator_epoch = ?, run_id = ?, recorded_at = CURRENT_TIMESTAMP(3)",
            (kind, epoch, run_id, epoch, run_id))


def forget_scatter_epoch(conn, kinds):
    """Removes the record of each of `kinds` (a clear of that scatter). Not committed."""
    kinds = (kinds,) if isinstance(kinds, str) else tuple(kinds)
    if kinds:
        conn.execute(f"DELETE FROM scatter_epochs WHERE kind IN ({', '.join('?' * len(kinds))})", kinds)


def scatter_epochs(conn):
    """`{kind: generator epoch or None}` of the scatters recorded, in the order of `SCATTER_KINDS`."""
    found = {row["kind"]: row["generator_epoch"] for row in conn.execute(
        "SELECT kind, generator_epoch FROM scatter_epochs").fetchall()}
    return {kind: found[kind] for kind in SCATTER_KINDS if kind in found}


def record_bright_star_scatter(conn, min_luminosity_sol, seed, mass_limit_sol=None):
    """Stores the threshold and seed a finished scatter used, so a fill
    reads them rather than today's constant, and the mass limit its mass
    pass placed every star from (GEN.185; `None` keeps the one already
    stored, as a staged band below the floor does)."""
    if mass_limit_sol is None:
        conn.execute("UPDATE galaxy_shape SET bright_star_min_luminosity_sol = ?, bright_star_seed = ? WHERE id = 1",
                     (min_luminosity_sol, seed))
        return
    conn.execute("UPDATE galaxy_shape SET bright_star_min_luminosity_sol = ?, bright_star_seed = ?,"
                 " bright_star_mass_limit_sol = ? WHERE id = 1", (min_luminosity_sol, seed, mass_limit_sol))


def bright_star_mass_limit(conn):
    """The mass limit (solar masses) the galaxy's star scatter placed every
    star from, or `None` when it had no mass pass (or none has run)."""
    row = conn.execute("SELECT bright_star_mass_limit_sol FROM galaxy_shape WHERE id = 1").fetchone()
    if row is None or row["bright_star_mass_limit_sol"] is None:
        return None
    return float(row["bright_star_mass_limit_sol"])


def scatter_limits_of(conn, skeleton=None):
    """
    The galaxy's scatter limits by place (`galaxy.scatter_limits.ScatterLimits`, GEN.333): the baseline the scatter
    recorded (`bright_star_scatter_settings`, `bright_star_mass_limit`) and the core and halo settings of
    `galaxy_shape`. Before any scatter ran (or with no galaxy) no tier has a limit.
    """
    skeleton = get_galaxy_shape(conn) if skeleton is None else skeleton
    if skeleton is None:
        return scatterLimits.ScatterLimits(None, None, None)
    settings = bright_star_scatter_settings(conn)
    return scatterLimits.ScatterLimits.of(skeleton, None if settings is None else settings[0],
                                          bright_star_mass_limit(conn))


def bright_star_scatter_settings(conn):
    """`(min_luminosity_sol, seed)` of the galaxy's scatter, or `None`
    when none has run (a fill then behaves exactly as before)."""
    row = conn.execute(
        "SELECT bright_star_min_luminosity_sol, bright_star_seed FROM galaxy_shape WHERE id = 1").fetchone()
    if row is None or row["bright_star_min_luminosity_sol"] is None:
        return None
    return row["bright_star_min_luminosity_sol"], row["bright_star_seed"]


def bright_star_marked_addresses(conn, limits):
    """
    The sectors the star scatter's mass pass marked (GEN.185): the
    `(shell, band, wedge)` addresses holding a star born at their band's mass limit
    or more that shines at their band's luminosity floor or more (`limits`, a
    `galaxy.scatter_limits.ScatterLimits`: the core, baseline and halo tiers each have their own, GEN.333).
    The luminosity pass skips them, since they already hold a star that
    bright. Read between the passes, when the table holds only the mass
    pass's stars.
    """
    tiers = [(found.mass_limit_sol, found.min_luminosity_sol * physical_constants.SOLAR_LUMINOSITY)
             for found in limits.tiers]
    if any(mass is None for mass, _watts in tiers):
        return set()
    flags = ", ".join(f"MAX(initial_mass_sol >= ? AND luminosity_w >= ?) AS t{index}" for index in range(3))
    values = [value for item in tiers for value in item]
    rows = conn.execute(
        f"SELECT shell_index, band_index, wedge_index, {flags} FROM bright_stars"
        " WHERE initial_mass_sol >= ? AND luminosity_w >= ?"
        " GROUP BY shell_index, band_index, wedge_index",
        (*values, min(mass for mass, _watts in tiers), min(watts for _mass, watts in tiers))).fetchall()
    cells = [(row["shell_index"], row["band_index"], row["wedge_index"]) for row in rows]
    found = _cell_tiers(limits, cells)
    return {cell for cell, row in zip(cells, rows) if row[f"t{found[cell]}"]}


def stamp_phenomenon_scatter_epoch(conn):
    """Records the database's orbit epoch (`get_orbit_epoch_unix`) as the time
    the just-drawn hypervelocity stars' positions hold at (GEN.137); NULL
    when no orbit update has run yet."""
    conn.execute("UPDATE phenomenon_scatter SET epoch_unix = ? WHERE kind = 'hypervelocity-star'",
                 (get_orbit_epoch_unix(conn),))


def estimated_rows(conn, table):
    """About how many rows `table` holds, from the server's own estimate (no
    count: that scans a table of millions). Used to predict the time of a
    statement that can't report progress (UX.84)."""
    row = conn.execute(
        "SELECT TABLE_ROWS AS n FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?",
        (table,)).fetchone()
    return int(row["n"] or 0) if row else 0


def clear_phenomenon_scatter(conn):
    """Empties `phenomenon_scatter` and forgets the scatter's seed -- a
    plan re-run or a new galaxy starts over. `TRUNCATE` (an implicit
    commit), since a real scatter leaves over a hundred million rows."""
    conn.execute("TRUNCATE TABLE phenomenon_scatter")
    scatter_systems.delete_all_member_trees(conn, "phenomenon_scatter")
    conn.execute("TRUNCATE TABLE phenomenon_scatter_classes")
    conn.execute("UPDATE galaxy_shape SET phenomenon_scatter_seed = NULL, phenomenon_min_mass_solar = NULL")
    forget_scatter_epoch(conn, "phenomena")
    conn.commit()


NEBULA_FIELD_COLUMNS = ("cell_i", "cell_j", "cell_k", "n", "nebula_class", "radius_ly", "center_x_pc", "center_y_pc",
                        "center_z_pc")
"""tuple: The columns `insert_nebula_field` writes (`nebulaField.table_rows` order)."""


def clear_nebula_field(conn):
    """Empties `nebula_field` (DB.24): a plan or a redo of the phenomena draws it again from the seed."""
    conn.execute("TRUNCATE TABLE nebula_field")
    forget_scatter_epoch(conn, "nebula-field")
    conn.commit()


def insert_nebula_field(conn, rows, batch_size=5000):
    """Writes `nebula_field` rows (`NEBULA_FIELD_COLUMNS` order), `batch_size` to a statement. Returns the rows
    written."""
    head = f"INSERT INTO nebula_field ({', '.join(NEBULA_FIELD_COLUMNS)}) VALUES "
    written = 0
    batch = []
    for row in rows:
        batch.append("(" + ",".join(map(_scatter_value, row)) + ")")
        if len(batch) >= batch_size:
            conn.execute_literal(head + ",".join(batch))
            written += len(batch)
            batch = []
    if batch:
        conn.execute_literal(head + ",".join(batch))
        written += len(batch)
    return written


def link_nebula_field_stored(conn):
    """Marks the `nebula_field` rows whose cloud is already a `nebulae` row (found by its centre, as
    `_insert_field_nebulae` finds one) as stored, with that nebula's object ID: after a plan redraws the table
    over a galaxy that already holds some of its clouds. Returns the rows linked."""
    cursor = conn.execute(
        "UPDATE nebula_field f JOIN nebulae n ON n.center_x_pc = f.center_x_pc AND n.center_y_pc = f.center_y_pc"
        " AND n.center_z_pc = f.center_z_pc SET f.is_stored = 1, f.uid = n.uid")
    conn.commit()
    return cursor.rowcount


def record_phenomenon_scatter(conn, seed, min_mass_solar):
    """Stores the seed and the mass cut (GEN.167) a finished phenomenon
    scatter used: their presence tells a sector fill to build its objects
    from `phenomenon_scatter` instead of rolling its own, and to draw the
    neutron stars and black holes below the cut itself (GEN.168)."""
    conn.execute("UPDATE galaxy_shape SET phenomenon_scatter_seed = ?, phenomenon_min_mass_solar = ? WHERE id = 1",
                 (seed, min_mass_solar))


def record_phenomenon_scatter_classes(conn, counts):
    """Stores how many phenomena of each class the scatter placed
    (`{(kind, subtype): count}`, `subtype` `""` for none), none built yet, so
    the Phenomena table counts them without scanning `phenomenon_scatter`."""
    conn.execute("DELETE FROM phenomenon_scatter_classes")
    conn.executemany(
        "INSERT INTO phenomenon_scatter_classes (kind, subtype, placed, built) VALUES (?, ?, ?, 0)",
        [(kind, subtype, int(count)) for (kind, subtype), count in sorted(counts.items()) if count])


def add_phenomenon_scatter_classes(conn, counts):
    """Adds to the classes `record_phenomenon_scatter_classes` stored (`{(kind, subtype): count}`): the special phenomena
    a scatter writes once its shells are done."""
    conn.executemany(
        "INSERT INTO phenomenon_scatter_classes (kind, subtype, placed, built) VALUES (?, ?, ?, 0) AS added"
        " ON DUPLICATE KEY UPDATE placed = phenomenon_scatter_classes.placed + added.placed",
        [(kind, subtype, int(count)) for (kind, subtype), count in sorted(counts.items()) if count])


def phenomenon_scatter_classes(conn):
    """`[{"kind", "subtype", "unbuilt"}]`: each scattered class with phenomena no sector has built yet."""
    return [{"kind": row["kind"], "subtype": row["subtype"], "unbuilt": int(row["unbuilt"])}
            for row in conn.execute(
                "SELECT kind, subtype, placed - built AS unbuilt FROM phenomenon_scatter_classes"
                " WHERE placed > built ORDER BY kind, subtype").fetchall()]


def phenomenon_scatter_settings(conn):
    """`(seed, min_mass_solar)` of the galaxy's phenomenon scatter, or
    `None` when none has run."""
    row = conn.execute(
        "SELECT phenomenon_scatter_seed, phenomenon_min_mass_solar FROM galaxy_shape WHERE id = 1").fetchone()
    if row is None or row["phenomenon_scatter_seed"] is None:
        return None
    cut = row["phenomenon_min_mass_solar"]
    return int(row["phenomenon_scatter_seed"]), (None if cut is None else float(cut))


def phenomenon_scatter_seed(conn):
    """The seed of the galaxy's phenomenon scatter, or `None` when none
    has run (a fill then rolls its own black holes, neutron stars,
    planetary nebulae and supernova remnants)."""
    row = conn.execute("SELECT phenomenon_scatter_seed FROM galaxy_shape WHERE id = 1").fetchone()
    if row is None or row["phenomenon_scatter_seed"] is None:
        return None
    return int(row["phenomenon_scatter_seed"])


PHENOMENON_SCATTER_COLUMNS = (
    "shell_index", "band_index", "wedge_index", "kind", "subtype",
    *point_columns.LAYOUTS["phenomenon_scatter"].named("position", "velocity"), "seed", "mass_solar",
    "map_luminosity_sol",
)
"""tuple: The `phenomenon_scatter` columns a scatter writes, in the order
`insert_phenomenon_scatter` expects each row's values."""


PHENOMENON_SCATTER_INDEXES = {
    "idx_phenomenon_scatter_address": "shell_index, band_index, wedge_index",
    "idx_phenomenon_scatter_class": "kind, subtype, mass_solar",
    "idx_phenomenon_scatter_cell": "cell_key, position_x_mpc, position_y_mpc, position_z_mpc",
    "idx_phenomenon_scatter_next_update_due": "next_update_due",
}
"""dict: The secondary indexes of `phenomenon_scatter` (as `schema.sql` has them) that a scatter drops while it
writes and builds again at the end, one at a time, each sorted (PERF.73, PERF.101). The cell index and the due-time
index are in it since PERF.167: a row lands at a random place in each, so once they outgrow the buffer pool every
insert is a page read. Over 4.2 million rows the writes took 266 s with them kept and 145 s without, and building both
afterwards takes 22 s; the gap grows with the table."""


GALAXY_COMPONENT_COLUMNS = (
    "kind", "name", "catalogue_id", "source", "parent_component_id", "center_x_pc", "center_y_pc", "center_z_pc",
    "velocity_x_kms", "velocity_y_kms", "velocity_z_kms", "epoch_unix", "mass_solar", "half_mass_radius_pc", "profile",
    "king_w0", "tidal_radius_pc", "age_gy", "metallicity_fe_h", "core_collapsed", "seed")
"""tuple: The columns `replace_galaxy_components` writes and `get_galaxy_components` reads (the table's, less `id`,
`cell_key` and `created_at`)."""


def replace_galaxy_components(conn, found, kind=None):
    """
    Replaces the stored components of one kind (every kind when `kind` is `None`) with `found` (GEN.159): a plan
    draws a kind's whole set again from the seed, so there is no partial update. Does not commit. The components keep
    no `id` here; a `parent_component_id` must name a row that is already stored.

    Args:
        conn (Connection): The connection to write through.
        found (iterable[galaxy.components.Component]): The components to write.
        kind (str, optional): The kind to replace; `found` may then hold only that kind.
    """
    found = list(found)
    if kind is not None and any(component.kind != kind for component in found):
        raise ValueError(f"components of a kind other than {kind!r} given to replace it")
    if kind is None:
        conn.execute("DELETE FROM galaxy_components")
    else:
        conn.execute("DELETE FROM galaxy_components WHERE kind = ?", (kind,))
    conn.executemany(
        f"INSERT INTO galaxy_components ({', '.join(GALAXY_COMPONENT_COLUMNS)})"
        f" VALUES ({', '.join('?' * len(GALAXY_COMPONENT_COLUMNS))})",
        [(c.kind, c.name, c.catalogue_id, c.source, c.parent_component_id, *c.center_pc, *c.velocity_kms, c.epoch_unix,
          c.mass_solar, c.half_mass_radius_pc, c.profile, c.king_w0, c.tidal_radius_pc, c.age_gy, c.metallicity_fe_h,
          int(c.core_collapsed), c.seed) for c in found])


def get_galaxy_components(conn, kind=None, profile=None):
    """
    The stored components, by id.

    Args:
        conn (Connection): The connection to read through.
        kind (str, optional): Only this kind.
        profile (str, optional): Only this density profile (the containers' read wants the 157 King clusters among a
            plan's tens of thousands of components, PERF.165).

    Returns:
        list[galaxy.components.Component]: Each with its `id`.
    """
    sql = f"SELECT id, {', '.join(GALAXY_COMPONENT_COLUMNS)} FROM galaxy_components"
    tests = [(column, value) for column, value in (("kind", kind), ("profile", profile)) if value]
    where = " WHERE " + " AND ".join(f"{column} = ?" for column, _ in tests) if tests else ""
    rows = conn.execute(sql + where + " ORDER BY id", tuple(value for _, value in tests)).fetchall()
    return [galaxyComponents.Component(
        id=row["id"], kind=row["kind"], name=row["name"], catalogue_id=row["catalogue_id"], source=row["source"],
        parent_component_id=row["parent_component_id"],
        center_pc=(row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]),
        velocity_kms=(row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"]),
        epoch_unix=row["epoch_unix"], mass_solar=row["mass_solar"], half_mass_radius_pc=row["half_mass_radius_pc"],
        profile=row["profile"], king_w0=row["king_w0"], tidal_radius_pc=row["tidal_radius_pc"], age_gy=row["age_gy"],
        metallicity_fe_h=row["metallicity_fe_h"], core_collapsed=bool(row["core_collapsed"]), seed=row["seed"])
        for row in rows]


def _scatter_value(value):
    """One value of a `phenomenon_scatter` row as SQL text: numbers (also numpy's) and `None` as they are, the
    kind and subtype names quoted."""
    if value is None:
        return "NULL"
    if isinstance(value, str):
        return "'" + value.replace("\\", "\\\\").replace("'", "''") + "'"
    return str(value)


LITERAL_STATEMENT_CHARS = 1_000_000
"""int: A multi-row `INSERT` built by `_insert_literal_rows` is sent once its text reaches this length, whatever
`batch_size` says: the driver's own `executemany` cuts its statements at the same size, and a server's
`max_allowed_packet` may be small."""


def _insert_literal_rows(conn, table, columns, rows, batch_size):
    """Bulk-writes `rows` (values in `columns` order) into `table`, `batch_size` rows to a statement whose text is
    built here (`_scatter_value`): numbers, `None` and a few fixed names only, so the driver's own per-value
    escaping, about a quarter of an insert over tens of millions of rows, is skipped (PERF.73, PERF.48). A value
    that is any other kind of string must not come through here. Returns the rows written."""
    head = f"INSERT INTO {table} ({', '.join(columns)}) VALUES "
    written = 0
    batch = []
    chars = 0
    for row in rows:
        text = "(" + ",".join(map(_scatter_value, row)) + ")"
        batch.append(text)
        chars += len(text) + 1
        if len(batch) >= batch_size or chars >= LITERAL_STATEMENT_CHARS:
            conn.execute_literal(head + ",".join(batch))
            written += len(batch)
            batch = []
            chars = 0
    if batch:
        conn.execute_literal(head + ",".join(batch))
        written += len(batch)
    return written


def insert_phenomenon_scatter(conn, rows, batch_size=5000):
    """Bulk-writes scattered phenomena (`PHENOMENON_SCATTER_COLUMNS` order),
    `batch_size` rows to a statement. The statement text is built here
    (`_insert_literal_rows`), which costs about half what the driver's own
    per-value escaping does over tens of millions of rows (PERF.73).
    Returns the rows written."""
    return _insert_literal_rows(conn, "phenomenon_scatter", PHENOMENON_SCATTER_COLUMNS, rows, batch_size)


def drop_phenomenon_scatter_indexes(conn):
    """Drops the secondary indexes of `phenomenon_scatter` that exist, so a scatter's many inserts do not keep them
    up to date row by row; `add_phenomenon_scatter_indexes` builds them again from the finished table."""
    present = {row["name"] for row in conn.execute(
        "SELECT DISTINCT INDEX_NAME AS name FROM information_schema.STATISTICS"
        " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'phenomenon_scatter'").fetchall()}
    gone = [name for name in PHENOMENON_SCATTER_INDEXES if name in present]
    if gone:
        conn.execute("ALTER TABLE phenomenon_scatter " + ", ".join(f"DROP INDEX {name}" for name in gone))


def missing_phenomenon_scatter_indexes(conn):
    """The secondary indexes of `phenomenon_scatter` not built yet, as `(name, columns)` in build order."""
    present = {row["name"] for row in conn.execute(
        "SELECT DISTINCT INDEX_NAME AS name FROM information_schema.STATISTICS"
        " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'phenomenon_scatter'").fetchall()}
    return [(name, columns) for name, columns in PHENOMENON_SCATTER_INDEXES.items() if name not in present]


def add_phenomenon_scatter_indexes(conn, on_index=None):
    """Builds the secondary indexes `drop_phenomenon_scatter_indexes` took off (those missing), one statement per
    index so there is something to count: `on_index(name)` is called after each is built (PERF.101)."""
    for name, columns in missing_phenomenon_scatter_indexes(conn):
        conn.execute(f"ALTER TABLE phenomenon_scatter ADD KEY {name} ({columns})")
        if on_index is not None:
            on_index(name)


def phenomena_for_sector(conn, shell_index, band_index, wedge_index):
    """The scattered phenomena in one sector cell not yet built, in
    placement order, as dicts of every `phenomenon_scatter` column."""
    return [dict(row) for row in conn.execute(
        "SELECT * FROM phenomenon_scatter WHERE shell_index = ? AND band_index = ? AND wedge_index = ?"
        " AND built_at IS NULL ORDER BY id", (shell_index, band_index, wedge_index)).fetchall()]


def mark_phenomena_built(conn, scatter_ids):
    """Stamps the scattered phenomena a sector's save has built."""
    ids = [int(value) for value in scatter_ids]
    if ids:
        marks = ", ".join("?" * len(ids))
        for row in conn.execute(
                f"SELECT kind, COALESCE(subtype, '') AS subtype, COUNT(*) AS n FROM phenomenon_scatter"
                f" WHERE id IN ({marks}) AND built_at IS NULL GROUP BY kind, COALESCE(subtype, '')", ids).fetchall():
            conn.execute("UPDATE phenomenon_scatter_classes SET built = built + ? WHERE kind = ? AND subtype = ?",
                         (int(row["n"]), row["kind"], row["subtype"]))
        conn.executemany("UPDATE phenomenon_scatter SET built_at = CURRENT_TIMESTAMP WHERE id = ?",
                         [(value,) for value in ids])


UNTOUCHED_LEVEL = -1.0
"""float: `sector_stats.bright_level_sol` of a sector no backfill has
drawn (GEN.44): it follows the galaxy scatter's level."""

FILLED_LEVEL = 0.0
"""float: `sector_stats.bright_level_sol` of a generated sector (GEN.44)."""


def _address_chunks(addresses, size=500):
    keys = sorted({(int(a[0]), int(a[1]), int(a[2])) for a in addresses})
    for start in range(0, len(keys), size):
        yield keys[start:start + size]


def sector_bright_levels(conn, addresses):
    """
    The stored bright-star level (GEN.44) of each of `addresses` that has
    a `sector_stats` row: -1 untouched, a positive L_sun for the dimmest a
    backfill drew it down to, 0 filled.

    Args:
        conn (Connection): An open connection.
        addresses (iterable): `(shell, band, wedge)` tuples.

    Returns:
        dict: `(shell, band, wedge)` -> level, for the sectors with a row.
    """
    levels = {}
    for chunk in _address_chunks(addresses):
        marks = ", ".join("(?, ?, ?)" for _ in chunk)
        for row in conn.execute(
            "SELECT shell_index, band_index, wedge_index, bright_level_sol FROM sector_stats"
            f" WHERE (shell_index, band_index, wedge_index) IN ({marks})",
            tuple(value for key in chunk for value in key),
        ).fetchall():
            levels[(row["shell_index"], row["band_index"], row["wedge_index"])] = row["bright_level_sol"]
    return levels


def sector_bright_masses(conn, addresses, current=False):
    """
    The initial mass (solar masses) each of `addresses` was backfilled down
    to (GEN.187, `sector_stats.bright_mass_sol`), for the sectors that
    were. `current` reads the latest committed value, not the transaction's
    snapshot (`FOR UPDATE`), as the draw does once it holds the rows' locks.

    Returns:
        dict: `(shell, band, wedge)` -> mass, only for those with one.
    """
    masses = {}
    for chunk in _address_chunks(addresses):
        marks = ", ".join("(?, ?, ?)" for _ in chunk)
        for row in conn.execute(
            "SELECT shell_index, band_index, wedge_index, bright_mass_sol FROM sector_stats FORCE INDEX (PRIMARY)"
            f" WHERE bright_mass_sol IS NOT NULL AND (shell_index, band_index, wedge_index) IN ({marks})"
            + (" FOR UPDATE" if current else ""),
            tuple(value for key in chunk for value in key),
        ).fetchall():
            masses[(row["shell_index"], row["band_index"], row["wedge_index"])] = row["bright_mass_sol"]
    return masses


def set_sector_bright_masses(conn, states):
    """Records how far down by mass some sectors' stars now go (GEN.187),
    inside the transaction `lock_sector_stats` started: `(shell, band,
    wedge) -> (mass_solar, level_sol)`, `level_sol` the luminosity floor to
    store with it (the galaxy scatter's, so the sector counts as having its
    own level; `None` leaves the stored level alone)."""
    for address, (mass, level) in sorted(states.items()):
        if level is None:
            conn.execute("UPDATE sector_stats SET bright_mass_sol = ? WHERE shell_index = ? AND band_index = ?"
                         " AND wedge_index = ?", (mass, *address))
        else:
            conn.execute("UPDATE sector_stats SET bright_mass_sol = ?, bright_level_sol = ? WHERE shell_index = ?"
                         " AND band_index = ? AND wedge_index = ?", (mass, level, *address))


def sector_bright_level_keys(conn):
    """
    Every sector a backfill took to its own level (`bright_level_sol` above
    0, GEN.44), as `(shell, band, wedge) -> level`: what a staged scatter
    (`planetgen scatter --band --down-to`) leaves out of its layers
    and tops up sector by sector instead.
    """
    rows = conn.execute("SELECT shell_index, band_index, wedge_index, bright_level_sol FROM sector_stats"
                        " WHERE bright_level_sol > 0").fetchall()
    return {(row["shell_index"], row["band_index"], row["wedge_index"]): row["bright_level_sol"] for row in rows}


def lock_sector_stats(conn, entries):
    """
    Takes the row locks on some sectors' `sector_stats` rows (GEN.44),
    making each row first (untouched, -1) with its expected density
    (PERF.11) if there is none, so two backfills never draw the same
    sector at once. Rows are locked in address order, so two callers
    wait rather than deadlock. Holds until the caller commits.

    Args:
        conn (Connection): An open connection.
        entries (iterable): `((shell, band, wedge), relative_density,
            expected_systems)` per sector.

    Returns:
        dict: `(shell, band, wedge)` -> the level now stored.
    """
    rows = sorted((int(address[0]), int(address[1]), int(address[2]), density, expected)
                  for address, density, expected in entries)
    for start in range(0, len(rows), 500):
        chunk = rows[start:start + 500]
        # ON DUPLICATE KEY UPDATE takes each row's exclusive lock at once;
        # INSERT IGNORE would take a shared one, and two workers both
        # upgrading it to exclusive for the SELECT below would deadlock. The
        # SELECT is held to the primary key: left to choose, the server reads
        # the (level, address) index, which locks the rows in level order.
        # A derived table rather than VALUES(col), deprecated on MySQL
        # (tests/test_sql_portability.py).
        selects = " UNION ALL ".join(
            ["SELECT ? AS shell_index, ? AS band_index, ? AS wedge_index, ? AS density, ? AS expected"]
            + ["SELECT ?, ?, ?, ?, ?"] * (len(chunk) - 1)
        )
        conn.execute(
            "INSERT INTO sector_stats (shell_index, band_index, wedge_index, relative_density, expected_systems)"
            f" SELECT * FROM ({selects}) AS incoming"
            " ON DUPLICATE KEY UPDATE relative_density = COALESCE(sector_stats.relative_density, incoming.density),"
            " expected_systems = COALESCE(sector_stats.expected_systems, incoming.expected)",
            tuple(value for row in chunk for value in row),
        )
    levels = {}
    for chunk in _address_chunks([row[:3] for row in rows]):
        marks = ", ".join("(?, ?, ?)" for _ in chunk)
        for row in conn.execute(
            "SELECT shell_index, band_index, wedge_index, bright_level_sol FROM sector_stats FORCE INDEX (PRIMARY)"
            f" WHERE (shell_index, band_index, wedge_index) IN ({marks}) FOR UPDATE",
            tuple(value for key in chunk for value in key),
        ).fetchall():
            levels[(row["shell_index"], row["band_index"], row["wedge_index"])] = row["bright_level_sol"]
    return levels


def set_sector_bright_levels(conn, levels):
    """Records how dim some sectors' stars now go (`(shell, band, wedge) ->
    L_sun`), inside the transaction `lock_sector_stats` started."""
    conn.executemany(
        "UPDATE sector_stats SET bright_level_sol = ? WHERE shell_index = ? AND band_index = ? AND wedge_index = ?",
        [(level, *address) for address, level in sorted(levels.items())],
    )


def bright_star_fill_level(conn, shell_index, band_index, wedge_index):
    """
    The luminosity one sector's pre-placed stars go down to: its own
    backfill level (GEN.44, `sector_stats`) when it has one, else the
    galaxy's scatter threshold, else `None` (no bright stars placed; a
    fill caps nothing).
    """
    if shell_index is not None and band_index is not None and wedge_index is not None:
        level = sector_bright_levels(conn, [(shell_index, band_index, wedge_index)]).get(
            (shell_index, band_index, wedge_index))
        if level is not None and level > 0:
            return level
    return scatter_limits_of(conn).level_at(shell_index, band_index)


def bright_star_fill_mass_limit(conn, shell_index, band_index, wedge_index):
    """
    The initial mass one sector's pre-placed stars go down to (GEN.187): its
    own backfill's (`sector_stats.bright_mass_sol`) when it has one, else the
    galaxy scatter's mass limit, else `None` (no star was placed by mass; a
    fill caps nothing by mass).
    """
    if shell_index is not None and band_index is not None and wedge_index is not None:
        mass = sector_bright_masses(conn, [(shell_index, band_index, wedge_index)]).get(
            (shell_index, band_index, wedge_index))
        if mass is not None:
            return mass
    return scatter_limits_of(conn).mass_limit_at(shell_index, band_index)


DENSITY_RATIO_WINDOW = 1000
"""int: The decaying average of actual against expected systems (PERF.11)
is a plain mean over the first fills, then weighs each new fill 1 in
this many, so it follows a galaxy whose model drifts."""


def sector_habitable_worlds(conn, sector_id):
    """
    MAP.143: how many planets and moons the systems of one sector hold that are potentially habitable: an overall PHI-6
    score that is Blue or Green (`habitability.POTENTIALLY_HABITABLE_MIN_SCORE`, the stored `phi4`; a body with no score
    is not counted).
    """
    floor = habitability.POTENTIALLY_HABITABLE_MIN_SCORE
    planets = conn.execute(
        "SELECT COUNT(*) AS n FROM planets p JOIN star_systems ss ON ss.id = p.star_system_id"
        " WHERE ss.sector_id = ? AND p.phi4 >= ?", (sector_id, floor)).fetchone()["n"]
    moons = conn.execute(
        "SELECT COUNT(*) AS n FROM moons m JOIN star_systems ss ON ss.id = m.star_system_id"
        " WHERE ss.sector_id = ? AND m.phi4 >= ?", (sector_id, floor)).fetchone()["n"]
    return int(planets) + int(moons)


def record_system_interests(conn, where, params):
    """
    MAP.305: puts the `system_interests` rows of the systems matching `where` (SQL on the table aliased `ss`, with `params`) right
    (`interests.record`): a row for each reason a system has a planet or moon for, none for the rest. The Galaxy Map's marker layer
    reads them.
    """
    from planetgen.db import interests

    interests.record(conn, where, params)


def recount_habitable_worlds(conn, system_ids):
    """
    MAP.143: puts `sector_stats.habitable_worlds` right again for the sectors of `system_ids` after their scores
    were worked out again (a re-score, a dose refresh), and (MAP.305) each of those systems' interest rows. Only a
    filled sector's row is touched. Returns the `(shell, band, wedge)` of each sector whose count changed, for
    `block_stats.refresh_sectors`.
    """
    system_ids = list(system_ids)
    for start in range(0, len(system_ids), 500):
        chunk = system_ids[start:start + 500]
        record_system_interests(conn, f"ss.id IN ({', '.join('?' * len(chunk))})", chunk)
    sector_ids = set()
    for system_id in system_ids:
        row = conn.execute("SELECT sector_id FROM star_systems WHERE id = ?", (system_id,)).fetchone()
        if row is not None and row["sector_id"] is not None:
            sector_ids.add(row["sector_id"])
    changed = []
    for sector_id in sorted(sector_ids):
        stats = conn.execute(
            "SELECT s.shell_index, s.band_index, s.wedge_index, st.habitable_worlds FROM sectors s"
            " JOIN sector_stats st ON st.shell_index = s.shell_index AND st.band_index = s.band_index"
            " AND st.wedge_index = s.wedge_index WHERE s.id = ? AND st.bright_level_sol = 0", (sector_id,)).fetchone()
        if stats is None:
            continue
        count = sector_habitable_worlds(conn, sector_id)
        if count != stats["habitable_worlds"]:
            address = (stats["shell_index"], stats["band_index"], stats["wedge_index"])
            conn.execute("UPDATE sector_stats SET habitable_worlds = ? WHERE shell_index = ? AND band_index = ?"
                         " AND wedge_index = ?", (count, *address))
            changed.append(address)
    return changed


def _content_figures(conn, sector_id):
    """What a sector's rows hold, as `sector_stats` states it: `(systems, stars, mean temperature K, mean luminosity in solar
    units, mean age Gy, summed luminosity in solar units, habitable worlds)`."""
    systems = conn.execute("SELECT COUNT(*) AS n FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["n"]
    stars = conn.execute(
        "SELECT COUNT(*) AS n, AVG(st.temperature_k) AS temperature, AVG(st.luminosity_w) AS luminosity,"
        " AVG(st.age_gy) AS age, SUM(st.luminosity_w) AS total_luminosity"
        " FROM stars st JOIN star_systems ss ON ss.id = st.star_system_id WHERE ss.sector_id = ?", (sector_id,),
    ).fetchone()
    luminosity = None if stars["luminosity"] is None else float(stars["luminosity"]) / physical_constants.SOLAR_LUMINOSITY
    temperature = None if stars["temperature"] is None else float(stars["temperature"])
    age = None if stars["age"] is None else float(stars["age"])
    total_luminosity = (None if stars["total_luminosity"] is None
                        else float(stars["total_luminosity"]) / physical_constants.SOLAR_LUMINOSITY)
    return (systems, stars["n"], temperature, luminosity, age, total_luminosity, sector_habitable_worlds(conn, sector_id))


def refresh_edited_sectors(conn, sector_ids):
    """
    A generated sector's rows were edited (DB.83): states what it holds again in its `sector_stats` row (the systems, stars,
    mean temperature, luminosity and age, summed luminosity and habitable worlds; its expected figures, density and bright-star
    level are the galaxy's and stay) and brings the totals of the blocks above it up to them (`block_stats.refresh_chain`, which
    also brings the region's). Call it inside the edit's transaction, last. A sector that is not generated at a grid address
    (a one-off system's sector, a sector with no stats row) is left alone. Returns the addresses refreshed.
    """
    from planetgen.db import block_stats

    done = []
    for sector_id in sorted({int(i) for i in sector_ids if i is not None}):
        row = conn.execute("SELECT shell_index, band_index, wedge_index FROM sectors WHERE id = ?", (sector_id,)).fetchone()
        if row is None or row["shell_index"] is None:
            continue
        address = (row["shell_index"], row["band_index"], row["wedge_index"])
        figures = _content_figures(conn, sector_id)
        found = conn.execute(
            "UPDATE sector_stats SET actual_systems = ?, actual_stars = ?, mean_temperature_k = ?, mean_luminosity_sol = ?,"
            " mean_age_gy = ?, total_luminosity_sol = ?, habitable_worlds = ?"
            " WHERE shell_index = ? AND band_index = ? AND wedge_index = ? AND bright_level_sol = 0", (*figures, *address))
        if found.rowcount:
            block_stats.refresh_chain(conn, address)
            done.append(address)
    return done


def record_sector_stats(conn, sector_id, address, center_pc):
    """
    Writes a sector's stats once it is generated (PERF.11, GEN.44): its
    expected density from the galaxy model, the systems and stars it got,
    the mean temperature, mean luminosity and age of those stars and
    their summed luminosity (DB.14; the map colors from these), the habitable
    worlds they hold (MAP.143, `sector_habitable_worlds`), the mass of everything in it
    (DB.108, `sector_mass.sector_mass_kg`), and
    level 0 (filled), keeping the level it had in
    `level_before_fill_sol`; then folds its actual-to-expected systems
    into `galaxy_shape`'s decaying average. Inside the sector's own save,
    last, so the galaxy row is locked only for the commit.

    Args:
        conn (Connection): An open connection, mid-transaction.
        sector_id (int): The new `sectors.id`.
        address (tuple): Its `(shell, band, wedge)`.
        center_pc (tuple): Its galaxy-frame center, parsecs.
    """
    from planetgen.galaxy import profile as galaxyProfile, resolution as galaxyResolution, spherical
    from planetgen.galaxy.density import relative_density

    skeleton = get_galaxy_shape(conn)
    density = expected = None
    if skeleton is not None:
        density = max(relative_density(center_pc, skeleton.shape), 0.0)
        expected = density * skeleton.expected_system_count_at_density_1
        cell_expected = expected * spherical.cell_volume(*address, skeleton.edge_pc) / skeleton.edge_pc ** 3
        if cell_expected > skeleton.cell_star_ceiling:
            # GEN.281: a capped cell aimed at the stars it resolves, so its ratio to what it got stays meaningful.
            found = galaxyResolution.cell_resolution(
                galaxyProfile.sector_profile(address, skeleton), skeleton.cell_star_ceiling,
                skeleton.fill_luminosity_cap_sol, scatter_limits_of(conn, skeleton).mass_limit_at(*address[:2]))
            if found.capped:
                expected = found.resolved_systems
    systems, star_count, temperature, luminosity, age, total_luminosity, habitable = _content_figures(conn, sector_id)
    mass = sector_mass.sector_mass_kg(conn, sector_id)
    record_system_interests(conn, "ss.sector_id = ?", (sector_id,))
    conn.execute(
        "INSERT INTO sector_stats (shell_index, band_index, wedge_index, bright_level_sol, level_before_fill_sol,"
        " relative_density, expected_systems, actual_systems, actual_stars, mean_temperature_k, mean_luminosity_sol,"
        " mean_age_gy, total_luminosity_sol, habitable_worlds, total_mass_kg, filled_at)"
        " VALUES (?, ?, ?, 0, -1, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP(3))"
        " ON DUPLICATE KEY UPDATE level_before_fill_sol = IF(bright_level_sol = 0, level_before_fill_sol,"
        " bright_level_sol), bright_level_sol = 0, relative_density = ?, expected_systems = ?, actual_systems = ?,"
        " actual_stars = ?, mean_temperature_k = ?, mean_luminosity_sol = ?, mean_age_gy = ?,"
        " total_luminosity_sol = ?, habitable_worlds = ?, total_mass_kg = ?, filled_at = CURRENT_TIMESTAMP(3)",
        (*address, *(2 * (density, expected, systems, star_count, temperature, luminosity, age, total_luminosity,
                          habitable, mass))),
    )
    if expected:
        ratio = systems / expected
        conn.execute(
            "UPDATE galaxy_shape SET density_ratio_avg = IF(density_ratio_avg IS NULL, ?,"
            " density_ratio_avg + (? - density_ratio_avg) / LEAST(density_ratio_samples + 1, ?)),"
            " density_ratio_samples = density_ratio_samples + 1 WHERE id = 1",
            (ratio, ratio, DENSITY_RATIO_WINDOW),
        )


def forget_sector_fill(conn, address):
    """
    A generated sector was deleted (GEN.44): its stats row goes back to
    the bright-star level it had before the fill, and its actual stats
    are cleared (the decaying average keeps what it learned). Nothing for
    a sector off the grid.
    """
    if address is None or any(value is None for value in address):
        return
    conn.execute(
        "UPDATE sector_stats SET bright_level_sol = COALESCE(level_before_fill_sol, -1), level_before_fill_sol = NULL,"
        " actual_systems = NULL, actual_stars = NULL, mean_temperature_k = NULL, mean_luminosity_sol = NULL,"
        " mean_age_gy = NULL, total_luminosity_sol = NULL, habitable_worlds = NULL, total_mass_kg = NULL,"
        " filled_at = NULL WHERE shell_index = ? AND band_index = ? AND wedge_index = ? AND bright_level_sol = 0",
        tuple(address),
    )
    from planetgen.db import block_lights, block_stats    # DB.25: the blocks above it no longer count it

    block_stats.refresh_chain(conn, address)
    block_lights.refresh_chain(conn, address)    # MAP.259: its stars that were built are pre-placed again


SECTOR_DELETIONS_COUNTER = "sdel"
"""str: The `id_counters.kind` of the galaxy's deletion epoch (PERF.38): one row, scope empty, whose number goes up by
one with every sector deleted. `query.galaxy_content_state` reads it where it used to count the placed sectors."""


def note_sector_deleted(conn):
    """
    PERF.38: a sector was deleted. Moves the deletion epoch up in the
    deleting transaction, so a cache that remembers the epoch knows a node
    may have lost a sector (a deletion leaves no row to find) without a
    linear count of the placed sectors on every check. The row survives a
    `planetgen reset`, like every `id_counters` row, so it only goes up.
    """
    conn.execute(
        "INSERT INTO id_counters (kind, scope, next_value) VALUES (?, ?, 1)"
        " ON DUPLICATE KEY UPDATE next_value = next_value + 1",
        (SECTOR_DELETIONS_COUNTER, b""),
    )


def sector_deletions(conn):
    """The deletion epoch (`note_sector_deleted`): 0 before any sector was deleted."""
    row = conn.execute("SELECT next_value FROM id_counters WHERE kind = ? AND scope = ?",
                       (SECTOR_DELETIONS_COUNTER, b"")).fetchone()
    return int(row["next_value"]) if row else 0


def get_sector_stats(conn, shell_index, band_index, wedge_index):
    """One sector's `sector_stats` row as a dict (PERF.11), or `None`."""
    return conn.execute(
        "SELECT * FROM sector_stats WHERE shell_index = ? AND band_index = ? AND wedge_index = ?",
        (shell_index, band_index, wedge_index),
    ).fetchone()


def galaxy_density_ratio(conn):
    """`(average, samples)`: the decaying average of the systems sector
    fills got against the systems expected (PERF.11); the average is
    `None` before any fill."""
    row = conn.execute("SELECT density_ratio_avg, density_ratio_samples FROM galaxy_shape WHERE id = 1").fetchone()
    if row is None:
        return None, 0
    return row["density_ratio_avg"], int(row["density_ratio_samples"])


def delete_unfinished_band(conn, limits, keep_addresses=(), batch_size=500):
    """
    Deletes the unbuilt bright stars below the galaxy's star-fill level outside `keep_addresses`: what a band run
    (`planetgen scatter --band --down-to`) that stopped part way left
    in the bands it got through, since the level only moves once a band
    is whole (GEN.32). `limits` (`galaxy.scatter_limits.ScatterLimits`, GEN.333) is the level before the run, and
    each band's level is its tier's. Every finished scatter or band is at or above the
    level, and the stars below it that belong there (a backfilled block's,
    a filled sector's) are in `keep_addresses` or built.

    Returns:
        int: Stars deleted.
    """
    keep = set(keep_addresses)
    floors = [found.min_luminosity_sol * physical_constants.SOLAR_LUMINOSITY for found in limits.tiers]
    rows = conn.execute(
        "SELECT shell_index, band_index, wedge_index FROM bright_stars"
        " WHERE luminosity_w < ? AND star_system_id IS NULL"
        " GROUP BY shell_index, band_index, wedge_index", (max(floors),)).fetchall()
    cells = [(row["shell_index"], row["band_index"], row["wedge_index"]) for row in rows]
    cells = [cell for cell in cells if cell not in keep]
    tiers = _cell_tiers(limits, cells)
    deleted = 0
    for tier in range(3):
        below_luminosity_w = floors[tier]
        mine = [cell for cell in cells if tiers[cell] == tier]
        for start in range(0, len(mine), batch_size):
            chunk = mine[start:start + batch_size]
            placeholders = ", ".join("(?, ?, ?)" for _ in chunk)
            scatter_systems.delete_member_trees_where(
                conn, "bright_stars", "t.luminosity_w < ? AND t.star_system_id IS NULL"
                f" AND (t.shell_index, t.band_index, t.wedge_index) IN ({placeholders})",
                (below_luminosity_w, *[value for cell in chunk for value in cell]))
            cur = conn.execute(
                "DELETE FROM bright_stars WHERE luminosity_w < ? AND star_system_id IS NULL"
                f" AND (shell_index, band_index, wedge_index) IN ({placeholders})",
                (below_luminosity_w, *[value for cell in chunk for value in cell]))
            deleted += cur.rowcount
    return deleted


def insert_bright_stars(conn, rows, batch_size=10000):
    """
    Bulk-writes scattered bright stars (`BRIGHT_STAR_COLUMNS` order) in
    batches of `batch_size`, which need no server setting (`LOAD DATA
    LOCAL INFILE` would need `local_infile` on both ends). Only numbers,
    `None` and the fixed population, spectral type and class names are in
    a row, so the statement text is built here (`_insert_literal_rows`,
    PERF.48) instead of through the driver's per-value escaping.

    Returns:
        int: Rows written.
    """
    return _insert_literal_rows(conn, "bright_stars", BRIGHT_STAR_COLUMNS, rows, batch_size)


def bright_stars_for_sector(conn, shell_index, band_index, wedge_index):
    """The pre-placed stars in one sector cell not yet built into a
    system, brightest first, as dicts of every `bright_stars` column."""
    return [dict(row) for row in conn.execute(
        "SELECT * FROM bright_stars WHERE shell_index = ? AND band_index = ? AND wedge_index = ?"
        " AND star_system_id IS NULL ORDER BY luminosity_w DESC, id",
        (shell_index, band_index, wedge_index),
    ).fetchall()]


def database_now(conn):
    """The database server's own clock (`NOW()`), to compare with
    `created_at` columns."""
    return conn.execute("SELECT NOW() AS now").fetchone()["now"]


def sector_centers_since(conn, since):
    """The `(x, y, z)` galaxy-frame centers, parsecs, of every galaxy
    sector created at or after `since` (`database_now`), in address order:
    not save order, which the worker count changes (GEN.39)."""
    return [
        (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
        for row in conn.execute(
            "SELECT center_x_pc, center_y_pc, center_z_pc FROM sectors WHERE created_at >= ?"
            " AND shell_index IS NOT NULL ORDER BY shell_index, band_index, wedge_index", (since,)
        ).fetchall()
    ]


def filled_sector_addresses(conn):
    """The `(shell, band, wedge)` of every galaxy sector already filled."""
    return {
        (row["shell_index"], row["band_index"], row["wedge_index"])
        for row in conn.execute(
            "SELECT shell_index, band_index, wedge_index FROM sectors WHERE shell_index IS NOT NULL"
        ).fetchall()
    }


def mark_bright_star_filled(conn, bright_star_id, star_system_id):
    """Links a pre-placed star to the system its sector's fill built
    around it."""
    conn.execute("UPDATE bright_stars SET star_system_id = ? WHERE id = ?", (star_system_id, bright_star_id))


def sectors_reached_by(conn, center_pc, radius_pc):
    """The ids of every placed sector a sphere of `radius_pc` around
    `center_pc` (galaxy-frame parsecs) overlaps. Padded by the farthest a cell
    reaches from its center (`spherical.CELL_REACH_EDGES` edges: a polar cap
    reaches about 4.4 pc at the 4 pc edge, past the cube's 3.46, NAV.53)."""
    edge_row = conn.execute("SELECT MAX(edge_mpc) AS edge FROM sectors").fetchone()
    if edge_row["edge"] is None:
        return []
    reach = radius_pc + CELL_REACH_EDGES * edge_row["edge"] / 1000.0
    rows = conn.execute(
        "SELECT id, center_x_pc, center_y_pc, center_z_pc FROM sectors"
        " WHERE center_x_pc BETWEEN ? AND ? AND center_y_pc BETWEEN ? AND ? AND center_z_pc BETWEEN ? AND ?",
        (center_pc[0] - reach, center_pc[0] + reach, center_pc[1] - reach, center_pc[1] + reach,
         center_pc[2] - reach, center_pc[2] + reach),
    ).fetchall()
    return [
        row["id"] for row in rows
        if detmath.dist(center_pc, (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])) <= reach
    ]


def _refresh_containment_around(conn, placement, radius_ly):
    """After placing a nebula or remnant: refresh every sector it reaches."""
    if placement is None or placement.get("center_x_pc") is None:
        return
    center = (placement["center_x_pc"], placement["center_y_pc"], placement["center_z_pc"])
    refresh_containment(conn, sectors_reached_by(conn, center, ly_to_pc(radius_ly)))


def insert_sector(conn, sector: SpaceSector, galaxy_position=None, link_neighbors=True, run_id=None) -> int:
    """
    Inserts a full `SpaceSector` -- the `sectors` row, every system it
    contains (with its placement), and every exotic phenomenon it contains
    (`sector.phenomena`, see `run_sector.generate_sector_phenomena`) -- into
    the database. A phenomenon's own sector-relative `(x, y, z)` position
    (in light-years, from `SpaceSector.add_phenomenon`) is converted to an
    absolute galaxy-frame center via `_galaxy_placement_from_sector_offset`
    when `galaxy_position` is given (`None` when this sector itself was
    placement) -- every phenomenon type since v28 (see `schema.sql`'s
    "v28" header note).

    Args:
        conn (Connection): An open, schema-initialized connection.
        sector (SpaceSector): The sector to persist.
        galaxy_position (dict, optional): This sector's galaxy-frame
            placement (see `schema.sql`'s "v4"/"v6/v7" header notes), or
            `None` (the default) for a sector never placed in a galaxy --
            `generation/run_sector.py`'s own standalone CLI keeps producing these.
            When given, must have keys `center_x_pc`, `center_y_pc`,
            `center_z_pc`, `galactic_radius_pc` (all required together --
            the schema's CHECK constraint enforces that), and optionally
            `shell_index`/`band_index`/`wedge_index`, the sector's
            cell in the spherical grid (`None`/omitted leaves them NULL
            -- a hand-placed position with no grid address).
        run_id (int, optional): The `generation_runs.id` of the run saving the sector, kept on the sector with the
            generator epoch (DB.16). `None` uses the run row for the running code's version that saves made outside a run
            share (`sector_run_id`).

    Rows are written in batches (PERF.13, `Connection.batched`): one
    multi-row INSERT per table and statement shape, with ids from
    `id_blocks`, instead of one INSERT per row.

    A grid-addressed sector's stats (`sector_stats`, PERF.11, GEN.44) are
    written last, from the rows just stored (`record_sector_stats`).

    Returns:
        int: The new `sectors.id`.
    """
    try:
        with conn.batched():
            sector_id = _insert_sector_rows(conn, sector, galaxy_position, link_neighbors, run_id)
    finally:
        issuer, conn.uid_issuer = conn.uid_issuer, None
    if issuer is None or not issuer.complete:
        assign_uids(conn, sector_id=sector_id)
    else:
        set_sector_uid(conn, sector_id)
    address = None if galaxy_position is None else tuple(
        galaxy_position.get(key) for key in ("shell_index", "band_index", "wedge_index"))
    if address is not None and None not in address:
        record_sector_stats(conn, sector_id, address, (galaxy_position["center_x_pc"], galaxy_position["center_y_pc"],
                                                       galaxy_position["center_z_pc"]))
    conn.parity_touch("id = ?", (sector_id,))
    return sector_id


def _set_reserved_name(obj, final, rename_bodies=False):
    """Gives a system or phenomenon the name `reserve_system_names` settled
    on, and its remnant's core with it. `rename_bodies`: also rename a
    system's stars, planets and moons (for a system already inserted)."""
    core = getattr(obj, "compact_remnant", None) if isinstance(obj, SupernovaRemnant) else None
    if core is not None and core.name == f"{obj.name} Core":
        core.name = f"{final} Core"
    obj.name = final
    if rename_bodies and hasattr(obj, "assign_names"):
        obj.assign_names(final)


def _save_progressed(conn):
    """One more system or phenomenon of a sector is saved (PERF.50)."""
    progress = getattr(conn, "save_progress", None)
    if progress is not None:
        progress.advance()


SAVE_RUN_COMMAND = "save"
"""str: The `generation_runs.command` of the row that sectors saved outside a run (a test, the admin pages) point at."""


def sector_run_id(conn, run_id=None):
    """
    The `generation_runs.id` a sector being saved is recorded under (DB.16): `run_id` when the run that saves it is
    known, otherwise the one row, per version of the code, that saves made outside a run share (command `save`, made
    the first time it is needed and kept on the connection).
    """
    if run_id is not None:
        return run_id
    cached = getattr(conn, "save_run_id", None)
    if cached is not None:
        return cached
    made_by = versionKey.current()
    values = (SAVE_RUN_COMMAND, made_by["version_key"], made_by["planetgen_version"], made_by["python_version"],
              made_by["platform"], generatorEpoch.generator_epoch())
    row = conn.execute(
        "SELECT id FROM generation_runs WHERE command = ? AND version_key = ? AND planetgen_version = ?"
        " AND python_version = ? AND platform = ? AND generator_epoch <=> ? ORDER BY id LIMIT 1", values).fetchone()
    if row is not None:
        conn.save_run_id = row["id"]
    else:
        conn.save_run_id = conn.execute(
            "INSERT INTO generation_runs (command, version_key, planetgen_version, python_version, platform,"
            " generator_epoch, finished_at, outcome) VALUES (?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP(3), 'ok')",
            values).lastrowid
    return conn.save_run_id


def _insert_sector_rows(conn, sector, galaxy_position, link_neighbors=True, run_id=None):
    run_id = sector_run_id(conn, run_id)  # DB.16: the run (and so the code) generating this sector
    epoch = generatorEpoch.generator_epoch()
    if galaxy_position is not None:
        cur = conn.execute(
            """
            INSERT INTO sectors (
                name, edge_mpc, center_x_pc, center_y_pc, center_z_pc,
                galactic_radius_pc, shell_index, band_index, wedge_index,
                generator_epoch, run_id
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                sector.name, ly_to_milliparsecs(sector.edge_ly),
                galaxy_position["center_x_pc"], galaxy_position["center_y_pc"],
                galaxy_position["center_z_pc"], galaxy_position["galactic_radius_pc"],
                galaxy_position.get("shell_index"), galaxy_position.get("band_index"),
                galaxy_position.get("wedge_index"), epoch, run_id,
            ),
        )
    else:
        cur = conn.execute(
            "INSERT INTO sectors (name, edge_mpc, generator_epoch, run_id) VALUES (?, ?, ?, ?)",
            (sector.name, ly_to_milliparsecs(sector.edge_ly), epoch, run_id),
        )
    sector_id = cur.lastrowid

    # Name-uniqueness (v24, names/uniqueness.py) -- reserved right after the
    # INSERT (the registry row points at it), mutating sector.name in
    # place, so every print/rendering call site that reads it afterward
    # (including this same sector's own systems, added below) sees the
    # final, collision-free name. See `reserve_sector_name`.
    final_name, _name_base = reserve_sector_name(conn, sector.name, sector_id)
    if final_name != sector.name:
        conn.execute("UPDATE sectors SET name = ?, parity_dirty = 1 WHERE id = ?", (final_name, sector_id))
    sector.name = final_name

    # Every system and uniquely named phenomenon in the sector reserves
    # its name in one go (PERF.14); the insert functions below take their
    # reservation from `prereserved_names` and leave their registry writes
    # in `deferred_name_confirmations` for one upsert at the end.
    conn.prereserved_names = {}
    conn.deferred_name_confirmations = []
    conn.uid_issuer = _uid_issuer_for_sector(conn, sector_id, galaxy_position)
    try:
        # Every placed phenomenon is named by its object ID (GEN.64),
        # claimed here for the whole sector in generation order; only the
        # rest go through the name registry. A bright-sweep system is named
        # by the registry like any system.
        by_id = _sector_object_ids(sector, galaxy_position)
        _claim_object_ids(conn, by_id, _address_of(galaxy_position))
        claimed = {id(obj) for obj, _kind, _center in by_id}
        named = [entry.star_system for entry in sector.entries if id(entry.star_system) not in claimed]
        named += [entry.phenomenon for entry in sector.phenomena
                  if _phenomenon_registers_name(entry.phenomenon) and id(entry.phenomenon) not in claimed]
        reservations = reserve_system_names(conn, [obj.name for obj in named])
        for obj, (final, base, diminutive_index) in zip(named, reservations):
            _set_reserved_name(obj, final)
            conn.prereserved_names[id(obj)] = (base, diminutive_index)

        bright_star_links = []
        built_tree_ids = []
        built_scatter_ids = []
        nearest = sector.all_nearest_neighbors(3)   # PERF.81: one spatial tree, not a distance to every system each
        for entry, neighbors in zip(sector.entries, nearest):
            star_system_id = insert_star_system(
                conn, entry.star_system, entry.system_config,
                sector_id=sector_id, position=entry.position,
                location=_location_for_entry(sector, entry, neighbors),
            )
            bright_star_id = getattr(entry, "bright_star_id", None)
            if bright_star_id is not None:
                bright_star_links.append((bright_star_id, star_system_id))
                if getattr(entry.star_system, "tree_built", False):
                    built_tree_ids.append(bright_star_id)
            if getattr(entry, "scatter_id", None) is not None:
                built_scatter_ids.append(entry.scatter_id)
            _save_progressed(conn)

        for entry in sector.phenomena:
            inserter = _PHENOMENON_INSERTERS.get(entry.phenomenon_type)
            if inserter is None:
                raise ValueError(f"Unknown phenomenon type: {entry.phenomenon_type!r}")
            placement = _galaxy_placement_from_sector_offset(galaxy_position, entry.position)
            inserter(conn, entry.phenomenon, sector_id=sector_id, placement=placement)
            if getattr(entry, "scatter_id", None) is not None:
                built_scatter_ids.append(entry.scatter_id)
            _save_progressed(conn)

        healed = confirm_system_names(conn, conn.deferred_name_confirmations)
        if healed:
            seen = set()
            for obj, (_final, base, _index) in zip(named, reservations):
                key = _name_key(base)
                if key not in seen and key in healed:
                    _set_reserved_name(obj, healed[key], rename_bodies=True)
                seen.add(key)
    finally:
        conn.prereserved_names = None
        conn.deferred_name_confirmations = None
    if bright_star_links:
        _update_by_id(conn, "bright_stars", ("star_system_id",), bright_star_links, touch=True)
    if built_tree_ids:  # GEN.315: a built row's members and nodes are its system's stars and orbit_nodes now
        scatter_systems.delete_member_trees(conn, "bright_stars", built_tree_ids)
    mark_phenomena_built(conn, built_scatter_ids)
    if _address_of(galaxy_position) is not None:
        fingerprintTree.store_leaf(conn, sector_id, _address_of(galaxy_position),
                                   fingerprintTree.leaf_digest(sector, _address_of(galaxy_position)))

    field_nebulae = getattr(sector, "field_nebulae", None)
    if galaxy_position is not None and (link_neighbors or field_nebulae):
        # Containment and nearest-neighbor lists read the sectors around
        # this one and rewrite theirs, so two neighbors saved at once
        # (parallel generation, PERF.7) would each miss the other. One
        # writer at a time does this last step, holding the lock until
        # commit; the slower part above still runs side by side. A run that
        # links its sectors later (PERF.45, `link_sector_neighbors`) only needs the
        # lock for a cloud's first sector.
        conn.lock_until_commit(_neighbor_lock_name(conn))
        _insert_field_nebulae(conn, sector, sector_id)
        if link_neighbors:
            refresh_containment(conn, [sector_id])
            _add_sector_to_nearest(conn, sector_id)

    return sector_id


FIELD_NEBULA_MATCH_PC = 1e-6
"""float: How close a stored nebula's center must be to a field cloud's to
be that cloud (the field draws the same center every time; this only
absorbs rounding)."""


def _insert_field_nebulae(conn, sector, sector_id):
    """
    Stores the galaxy's molecular clouds that reach `sector`
    (`sector.field_nebulae`, GEN.47) that no earlier sector stored: the
    first sector saved that a cloud reaches becomes its home sector. Runs
    under the neighbor lock, so two sectors saved at once can't both store
    the same cloud. `insert_nebula` names each by its object ID and
    refreshes containment in every sector it reaches.
    """
    cloud_field = getattr(sector, "cloud_field", None)
    births = {}
    for nebula, center in getattr(sector, "field_nebulae", None) or ():
        x, y, z = center
        pad = FIELD_NEBULA_MATCH_PC
        stored = conn.execute(
            "SELECT id FROM nebulae WHERE center_x_pc BETWEEN ? AND ? AND center_y_pc BETWEEN ? AND ?"
            " AND center_z_pc BETWEEN ? AND ? LIMIT 1",
            (x - pad, x + pad, y - pad, y + pad, z - pad, z + pad),
        ).fetchone()
        if stored is not None:
            continue
        placement = {"center_x_pc": x, "center_y_pc": y, "center_z_pc": z,
                     "galactic_radius_pc": math.sqrt(x * x + y * y + z * z)}
        if conn.uid_issuer is not None and cloud_field is not None:
            # GEN.176: born in the sector holding its centroid, with a serial worked out from the seed alone, so the
            # ID does not depend on which sector reached the cloud first.
            edge_pc = sector.sector_edge_pc or tuning.DEFAULT_SECTOR_EDGE_PC
            conn.uid_issuer.given_nebula_uid = objectUid.to_bytes(nebulaField.cloud_object_id(
                cloud_field[0], cloud_field[1], nebula, center, edge_pc, births))
        nebula_id = insert_nebula(conn, nebula, sector_id=sector_id, placement=placement)
        if conn.uid_issuer is not None:
            conn.uid_issuer.given_nebula_uid = None
        _mark_nebula_field_stored(conn, nebula_id, center)


def _mark_nebula_field_stored(conn, nebula_id, center):
    """Sets `stored` and the object ID on the `nebula_field` row of the field cloud just stored as nebula `nebula_id`
    (DB.24), found by its cell and exact centre; nothing if the table has no such row (a class P or Q cloud, or a
    galaxy planned before v83)."""
    x, y, z = center
    i, j, k = nebulaField.cell_index(center)
    conn.execute(
        "UPDATE nebula_field SET is_stored = 1, uid = (SELECT uid FROM nebulae WHERE id = ?)"
        " WHERE cell_i = ? AND cell_j = ? AND cell_k = ? AND center_x_pc = ? AND center_y_pc = ? AND center_z_pc = ?",
        (nebula_id, i, j, k, x, y, z))


def _sector_object_ids(sector, galaxy_position):
    """
    `(obj, kind, center_pc)` for everything in a galaxy-placed sector that
    is named by its object ID (GEN.64): each phenomenon, then each
    bright-sweep system (`bright_star_id`). Empty for a sector never
    placed in the galaxy, whose objects keep their generated names. A
    registry-named phenomenon whose name was given by hand (`name_given`)
    keeps it.
    """
    if galaxy_position is None:
        return []
    items = []
    for entry in sector.phenomena:
        if getattr(entry.phenomenon, "name_given", False) and _phenomenon_registers_name(entry.phenomenon):
            continue  # a name given by hand stays
        if entry.phenomenon_type == "quasar":
            center = _placement_center(GALACTIC_CENTER_PLACEMENT)
        else:
            center = _placement_center(_galaxy_placement_from_sector_offset(galaxy_position, entry.position))
        items.append((entry.phenomenon, entry.phenomenon_type, center))
        core_item = _remnant_core_item(entry.phenomenon, center)
        if core_item is not None:
            items.append(core_item)
    return items


def _remnant_core_item(remnant, remnant_center_pc):
    """`(core, kind, center_pc)` for a supernova remnant's collapsed core
    (GEN.64: it has its own object ID), at the remnant's center moved by
    the core's birth kick, or `None` when there's no core."""
    core = getattr(remnant, "compact_remnant", None) if isinstance(remnant, SupernovaRemnant) else None
    if core is None:
        return None
    kind = "black-hole-core" if isinstance(core, BlackHole) else "neutron-star-core"
    offset_ly = getattr(remnant, "compact_offset_ly", None) or (0.0, 0.0, 0.0)
    return core, kind, tuple(remnant_center_pc[i] + ly_to_pc(offset_ly[i]) for i in range(3))


def _phenomenon_registers_name(phenomenon):
    """Whether a sector phenomenon's name goes through
    `system_name_registry` (`NAMED_PHENOMENON_TABLES`) -- every
    standalone black hole, neutron star, nebula, supernova remnant, rogue
    planet and quasar; not interstellar comets or asteroid fields, which
    get designations."""
    return isinstance(phenomenon, (BlackHole, NeutronStar, Nebula, SupernovaRemnant, RoguePlanet, Quasar))


def sector_for_placement(conn, sector_id):
    """
    A `SpaceSector` holding just enough of a stored sector to place one
    more system in it (`SpaceSector.add_system`): its edge and grid cell,
    and every stored system as a light stand-in carrying its name,
    position and Hill-sphere radius (`system_perimeter`), not the full
    object graph `load_sector` rebuilds.

    Args:
        conn (Connection): An open connection.
        sector_id (int): The `sectors.id`.

    Returns:
        SpaceSector: The stand-in sector (`name` the stored name).

    Raises:
        ValueError: If no such sector exists.
    """
    row = conn.execute("SELECT name, edge_mpc, shell_index, band_index, wedge_index FROM sectors WHERE id = ?",
                       (sector_id,)).fetchone()
    if row is None:
        raise ValueError(f"no sectors row with id {sector_id}")
    edge_ly = milliparsecs_to_ly(row["edge_mpc"])
    cell = _sector_cell(row, edge_ly)
    sector = SpaceSector(row["name"], edge_ly=edge_ly, cell=cell)
    systems = conn.execute(
        "SELECT ss.id, ss.name, ss.position_x_mpc, ss.position_y_mpc, ss.position_z_mpc,"
        " ss.binary_configuration, ss.binary_system_perimeter_km,"
        " (SELECT MAX(st.system_perimeter_km) FROM stars st WHERE st.star_system_id = ss.id) AS star_perimeter_km"
        " FROM star_systems ss WHERE ss.sector_id = ? AND ss.position_x_mpc IS NOT NULL ORDER BY ss.id",
        (sector_id,),
    ).fetchall()
    for system in systems:
        perimeter_km = (system["binary_system_perimeter_km"] if system["binary_configuration"] == "close"
                        else system["star_perimeter_km"]) or 0.0
        stand_in = SimpleNamespace(
            name=system["name"], system_config=None,
            star=SimpleNamespace(system_perimeter=_km_to_au(perimeter_km)))
        position = tuple(milliparsecs_to_ly(system[f"position_{axis}_mpc"]) for axis in "xyz")
        sector.entries.append(SectorSystemEntry(stand_in, position, sector_center_ly=sector.center_galactic_ly))
    return sector


def add_system_to_sector(conn, sector_id, star_system, system_config, position=None):
    """
    Saves one generated system into an existing sector: placed clear of
    every stored system's Hill sphere (`SpaceSector.add_system`), or at
    `position` when given, with its location text, and for a galaxy-placed
    sector its containment and nearest-systems rows (and its neighbors').

    Args:
        conn (Connection): An open connection.
        sector_id (int): The `sectors.id`.
        star_system (StarSystem): The generated system.
        system_config (SystemConfig): Its recipe.
        position (tuple, optional): Sector-local `(x, y, z)` in light-years;
            it must lie inside the sector.

    Returns:
        tuple: `(star_systems.id, (x, y, z))`.

    Raises:
        ValueError: No such sector, a position outside it, or no room left.
    """
    sector = sector_for_placement(conn, sector_id)
    if position is not None and not sector.contains(tuple(position)):
        raise ValueError(f"position {tuple(position)!r} is outside sector {sector_id}")
    entry = sector.add_system(star_system, position=position, system_config=system_config)
    system_id = insert_star_system(conn, star_system, system_config, sector_id=sector_id,
                                   position=entry.position, location=_location_for_entry(sector, entry))
    assign_uids(conn, sector_id=sector_id)
    if get_sector_galaxy_position(conn, sector_id) is not None:
        refresh_containment(conn, [sector_id])
        _add_sector_to_nearest(conn, sector_id)
    return system_id, entry.position


_SYSTEM_CONTENT_COLUMNS = (
    "system_config_id", "is_binary", "binary_configuration",
    "binary_separation_km", "binary_eccentricity", "binary_periapsis_km", "binary_apoapsis_km",
    "binary_type", "binary_temperature_k", "binary_radius_km",
    "binary_effective_mass_kg", "binary_effective_luminosity_w", "binary_age_gy", "binary_lifespan_gy",
    "binary_habitable_zone_inner_km", "binary_habitable_zone_outer_km",
    "binary_system_perimeter_km", "binary_heliosphere_radius_km",
    "binary_galactic_orbital_period_gy",
    "binary_galactic_orbital_phase_deg", "binary_galactic_min_update_interval_years",
    "binary_mutual_orbital_period_years", "binary_mutual_orbital_speed_kms",
    "binary_mutual_orbital_inclination_deg", "binary_mutual_orbital_ascending_node_deg",
    "binary_mutual_orbital_phase_deg", "binary_epoch_unix", "binary_next_update_due",
    "binary_mutual_min_update_interval_years", "binary_mutual_position_x_km", "binary_mutual_position_y_km", "binary_mutual_position_z_km",
    "binary_primary_position_x_km", "binary_primary_position_y_km", "binary_primary_position_z_km",
    "binary_secondary_position_x_km", "binary_secondary_position_y_km", "binary_secondary_position_z_km",
    "binary_secondary_mass_fraction",
    "binary_planetary_wobble_x_km", "binary_planetary_wobble_y_km", "binary_planetary_wobble_z_km",
    "system_flavor_text", "runaway_class", "runaway_speed_kms", "schema_version",
)
"""tuple: The `star_systems` columns that describe a system's generated
content (everything `insert_star_system` writes except its name and
placement), which `replace_star_system_content` swaps."""

_SYSTEM_CONTENT_TABLES = ("stars", "planets", "moons", "asteroid_belts", "comets")
"""tuple: The tables holding a system's bodies, keyed by `star_system_id`."""


def system_content_blockers(conn, star_system_id):
    """
    What stops `replace_star_system_content` from swapping a system's
    bodies without losing something: `{"facilities": n}` hosted on the
    system or its bodies (they would be deleted with the bodies), and
    `"bright_star": True` when the system was built around a pre-placed
    bright star (its star is fixed by the galaxy scatter).
    """
    facilities = conn.execute("SELECT COUNT(*) AS n FROM facilities WHERE star_system_id = ?",
                              (star_system_id,)).fetchone()["n"]
    bright = conn.execute("SELECT 1 FROM bright_stars WHERE star_system_id = ? LIMIT 1",
                          (star_system_id,)).fetchone() is not None
    return {"facilities": facilities, "bright_star": bright}


def _body_slots(conn, system_id):
    """
    `{(table, row id): slot}` for the bodies of a system: where each one sits in the system, which is what its ID
    follows when the system's content is drawn again (DB.75).

    A slot is `(kind, host, index)`: a star's is its role (`primary`, `secondary`, `single`, or `member` with its place
    among the members); a planet's or a belt's
    is the role of the star it orbits (`system` for a circumbinary one) and its place in that star's list of orbits; a
    moon's is its planet's slot and its place among the planet's moons; a comet's is the star it belongs to and its
    rank among the system's comets.
    """
    roles = {row["id"]: row["role"] for row in conn.execute(
        "SELECT id, role FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()}
    seen = collections.Counter()
    slots = {}
    for star_id, role in roles.items():    # ordered by id: a role held by more than one star (`member`) numbers them
        slots[("stars", star_id)] = ("star", role, seen[role])
        seen[role] += 1
    planet_slots = {}
    for table, kind in (("planets", "planet"), ("asteroid_belts", "belt")):
        for row in conn.execute(f"SELECT id, star_id, orbital_index FROM {table} WHERE star_system_id = ?",
                                (system_id,)).fetchall():
            slot = (kind, roles.get(row["star_id"], "system"), row["orbital_index"])
            slots[table, row["id"]] = slot
            if table == "planets":
                planet_slots[row["id"]] = slot
    for row in conn.execute("SELECT id, planet_id, orbital_index FROM moons WHERE star_system_id = ?",
                            (system_id,)).fetchall():
        slots["moons", row["id"]] = ("moon", planet_slots[row["planet_id"]], row["orbital_index"])
    comet_rank = {}
    for row in conn.execute("SELECT id, star_id FROM comets WHERE star_system_id = ? ORDER BY id",
                            (system_id,)).fetchall():
        host = roles.get(row["star_id"], "system")
        slots["comets", row["id"]] = ("comet", host, comet_rank.get(host, 0))
        comet_rank[host] = comet_rank.get(host, 0) + 1
    return slots


def _slot_uids(conn, system_id):
    """`{slot: uid bytes}` of the bodies of a system that have an ID."""
    out = {}
    for (table, row_id), slot in _body_slots(conn, system_id).items():
        row = conn.execute(f"SELECT uid FROM {table} WHERE id = ?", (row_id,)).fetchone()
        if row is not None and row["uid"] is not None:
            out[slot] = row["uid"]
    return out


def _carry_uids(conn, system_id, kept):
    """Gives the bodies of a system that have no ID yet the ID of the body that held their slot (`kept`, from
    `_slot_uids`): a body of the same kind in the same place is the same body to the outside."""
    for (table, row_id), slot in _body_slots(conn, system_id).items():
        if slot in kept:
            conn.execute(f"UPDATE {table} SET uid = ? WHERE id = ? AND uid IS NULL", (kept[slot], row_id))


def replace_star_system_content(conn, star_system_id, star_system, system_config):
    """
    Swaps a stored system's generated content -- its stars, planets,
    moons, asteroid belts, comets and the system-level star/binary fields
    -- for a newly generated system's, keeping the system's id, name,
    sector, position, location, and containment. The new
    bodies are renamed from the system's name. Facilities hosted on the
    old system or its bodies are deleted with them (check
    `system_content_blockers` first).

    Args:
        conn (Connection): An open connection, inside a transaction.
        star_system_id (int): The system to change.
        star_system (StarSystem): The new content.
        system_config (SystemConfig): Its recipe.

    Returns:
        bool: `False` if no such system exists.
    """
    row = conn.execute("SELECT name FROM star_systems WHERE id = ?", (star_system_id,)).fetchone()
    if row is None:
        return False
    storage_checksum.mark_system(conn, star_system_id)
    # Saved standalone first (the registry gives it a temporary name),
    # then its bodies are moved under the kept row and the stand-in row
    # dropped, which also drops its name reservation.
    temp_id = insert_star_system(conn, star_system, system_config)
    temp_name = conn.execute("SELECT name FROM star_systems WHERE id = ?", (temp_id,)).fetchone()["name"]
    keep_body_numbers(conn, star_system_id)  # a number the old bodies held is not given to another body
    kept_uids = _slot_uids(conn, star_system_id)
    conn.execute("DELETE FROM facilities WHERE star_system_id = ?", (star_system_id,))
    for table in ("moons", "planets", "asteroid_belts", "comets", "stars"):
        conn.execute(f"DELETE FROM {table} WHERE star_system_id = ?", (star_system_id,))
    for table in _SYSTEM_CONTENT_TABLES:
        conn.execute(f"UPDATE {table} SET star_system_id = ? WHERE star_system_id = ?", (star_system_id, temp_id))
    assignments = ", ".join(f"kept.{column} = made.{column}" for column in _SYSTEM_CONTENT_COLUMNS)
    conn.execute(
        f"UPDATE star_systems kept JOIN star_systems made ON made.id = ? SET {assignments},"
        " kept.next_update_due = NULL, kept.modified_at = CURRENT_TIMESTAMP(3) WHERE kept.id = ?",
        (temp_id, star_system_id),
    )
    conn.execute("DELETE FROM star_systems WHERE id = ?", (temp_id,))
    kept = conn.execute("SELECT sector_id, position_x_mpc, position_y_mpc, position_z_mpc FROM star_systems"
                        " WHERE id = ?", (star_system_id,)).fetchone()
    if kept["position_x_mpc"] is not None:  # the new star moves at its own speed from the kept place
        velocity = _system_velocity_kms(conn, star_system, kept["sector_id"], tuple(
            mpc_to_pc(kept[f"position_{axis}_mpc"]) for axis in "xyz"))
        conn.execute("UPDATE star_systems SET velocity_x_kms = ?, velocity_y_kms = ?, velocity_z_kms = ?,"
                     " modified_at = modified_at WHERE id = ?", (*velocity, star_system_id))
    _rename_bodies_with_prefix(conn, star_system_id, temp_name, row["name"])
    _carry_uids(conn, star_system_id, kept_uids)  # a body keeps the ID of the one that held its place (DB.75)
    assign_uids(conn, system_ids=[star_system_id])  # the others take the system's next body numbers
    return True


def get_sector_galaxy_position(conn, sector_id):
    """
    Reads back a sector's stored galaxy-frame placement (see `schema.sql`'s
    "v4" header note) -- used by `generation/run_galaxy.py`'s local-neighborhood mode
    to look up an existing sector's own center before enumerating its
    neighbors (`spherical.enumerate_sectors_within_radius`).

    Args:
        conn (Connection): An open, schema-initialized connection.
        sector_id (int): The `sectors.id` to look up.

    Returns:
        dict or None: A dict with keys `center_x_pc`, `center_y_pc`,
            `center_z_pc`, `galactic_radius_pc`, `shell_index`,
            `band_index`, `wedge_index`, or `None` if
            this sector has never been placed in a galaxy (the four
            galaxy-position columns NULL).

    Raises:
        ValueError: If no such `sectors` row exists.
    """
    row = conn.execute(
        """
        SELECT center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc,
               shell_index, band_index, wedge_index
        FROM sectors WHERE id = ?
        """,
        (sector_id,),
    ).fetchone()
    if row is None:
        raise ValueError(f"no sectors row with id {sector_id}")
    if row["center_x_pc"] is None:
        return None

    return dict(row)


def _galaxy_placement_from_sector_offset(galaxy_position, offset_ly):
    """
    Converts a phenomenon's sector-relative `(x, y, z)` offset in
    light-years (as `SpaceSector.add_phenomenon`/`SectorPhenomenonEntry`
    placed it -- the exact Hill-sphere-checked position a standalone black
    hole/neutron star was actually generated at) into an absolute
    galaxy-frame center in parsecs, given the owning sector's own stored
    `galaxy_position` -- see `schema.sql`'s "v21" header note.

    `offset_ly` is expressed along the sector's own local axes (local +X
    radially outward from the galactic axis, +Y along the longitude, +Z
    galactic north -- `frame.sector_orientation`), the same frame
    `SpaceSector.add_system` places a star system's position in, so it is
    rotated into the galaxy frame before being added to the sector's
    center -- the identical transform `planetgen/web/maps/starmap.py` applies to a
    star system's position at render time.

    Unlike `compute_phenomenon_placement` (an independent random jitter,
    for `phenomenonGen.py`'s own standalone `--sector-id` use, where no
    specific in-sector position was ever computed), this is a plain
    coordinate conversion: `insert_sector` already knows exactly where
    within the sector's cube each phenomenon landed, so it reuses that
    real position instead of re-randomizing one.

    Args:
        galaxy_position (dict or None): The owning sector's own galaxy-frame
            placement (`center_x_pc`/`center_y_pc`/`center_z_pc`), or `None`
            if that sector itself was never placed in the galaxy
            (`generation/run_sector.py`'s own standalone CLI).
        offset_ly (tuple): The phenomenon's `(x, y, z)` position in
            light-years, relative to the sector's own center, along that
            sector's own local cube axes.

    Returns:
        dict or None: `center_x_pc`/`center_y_pc`/`center_z_pc`/
            `galactic_radius_pc`, or `None` if `galaxy_position` is `None`.
    """
    if galaxy_position is None:
        return None

    center_pc = (
        galaxy_position["center_x_pc"], galaxy_position["center_y_pc"], galaxy_position["center_z_pc"],
    )
    x, y, z = local_to_galaxy_pc(center_pc, tuple(ly_to_pc(coordinate) for coordinate in offset_ly))
    return {
        "center_x_pc": x, "center_y_pc": y, "center_z_pc": z,
        "galactic_radius_pc": math.sqrt(x * x + y * y + z * z),
    }


def compute_phenomenon_placement(conn, sector_id):
    """
    Picks a galaxy-frame center for a nebula/asteroid field placed "at"
    `sector_id` (`phenomenonGen.py --sector-id`) -- see `schema.sql`'s
    "v18" header note for why these phenomena get their own galaxy-frame
    sphere rather than a sector-relative offset the way `star_systems`
    does.

    The center is a uniformly random point inside `sector_id`'s own grid
    cell (`spherical.SectorCell`; a sector with no grid address falls
    back to a `+/- edge_pc / 2` cube around its center) -- so the
    phenomenon's center always falls inside that sector, regardless of the phenomenon's own `radius_ly` (which may be
    far larger than the sector itself, and is free to spill into
    neighboring sectors -- see `db.query.phenomena_near_sector`, which finds
    those by real distance, not by this row's `sector_id`).

    Args:
        conn (Connection): An open, schema-initialized connection.
        sector_id (int): The `sectors.id` to place this phenomenon at.

    Returns:
        dict: `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc`.

    Raises:
        ValueError: If `sector_id` doesn't exist, or exists but has never
            been placed in the galaxy itself (`generation/run_galaxy.py`) -- there is
            no galaxy-frame position to jitter around in that case.
    """
    position = get_sector_galaxy_position(conn, sector_id)
    if position is None:
        raise ValueError(
            f"sector {sector_id} has no galaxy placement of its own -- see generation/run_galaxy.py"
        )
    edge_row = conn.execute("SELECT edge_mpc FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    edge_pc = mpc_to_pc(edge_row["edge_mpc"])

    center_pc = (position["center_x_pc"], position["center_y_pc"], position["center_z_pc"])
    if position["shell_index"] is not None:
        offset_pc = SectorCell.for_address(
            position["shell_index"], position["band_index"], position["wedge_index"], edge_pc).sample(random)
        center_x, center_y, center_z = local_to_galaxy_pc(center_pc, offset_pc)
    else:
        center_x, center_y, center_z = (c + draw.uniform(-edge_pc / 2, edge_pc / 2) for c in center_pc)
    return {
        "center_x_pc": center_x, "center_y_pc": center_y, "center_z_pc": center_z,
        "galactic_radius_pc": math.sqrt(center_x ** 2 + center_y ** 2 + center_z ** 2),
    }


def get_occupied_addresses(conn, shell_indices):
    """
    Returns every already-occupied `(shell_index, band_index,
    wedge_index)` address among the given shells -- used by
    `planetgen galaxy`'s batch and neighborhood modes to skip addresses a
    sector already exists at, in one query rather than one per candidate.

    Args:
        conn (Connection): An open, schema-initialized connection.
        shell_indices (iterable): Shell indices to check.

    Returns:
        set: Address tuples already present in `sectors`. Empty if
            `shell_indices` is empty.
    """
    shell_indices = sorted(set(shell_indices))
    if not shell_indices:
        return set()

    placeholders = ", ".join("?" for _ in shell_indices)
    rows = conn.execute(
        f"SELECT shell_index, band_index, wedge_index FROM sectors "
        f"WHERE shell_index IN ({placeholders}) AND wedge_index IS NOT NULL",
        tuple(shell_indices),
    ).fetchall()
    return {(row["shell_index"], row["band_index"], row["wedge_index"]) for row in rows}


def _sector_cell(row, edge_ly):
    """The `spherical.SectorCell` of a `sectors` row with a grid address (it
    needs `shell_index`, `band_index` and `wedge_index`), in light-years, else
    `None`."""
    address = (row["shell_index"], row["band_index"], row["wedge_index"])
    return None if None in address else SectorCell.for_address(*address, edge_ly)


def get_sector_id_at(conn, shell_index, band_index, wedge_index):
    """
    Looks up the `sectors.id` already generated at one grid address, if
    any -- used to check one address at a time.

    Returns:
        int or None: The existing `sectors.id`, or `None`.
    """
    row = conn.execute(
        "SELECT id FROM sectors WHERE shell_index = ? AND band_index = ? AND wedge_index = ?",
        (shell_index, band_index, wedge_index),
    ).fetchone()
    return row["id"] if row is not None else None


GalaxySkeletonInfo = namedtuple(
    "GalaxySkeletonInfo",
    ["shape", "edge_pc", "outer_shell_index", "expected_system_count_at_density_1", "galaxy_seed",
     "reach_shell_index", "cell_star_ceiling", "fill_luminosity_cap_sol", "core_min_luminosity_sol",
     "core_mass_limit_sol", "halo_min_luminosity_sol", "pinned_stellar_mass_solar", "pins"],
    defaults=(None, 65535, tuning.CELL_STAR_CEILING, tuning.FILL_LUMINOSITY_CAP_SOL, None, None,
              tuning.HALO_MIN_LUMINOSITY_SOL, None, None),
)
"""The galaxy's stored skeleton -- everything needed to recompute any
sector's exact position/density on demand (see schema.sql's "v8" header
note). `shape` is a `galaxy.density.GalaxyShape`; the other fields are
`galaxy_shape`'s own remaining columns, `galaxy_seed` (v51, GEN.39) the
16-byte seed or `None` for a galaxy planned before it, `reach_shell_index`
(v93, GEN.233) the outermost shell the grid has addresses for, `cell_star_ceiling` and `fill_luminosity_cap_sol` (v98,
GEN.281) the systems a cell may expect before it is filled only down to the cap, and that cap in solar luminosities,
`core_min_luminosity_sol` and `core_mass_limit_sol` (v101, DB.54) the core's scatter limits (`None`: follow the
baseline times `tuning.CORE_TIER_FACTOR`, see `galaxy.scatter_limits`) and `halo_min_luminosity_sol` the halo's,
`pinned_stellar_mass_solar` (v114, DB.81) the stellar mass the plan was pinned to (`None`: the seed drew it) and `pins`
the shape levers the plan was told, `{lever name: value}` (`galaxy_shape_pins`; `None` for a skeleton not read from the
database)."""


def save_galaxy_shape(shape: GalaxyShape, edge_pc, outer_shell_index,
                       expected_system_count_at_density_1, config=None, galaxy_seed=None,
                       reach_shell_index=65535, cell_star_ceiling=None, fill_luminosity_cap_sol=None,
                       core_min_luminosity_sol=None, core_mass_limit_sol=None, halo_min_luminosity_sol=None,
                       pinned_stellar_mass_solar=None, pins=None):
    """
    Replaces the galaxy's singleton `galaxy_shape` row -- there is exactly
    one galaxy, so this always overwrites whatever was there before rather
    than inserting a second row (`galaxyPlan.py` calls this once per full
    skeleton (re)build).

    Args:
        shape (galaxy.density.GalaxyShape): The galaxy's shape parameters.
        edge_pc (float): The sector edge length this skeleton was built
                         at, parsecs.
        outer_shell_index (int): The last shell with any qualifying
            content (`planetgen plan`'s own discovered galaxy edge).
        reach_shell_index (int): The outermost shell the grid has
            addresses for (GEN.233); the plan has checked the object ID
            can name it (`spherical.check_reach`).
        cell_star_ceiling (float, optional): The systems a cell may expect before it is filled only down to the
            fill luminosity cap (GEN.281); `None` is `tuning.CELL_STAR_CEILING`.
        fill_luminosity_cap_sol (float, optional): That cap, solar luminosities; `None` is
            `tuning.FILL_LUMINOSITY_CAP_SOL`.
        core_min_luminosity_sol (float, optional): The core's luminosity floor (DB.54); `None` follows the baseline.
        core_mass_limit_sol (float, optional): The core's mass limit (DB.54); `None` follows the baseline.
        halo_min_luminosity_sol (float, optional): The halo's luminosity floor (DB.54); `None` is
            `tuning.HALO_MIN_LUMINOSITY_SOL`.
        pinned_stellar_mass_solar (float, optional): The stellar mass the plan pinned (DB.81), solar masses, which a
            loaded shape draws its mass chain and structure at; `None` leaves the draw to the seed.
        pins (dict, optional): The shape levers the plan was told, `{lever name: value}` in the units of
            `galaxy.dials` (`shape_draw.pins_of`), kept in `galaxy_shape_pins`; replaces the ones stored.
        expected_system_count_at_density_1 (float): See
            `galaxy.skeleton.expected_system_count_at_density_1`.
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.
        galaxy_seed (bytes, optional): The galaxy's 16-byte seed (GEN.39).
            `None` keeps the one stored, or draws one
            (`galaxySeed.new_seed`) when there is none. A seed that isn't
            the stored one also records the running code's version key
            and versions (DB.6).

    Returns:
        bytes: The galaxy seed now stored.
    """
    conn = get_connection(config)
    try:
        with conn:
            stored_seed = get_galaxy_seed(conn)
            if galaxy_seed is None:
                galaxy_seed = stored_seed or galaxySeed.new_seed()
            if cell_star_ceiling is None:
                cell_star_ceiling = tuning.CELL_STAR_CEILING
            if fill_luminosity_cap_sol is None:
                fill_luminosity_cap_sol = tuning.FILL_LUMINOSITY_CAP_SOL
            if halo_min_luminosity_sol is None:
                halo_min_luminosity_sol = tuning.HALO_MIN_LUMINOSITY_SOL
            family_values = (shape.family, None if shape.bar_present is None else int(shape.bar_present),
                             shape.bar_axis_ratio_b, shape.bar_axis_ratio_c, shape.bar_angle_rad,
                             shape.thick_disk_amplitude, shape.spheroid_amplitude, shape.spheroid_sersic_n,
                             shape.spheroid_effective_radius_pc, shape.spheroid_axis_ratio_b,
                             shape.spheroid_axis_ratio_c, *shape.spheroid_orientation_rad, shape.body_clump_share)
            conn.execute(
                """
                INSERT INTO galaxy_shape (
                    id, disk_scale_length_pc, disk_scale_height_pc,
                    bulge_scale_radius_pc, bulge_amplitude, arm_count,
                    pitch_angle_rad, arm_amplitude, arm_level, nuclear_amplitude, spiral_reference_radius_pc,
                    spiral_reference_angle_rad, k_norm, edge_pc,
                    expected_system_count_at_density_1, outer_shell_index, galaxy_seed, reach_shell_index,
                    cell_star_ceiling, fill_luminosity_cap_sol, core_min_luminosity_sol, core_mass_limit_sol,
                    halo_min_luminosity_sol, pinned_stellar_mass_solar,
                    family, bar_present, bar_axis_ratio_b, bar_axis_ratio_c, bar_angle_rad, thick_disk_amplitude,
                    spheroid_amplitude, spheroid_sersic_n, spheroid_effective_radius_pc, spheroid_axis_ratio_b,
                    spheroid_axis_ratio_c, spheroid_orientation_1_rad, spheroid_orientation_2_rad,
                    spheroid_orientation_3_rad, body_clump_share
                ) VALUES (1, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                          ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON DUPLICATE KEY UPDATE
                    disk_scale_length_pc = VALUES(disk_scale_length_pc),
                    disk_scale_height_pc = VALUES(disk_scale_height_pc),
                    bulge_scale_radius_pc = VALUES(bulge_scale_radius_pc),
                    bulge_amplitude = VALUES(bulge_amplitude),
                    arm_count = VALUES(arm_count),
                    pitch_angle_rad = VALUES(pitch_angle_rad),
                    arm_amplitude = VALUES(arm_amplitude),
                    arm_level = VALUES(arm_level),
                    nuclear_amplitude = VALUES(nuclear_amplitude),
                    spiral_reference_radius_pc = VALUES(spiral_reference_radius_pc),
                    spiral_reference_angle_rad = VALUES(spiral_reference_angle_rad),
                    k_norm = VALUES(k_norm),
                    edge_pc = VALUES(edge_pc),
                    expected_system_count_at_density_1 = VALUES(expected_system_count_at_density_1),
                    outer_shell_index = VALUES(outer_shell_index),
                    reach_shell_index = ?,
                    cell_star_ceiling = ?,
                    fill_luminosity_cap_sol = ?,
                    core_min_luminosity_sol = ?,
                    core_mass_limit_sol = ?,
                    halo_min_luminosity_sol = ?,
                    pinned_stellar_mass_solar = ?,
                    bright_star_min_luminosity_sol = NULL,
                    bright_star_seed = NULL,
                    galaxy_seed = ?,
                    family = ?, bar_present = ?, bar_axis_ratio_b = ?, bar_axis_ratio_c = ?, bar_angle_rad = ?,
                    thick_disk_amplitude = ?, spheroid_amplitude = ?, spheroid_sersic_n = ?,
                    spheroid_effective_radius_pc = ?, spheroid_axis_ratio_b = ?, spheroid_axis_ratio_c = ?,
                    spheroid_orientation_1_rad = ?, spheroid_orientation_2_rad = ?, spheroid_orientation_3_rad = ?,
                    body_clump_share = ?
                """,
                (
                    shape.disk_scale_length_pc, shape.disk_scale_height_pc,
                    shape.bulge_scale_radius_pc, shape.bulge_amplitude, shape.arm_count,
                    shape.pitch_angle_rad, shape.arm_amplitude, shape.arm_level, shape.nuclear_amplitude,
                    shape.spiral_reference_radius_pc,
                    shape.spiral_reference_angle_rad, shape.k_norm, edge_pc,
                    expected_system_count_at_density_1, outer_shell_index, bytes(galaxy_seed), reach_shell_index,
                    cell_star_ceiling, fill_luminosity_cap_sol, core_min_luminosity_sol, core_mass_limit_sol,
                    halo_min_luminosity_sol, pinned_stellar_mass_solar, *family_values,
                    reach_shell_index, cell_star_ceiling, fill_luminosity_cap_sol, core_min_luminosity_sol,
                    core_mass_limit_sol, halo_min_luminosity_sol, pinned_stellar_mass_solar, bytes(galaxy_seed),
                    *family_values,
                ),
            )
            conn.execute("DELETE FROM galaxy_shape_pins")
            if pins:
                conn.executemany("INSERT INTO galaxy_shape_pins (name, value) VALUES (?, ?)",
                                 [(name, float(value)) for name, value in sorted(pins.items())])
            if stored_seed != bytes(galaxy_seed):
                # DB.6: a new seed is a new galaxy -- record the code making it.
                made_by = versionKey.current()
                battery = generatorBattery.running_digest()
                conn.execute(
                    "UPDATE galaxy_shape SET version_key = ?, planetgen_version = ?, python_version = ?, platform = ?,"
                    " generator_epoch = ?, fp_spec = ?, battery_digest = ?, battery_cases = ? WHERE id = 1",
                    (made_by["version_key"], made_by["planetgen_version"], made_by["python_version"],
                     made_by["platform"], generatorEpoch.generator_epoch(), fingerprintEncoding.FP_SPEC,
                     battery.overall, json.dumps(generatorBattery.cases_of(battery), sort_keys=True)),
                )
    finally:
        conn.close()
    return bytes(galaxy_seed)


def has_galaxy_sectors(conn):
    """Whether any sector has been placed in the galaxy grid."""
    return conn.execute("SELECT 1 FROM sectors WHERE shell_index IS NOT NULL LIMIT 1").fetchone() is not None


GalaxyMaker = namedtuple(
    "GalaxyMaker", ["version_key", "planetgen_version", "python_version", "platform", "generator_epoch", "fp_spec",
                    "battery_digest", "battery_cases"], defaults=(None, None, None, None))
"""What made the galaxy: `galaxy_shape`'s version columns (v52, DB.6), the generator epoch, the fingerprint format and
the generator battery's overall digest with its `{case: full digest}` (v131, OPS.58); the last four are `None` for a
galaxy planned before v131."""


def get_galaxy_maker(conn):
    """The version key and versions that made the galaxy, or `None` when it
    has never been planned (or was planned before v52)."""
    row = conn.execute("SELECT * FROM galaxy_shape WHERE id = 1").fetchone()
    if row is None or row.get("version_key") is None:
        return None
    cases = row["battery_cases"]
    return GalaxyMaker(row["version_key"], row["planetgen_version"], row["python_version"], row["platform"],
                       row["generator_epoch"], row["fp_spec"], row["battery_digest"],
                       json.loads(cases) if isinstance(cases, (str, bytes)) else cases)


def sector_versions(conn):
    """
    The code that generated the galaxy's sectors (v67, DB.7; from the run each sector points at since v121, DB.16): one
    dict per distinct `version_key` (`None` for sectors with no recorded version) with `planetgen_version`,
    `python_version`, `platform`, `generator_epoch` (OPS.58; `None` for runs from before v131) and the `sectors` count,
    the largest group first.
    """
    rows = conn.execute(
        "SELECT r.version_key AS version_key, MIN(r.planetgen_version) AS planetgen_version,"
        " MIN(r.python_version) AS python_version, MIN(r.platform) AS platform,"
        " MIN(r.generator_epoch) AS generator_epoch, COUNT(*) AS sectors"
        " FROM sectors s LEFT JOIN generation_runs r ON r.id = s.run_id GROUP BY r.version_key"
        " ORDER BY sectors DESC, version_key").fetchall()
    return [dict(row) for row in rows]


GENERATION_ARGUMENT_LENGTH = 1024
"""int: The longest command-line argument `generation_run_arguments` keeps."""


def start_generation_run(conn, command, arguments, run_seed=None):
    """
    Records a run that changes the galaxy (v52, DB.6) as started, with the
    running code's version key and versions and the galaxy seed it runs
    against, and commits.

    Args:
        conn (Connection): An open connection to the galaxy database.
        command (str): The `planetgen` subcommand.
        arguments (list): Its command line, without the --mysql-* and
            --debug options (one `generation_run_arguments` row each).
        run_seed (int, optional): The run's own 128-bit seed.

    Returns:
        int: The run's `generation_runs.id`.
    """
    made_by = versionKey.current()
    seed_bytes = None if run_seed is None else int(run_seed).to_bytes(galaxySeed.SEED_BYTES, "big")
    run_id = conn.execute(
        "INSERT INTO generation_runs (command, run_seed, galaxy_seed, version_key, planetgen_version,"
        " python_version, platform, generator_epoch) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        (command, seed_bytes, get_galaxy_seed(conn), made_by["version_key"],
         made_by["planetgen_version"], made_by["python_version"], made_by["platform"],
         generatorEpoch.generator_epoch()),
    ).lastrowid
    for position, argument in enumerate(arguments):
        conn.execute("INSERT INTO generation_run_arguments (run_id, position, value) VALUES (?, ?, ?)",
                     (run_id, position, str(argument)[:GENERATION_ARGUMENT_LENGTH]))
    conn.commit()
    return run_id


def newest_sectors_run_start(conn):
    """When the newest `sectors` or `galaxy` run started (`generation_runs.started_at`, the database's clock), the
    default `--since` of `planetgen link`, `backfill` and `settle`; `None` when no such run is recorded."""
    row = conn.execute("SELECT started_at FROM generation_runs WHERE command IN ('sectors', 'galaxy') "
                       "ORDER BY id DESC LIMIT 1").fetchone()
    return None if row is None else row["started_at"]


def finish_generation_run(conn, run_id, outcome):
    """Records how run `run_id` ended (`ok`, `failed` or `interrupted`) and
    when, and commits."""
    conn.execute("UPDATE generation_runs SET finished_at = CURRENT_TIMESTAMP(3), outcome = ? WHERE id = ?",
                 (outcome, run_id))
    conn.commit()


def get_galaxy_seed(conn):
    """The galaxy's 16-byte seed (v51, GEN.39), or `None` when it has
    never been planned (or was planned before v51)."""
    row = conn.execute("SELECT galaxy_seed FROM galaxy_shape WHERE id = 1").fetchone()
    if row is None or row["galaxy_seed"] is None:
        return None
    return bytes(row["galaxy_seed"])


@functools.lru_cache(maxsize=8)
def _field_of(found, edge_pc, expected_at_density_1):
    return componentField.ComponentField(found, edge_pc, expected_at_density_1)


_FIELD_BY_SIGNATURE = collections.OrderedDict()
"""OrderedDict: The last few `ComponentField`s by the signature of the rows they were made from, so a read that finds the
table unchanged does not load it again (GEN.244: a galaxy holds up to `tuning.OPEN_CLUSTER_MAX_COUNT` open clusters)."""
_FIELD_BY_SIGNATURE_LIMIT = 8


def _component_field(conn, edge_pc, expected_at_density_1):
    """The galaxy's `ComponentField` (GEN.160, GEN.237) from its stored clusters and satellites, or `None` when it has
    none (or the table is not there yet: a migration reads the skeleton before v97). Fields made from equal rows are one
    object, so the density model's caches by shape keep hitting. The rows are read only when the table's signature (the
    database, the row count, the first and last id and the mass sum) is not the one of a field already made."""
    if not conn.execute(
            "SELECT 1 FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'galaxy_components'"
    ).fetchone():
        return None
    row = conn.execute(
        "SELECT DATABASE() AS db, COUNT(*) AS n, MIN(id) AS first_id, MAX(id) AS last_id, SUM(mass_solar) AS mass"
        " FROM galaxy_components").fetchone()
    if not row["n"]:
        return None
    key = (row["db"], row["n"], row["first_id"], row["last_id"], row["mass"], float(edge_pc),
           float(expected_at_density_1))
    field = _FIELD_BY_SIGNATURE.get(key)
    if field is not None:
        with contextlib.suppress(KeyError):       # another thread dropped it meanwhile
            _FIELD_BY_SIGNATURE.move_to_end(key)
        log.trace(f"component field: cache hit for {row['n']} stored components (database {row['db']})", area="cache")
        return field
    log.trace(f"component field: cache miss, reading {row['n']} stored components (database {row['db']})", area="cache")
    found = tuple(c for c in get_galaxy_components(conn) if c.kind in componentField.KINDS_IN_THE_GATE)
    if not found:
        return None
    field = _field_of(tuple(dataclasses.replace(c, id=None) for c in found), float(edge_pc),
                      float(expected_at_density_1))
    _FIELD_BY_SIGNATURE[key] = field
    while len(_FIELD_BY_SIGNATURE) > _FIELD_BY_SIGNATURE_LIMIT:
        _FIELD_BY_SIGNATURE.popitem(last=False)
    return field


def get_galaxy_shape(conn):
    """
    Reads back the galaxy's stored skeleton parameters.

    Args:
        conn (Connection): An open, schema-initialized connection.

    Returns:
        GalaxySkeletonInfo or None: `None` if `galaxyPlan.py` has never
            been run against this database (the `galaxy_shape` singleton
            row doesn't exist yet).
    """
    # Every column, not a list: older migration steps read the skeleton
    # before v51 adds `galaxy_seed`.
    row = conn.execute("SELECT * FROM galaxy_shape WHERE id = 1").fetchone()
    if row is None:
        return None

    shape = GalaxyShape(
        disk_scale_length_pc=row["disk_scale_length_pc"],
        disk_scale_height_pc=row["disk_scale_height_pc"],
        bulge_scale_radius_pc=row["bulge_scale_radius_pc"],
        bulge_amplitude=row["bulge_amplitude"],
        arm_count=row["arm_count"],
        pitch_angle_rad=row["pitch_angle_rad"],
        arm_amplitude=row["arm_amplitude"],
        arm_level=row["arm_level"],
        nuclear_amplitude=row["nuclear_amplitude"],
        spiral_reference_radius_pc=row["spiral_reference_radius_pc"],
        spiral_reference_angle_rad=row["spiral_reference_angle_rad"],
        k_norm=row["k_norm"],
    )
    seed = None if row.get("galaxy_seed") is None else bytes(row["galaxy_seed"])
    pinned_mass = row.get("pinned_stellar_mass_solar")
    pins = get_galaxy_shape_pins(conn)
    # The thick disk, bar, stellar halo and motions are not columns: the seed draws them again (GEN.532, GEN.533), at
    # the mass and halo the plan pinned when it did (DB.81, GEN.606).
    shape = shape._replace(components=_component_field(conn, row["edge_pc"], row["expected_system_count_at_density_1"]),
                           structure=shape_draw.structure_for(seed, pinned_mass, pins.get("stellar_halo_share")))
    if "family" in row:      # DB.66 (GEN.379): the family and component numbers; a NULL bar or thick disk number is the draw's
        shape = with_family_columns(shape, row)
    if seed is not None:      # a galaxy with no seed has no masses: its disk's size stands for them
        chain = massChain.draw_mass_chain(seed, pinned_mass, halo_mass_pinned(seed, pinned_mass, pins))
        shape = shape._replace(masses=(chain.stellar_mass_solar, chain.halo_mass_solar))
    if pins.get("richness") is not None:     # the Census dial (GEN.607): the mean system mass follows the pinned richness
        shape = with_system_richness(shape, pins["richness"])
    return GalaxySkeletonInfo(
        shape=shape, edge_pc=row["edge_pc"], outer_shell_index=row["outer_shell_index"],
        expected_system_count_at_density_1=row["expected_system_count_at_density_1"],
        galaxy_seed=seed,
        reach_shell_index=row.get("reach_shell_index", 65535),
        cell_star_ceiling=row.get("cell_star_ceiling", tuning.CELL_STAR_CEILING),
        fill_luminosity_cap_sol=row.get("fill_luminosity_cap_sol", tuning.FILL_LUMINOSITY_CAP_SOL),
        core_min_luminosity_sol=row.get("core_min_luminosity_sol"),
        core_mass_limit_sol=row.get("core_mass_limit_sol"),
        halo_min_luminosity_sol=row.get("halo_min_luminosity_sol", tuning.HALO_MIN_LUMINOSITY_SOL),
        pinned_stellar_mass_solar=pinned_mass, pins=pins,
    )


def halo_mass_pinned(seed, pinned_stellar_mass_solar, pins):
    """The dark halo mass, solar masses, the stored `pins` fix (the Halo dial's `halo_ratio` times the stellar mass the
    galaxy has), or `None` when the seed draws it."""
    ratio = pins.get("halo_ratio")
    if ratio is None:
        return None
    return shape_draw.halo_mass_of(massChain.draw_mass_chain(seed, pinned_stellar_mass_solar).stellar_mass_solar, ratio)


def get_galaxy_shape_pins(conn):
    """The shape levers the stored plan was told (`galaxy_shape_pins`, DB.81): `{lever name: value}`, the values in the
    units of `galaxy.dials`; empty when everything was drawn from the seed."""
    return {row["name"]: row["value"] for row in conn.execute("SELECT name, value FROM galaxy_shape_pins").fetchall()}


def replace_galaxy_shells(shell_extents, config=None, conn=None):
    """
    Replaces the galaxy's stored outline wholesale: every `galaxy_shell` row
    (each shell's highest band and its cell count) -- a full skeleton build
    always produces the whole outline in one pass, so there is no partial
    update (matches `save_galaxy_shape`).

    Args:
        shell_extents (iterable): `(shell_index, high_band_index)` pairs
            (`galaxy.skeleton.build_shell_extents`), any order.
        config (MySQLConfig, optional): Connection parameters. Defaults
            to `DEFAULT_MYSQL_CONFIG`. Ignored when `conn` is given.
        conn (Connection, optional): Write through this connection, inside
            the caller's own transaction (the v91 migration does), instead
            of opening and committing a new one.
    """
    from planetgen.galaxy.skeleton import shell_candidate_counts

    rows = shell_candidate_counts(sorted(shell_extents))

    def _write(c):
        c.execute("DELETE FROM galaxy_shell")
        c.executemany(
            "INSERT INTO galaxy_shell (shell_index, high_band_index, cell_count) VALUES (?, ?, ?)", rows,
        )

    if conn is not None:
        _write(conn)
        return
    own = get_connection(config)
    try:
        with own:
            _write(own)
    finally:
        own.close()


def get_galaxy_bounds(conn):
    """
    The galaxy's stored outline as a `galaxy.skeleton.GalaxyBounds`, the
    object every generation path checks an address against before
    generating anything there.

    Returns:
        GalaxyBounds or None: `None` if the skeleton was never built (no
            `galaxy_shape` row).
    """
    from planetgen.galaxy.skeleton import GalaxyBounds

    skeleton = get_galaxy_shape(conn)
    if skeleton is None:
        return None
    return GalaxyBounds(get_galaxy_shells(conn), skeleton.edge_pc, skeleton.reach_shell_index)


def get_galaxy_shells(conn):
    """
    The galaxy's whole stored outline, innermost shell first.

    Returns:
        list[tuple]: `(shell_index, high_band_index, cell_count)` triples;
            empty if the skeleton was never built.
    """
    rows = conn.execute(
        "SELECT shell_index, high_band_index, cell_count FROM galaxy_shell ORDER BY shell_index"
    ).fetchall()
    return [(row["shell_index"], row["high_band_index"], row["cell_count"]) for row in rows]


PHENOMENON_UID_TABLES = kind_tables.tables("phenomenon")
"""tuple: The phenomenon tables a sector holds, each with a `uid` (v58)."""

def _sector_address_of(galaxy_position):
    address = None if galaxy_position is None else tuple(
        galaxy_position.get(key) for key in ("shell_index", "band_index", "wedge_index"))
    return address if address is not None and None not in address else None


def _uid_issuer_for_sector(conn, sector_id, galaxy_position):
    """The `_UidIssuer` for a sector just inserted (its `uid` is set by
    `set_sector_uid`), or `None` when the sector has no grid address to be
    born in: its rows then take run-time IDs from `assign_uids`."""
    address = _sector_address_of(galaxy_position)
    return None if address is None else _UidIssuer(sector_id, address)


def set_sector_uid(conn, sector_id):
    """Writes a sector's own `uid` when it has none: its designation, or for a sector with no grid address
    `galaxyUid.unplaced_sector_uid` (API.23: every sector has a public ID)."""
    row = conn.execute("SELECT shell_index, band_index, wedge_index, uid FROM sectors WHERE id = ?",
                       (sector_id,)).fetchone()
    if row is None or row["uid"] is not None:
        return
    address = (row["shell_index"], row["band_index"], row["wedge_index"])
    uid = galaxyUid.sector_uid(*address) if None not in address else galaxyUid.unplaced_sector_uid(sector_id)
    conn.execute("UPDATE sectors SET uid = ?, parity_dirty = 1, modified_at = modified_at WHERE id = ?", (uid, sector_id))


_UID_ISSUED_TABLES = frozenset(("star_systems", "stars", "planets", "moons", "asteroid_belts", "comets",
                                *PHENOMENON_UID_TABLES))
"""frozenset: Tables whose new rows `_UidIssuer` gives a `uid` as they are inserted."""

BODY_UID_TABLES = ("stars", "planets", "moons", "asteroid_belts", "comets")
"""tuple: The tables whose rows are bodies born in a system (body numbers 1 and up, one counter across them)."""

NO_BIRTH_SECTOR = (0, 0, 0)
"""tuple: The birth sector `(shell, band, wedge)` of an object saved with no sector of its own. Its IDs are run-time
births, drawn from the same `id_counters` row as any other run-time birth in that sector."""

RUNTIME_COUNTER = "sector"
BODY_COUNTER = "system"
"""str: The `id_counters.kind`s: run-time serials per birth sector, body numbers per system."""


class _UidIssuer:
    """
    GEN.171: the object IDs (`galaxy/object_uid.py`) of a sector's rows,
    worked out as the rows are inserted. The sector's systems and
    phenomena take generated serials 0, 1, 2 ... in the order they are
    inserted, which is the order the fill generates them in, so the same
    galaxy seed gives the same IDs. A system's stars, planets, moons, belts
    and comets take body numbers 1, 2, 3 ... in insertion order, one counter
    across the five kinds. `Connection.execute` asks for one per INSERT into
    `_UID_ISSUED_TABLES` while `conn.uid_issuer` is set (`insert_sector`).

    `complete` turns false when a row cannot be given one here (its system
    was not issued by this sector, or it would be the 4,096th body of a
    system); `insert_sector` then runs `assign_uids` for the rows with none.
    """

    def __init__(self, sector_id, address):
        self.sector_id = sector_id
        self.address = tuple(address)
        self.complete = True
        self._next_serial = 0
        self.given_nebula_uid = None
        """bytes or `None`: The next nebula's own ID, when it is a field cloud born in another sector (GEN.176)."""
        self._system_uid = {}
        self._rogue_uid = {}
        self._body_count = {}

    def _top_level(self):
        serial = self._next_serial
        if serial >= (1 << objectUid.LAYOUT.serial_count_bits):
            self.complete = False
            return None
        self._next_serial = serial + 1
        return objectUid.pack(*self.address, objectUid.SERIAL_GENERATED, serial)

    def issue(self, table, columns, params):
        """The `uid` bytes for a row about to be inserted, or `None` when
        the row gets none here."""
        row = dict(zip(columns, params))
        if table == "nebulae" and self.given_nebula_uid is not None:
            uid, self.given_nebula_uid = self.given_nebula_uid, None
            return uid
        if table == "star_systems" or table in PHENOMENON_UID_TABLES:
            if row.get("sector_id") != self.sector_id:
                return None
            value = self._top_level()
            if value is None:
                return None
            if table == "star_systems":
                self._system_uid[row["id"]] = value
            elif table == "rogue_planets":
                self._rogue_uid[row["id"]] = value
            return objectUid.to_bytes(value)
        system_id = row.get("star_system_id")
        parent = self._system_uid.get(system_id)
        if system_id is None and row.get("rogue_planet_id") is not None:
            # DB.73: a moon of a free planet numbers its bodies under the rogue planet's own ID.
            system_id = ("rogue_planets", row["rogue_planet_id"])
            parent = self._rogue_uid.get(row["rogue_planet_id"])
        number = self._body_count.get(system_id, 0) + 1
        if parent is None or number > objectUid.max_body(objectUid.LAYOUT):
            self.complete = False
            return None
        self._body_count[system_id] = number
        return objectUid.to_bytes(objectUid.with_body(parent, number))


def _reserve_counter(config, kind, scope, size, floor):
    """Reserves `size` numbers of the `id_counters` row (`kind`, `scope`) and
    returns the first. Like `_reserve_id_block`: its own autocommitted
    connection, so writers never queue behind each other's transactions and a
    rolled-back save only leaves a gap. `floor` is the lowest number allowed:
    the row only moves up, so a number is never handed out twice, whatever
    happens to the object that took it (DB.20)."""
    raw = _get_pool(config).connect()
    try:
        cur = raw.cursor()
        cur.execute(
            "INSERT INTO id_counters (kind, scope, next_value) VALUES (%s, %s, LAST_INSERT_ID(%s + %s))"
            " ON DUPLICATE KEY UPDATE next_value = LAST_INSERT_ID(GREATEST(next_value, %s) + %s)",
            (kind, scope, floor, size, floor, size),
        )
        cur.execute("SELECT LAST_INSERT_ID() AS end_id")
        end = cur.fetchone()["end_id"]
        raw.commit()
    except Exception:
        raw.rollback()
        raise
    finally:
        raw.close()
    return end - size


def _runtime_uids(conn, address, count):
    """`count` new top-level IDs born at run time in the sector at `address`
    (`None`: `NO_BIRTH_SECTOR`): serials of the run-time kind from the sector's counter."""
    address = NO_BIRTH_SECTOR if address is None else tuple(address)
    scope = objectUid.sector_bytes(*address)
    first = _reserve_counter(conn._config, RUNTIME_COUNTER, scope, count, 0)
    limit = 1 << objectUid.LAYOUT.serial_count_bits
    if first + count > limit:
        raise ValueError(f"sector {address} has used all of its run-time object IDs")
    return [objectUid.to_bytes(objectUid.pack(*address, objectUid.SERIAL_RUNTIME, first + n)) for n in range(count)]


def _next_body_numbers(conn, system_uid, system_id, count):
    """`count` new body numbers for system `system_id`, whose own ID is
    `system_uid` (bytes): the counter, started past every body the system
    holds (rows written before it had a counter row)."""
    highest = 0
    for table in BODY_UID_TABLES:
        for row in conn.execute(f"SELECT uid FROM {table} WHERE star_system_id = ?", (system_id,)).fetchall():
            if row["uid"] is not None:
                highest = max(highest, objectUid.body_of(objectUid.from_bytes(row["uid"])))
    first = _reserve_counter(conn._config, BODY_COUNTER, bytes(system_uid), count, highest + 1)
    if first + count - 1 > objectUid.max_body(objectUid.LAYOUT):
        raise ValueError(f"system {objectUid.format_id(objectUid.from_bytes(system_uid))} has used all of its body numbers")
    return first


def keep_body_numbers(conn, system_id):
    """Makes sure the system's body counter is past every body it holds now, so a body deleted afterwards
    never has its number given to another (the counter row is made on first use otherwise)."""
    row = conn.execute("SELECT uid FROM star_systems WHERE id = ?", (system_id,)).fetchone()
    if row is not None and row["uid"] is not None:
        _next_body_numbers(conn, row["uid"], system_id, 0)


def _sector_birth_address(conn, sector_id):
    """The `(shell, band, wedge)` a sector's objects are born at, or `None`."""
    if sector_id is None:
        return None
    row = conn.execute("SELECT shell_index, band_index, wedge_index FROM sectors WHERE id = ?",
                       (sector_id,)).fetchone()
    if row is None:
        return None
    address = (row["shell_index"], row["band_index"], row["wedge_index"])
    return None if None in address else address


def assign_uids(conn, sector_id=None, system_ids=(), phenomenon=None):
    """
    Writes the object ID (`galaxy/object_uid.py`) of every row below that
    has none yet: a row that has one keeps it, so saving again never
    changes an ID. The rows of a sector saved by `insert_sector` already
    have theirs (`_UidIssuer`); this is for rows that arrive later or that
    the issuer could not number, and they are run-time births: a serial
    from the sector's counter, or for a body the system's next body number.

    Args:
        conn (Connection): Part of the caller's transaction, with the
            rows already written (not mid-`batched`).
        sector_id (int, optional): A sector: its own ID (its designation),
            then its systems and their contents, then its phenomena.
        system_ids (iterable of int): Systems and everything in them.
        phenomenon (tuple, optional): `(table, id)` of one phenomenon.
    """
    if sector_id is not None:
        _assign_sector_uids(conn, sector_id)
    system_ids = list(system_ids)
    if system_ids:
        _assign_system_uids(conn, system_ids)
    if phenomenon is not None:
        table, row_id = phenomenon
        _assign_top_level_uids(conn, table, "id", [row_id])
        if table == "rogue_planets":
            _assign_rogue_moon_uids(conn, [row_id])


def _assign_sector_uids(conn, sector_id):
    set_sector_uid(conn, sector_id)
    ids = [row["id"] for row in conn.execute("SELECT id FROM star_systems WHERE sector_id = ? ORDER BY id",
                                             (sector_id,)).fetchall()]
    for first in range(0, len(ids), _UID_SYSTEM_BATCH):
        _assign_system_uids(conn, ids[first:first + _UID_SYSTEM_BATCH])
    for table in PHENOMENON_UID_TABLES:
        _assign_top_level_uids(conn, table, "sector_id", [sector_id])
    rogues = [row["id"] for row in conn.execute("SELECT id FROM rogue_planets WHERE sector_id = ? ORDER BY id",
                                                (sector_id,)).fetchall()]
    for first in range(0, len(rogues), _UID_SYSTEM_BATCH):
        _assign_rogue_moon_uids(conn, rogues[first:first + _UID_SYSTEM_BATCH])


_UID_SYSTEM_BATCH = 200


def _assign_rogue_moon_uids(conn, rogue_planet_ids):
    """Fills the NULL `uid`s of the moons of these rogue planets (DB.73) with the next body numbers under each rogue
    planet's own ID, as a system's bodies are numbered under the system's."""
    for rogue_id in rogue_planet_ids:
        rogue = conn.execute("SELECT uid FROM rogue_planets WHERE id = ?", (rogue_id,)).fetchone()
        if rogue is None or rogue["uid"] is None:
            continue
        rows = conn.execute("SELECT id FROM moons WHERE rogue_planet_id = ? AND uid IS NULL ORDER BY id",
                            (rogue_id,)).fetchall()
        if not rows:
            continue
        highest = 0
        for taken in conn.execute("SELECT uid FROM moons WHERE rogue_planet_id = ? AND uid IS NOT NULL",
                                  (rogue_id,)).fetchall():
            highest = max(highest, objectUid.body_of(objectUid.from_bytes(taken["uid"])))
        first = _reserve_counter(conn._config, BODY_COUNTER, bytes(rogue["uid"]), len(rows), highest + 1)
        if first + len(rows) - 1 > objectUid.max_body(objectUid.LAYOUT):
            raise ValueError("a rogue planet has used all of its body numbers")
        parent = objectUid.from_bytes(rogue["uid"])
        _update_by_id(conn, "moons", ("uid",),
                      [(row["id"], objectUid.to_bytes(objectUid.with_body(parent, first + n)))
                       for n, row in enumerate(rows)], touch=True)


def _assign_top_level_uids(conn, table, column, values):
    """Run-time IDs for the rows of `table` with none, selected by `column` in `values`
    (a sector's phenomena, or one phenomenon by its row id), in row id order."""
    marks = ", ".join("?" * len(values))
    rows = conn.execute(f"SELECT id, sector_id FROM {table} WHERE {column} IN ({marks}) AND uid IS NULL ORDER BY id",
                        tuple(values)).fetchall()
    by_sector = {}
    for row in rows:
        by_sector.setdefault(row["sector_id"], []).append(row["id"])
    for sector_id, row_ids in by_sector.items():
        uids = _runtime_uids(conn, _sector_birth_address(conn, sector_id), len(row_ids))
        _update_by_id(conn, table, ("uid",), list(zip(row_ids, uids)))


def _assign_system_uids(conn, system_ids):
    """Fills the NULL `uid`s of these systems, then of the stars, planets, moons, belts and comets in
    them: a system as a run-time birth in its sector, a body with the system's next body number."""
    ids = tuple(system_ids)
    marks = ", ".join("?" * len(ids))
    _assign_top_level_uids(conn, "star_systems", "id", ids)
    for system in conn.execute(f"SELECT id, uid FROM star_systems WHERE id IN ({marks}) ORDER BY id", ids).fetchall():
        for table in BODY_UID_TABLES:
            rows = conn.execute(f"SELECT id FROM {table} WHERE star_system_id = ? AND uid IS NULL ORDER BY id",
                                (system["id"],)).fetchall()
            if not rows:
                continue
            first = _next_body_numbers(conn, system["uid"], system["id"], len(rows))
            parent = objectUid.from_bytes(system["uid"])
            _update_by_id(conn, table, ("uid",),
                          [(row["id"], objectUid.to_bytes(objectUid.with_body(parent, first + n)))
                           for n, row in enumerate(rows)], touch=True)  # these tables have no `modified_at`


def save_system(star_system: StarSystem, system_config: SystemConfig, config=None) -> int:
    """
    Opens the database and persists a single, standalone `StarSystem` (no
    sector -- `sector_id`/`position` are left `None`) in one transaction.
    The single-system counterpart to `save_sector`, for `generation/run_system.py`
    (which, unlike `generation/run_sector.py`, generates one system with no natural
    sector placement of its own). Retried like `save_sector` on a
    deadlock or lock wait timeout (several one-off systems saved at once).

    Args:
        star_system (StarSystem): The generated system to persist.
        system_config (SystemConfig): The config it was generated from.
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.

    Returns:
        int: The new `star_systems.id`.
    """
    names = [(star_system, star_system.name)]
    def insert(conn):
        system_id = insert_star_system(conn, star_system, system_config)
        assign_uids(conn, system_ids=[system_id])
        return system_id

    return _save_with_retries(config, names, insert)


RETRYABLE_ERRORS = (1213, 1205)
"""tuple: MySQL error codes a whole sector save is retried on (PERF.14):
a deadlock, and a lock wait timeout."""

SECTOR_SAVE_ATTEMPTS = 8


def _sector_names(sector):
    """Every name `insert_sector` may change, so a retried save can start
    over from the generated names."""
    objects = [sector] + [entry.star_system for entry in sector.entries]
    for entry in sector.phenomena:
        objects.append(entry.phenomenon)
        core = getattr(entry.phenomenon, "compact_remnant", None)
        if core is not None:
            objects.append(core)
    return [(obj, obj.name) for obj in objects]


def _neighbor_lock_name(conn):
    """The named lock `insert_sector` holds while it links a placed
    sector to its neighbors: one per database (names are server-wide,
    at most 64 characters)."""
    database = conn._config.database if conn._config is not None else ""
    return f"planetgen.neighbors.{database}"[:64]


def save_sector(sector: SpaceSector, config=None, galaxy_position=None, link_neighbors=True, progress=None,
                run_id=None) -> int:
    """
    Opens the database and persists a full `SpaceSector` to it in one
    transaction.

    The transaction runs at READ COMMITTED, so its reads take no gap
    locks, and is retried from the start, with the sector's generated
    names restored, when MySQL reports a deadlock or lock wait timeout
    (`RETRYABLE_ERRORS`) -- the case several sector writers at once
    (PERF.8) can still meet on shared rows such as a neighbor's nearest
    systems (PERF.14).

    Args:
        sector (SpaceSector): The sector to persist.
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.
        galaxy_position (dict, optional): This sector's galaxy-frame
            placement -- see `insert_sector`'s docstring. `None` (the
            default) for a sector never placed in a galaxy.
        link_neighbors (bool): Link the sector to its neighbours now
            (containment, nearest systems). `False` leaves that to the
            caller's later `link_sector_neighbors` (a `galaxy` run, PERF.45).
        progress (optional): A step (`generation.steps`) that is told, as
            each system and phenomenon is inserted, how far the save has
            got (PERF.50); it restarts at 0 when a save is retried.
        run_id (int, optional): The `generation_runs.id` of the run saving the sector (see `insert_sector`).

    Returns:
        int: The new `sectors.id`.
    """
    def insert(conn):
        conn.save_progress = progress
        if progress is not None:
            progress.update(completed=0)
        return insert_sector(conn, sector, galaxy_position=galaxy_position, link_neighbors=link_neighbors,
                             run_id=run_id)

    return _save_with_retries(config, _sector_names(sector), insert)


_RETRY_JITTER = random.Random()
"""The retry pause's own stream (not a generation draw), so waiting never
moves the unit's draw stream (GEN.39, GEN.56)."""


def _save_with_retries(config, names, insert):
    """
    Runs `insert(conn)` in one READ COMMITTED transaction and returns its
    result, starting over (up to `SECTOR_SAVE_ATTEMPTS` times, with
    `names` -- `(object, generated name)` pairs -- restored) on a deadlock
    or lock wait timeout. Shared by `save_sector` and `save_system`. When
    it finally fails, the names are restored too: the rollback dropped the
    reservations behind any renaming, and a later save of the same objects
    must start from their generated names (TEST.15).

    READ COMMITTED matters for the name registry too, not only for gap
    locks: a writer that waited on another's registry row must then see
    that writer's new system to rename it "Alpha ..." (a REPEATABLE READ
    snapshot taken before the wait doesn't).
    """
    # A retry draws what the first try drew (GEN.39): the sector's numbers
    # can't depend on whether another worker's save got in its way.
    state = draw.getstate()
    for attempt in range(1, SECTOR_SAVE_ATTEMPTS + 1):
        conn = get_connection(config)
        claim = None
        try:
            conn.commit()       # a connection that just checked the schema has a transaction open, and MySQL refuses the next line inside one
            conn.execute("SET TRANSACTION ISOLATION LEVEL READ COMMITTED")
            # PERF.49: names are claimed on a connection of their own that
            # commits at once, so other writers wait for a claim, not for
            # this whole save.
            claim = get_connection(config)
            claim.commit()
            claim.execute("SET TRANSACTION ISOLATION LEVEL READ COMMITTED")
            conn.claim_connection = claim
            draw.setstate(state)
            with conn:
                return insert(conn)
        except Exception as exc:
            if claim is not None and conn.name_claims:
                _release_name_claims(claim, conn.name_claims)
            retry = (isinstance(exc, pymysql.err.OperationalError) and exc.args
                     and exc.args[0] in RETRYABLE_ERRORS and attempt < SECTOR_SAVE_ATTEMPTS)
            for obj, name in names:
                obj.name = name
            if not retry:
                raise
            log.debug(f"Save hit MySQL error {exc.args[0]} ({exc.args[1] if len(exc.args) > 1 else ''}); "
                      f"retrying ({attempt}/{SECTOR_SAVE_ATTEMPTS - 1}).")
            time.sleep(_RETRY_JITTER.uniform(0.05, 0.25) * attempt)
        finally:
            conn.claim_connection = None
            if claim is not None:
                claim.close()
            conn.close()


# ---------------------------------------------------------------------------
# Read path -- reconstructs live objects from database rows. See the module
# docstring for why this maps rows directly to each leaf class's own
# `from_dict` shape rather than routing everything through a single nested
# `StarSystem.from_dict` call.
# ---------------------------------------------------------------------------

def _tristate_from_db(value):
    """
    Inverts `_tristate`: converts the schema's nullable `INTEGER` `0`/`1`/
    `NULL` representation back to a `SystemConfig` tri-state value.

    Args:
        value (int or None): The stored `0`/`1`/`NULL` value.

    Returns:
        bool or None: The tri-state value.
    """
    return None if value is None else bool(value)


def load_system_config(conn, config_id) -> SystemConfig:
    """
    Reconstructs a `SystemConfig` from a `system_configs` row (plus its
    `system_config_slots` child rows, if any).

    Args:
        conn (Connection): An open, schema-initialized connection.
        config_id (int): The `system_configs.id` to load.

    Returns:
        SystemConfig: The reconstructed config.

    Raises:
        ValueError: If no such row exists.
    """
    row = conn.execute("SELECT * FROM system_configs WHERE id = ?", (config_id,)).fetchone()
    if row is None:
        raise ValueError(f"no system_configs row with id {config_id}")

    slot_rows = conn.execute(
        """
        SELECT orbit_index, type, planet_class, moons
        FROM system_config_slots WHERE config_id = ? ORDER BY orbit_index
        """,
        (config_id,),
    ).fetchall()

    slots = None
    if slot_rows:
        slots = [None] * (max(r["orbit_index"] for r in slot_rows) + 1)
        for r in slot_rows:
            slots[r["orbit_index"]] = {
                "type": r["type"], "planet_class": r["planet_class"], "moons": r["moons"],
            }

    return SystemConfig.from_dict({
        "habitable_world": _tristate_from_db(row["habitable_world"]),
        "asteroid_belt": _tristate_from_db(row["asteroid_belt"]),
        "comets": _tristate_from_db(row["comets"]),
        "large_star": _tristate_from_db(row["large_star"]),
        "moons": _tristate_from_db(row["moons"]),
        "max_planets": _tristate_from_db(row["max_planets"]),
        "planets": _tristate_from_db(row["planets"]),
        "star_type": row["star_type"],
        "name": row["name"],
        "age": row["age"],
        "intelligent_life": _tristate_from_db(row["intelligent_life"]),
        "binary_system": _tristate_from_db(row["binary_system"]),
        "wide_binary": _tristate_from_db(row["wide_binary"]),
        "num_orbits": row["num_orbits"],
        "stars": row["stars"],
        "max_stars": row["max_stars"],
        "slots": slots,
        "prevalence": {feature: row[f"prevalence_{feature}"] for feature in PREVALENCE_FEATURES
                       if row[f"prevalence_{feature}"]},
    })


def _star_row_to_dict(row, speed_kms):
    """Maps a `stars` row to `Star.from_dict`'s expected dict shape: `STAR_FIELDS` inverts every conversion
    `insert_star` applies; the galactic orbital speed is not stored (the caller derives it)."""
    data = object_rows.load_dict(STAR_FIELDS, row)
    data["galactic_orbital_speed_kms"] = speed_kms
    return data


def _black_hole_row_to_dict(star_row, black_hole_row, speed_kms):
    """
    Maps a `stars` row plus its owning `black_holes` satellite row to
    `BlackHole.from_dict`'s expected dict shape -- `_star_row_to_dict`'s
    base fields (every `Star.SERIALIZABLE_FIELDS` entry) layered with the
    black hole's own extra fields (`mass_solar`, `event_horizon_radius_km`,
    `spin`, `has_accretion_disk`), which need no unit conversion (already
    stored in the same units `BlackHole`'s own attributes use).

    Args:
        star_row: The owning `stars` row (`role = 'single'`).
        black_hole_row: The `black_holes` row (`star_id` -> `star_row["id"]`).

    Returns:
        dict: In the shape `BlackHole.to_dict()` produces.
    """
    data = _star_row_to_dict(star_row, speed_kms)
    data.update(object_rows.load_dict(BLACK_HOLE_FIELDS, black_hole_row, only=_BLACK_HOLE_OWN_COLUMNS))
    return data


def _neutron_star_row_to_dict(star_row, neutron_star_row, speed_kms):
    """
    Maps a `stars` row plus its owning `neutron_stars` satellite row to
    `NeutronStar.from_dict`'s expected dict shape -- see
    `_black_hole_row_to_dict`'s identical role/reasoning.

    Args:
        star_row: The owning `stars` row (`role = 'single'`).
        neutron_star_row: The `neutron_stars` row (`star_id` -> `star_row["id"]`).

    Returns:
        dict: In the shape `NeutronStar.to_dict()` produces.
    """
    data = _star_row_to_dict(star_row, speed_kms)
    data.update(object_rows.load_dict(NEUTRON_STAR_FIELDS, neutron_star_row, only=_NEUTRON_STAR_OWN_COLUMNS))
    return data


def _load_single_star(conn, star_system_id, system_config, speed_kms):
    """
    Loads the `stars` row for a single (non-binary) system and reconstructs
    the correct Python class from it: a `BlackHole`/`NeutronStar` when a
    matching `black_holes`/`neutron_stars` satellite row exists for it
    (`phenomenonGen.py --anchor-system`, see `generation/phenomena/compact_remnant.py`'s module
    docstring), or a plain `Star` otherwise.

    Without this dispatch, `Star.from_dict` alone would silently reconstruct
    an anchored compact remnant as an ordinary `Star` carrying the literal
    `yerkes_class` marker `'BH'`/`'NS'` -- losing every remnant-specific
    field (`mass_solar`, `event_horizon_radius_km`, `spin`, ... /
    `spin_period_ms`, `magnetic_field_gauss`, `pulsar_type`, ...) and
    rendering via `Star`'s own paragraph methods instead of the remnant's.

    Args:
        conn (Connection): An open, schema-initialized connection.
        star_system_id (int): The owning `star_systems.id`.
        system_config (SystemConfig): The system's shared config.
        speed_kms (float): The system's galactic speed, the length of its stored vector (`_speed_kms`).

    Returns:
        Star, BlackHole, or NeutronStar: The reconstructed single star.
    """
    star_row = conn.execute(
        "SELECT * FROM stars WHERE star_system_id = ? AND role = 'single'", (star_system_id,)
    ).fetchone()

    black_hole_row = conn.execute("SELECT * FROM black_holes WHERE star_id = ?", (star_row["id"],)).fetchone()
    neutron_star_row = None
    if black_hole_row is None:
        neutron_star_row = conn.execute("SELECT * FROM neutron_stars WHERE star_id = ?", (star_row["id"],)).fetchone()

    if black_hole_row is not None:
        star = BlackHole.from_dict(_black_hole_row_to_dict(star_row, black_hole_row, speed_kms), system_config)
    elif neutron_star_row is not None:
        star = NeutronStar.from_dict(_neutron_star_row_to_dict(star_row, neutron_star_row, speed_kms), system_config)
    else:
        star = Star.from_dict(_star_row_to_dict(star_row, speed_kms), system_config)
    star.db_id = star_row["id"]
    return star


def _binary_proxy_row_to_dict(star_system_row, primary_dict, secondary_dict):
    """Maps a `star_systems` row's `binary_*` columns to
    `BinaryStarProxy.from_dict`'s expected dict shape, inverting every unit
    conversion `_proxy_only_binary_fields`/`_mutual_orbit_fields_from_proxy`
    apply."""
    row = star_system_row
    return {
        "name": row["name"],
        "type": row["binary_type"],
        "temperature": row["binary_temperature_k"],
        "radius": row["binary_radius_km"],
        "age": row["binary_age_gy"],
        "lifespan": row["binary_lifespan_gy"],
        "habitable_zone": [
            _km_to_au(row["binary_habitable_zone_inner_km"]),
            _km_to_au(row["binary_habitable_zone_outer_km"]),
        ],
        "system_perimeter": _km_to_au(row["binary_system_perimeter_km"]),
        "heliosphere_radius": _km_to_au(row["binary_heliosphere_radius_km"]),
        "galactic_orbital_speed_kms": _speed_kms(row),
        "galactic_orbital_period_gy": row["binary_galactic_orbital_period_gy"],
        "galactic_orbital_phase_deg": row["binary_galactic_orbital_phase_deg"],
        "galactic_min_update_interval_years": row["binary_galactic_min_update_interval_years"],
        "binary_mutual_orbital_period_years": row["binary_mutual_orbital_period_years"],
        "binary_mutual_orbital_speed_kms": row["binary_mutual_orbital_speed_kms"],
        "binary_mutual_orbital_inclination_deg": row["binary_mutual_orbital_inclination_deg"],
        "binary_mutual_orbital_ascending_node_deg": row["binary_mutual_orbital_ascending_node_deg"],
        "binary_mutual_orbital_phase_deg": row["binary_mutual_orbital_phase_deg"],
        "binary_mutual_min_update_interval_years": row["binary_mutual_min_update_interval_years"],
        "binary_mutual_position_x": _km_to_au(row["binary_mutual_position_x_km"]),
        "binary_mutual_position_y": _km_to_au(row["binary_mutual_position_y_km"]),
        "binary_mutual_position_z": _km_to_au(row["binary_mutual_position_z_km"]),
        "binary_primary_position_x": _km_to_au(row["binary_primary_position_x_km"]),
        "binary_primary_position_y": _km_to_au(row["binary_primary_position_y_km"]),
        "binary_primary_position_z": _km_to_au(row["binary_primary_position_z_km"]),
        "binary_secondary_position_x": _km_to_au(row["binary_secondary_position_x_km"]),
        "binary_secondary_position_y": _km_to_au(row["binary_secondary_position_y_km"]),
        "binary_secondary_position_z": _km_to_au(row["binary_secondary_position_z_km"]),
        "binary_secondary_mass_fraction": row["binary_secondary_mass_fraction"],
        "_binary_separation_au": _km_to_au(row["binary_separation_km"]),
        "_effective_mass": row["binary_effective_mass_kg"],
        "_effective_luminosity": row["binary_effective_luminosity_w"],
        "primary": primary_dict,
        "secondary": secondary_dict,
    }


def _wide_binary_row_to_dict(star_system_row):
    """
    Maps a `star_systems` row's reused `binary_separation_km`/
    `binary_mutual_*` columns plus the new `binary_eccentricity`/
    `binary_periapsis_km`/`binary_apoapsis_km` columns to
    `WideBinaryPair.from_dict`'s expected dict shape (its
    `SERIALIZABLE_FIELDS`), inverting every unit conversion
    `_mutual_orbit_fields_from_wide_binary` applies. Only ever called for a
    `binary_configuration == 'wide'` row -- see `load_star_system`.
    """
    row = star_system_row
    return {
        "separation_au": _km_to_au(row["binary_separation_km"]),
        "eccentricity": row["binary_eccentricity"],
        "period_years": row["binary_mutual_orbital_period_years"],
        "speed_kms": row["binary_mutual_orbital_speed_kms"],
        "periapsis_au": _km_to_au(row["binary_periapsis_km"]),
        "apoapsis_au": _km_to_au(row["binary_apoapsis_km"]),
        "inclination_deg": row["binary_mutual_orbital_inclination_deg"],
        "ascending_node_deg": row["binary_mutual_orbital_ascending_node_deg"],
        "phase_deg": row["binary_mutual_orbital_phase_deg"],
        "min_update_interval_years": row["binary_mutual_min_update_interval_years"],
        "position_x_au": _km_to_au(row["binary_mutual_position_x_km"]),
        "position_y_au": _km_to_au(row["binary_mutual_position_y_km"]),
        "position_z_au": _km_to_au(row["binary_mutual_position_z_km"]),
        "primary_position_x_au": _km_to_au(row["binary_primary_position_x_km"]),
        "primary_position_y_au": _km_to_au(row["binary_primary_position_y_km"]),
        "primary_position_z_au": _km_to_au(row["binary_primary_position_z_km"]),
        "secondary_position_x_au": _km_to_au(row["binary_secondary_position_x_km"]),
        "secondary_position_y_au": _km_to_au(row["binary_secondary_position_y_km"]),
        "secondary_position_z_au": _km_to_au(row["binary_secondary_position_z_km"]),
        "secondary_mass_fraction": row["binary_secondary_mass_fraction"],
    }


def _planet_or_moon_row_to_dict(conn, row, is_moon):
    """
    Maps a `planets` or `moons` row (identical column shape -- see
    `schema.sql`'s "v2" header note) to `Planet.from_dict`'s expected dict
    shape, inverting every unit conversion `insert_planet`/`insert_moon`
    apply and pulling in the row's reflection-spectrum child rows.

    Args:
        conn (Connection): An open, schema-initialized connection.
        row (dict): The `planets` or `moons` row.
        is_moon (bool): Which pair of child tables to query
                        (`planet_*`/`moon_*`) and id column to filter by.

    Returns:
        dict: In the shape `Planet.to_dict()` produces (`moons` left as
             `[]` -- the caller fills it in for a top-level planet).
    """
    table_prefix = "moon" if is_moon else "planet"
    id_column = "moon_id" if is_moon else "planet_id"

    spectrum_rows = conn.execute(
        f"SELECT spectrum_type, value FROM {table_prefix}_reflection_spectrum "
        f"WHERE {id_column} = ? ORDER BY position",
        (row["id"],),
    ).fetchall()
    visible = [r["value"] for r in spectrum_rows if r["spectrum_type"] == "visible"]
    non_visible = [r["value"] for r in spectrum_rows if r["spectrum_type"] == "non_visible"]

    layer_rows = conn.execute(
        f"SELECT {', '.join(INTERIOR_LAYER_COLUMNS)} FROM {table_prefix}_interior_layers "
        f"WHERE {id_column} = ? ORDER BY layer_index",
        (row["id"],),
    ).fetchall()

    data = object_rows.load_dict(PLANET_BODY_FIELDS if not is_moon else BODY_FIELDS, row)
    data.update({
        atmosphere.GAS_FIELD: body_gases(conn, [row["id"]], is_moon).get(row["id"], {}),
        life_classes.SCORES_FIELD: body_class_scores(conn, [row["id"]], is_moon).get(row["id"]),
        life_classes.LIVING_FIELD: dominant_first(body_living(conn, [row["id"]], is_moon).get(row["id"]), row["life_class"], row["life_class_variant"]),
        "is_moon": is_moon,
        "reflection_spectrum_visible": visible or None,
        "reflection_spectrum_non_visible": non_visible or None,
        giant_interior.LAYERS_ATTRIBUTE: [{name: layer[name] for name in INTERIOR_LAYER_COLUMNS} for layer in layer_rows] or None,
        "moons": [],
    })
    return data


def _belt_row_to_dict(row, composition_pairs):
    """Maps an `asteroid_belts` row (plus its `asteroid_belt_composition`
    child rows) to `AsteroidBelt.from_dict`'s expected dict shape."""
    data = object_rows.load_dict(BELT_FIELDS, row)
    data["body_type"] = "a"
    data["composition"] = [[pair["component"], pair["concentration"]] for pair in composition_pairs]
    return data


def _comet_row_to_dict(row, composition_rows):
    """Maps a `comets` row (plus its `comet_composition` child rows) to
    `Comet.from_dict`'s expected dict shape."""
    data = object_rows.load_dict(COMET_FIELDS, row)
    data["composition"] = [r["component"] for r in composition_rows]
    return data


def _runaway_origin(conn, star_system_id):
    """How the runaway star of a system came to move so fast (`supernova` or `dynamical`, GEN.329), from the flags of its
    first star row; `None` for a runaway saved without one."""
    row = conn.execute("SELECT flags FROM stars WHERE star_system_id = ? AND role IN ('primary', 'single')"
                       " ORDER BY id LIMIT 1", (star_system_id,)).fetchone()
    flags = 0 if row is None or row["flags"] is None else row["flags"]
    for origin, bit in multistar.RUNAWAY_ORIGIN_FLAGS.items():
        if flags & bit:
            return origin
    return None


def load_star_system(conn, star_system_id) -> StarSystem:
    """
    Reconstructs a full `StarSystem` -- config, star(s), every planet/
    moon/asteroid belt it contains (in original orbital order), and every
    comet it contains -- from a `star_systems` row and its related rows.

    Each reconstructed star/planet/moon/belt/comet also carries `db_id`,
    its own row's `id`, so a caller rendering per-body text
    (`describe.model.split_by_body`) can match it back to the
    rows `db.query.system_detail` returns.

    Args:
        conn (Connection): An open, schema-initialized connection.
        star_system_id (int): The `star_systems.id` to load.

    Returns:
        StarSystem: The reconstructed system.

    Raises:
        ValueError: If no such row exists.
    """
    row = conn.execute("SELECT * FROM star_systems WHERE id = ?", (star_system_id,)).fetchone()
    if row is None:
        raise ValueError(f"no star_systems row with id {star_system_id}")

    system_config = load_system_config(conn, row["system_config_id"])
    speed_kms = _speed_kms(row)

    # binary_configuration is the authoritative discriminator (see
    # schema.sql's "v15" header note); a row written before that column
    # existed has it NULL, so fall back to the pre-v15 meaning of
    # is_binary (which only ever meant a 'close'/P-type pair back then).
    binary_configuration = row["binary_configuration"]
    if binary_configuration is None and row["is_binary"]:
        binary_configuration = "close"

    secondary_star = None
    wide_binary = None

    if binary_configuration == "close":
        primary_row = conn.execute(
            "SELECT * FROM stars WHERE star_system_id = ? AND role = 'primary'", (star_system_id,)
        ).fetchone()
        secondary_row = conn.execute(
            "SELECT * FROM stars WHERE star_system_id = ? AND role = 'secondary'", (star_system_id,)
        ).fetchone()
        proxy_data = _binary_proxy_row_to_dict(
            row, _star_row_to_dict(primary_row, speed_kms), _star_row_to_dict(secondary_row, speed_kms)
        )
        star = BinaryStarProxy.from_dict(proxy_data, system_config)
    elif binary_configuration == "wide":
        primary_row = conn.execute(
            "SELECT * FROM stars WHERE star_system_id = ? AND role = 'primary'", (star_system_id,)
        ).fetchone()
        secondary_row = conn.execute(
            "SELECT * FROM stars WHERE star_system_id = ? AND role = 'secondary'", (star_system_id,)
        ).fetchone()
        star = Star.from_dict(_star_row_to_dict(primary_row, speed_kms), system_config)
        secondary_star = Star.from_dict(_star_row_to_dict(secondary_row, speed_kms), system_config)
        wide_binary = WideBinaryPair.from_dict(
            _wide_binary_row_to_dict(row), system_config, star, secondary_star
        )
    else:
        star = _load_single_star(conn, star_system_id, system_config, speed_kms)

    # The page is titled with the system's name (`star_systems.name`). A
    # single star, or a close pair's proxy, shares it -- a close pair's
    # proxy already takes it (`_binary_proxy_row_to_dict`). A wide pair's
    # `star` is its primary, which keeps its own `stars.name` (see
    # names/bodies.py).
    if binary_configuration != "wide":
        star.name = row["name"]

    system = object.__new__(StarSystem)
    system.surrounding_cloud = surrounding_cloud(conn, row)
    system._name = row["name"]
    system.star_words = None
    system.system_config = system_config
    system.star = star
    system.binary_type = binary_configuration
    system.wide_binary = wide_binary
    system.orbit_nodes = orbit_trees.read_nodes(conn, star_system_id)    # GEN.311
    system.dark_members = orbit_trees.read_dark_members(conn, star_system_id)    # GEN.1013

    # v20: only meaningful (non-NULL) for a 'close' pair -- see
    # schema.sql's "v20" header note and StarSystem.__init__'s own comment.
    system.binary_planetary_wobble_x = (
        _km_to_au(row["binary_planetary_wobble_x_km"])
        if row["binary_planetary_wobble_x_km"] is not None else 0.0
    )
    system.binary_planetary_wobble_y = (
        _km_to_au(row["binary_planetary_wobble_y_km"])
        if row["binary_planetary_wobble_y_km"] is not None else 0.0
    )
    system.binary_planetary_wobble_z = (
        _km_to_au(row["binary_planetary_wobble_z_km"])
        if row["binary_planetary_wobble_z_km"] is not None else 0.0
    )

    if binary_configuration == "close":
        # `primary_star`/`secondary_star` keep generation order (the
        # `stars.role` rows), which can differ from the proxy's own
        # heavier-first `_primary`/`_secondary` -- `StarSystem.__str__`
        # renders the per-star sections in this order.
        if secondary_row["mass_kg"] > primary_row["mass_kg"]:
            system.primary_star, system.secondary_star = star._secondary, star._primary
        else:
            system.primary_star, system.secondary_star = star._primary, star._secondary
        system.primary_star.db_id = primary_row["id"]
        system.secondary_star.db_id = secondary_row["id"]
        system.stars = [system.primary_star, system.secondary_star]
    elif binary_configuration == "wide":
        system.primary_star = star
        system.secondary_star = secondary_star
        star.db_id = primary_row["id"]
        secondary_star.db_id = secondary_row["id"]
        system.stars = [star, secondary_star]
    else:
        system.primary_star = star
        system.stars = [star]

    member_rows = conn.execute(
        "SELECT * FROM stars WHERE star_system_id = ? AND role = 'member' ORDER BY id", (star_system_id,)).fetchall()
    if member_rows:    # GEN.630: the stars beyond the pair, and where each hangs in the tree
        members = [Star.from_dict(_star_row_to_dict(member_row, speed_kms), system_config) for member_row in member_rows]
        for member, member_row in zip(members, member_rows):
            member.db_id = member_row["id"]
        system.member_stars = tuple(members)
        system.stars = system.stars + members
    if member_rows or system.dark_members:    # (GEN.1013: a tree with a dark member says where its stars hang too)
        places = (primary_row, secondary_row, *member_rows) if binary_configuration in ("close", "wide") else (
            conn.execute("SELECT parent_node_id, slot FROM stars WHERE star_system_id = ? AND role = 'single'",
                         (star_system_id,)).fetchone(),)
        system.star_slots = tuple((place["parent_node_id"], place["slot"]) for place in places)

    planet_rows = conn.execute(
        "SELECT * FROM planets WHERE star_system_id = ? ORDER BY orbital_index", (star_system_id,)
    ).fetchall()
    belt_rows = conn.execute(
        "SELECT * FROM asteroid_belts WHERE star_system_id = ? ORDER BY orbital_index", (star_system_id,)
    ).fetchall()
    comet_rows = conn.execute(
        "SELECT * FROM comets WHERE star_system_id = ? ORDER BY id", (star_system_id,)
    ).fetchall()

    def _build_comet_list(comet_rows_subset):
        # Comets have no orbital_index (see insert_comet's own docstring),
        # so unlike _build_object_list there's no interleaved order to
        # restore -- just reconstruct each one.
        comets = []
        for r in comet_rows_subset:
            comp_rows = conn.execute(
                "SELECT component FROM comet_composition WHERE comet_id = ? ORDER BY position",
                (r["id"],),
            ).fetchall()
            comet = Comet.from_dict(_comet_row_to_dict(r, comp_rows), system_config)
            comet.db_id = r["id"]
            comets.append(comet)
        return comets

    def _build_object_list(planet_rows_subset, belt_rows_subset, owning_star):
        # planets/belts owned by the same star share one orbital_index
        # space (see insert_star_system's enumerate over
        # star_system.planets/secondary_planets) -- merge and re-sort by
        # it to restore that original interleaved order.
        combined = [("p", r) for r in planet_rows_subset] + [("b", r) for r in belt_rows_subset]
        combined.sort(key=lambda item: item[1]["orbital_index"])

        objects = []
        for kind, r in combined:
            if kind == "b":
                comp_rows = conn.execute(
                    "SELECT component, concentration FROM asteroid_belt_composition "
                    "WHERE belt_id = ? ORDER BY position",
                    (r["id"],),
                ).fetchall()
                belt = AsteroidBelt.from_dict(_belt_row_to_dict(r, comp_rows), system_config)
                belt.db_id = r["id"]
                objects.append(belt)
            else:
                planet_data = _planet_or_moon_row_to_dict(conn, r, is_moon=False)
                moon_rows = conn.execute(
                    "SELECT * FROM moons WHERE planet_id = ? ORDER BY orbital_index", (r["id"],)
                ).fetchall()
                planet_data["moons"] = [_planet_or_moon_row_to_dict(conn, mr, is_moon=True) for mr in moon_rows]
                planet = Planet.from_dict(planet_data, owning_star, system_config)
                planet.db_id = r["id"]
                for moon, moon_row in zip(planet.moons, moon_rows):
                    moon.db_id = moon_row["id"]
                objects.append(planet)
        return objects

    if binary_configuration == "wide":
        # star_id disambiguates the two stars' independent orbital-index
        # spaces (both restart at 0 -- see insert_star_system) -- group by
        # it before restoring each star's own orbital order.
        primary_db_id, secondary_db_id = primary_row["id"], secondary_row["id"]
        system.planets = _build_object_list(
            [r for r in planet_rows if r["star_id"] == primary_db_id],
            [r for r in belt_rows if r["star_id"] == primary_db_id],
            star,
        )
        system.secondary_planets = _build_object_list(
            [r for r in planet_rows if r["star_id"] == secondary_db_id],
            [r for r in belt_rows if r["star_id"] == secondary_db_id],
            secondary_star,
        )
        system.comets = _build_comet_list([r for r in comet_rows if r["star_id"] == primary_db_id])
        system.secondary_comets = _build_comet_list([r for r in comet_rows if r["star_id"] == secondary_db_id])
    else:
        member_star_ids = {member.db_id for member in system.member_stars}
        system.planets = _build_object_list([r for r in planet_rows if r["star_id"] not in member_star_ids],
                                            [r for r in belt_rows if r["star_id"] not in member_star_ids], star)
        system.secondary_planets = []
        system.comets = _build_comet_list([r for r in comet_rows if r["star_id"] not in member_star_ids])
        system.secondary_comets = []
    # GEN.332: a member star's planets, each its own orbital-index space
    system.member_planets = [
        _build_object_list([r for r in planet_rows if r["star_id"] == member.db_id],
                           [r for r in belt_rows if r["star_id"] == member.db_id], member)
        for member in system.member_stars]
    system.set_stability_limits()

    system.system_flavor_text = row["system_flavor_text"]
    system.runaway_class = row["runaway_class"]
    system.runaway_speed_kms = row["runaway_speed_kms"]
    system.runaway_origin = _runaway_origin(conn, star_system_id) if row["runaway_class"] == "runaway" else None
    system.planet_count, system.belt_count, system.moon_count = system.count_objects()
    system.hab_count, system.m_count = system.count_habitable()
    system.habitable_class_count = system.count_habitable_class_worlds()
    system.habitable_request = system.habitable_request_status()
    system.comet_count = system.count_comets()

    return system


def load_sector(conn, sector_id) -> SpaceSector:
    """
    Reconstructs a full `SpaceSector` -- every system it contains, with its
    placement -- from a `sectors` row and its related rows.

    Args:
        conn (Connection): An open, schema-initialized connection.
        sector_id (int): The `sectors.id` to load.

    Returns:
        SpaceSector: The reconstructed sector.

    Raises:
        ValueError: If no such row exists.
    """
    row = conn.execute("SELECT * FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    if row is None:
        raise ValueError(f"no sectors row with id {sector_id}")

    edge_ly = milliparsecs_to_ly(row["edge_mpc"])
    # The cell isn't stored, but a galaxy-placed sector's is fully
    # determined by its address and edge (as planetgen builds it).
    cell = _sector_cell(row, edge_ly)
    sector = SpaceSector(row["name"], edge_ly=edge_ly, cell=cell)
    sector.epoch_unix = get_orbit_epoch_unix(conn)
    if row["center_x_pc"] is not None:
        sector.place_in_galaxy(tuple(pc_to_ly(row[f"center_{axis}_pc"]) for axis in "xyz"))

    system_rows = conn.execute(
        f"SELECT id, {_SYSTEM_POINT.select('position', 'velocity')} FROM star_systems WHERE sector_id = ? ORDER BY id",
        (sector_id,),
    ).fetchall()

    for r in system_rows:
        star_system = load_star_system(conn, r["id"])
        position = (
            milliparsecs_to_ly(r["position_x_mpc"]),
            milliparsecs_to_ly(r["position_y_mpc"]),
            milliparsecs_to_ly(r["position_z_mpc"]),
        )
        sector.entries.append(SectorSystemEntry(star_system, position, system_config=star_system.system_config,
                                                sector_center_ly=sector.center_galactic_ly,
                                                sector_edge_pc=sector.sector_edge_pc,
                                                epoch_unix=sector.epoch_unix,
                                                velocity_ms=tuple(r[f"velocity_{axis}_kms"] * 1000.0 for axis in "xyz")
                                                if row["center_x_pc"] is not None else None))

    return sector


TIMESTAMPED_TABLES = kind_tables.tables("timestamped", order="timestamped")
"""tuple: The top-level tables carrying v27's `created_at`/`modified_at`
row timestamps -- see `schema.sql`'s "v27" header note. Child rows
(stars, planets, moons, belts, comets, ...) have no timestamps of their
own; a change to one bumps its parent system's `modified_at` instead
(`touch_star_system`)."""


ROGUE_SURFACE_COLUMNS = (
    ("age_gy", "DOUBLE"), ("internal_heat_flux_w_m2", "DOUBLE"), ("effective_temperature_k", "DOUBLE"),
    ("surface_regime", "VARCHAR(24)"), ("surface_temperature_k", "DOUBLE"), ("surface_pressure_pa", "DOUBLE"),
    ("ice_shell_thickness_km", "DOUBLE"), ("ocean_depth_km", "DOUBLE"), ("has_liquid_water", "TINYINT(1)"),
    ("hp_ice_km", "DOUBLE"),
)
"""tuple: `(column, type)` of each v48 `rogue_planets` column, in
`ROGUE_SURFACE_FIELDS` order."""


STATE_SPEED_SQL = "SQRT(velocity_x_kms * velocity_x_kms + velocity_y_kms * velocity_y_kms + velocity_z_kms * velocity_z_kms)"
"""str: SQL for the speed, km/s, of a planet or moon in its stored state (GEN.727: on an ellipse it changes along the
orbit, fastest at periapsis, so the next update is due when it has moved its threshold at the speed it has)."""


def touch_star_system(conn, star_system_id):
    """
    Bumps one `star_systems` row's `modified_at` to now -- how a change to
    one of its child rows (a planet's or moon's rename, say) shows up as a
    change to the system, since child tables carry no timestamps of their
    own (see `schema.sql`'s "v27" header note). A no-op for a
    `star_system_id` that no longer exists.

    Args:
        conn (Connection): Part of the same transaction as the child
            row's own change.
        star_system_id (int): The parent system's `id`.
    """
    conn.execute(
        "UPDATE star_systems SET modified_at = CURRENT_TIMESTAMP(3) WHERE id = ?", (star_system_id,),
    )
    storage_checksum.mark_system(conn, star_system_id)


def touch_sector(conn, sector_id):
    """
    Bumps one `sectors` row's `modified_at` to now -- `touch_star_system`'s
    counterpart one level up, for a change `sectors`' own columns can't
    show (a system in it being deleted, say). A no-op for a `None` or
    no-longer-existing `sector_id`.

    Args:
        conn (Connection): Part of the same transaction as the change.
        sector_id (int or None): The sector's `id`.
    """
    if sector_id is None:
        return
    conn.execute("UPDATE sectors SET modified_at = CURRENT_TIMESTAMP(3), parity_dirty = 1 WHERE id = ?", (sector_id,))


def _schema_version(conn):
    row = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()
    return row["version"]


def schema_status(config=None):
    """
    `(current version, number of migration steps pending)` for a
    database, without migrating it -- `planetgen.cli.migrate --status`, which
    update.sh reads to decide whether to ask about the database at all.
    Raises `SchemaTooNewError` for a database past `SCHEMA_VERSION`.
    """
    # No DDL here (TEST.62): it "changes nothing", and must work for an
    # account that can only read. A database with no version yet is
    # created at the current schema by the migration, so nothing is
    # pending -- unless it already has tables, which is refused (DB.97).
    conn = get_connection(config, ensure_schema=False)
    try:
        version = _stored_version(conn, "schema_migrations")
        _refuse_unversioned(conn, version)
        if version is None:
            version = SCHEMA_VERSION
        _refuse_newer(conn, version, SCHEMA_VERSION)
        _refuse_older(conn, version)
        return version, len(alembic_runner.pending(version))
    finally:
        conn.close()


def migrate_database(config=None, on_step=None, on_progress=None):
    """
    Brings a database's `schema_migrations` bookkeeping up to
    `SCHEMA_VERSION`, applying any migration step in between.

    `get_connection`/`_ensure_schema` always create a brand-new database
    already at `SCHEMA_VERSION` (every `CREATE TABLE IF NOT EXISTS` in
    `schema.sql` reflects the current shape directly), so this function
    only has real work to do against a database created by an older
    version of this project: it runs each Alembic revision
    (`planetgen/db/migrations/`) after the database's version. Databases
    older than the baseline (v130) are refused with `SchemaTooOldError`.
    See `planetgen.cli.migrate` for the CLI wrapper around this.

    Args:
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.
        on_step (callable, optional): Called as `on_step(number, total,
            from_version, to_version)` just before each pending step
            runs (`number` counts from 1 to `total`), so a caller can
            show progress (`planetgen.cli.migrate`'s progress bar).
        on_progress (callable, optional): Called as `on_progress(label,
            done, total)` by a revision that works in batches, for the
            bar of that one step (DB.15).

    Returns:
        int: The database's `schema_migrations` version (always
            `SCHEMA_VERSION` after this call).

    Raises:
        SchemaTooNewError: The database is past `SCHEMA_VERSION` already
            (a newer planetGen migrated it); nothing is changed.
        SchemaTooOldError: The database is older than the baseline.
    """
    conn = get_connection(config, ensure_schema=False)
    try:
        # The whole migration under the schema lock (DB.5), so a first
        # connection elsewhere waits for it instead of stamping the
        # database current halfway through.
        with _schema_lock(conn):
            # Always, not once per process (PERF.12): a migration is when a
            # database's missing tables and its baseline version row appear.
            _apply_schema(conn)
            _schema_ensured.add((config or DEFAULT_MYSQL_CONFIG)._key())
            version = _schema_version(conn)
            _refuse_newer(conn, version, SCHEMA_VERSION)
            _refuse_older(conn, version)

            # Each Alembic revision reached is mirrored into
            # `schema_migrations`.
            def record(reached):
                nonlocal version
                conn.execute("INSERT INTO schema_migrations (version) VALUES (?)", (reached,))
                conn.commit()
                activity_log.event("DB", "migrate", db=(config or DEFAULT_MYSQL_CONFIG).database,
                                  from_version=version, to_version=reached)
                version = reached

            version = alembic_runner.upgrade(config or DEFAULT_MYSQL_CONFIG, version, record, on_step=on_step,
                                             on_progress=on_progress)
        return version
    finally:
        conn.close()


def get_orbit_update_elapsed_years(conn):
    """
    Returns how many years have elapsed since `planetgen.cli.orbits` last
    advanced this database's orbital phases, or `None` if it has never run
    against this database before (nothing to measure elapsed time from
    yet).

    Computed server-side via `TIMESTAMPDIFF` rather than comparing
    `orbit_simulation_state.last_updated_at` against this process' own
    clock (`datetime.now()`) -- correct even when the script runs on a
    different host than the database server, whose clocks aren't
    guaranteed to agree.

    Args:
        conn (Connection): An open, schema-initialized connection.

    Returns:
        float or None.
    """
    row = conn.execute(
        "SELECT TIMESTAMPDIFF(SECOND, last_updated_at, NOW()) AS elapsed_seconds "
        "FROM orbit_simulation_state WHERE id = 1"
    ).fetchone()
    if row is None:
        return None
    return _s_to_yr(row["elapsed_seconds"])


def get_orbit_epoch_unix(conn):
    """
    When `planetgen.cli.orbits` last advanced this database's orbits, in Unix
    seconds: the time every stored position and velocity holds at. `None`
    when it has never run (`systemscene._epoch` reads the same row for the
    scene).

    Args:
        conn (Connection): An open, schema-initialized connection.

    Returns:
        int or None.
    """
    row = conn.execute("SELECT UNIX_TIMESTAMP(last_updated_at) AS unix FROM orbit_simulation_state WHERE id = 1").fetchone()
    return None if row is None or row["unix"] is None else int(row["unix"])


OrbitClock = namedtuple("OrbitClock", ("start_unix", "end_unix"))
OrbitClock.__doc__ = """The span one orbit update moves the database over, in Unix seconds by
the database server's clock: from `start_unix`, when every position without
its own `epoch_unix` holds (the last update), to `end_unix`."""


def orbit_clock(conn, elapsed_years=None):
    """
    The `OrbitClock` for an orbit update: from the last update (or now, when
    none has run) to now, or to `elapsed_years` after the last update when
    given (a simulated step). A server clock that went back since the last
    update gives an empty span at now, so nothing moves backwards.

    Args:
        conn (Connection): An open, schema-initialized connection.
        elapsed_years (float, optional): How far to step instead of to now.

    Returns:
        OrbitClock.

    Raises:
        ValueError: If `elapsed_years` is negative.
    """
    if elapsed_years is not None and elapsed_years < 0:
        raise ValueError(f"elapsed_years must be >= 0, got {elapsed_years}")
    now = float(conn.execute("SELECT UNIX_TIMESTAMP(NOW(6)) AS now").fetchone()["now"])
    last = get_orbit_epoch_unix(conn)
    start = now if last is None else float(last)
    if elapsed_years is not None:
        return OrbitClock(start, start + _yr_to_s(elapsed_years))
    return OrbitClock(min(start, now), now)


def finish_orbit_update(conn, clock):
    """Records `clock.end_unix` as the time every stored position without its
    own epoch holds (`orbit_simulation_state.last_updated_at`): the start of
    the next update. `planetgen.cli.orbits` calls it after the last step.

    The column holds whole seconds. MySQL 8.4 rounds a fraction on insert, which
    could store a second after "now" and make the next run read a clock that
    went backwards (TEST.129), so the fraction is dropped here."""
    whole_seconds = math.floor(clock.end_unix)
    conn.execute(
        "INSERT INTO orbit_simulation_state (id, last_updated_at) VALUES (1, FROM_UNIXTIME(?)) "
        "ON DUPLICATE KEY UPDATE last_updated_at = FROM_UNIXTIME(?)",
        (whole_seconds, whole_seconds),
    )
    conn.commit()


def _budget_tables():
    """The tables an orbital-update run moves, as `{key: (table, due column, where)}` (PERF.103): what the run counts
    to plan its ceiling and what each step asks its share of. `where` selects the rows the step moves."""
    tables = {
        "planets": ("planets", "next_update_due", "period_years > 0"),
        "moons": ("moons", "next_update_due", "period_years > 0"),
        "binary_mutual_orbits": ("star_systems", "binary_next_update_due",
                                 "is_binary = 1 AND binary_mutual_orbital_period_years > 0"),
        "comets": ("comets", "next_update_due", "1 = 1"),
        "facilities": ("facilities", "next_update_due", "placement IN ('orbital', 'asteroid') AND orbit_period_years > 0"),
        "star_systems": ("star_systems", "next_update_due", "1 = 1"),
    }
    for table in PLACED_PHENOMENON_TABLES:
        if table != "quasars":
            anchored = "star_id IS NULL" if table in ("black_holes", "neutron_stars") else "1 = 1"
            tables[table] = (table, "next_update_due", anchored)
    tables["space_facilities"] = ("facilities", "next_update_due", "host_type = 'space' AND center_x_pc IS NOT NULL")
    for table in SCATTER_MOTION_TABLES:
        tables[table] = (table, "next_update_due", _SCATTER_ROWS[table])
    return tables


class RowBudget:
    """
    The most rows one orbital-update run moves (PERF.103).

    `plan` counts what is due in every table the run moves. When that is no more than the ceiling every table is
    unlimited and the run is as it always was. When it is more, each table gets a limit: a floor so a small table is
    never starved (a quarter of the ceiling shared out equally), and the rest of the ceiling in proportion to
    what is left due in each table. A step moves its table's oldest-due rows up to its limit (`limit` for the
    set-based steps, `take` for the Python loops); the rest stay due, are counted by `defer`, and are first in line
    the next night. The limits add up to at most the ceiling, so a run never moves more than it.

    Args:
        ceiling (int or None): The ceiling; None or 0 for no ceiling.
    """

    def __init__(self, ceiling=None):
        self.ceiling = ceiling or None
        self.candidates = {}
        self.limits = {}
        self.used = {}
        self.deferred = {}

    def plan(self, conn, clock):
        """Counts the rows each table could move, sets the limits, and returns the candidate count of every table."""
        tables = _budget_tables()
        self.candidates = {}
        for key, (table, due, where) in tables.items():
            self.candidates[key] = conn.execute(
                f"SELECT COUNT(*) AS n FROM {table} WHERE {where} AND ({due} IS NULL OR {due} <= ?)",
                (clock.end_unix,)).fetchone()["n"]
        self.allocate()
        return self.candidates

    def allocate(self):
        """Sets the limits from `candidates`: none when everything fits under the ceiling."""
        self.limits = {}
        if self.ceiling is None or sum(self.candidates.values()) <= self.ceiling:
            return
        floor = self.ceiling // (4 * max(len(self.candidates), 1))
        self.limits = {key: min(count, floor) for key, count in self.candidates.items()}
        spare = self.ceiling - sum(self.limits.values())
        left = {key: count - self.limits[key] for key, count in self.candidates.items()}
        left_total = sum(left.values())
        if left_total:
            shares = {key: spare * count // left_total for key, count in left.items()}
            for key in sorted(left, key=lambda k: (spare * left[k]) % left_total, reverse=True):
                if sum(shares.values()) >= spare:
                    break
                if shares[key] < left[key]:
                    shares[key] += 1
            for key, share in shares.items():
                self.limits[key] += min(share, left[key])

    def limit(self, key):
        """How many rows the set-based step for `key` may move, or None for all that are due."""
        return self.limits.get(key)

    def take(self, key):
        """For a Python loop: whether it may move one more row of `key`. A refusal is counted as deferred."""
        limit = self.limits.get(key)
        if limit is not None and self.used.get(key, 0) >= limit:
            self.defer(key)
            return False
        self.used[key] = self.used.get(key, 0) + 1
        return True

    def defer(self, key, count=1):
        """Counts `count` rows of `key` that were due and left for the next run."""
        self.deferred[key] = self.deferred.get(key, 0) + count

    def defer_remaining(self, conn, clock, key):
        """After a set-based step that ran with a limit: counts what is still due in `key`'s table."""
        table, due, where = _budget_tables()[key]
        left = conn.execute(f"SELECT COUNT(*) AS n FROM {table} WHERE {where} AND {due} <= ?",
                            (clock.end_unix,)).fetchone()["n"]
        if left:
            self.defer(key, left)

    @property
    def backlog(self):
        """How many due rows the run left for the next one."""
        return sum(self.deferred.values())

    def warning(self):
        """The warning text when rows were left (None when none were): the backlog and where it sits."""
        if not self.backlog:
            return None
        where = ", ".join(f"{key} {count}" for key, count in sorted(self.deferred.items(), key=lambda kv: -kv[1]) if count)
        return (f"{self.backlog} due row(s) are over the ceiling of {self.ceiling} rows a run and stay due for the next "
                f"run, oldest first ({where}). Raise the ceiling (--row-ceiling) if the backlog keeps growing.")


def _ordered_limit(budget, key):
    """The SQL tail of a set-based step that has a limit (`ORDER BY` the due time, `LIMIT`) and that limit, or empty."""
    limit = budget.limit(key) if budget is not None else None
    return limit, ("" if limit is None else f" ORDER BY {_budget_tables()[key][1]}, id LIMIT {int(limit)}")


def _update_interval_sql(speed_kms, scale):
    """SQL for `position.update_interval_s` of a speed in km/s (`speed_kms` is
    an SQL expression): seconds until the object has moved its scale's
    threshold, capped."""
    threshold_km = THRESHOLDS_M[scale] / 1000.0
    cap = MAX_UPDATE_INTERVAL_S
    return f"(CASE WHEN {speed_kms} > 0 THEN LEAST({cap!r}, {threshold_km!r} / ({speed_kms})) ELSE {cap!r} END)"


_FACILITY_SPEED_SQL = (f"(orbit_distance_km * 2 * PI() / (orbit_period_years"
                       f" * {_SECONDS_PER_YEAR!r}))")
"""str: SQL for an orbiting facility's speed, km/s, on its circular orbit."""

_FACILITY_INTERVAL_SQL = (
    f"(CASE WHEN host_type IN ('planet', 'moon') THEN {_update_interval_sql(_FACILITY_SPEED_SQL, 'planetary')}"
    f" ELSE {_update_interval_sql(_FACILITY_SPEED_SQL, 'system')} END)")
"""str: A facility orbiting a planet or moon is on the planetary scale, one
round a star or in a belt on the system scale."""


def _schedule_updates(conn, clock, table, interval_sql, where, prefix="", keep_modified=False):
    """Works out `next_update_due` for `table`'s rows that don't have one yet
    (new or edited since the last update): their epoch (the last update's,
    `clock.start_unix`, when they have none) plus `interval_sql`. Rows
    matching `where` only."""
    keep = ", modified_at = modified_at" if keep_modified else ""
    conn.execute(
        f"UPDATE {table} SET {prefix}epoch_unix = COALESCE({prefix}epoch_unix, ?),"
        f" {prefix}next_update_due = COALESCE({prefix}epoch_unix, ?) + {interval_sql}{keep}"
        f" WHERE {prefix}next_update_due IS NULL AND {where}",
        (clock.start_unix, clock.start_unix),
    )


def _elapsed_years_sql(prefix=""):
    """SQL for the years from a row's epoch to the update's end (the `?`)."""
    return f"((? - {prefix}epoch_unix) / {_SECONDS_PER_YEAR!r})"


def _phase_turn_sql(period_column, prefix=""):
    """SQL for the degrees an orbit of `period_column` (years) turns from a
    row's epoch to the update's end (the `?`), whole orbits dropped first
    (GEN.108): a step many orbits long wraps exactly, where elapsed / period
    would leave few digits for the fraction."""
    return f"(MOD({_elapsed_years_sql(prefix)}, {period_column}) / {period_column} * 360)"


def _due_interval_s(row, clock, speed_kms, scale, prefix=""):
    """For the rows the Python steps move: `(epoch, interval, due)` -- the
    row's epoch (the last update's when it has none), the seconds it takes
    to move its threshold at `speed_kms`, and when it is due (worked out
    from those when the row has no `next_update_due` yet)."""
    epoch = row[f"{prefix}epoch_unix"]
    epoch = clock.start_unix if epoch is None else epoch
    interval = update_interval_s((speed_kms or 0.0) * 1000.0, scale)
    due = row[f"{prefix}next_update_due"]
    return epoch, interval, (epoch + interval if due is None else due)


ORBITAL_PHASE_UPDATES = 6
"""int: How many table updates `advance_orbital_phases` runs (the `total` it
reports progress against)."""


def _report(on_progress, label, done, total):
    """Calls `on_progress(label, done, total)` when the caller gave one --
    the hook a command-line run draws its progress bar from."""
    if on_progress is not None:
        on_progress(label, done, total)


STATE_BATCH_ROWS = 20_000
"""int: How many due planets or moons `_advance_body_states` reads, steps and writes at a time."""

STATE_WRITE_ROWS = 2_000
"""int: How many stepped rows go into one `UPDATE` (a `VALUES ROW` table joined to the table)."""

_STATE_COLUMNS = ("position_x_km", "position_y_km", "position_z_km", "velocity_x_kms", "velocity_y_kms", "velocity_z_kms")


def _advance_body_states(conn, clock, table, scale, budget=None):
    """
    Moves every due planet or moon of `table` (`planets` or `moons`) from its own epoch to `clock.end_unix` along its
    orbit (GEN.289): its stored position and velocity are stepped by `kepler.universal_step`
    (`physics/state_kernel.py`, compiled, the same arithmetic) on the pull `mu = (2 pi / period)^2 distance^3` of what it
    goes round (`physics.body_state.mu_km3_s2`), and the epoch becomes the update's end and the next due time that plus
    the interval its new speed gives (`position.update_interval_s`). `batch` rows are read, stepped and written at a
    time, oldest due first, up to the run's limit for the table when `budget` has one.

    A row the step cannot move (no convergence, a body on its primary) is logged and keeps its state and epoch; its
    next due time still moves on, so it is tried again in a later run and a bad row never holds the others up.

    Returns:
        int: The rows moved.
    """
    from planetgen.physics import state_kernel  # compiled: kept off the web app's import path (OPS.42)
    import numpy as np

    limit = budget.limit(table) if budget is not None else None
    threshold_km = THRESHOLDS_M[scale] / 1000.0
    seconds_per_year = _yr_to_s(1.0)
    columns = ", ".join(_STATE_COLUMNS)
    moved = seen = 0
    while limit is None or seen < limit:
        take = STATE_BATCH_ROWS if limit is None else min(STATE_BATCH_ROWS, limit - seen)
        rows = conn.execute(
            f"SELECT id, {columns}, distance_km, period_years, epoch_unix FROM {table}"
            f" WHERE period_years > 0 AND next_update_due <= ? ORDER BY next_update_due, id LIMIT {int(take)}",
            (clock.end_unix,)).fetchall()
        if not rows:
            break
        seen += len(rows)
        ids = [row["id"] for row in rows]
        states = np.array([[row[name] for name in _STATE_COLUMNS] for row in rows], dtype=np.float64)
        distance = np.array([row["distance_km"] for row in rows], dtype=np.float64)
        period_s = np.array([row["period_years"] for row in rows], dtype=np.float64) * seconds_per_year
        epoch = np.array([row["epoch_unix"] for row in rows], dtype=np.float64)
        mus = (2.0 * math.pi / period_s) ** 2 * distance ** 3
        stepped, status = state_kernel.step_many(states, mus, clock.end_unix - epoch)
        failed = status != state_kernel.OK
        for index in np.flatnonzero(failed):
            log.error(f"{table} {ids[index]}: the orbital step did not move it (status {int(status[index])}); "
                      f"its state and epoch are kept")
        new_epoch = np.where(failed, epoch, clock.end_unix)
        speed = np.sqrt(stepped[:, 3] ** 2 + stepped[:, 4] ** 2 + stepped[:, 5] ** 2)
        with np.errstate(divide="ignore"):
            interval = np.where(speed > 0.0, np.minimum(threshold_km / np.where(speed > 0.0, speed, 1.0),
                                                        MAX_UPDATE_INTERVAL_S), MAX_UPDATE_INTERVAL_S)
        due = clock.end_unix + interval
        for start in range(0, len(ids), STATE_WRITE_ROWS):
            stop = start + STATE_WRITE_ROWS
            values = ",".join(
                "ROW(" + ",".join(map(_scatter_value, (ids[i], *stepped[i].tolist(), float(new_epoch[i]), float(due[i]))))
                + ")" for i in range(start, min(stop, len(ids))))
            conn.execute_literal(
                f"UPDATE {table} t JOIN (VALUES {values}) AS v(id, px, py, pz, vx, vy, vz, epoch, due) ON t.id = v.id"
                " SET t.position_x_km = v.px, t.position_y_km = v.py, t.position_z_km = v.pz,"
                " t.velocity_x_kms = v.vx, t.velocity_y_kms = v.vy, t.velocity_z_kms = v.vz,"
                " t.epoch_unix = v.epoch, t.next_update_due = v.due")
        moved += int(len(ids) - failed.sum())
    return moved


def advance_orbital_phases(conn, clock, on_progress=None, budget=None):
    """
    Advances the orbits inside star systems that are due (GEN.106): every
    planet's and moon's state vector (GEN.289, `_advance_body_states`: its
    stored position and velocity stepped on its orbit by
    `kepler.universal_step`), and a binary pair's
    `binary_mutual_orbital_phase_deg`, each by the time from its own epoch
    to `clock.end_unix`, given its own stored period -- the bodies in
    batches through a compiled step, the pairs by one set-based `UPDATE`
    (see `planetgen.cli.orbits`). Galactic orbits are
    `advance_galactic_positions`'.

    Due: a row's `next_update_due` (indexed) is at or before the update's
    end, so it has moved at least its threshold since its epoch (0.01 AU
    for a planet or a pair's mutual orbit, 100,000 km for a moon; see
    `physics.position.THRESHOLDS_M`). Rows without one (new or edited since
    the last update) get one first, from their epoch and orbital speed.
    A row that isn't due is neither moved nor counted. A moved row's epoch
    becomes the update's end and its next due time that plus its interval.

    A planet or moon has no phase: the step reads its stored position and
    velocity, so an orbit another body's pull has changed is carried on as
    it now is. A quiet circular orbit comes out where the circular phase
    put it, to rounding (`tests/test_state_kernel.py`).

    `binary_mutual_position_x/y/z_km` (the secondary's position relative to
    the primary) and each star's offset from the pair's barycenter
    (`binary_primary/secondary_position_*_km`, from the constant
    `binary_secondary_mass_fraction`) move with the mutual phase the same
    way, for both binary configurations.

    v20 "reflex offset"/"wobble" values (`stars.reflex_offset_*_km` from
    each star's planets, `planets.reflex_offset_*_km` from each planet's
    moons, `star_systems.binary_planetary_wobble_*_km` from a close pair's
    circumbinary planets) are derived from the children's positions, so
    they're recomputed on every call with a correlated subquery
    (`SELECT SUM(...) FROM <children> WHERE <parent link>`) -- see
    `orbits.calculate_reflex_offset`.

    v27: `UPDATE`s on `star_systems` set `modified_at = modified_at`, which
    stops MySQL's `ON UPDATE CURRENT_TIMESTAMP` from firing -- an orbit
    tick is the simulation clock moving, not the row being edited.

    Args:
        conn (Connection): An open, schema-initialized, read-write
                           connection.
        clock (OrbitClock): The span to move over (`orbit_clock`).
        on_progress (callable, optional): `on_progress(label, done, total)`
                               after each of the `total` table updates
                               (see `_report`); `planetgen.cli.orbits`
                               draws its progress bar from it.
        budget (RowBudget, optional): The run's row ceiling (PERF.103); a
                               table over its limit moves its oldest-due
                               rows and leaves the rest due.

    Returns:
        dict: `{name: rows_updated}` -- `"planets"`, `"moons"`,
            `"star_reflex_offsets"`, `"planet_reflex_offsets"`,
            `"binary_mutual_orbits"`, `"binary_planetary_wobbles"`.
    """
    counts = {}
    for table, scale in (("planets", "system"), ("moons", "planetary")):
        interval = _update_interval_sql(STATE_SPEED_SQL, scale)
        _schedule_updates(conn, clock, table, interval, "period_years > 0")
        counts[table] = _advance_body_states(conn, clock, table, scale, budget)
        if budget is not None and budget.limit(table) is not None:
            budget.defer_remaining(conn, clock, table)
        _report(on_progress, table, len(counts), ORBITAL_PHASE_UPDATES)

    # v20: each star's own reflex-offset "wobble" from the planets it
    # hosts (planets.star_id) -- a correlated subquery summing every
    # hosted planet's individual pairwise pull, the same
    # orbits.calculate_reflex_offset formula generation time uses (see its
    # docstring). Recomputed on every run, using each planet's
    # position_x/y/z_km as it now stands -- a cheap derived value, not an
    # independently advancing phase, so it has no due time of its own. NULL/0 rows (no planets) are simply left alone by
    # the WHERE EXISTS guard, matching those rows' already-NULL default.
    cur = conn.execute(
        """
        UPDATE stars s
        SET reflex_offset_x_km = -(
                SELECT COALESCE(SUM((p.mass_kg / (s.mass_kg + p.mass_kg)) * p.position_x_km), 0)
                FROM planets p WHERE p.star_id = s.id
            ),
            reflex_offset_y_km = -(
                SELECT COALESCE(SUM((p.mass_kg / (s.mass_kg + p.mass_kg)) * p.position_y_km), 0)
                FROM planets p WHERE p.star_id = s.id
            ),
            reflex_offset_z_km = -(
                SELECT COALESCE(SUM((p.mass_kg / (s.mass_kg + p.mass_kg)) * p.position_z_km), 0)
                FROM planets p WHERE p.star_id = s.id
            )
        WHERE EXISTS (SELECT 1 FROM planets p WHERE p.star_id = s.id)
        """
    )
    counts["star_reflex_offsets"] = cur.rowcount
    _report(on_progress, "star_reflex_offsets", len(counts), ORBITAL_PHASE_UPDATES)

    # v20: each planet's own reflex-offset "wobble" from the moons it
    # hosts -- identical shape/reasoning to the stars UPDATE above, one
    # level down (moons.planet_id).
    cur = conn.execute(
        """
        UPDATE planets pl
        SET reflex_offset_x_km = -(
                SELECT COALESCE(SUM((m.mass_kg / (pl.mass_kg + m.mass_kg)) * m.position_x_km), 0)
                FROM moons m WHERE m.planet_id = pl.id
            ),
            reflex_offset_y_km = -(
                SELECT COALESCE(SUM((m.mass_kg / (pl.mass_kg + m.mass_kg)) * m.position_y_km), 0)
                FROM moons m WHERE m.planet_id = pl.id
            ),
            reflex_offset_z_km = -(
                SELECT COALESCE(SUM((m.mass_kg / (pl.mass_kg + m.mass_kg)) * m.position_z_km), 0)
                FROM moons m WHERE m.planet_id = pl.id
            )
        WHERE EXISTS (SELECT 1 FROM moons m WHERE m.planet_id = pl.id)
        """
    )
    counts["planet_reflex_offsets"] = cur.rowcount
    _report(on_progress, "planet_reflex_offsets", len(counts), ORBITAL_PHASE_UPDATES)

    # Mutual orbit: shared by both binary configurations, with its own
    # clock (binary_epoch_unix/binary_next_update_due). v20: also
    # recomputes each star's own offset from the pair's barycenter
    # (binary_primary/secondary_position_*_km) from the freshly-advanced
    # binary_mutual_position_*_km above and the constant
    # binary_secondary_mass_fraction, using the same left-to-right
    # single-table SET evaluation trick binary_mutual_position_*_km's own
    # computation already relies on (each position expression here reads
    # binary_mutual_position_*_km's *new* value, assigned earlier in this
    # same SET list).
    binary_interval = _update_interval_sql(
        f"(binary_separation_km * 2 * PI() / (binary_mutual_orbital_period_years"
        f" * {_SECONDS_PER_YEAR!r}))", "system")
    _schedule_updates(conn, clock, "star_systems", binary_interval,
                      "is_binary = 1 AND binary_mutual_orbital_period_years > 0", prefix="binary_", keep_modified=True)
    limited, tail = _ordered_limit(budget, "binary_mutual_orbits")
    cur = conn.execute(
        f"""
        UPDATE star_systems
        SET binary_mutual_orbital_phase_deg = MOD(binary_mutual_orbital_phase_deg
                + {_phase_turn_sql("binary_mutual_orbital_period_years", "binary_")}, 360),
            binary_mutual_position_x_km = binary_separation_km * (
                COS(RADIANS(binary_mutual_orbital_ascending_node_deg)) * COS(RADIANS(binary_mutual_orbital_phase_deg))
                - SIN(RADIANS(binary_mutual_orbital_ascending_node_deg)) * SIN(RADIANS(binary_mutual_orbital_phase_deg))
                  * COS(RADIANS(binary_mutual_orbital_inclination_deg))
            ),
            binary_mutual_position_y_km = binary_separation_km * (
                SIN(RADIANS(binary_mutual_orbital_ascending_node_deg)) * COS(RADIANS(binary_mutual_orbital_phase_deg))
                + COS(RADIANS(binary_mutual_orbital_ascending_node_deg)) * SIN(RADIANS(binary_mutual_orbital_phase_deg))
                  * COS(RADIANS(binary_mutual_orbital_inclination_deg))
            ),
            binary_mutual_position_z_km = binary_separation_km
                * SIN(RADIANS(binary_mutual_orbital_phase_deg)) * SIN(RADIANS(binary_mutual_orbital_inclination_deg)),
            binary_primary_position_x_km = -binary_secondary_mass_fraction * binary_mutual_position_x_km,
            binary_primary_position_y_km = -binary_secondary_mass_fraction * binary_mutual_position_y_km,
            binary_primary_position_z_km = -binary_secondary_mass_fraction * binary_mutual_position_z_km,
            binary_secondary_position_x_km = (1 - binary_secondary_mass_fraction) * binary_mutual_position_x_km,
            binary_secondary_position_y_km = (1 - binary_secondary_mass_fraction) * binary_mutual_position_y_km,
            binary_secondary_position_z_km = (1 - binary_secondary_mass_fraction) * binary_mutual_position_z_km,
            binary_epoch_unix = ?,
            binary_next_update_due = ? + {binary_interval},
            modified_at = modified_at
        WHERE is_binary = 1
          AND binary_mutual_orbital_period_years > 0
          AND binary_next_update_due <= ?{tail}
        """,
        (clock.end_unix, clock.end_unix, clock.end_unix, clock.end_unix),
    )
    if limited is not None:
        budget.defer_remaining(conn, clock, "binary_mutual_orbits")
    counts["binary_mutual_orbits"] = cur.rowcount
    _report(on_progress, "binary_mutual_orbits", len(counts), ORBITAL_PHASE_UPDATES)

    # v20: circumbinary (P-type) planets' combined pull on the whole pair
    # -- same correlated-subquery shape as the stars/planets reflex-offset
    # UPDATEs above, grouped by star_system_id instead (circumbinary
    # planets have star_id IS NULL, so there's no stars row to attach this
    # to -- see schema.sql's "v20" header note on why it's modeled as one
    # shared wobble rather than split between primary/secondary). Also
    # recomputed on every run, no due time of its own.
    cur = conn.execute(
        """
        UPDATE star_systems ss
        SET binary_planetary_wobble_x_km = -(
                SELECT COALESCE(SUM((p.mass_kg / (ss.binary_effective_mass_kg + p.mass_kg)) * p.position_x_km), 0)
                FROM planets p WHERE p.star_system_id = ss.id AND p.star_id IS NULL
            ),
            binary_planetary_wobble_y_km = -(
                SELECT COALESCE(SUM((p.mass_kg / (ss.binary_effective_mass_kg + p.mass_kg)) * p.position_y_km), 0)
                FROM planets p WHERE p.star_system_id = ss.id AND p.star_id IS NULL
            ),
            binary_planetary_wobble_z_km = -(
                SELECT COALESCE(SUM((p.mass_kg / (ss.binary_effective_mass_kg + p.mass_kg)) * p.position_z_km), 0)
                FROM planets p WHERE p.star_system_id = ss.id AND p.star_id IS NULL
            ),
            modified_at = modified_at
        WHERE ss.binary_configuration = 'close'
        """
    )
    counts["binary_planetary_wobbles"] = cur.rowcount
    _report(on_progress, "binary_planetary_wobbles", len(counts), ORBITAL_PHASE_UPDATES)

    return counts


def advance_comet_orbits(conn, clock, budget=None):
    """
    Advances every due comet's own orbital anomaly (`mean_anomaly_deg` for
    an elliptical comet, `parabolic_mean_anomaly` for a parabolic one) from
    its epoch to `clock.end_unix`, and recomputes `distance_km`/
    `position_x/y/z_km`/`orbital_speed_kms` and the velocity from the new
    anomaly via `kepler.comet_orbital_state` -- the Kepler/Barker-equation
    analog of `advance_orbital_phases`'s planet/moon handling, called
    separately by `planetgen.cli.orbits` alongside it.

    A comet's position is NOT a linear SQL expression of time: turning an
    advanced anomaly into a distance/position requires solving Kepler's
    equation (Newton-Raphson, elliptical) or Barker's equation (a
    real-cube-root closed form, parabolic). So this fetches the due rows
    and does that in Python (a small table -- see
    `tuning.SYSTEM_COMET_COUNT_RANGE`), then writes them back in one
    batched `executemany`.

    Due (GEN.106): like a planet's, a comet's `next_update_due` is when it
    will have moved 0.01 AU at its current speed; a row without one gets
    one first. A comet that isn't due is neither moved nor counted. An
    elliptical comet's `mean_anomaly_deg` advances
    `MOD(current + (elapsed_years / period_years) * 360, 360)` like
    `orbital_phase_deg`; a parabolic comet's `parabolic_mean_anomaly`
    advances LINEARLY (via `kepler.parabolic_mean_anomaly`) and does NOT
    wrap (a parabolic pass is a one-shot event, not periodic).

    Args:
        conn (Connection): An open, schema-initialized, read-write
                           connection.
        clock (OrbitClock): The span to move over (`orbit_clock`).
        budget (RowBudget, optional): The run's row ceiling (PERF.103); the
            comets over their limit stay due.

    Returns:
        int: The number of `comets` rows moved.
    """
    rows = conn.execute(
        "SELECT id, orbit_type, perihelion_distance_km, eccentricity, inclination_deg, "
        "arg_periapsis_deg, ascending_node_deg, orbital_period_years, mean_anomaly_deg, "
        "parabolic_mean_anomaly, primary_mass_solar, orbital_speed_kms, epoch_unix, next_update_due "
        "FROM comets WHERE next_update_due IS NULL OR next_update_due <= ? ORDER BY next_update_due, id",
        (clock.end_unix,),
    ).fetchall()

    elapsed_by_id, scheduled = {}, []
    for row in rows:
        epoch, _interval, due = _due_interval_s(row, clock, row["orbital_speed_kms"], "system")
        if due > clock.end_unix or (budget is not None and not budget.take("comets")):
            scheduled.append((epoch, due, row["id"]))
        else:
            elapsed_by_id[row["id"]] = _s_to_yr(clock.end_unix - epoch)
    if scheduled:
        conn.executemany("UPDATE comets SET epoch_unix = ?, next_update_due = ? WHERE id = ?", scheduled)
    rows = [row for row in rows if row["id"] in elapsed_by_id]

    update_params = []
    to_kms = _au_yr_to_kms(1.0)
    # Kepler's equation for every elliptical comet that moves, in one vectorised solve.
    moving = [row for row in rows if row["orbit_type"] == "elliptical"]
    new_mean_anomalies_deg = {
        row["id"]: (row["mean_anomaly_deg"]
                    + math.fmod(elapsed_by_id[row["id"]], row["orbital_period_years"]) / row["orbital_period_years"] * 360)
        % 360
        for row in moving}
    eccentric_anomalies = dict(zip(
        (row["id"] for row in moving),
        kepler.solve_eccentric_anomalies(
            [math.radians(new_mean_anomalies_deg[row["id"]]) for row in moving],
            [row["eccentricity"] for row in moving])))
    for row in rows:
        perihelion_distance_au = _km_to_au(row["perihelion_distance_km"])

        eccentric_anomaly_rad = None
        if row["orbit_type"] == "elliptical":
            new_mean_anomaly_deg = new_mean_anomalies_deg[row["id"]]
            eccentric_anomaly_rad = float(eccentric_anomalies[row["id"]])
            mean_anomaly_rad = math.radians(new_mean_anomaly_deg)
            new_parabolic_mean_anomaly = None
            parabolic_mean_anomaly_value = None
        else:
            new_mean_anomaly_deg = None
            mean_anomaly_rad = None
            new_parabolic_mean_anomaly = row["parabolic_mean_anomaly"] + kepler.parabolic_mean_anomaly(
                elapsed_by_id[row["id"]], perihelion_distance_au, row["primary_mass_solar"]
            )
            parabolic_mean_anomaly_value = new_parabolic_mean_anomaly

        state = kepler.comet_orbital_state(
            row["orbit_type"], perihelion_distance_au, row["eccentricity"],
            row["inclination_deg"], row["arg_periapsis_deg"], row["ascending_node_deg"],
            row["primary_mass_solar"],
            mean_anomaly_rad=mean_anomaly_rad,
            parabolic_mean_anomaly_value=parabolic_mean_anomaly_value,
            orbital_period_years=row["orbital_period_years"],
            eccentric_anomaly_rad=eccentric_anomaly_rad,
        )

        update_params.append((
            new_mean_anomaly_deg, new_parabolic_mean_anomaly,
            _au_to_km(state["distance_au"]),
            _au_to_km(state["position_x_au"]),
            _au_to_km(state["position_y_au"]),
            _au_to_km(state["position_z_au"]),
            state["orbital_speed_kms"],
            state["velocity_x_au_per_year"] * to_kms,
            state["velocity_y_au_per_year"] * to_kms,
            state["velocity_z_au_per_year"] * to_kms,
            clock.end_unix,
            clock.end_unix + update_interval_s(state["orbital_speed_kms"] * 1000.0, "system"),
            row["id"],
        ))

    # One batched round trip for every row that needs writing, not one
    # `execute` per comet -- see this function's own docstring.
    if update_params:
        conn.executemany(
            """
            UPDATE comets
            SET mean_anomaly_deg = ?, parabolic_mean_anomaly = ?,
                distance_km = ?, position_x_km = ?, position_y_km = ?, position_z_km = ?,
                orbital_speed_kms = ?, velocity_x_kms = ?, velocity_y_kms = ?, velocity_z_kms = ?,
                epoch_unix = ?, next_update_due = ?
            WHERE id = ?
            """,
            update_params,
        )

    conn.commit()
    return len(update_params)
