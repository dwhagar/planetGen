# planetgen/db/store.py

"""
Database Persistence (private)
===============================

Writes already-generated `StarSystem`/`SpaceSector` objects into the
MySQL database described by `planetgen/db/schema.sql` and
`docs/database-schema.md`. Leading underscore -- this module is an internal
implementation detail of the persistence boundary, not part of the
package's public generation API (`StarSystem`, `SpaceSector`, `Planet`,
etc. remain the public surface).

This module owns every unit conversion at the point of writing a value
into the database -- generation/physics code elsewhere in the package
keeps its own native units (km/AU/ly) throughout and never imports this
module. See `docs/database-schema.md` for the two-tier distance-unit convention
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
connection pool (`DBUtils.PooledDB`) rather than opening a fresh TCP
connection per call, per TODO.md's "add real connection pooling" note.
"""

import contextlib
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
from dbutils.pooled_db import PooledDB

from planetgen.admin import activity_log
from planetgen.names import object_id as objectId
from planetgen.galaxy import seed as galaxySeed, version_key as versionKey
from planetgen.physics import constants as physical_constants, kepler
from planetgen import tuning
from planetgen.util import log
from planetgen.util.appconfig import load_config
from planetgen.generation.belt import AsteroidBelt
from planetgen.generation.phenomena.asteroid_field import AsteroidField, asteroid_field_designation
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.population import facilities as facility_rules
from planetgen.generation.comet import Comet, comet_designation, rename_comet_designation
from planetgen.generation.config import SystemConfig
from planetgen.generation.binary import BinaryStarProxy
from planetgen.galaxy.density import GalaxyShape
from planetgen.galaxy.drill import DrillBlock
from planetgen.galaxy.geometry import (
    SectorCell, galaxy_to_local_pc, local_to_galaxy_pc, provisional_sector_designation, sector_address_at,
)
from planetgen.names.bodies import rename_prefix, wide_pair_first_word
from planetgen.names.wordlists import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from planetgen.names.uniqueness import (
    MAX_SYSTEM_NAME_WORDS, fits_word_limit, has_diminutive, resolve_diminutive, resolve_greek_roman_collision,
    word_limit_for,
)
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.planet import Planet
from planetgen.physics.rogue_surface import ROGUE_SURFACE_FIELDS, rogue_surface_conditions
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet, default_rogue_planet_class, interstellar_comet_designation
from planetgen.galaxy.sector import SectorSystemEntry, SpaceSector, classify_octant, distance_between
from planetgen.generation.star import Star
from planetgen.generation.phenomena.quasar import Quasar
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.generation.system import StarSystem
from stellarObjects.utils import (
    calculate_galactic_orbit, generate_phoneme_salad_name, generate_sector_name, ly_to_milliparsecs, ly_to_pc,
    milliparsecs_to_ly, mpc_to_pc, pc_to_ly,
)
from planetgen.generation.wide_binary import WideBinaryPair

SCHEMA_VERSION = 53
"""int: Matches `star_systems.schema_version` and the highest row in the
`schema_migrations` table (see `planetgen/db/schema.sql`'s header
comment). Also the target version `migrate_database` brings a database's
`schema_migrations` bookkeeping up to."""



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

CONTROL_SCHEMA_VERSION = 7
"""int: Version counter for `control_schema.sql`, independent of
`SCHEMA_VERSION` above -- see that file's header comment for why the
control plane (admin identities/sessions/API keys/audit log) is a
separate schema with its own versioning. v2 added `login_throttle`
(SEC.1, SEC.21), v5 the work queue's `work_jobs`/`work_tasks`/
`work_lease` (PERF.8), v6 `generation_stats`/`generation_size` (PERF.3,
PERF.10), v7 the job tree's columns on `work_jobs` and the queue pause
on `work_lease` (ADM.12, ADM.10). New tables need nothing more than
`CREATE TABLE IF NOT EXISTS`; new columns on an existing table are
added by `_add_control_columns`."""

_CONTROL_COLUMNS = {
    "work_jobs": [
        ("parent_id", "VARCHAR(32) NULL"),
        ("root_id", "VARCHAR(32) NULL"),
        ("kind", "VARCHAR(32) NOT NULL DEFAULT 'queue'"),
        ("seconds", "DOUBLE NULL"),
        ("tasks_total", "INT UNSIGNED NULL"),
        ("web_job_id", "VARCHAR(32) NULL"),
        ("database_name", "VARCHAR(64) NULL"),
        ("argv", "TEXT NULL"),
        ("control", "VARCHAR(16) NULL"),
    ],
    "work_lease": [
        ("paused", "TINYINT(1) NOT NULL DEFAULT 0"),
        ("paused_by", "VARCHAR(64) NULL"),
        ("paused_at", "DATETIME(6) NULL"),
    ],
}
"""dict: Columns added to existing control tables after they were first
created (v7), as `table -> [(column, definition)]`, in order. Must match
`control_schema.sql`'s `CREATE TABLE` statements."""

CONTROL_SCHEMA_PATH = os.path.join(_PACKAGE_DIR, "control_schema.sql")
"""str: Path to the DDL file applied by `_ensure_control_schema`."""

CONTROL_DB_ENV_VAR = "PLANETGEN_CONTROL_DATABASE"
"""str: Env var naming the one MySQL schema the control plane lives in
(admin identities are global to a deployment, not per-galaxy -- see
`control_schema.sql`'s header comment). Falls back to `config.json`'s
`control_database` (see `planetgen.util.appconfig`), then
`DEFAULT_CONTROL_DATABASE`, when unset."""

DEFAULT_CONTROL_DATABASE = "planetgen_control"
"""str: Default control-schema name when neither `CONTROL_DB_ENV_VAR` nor
`config.json`'s `control_database` is set."""


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
    mirrors every other entry point in this project (`sectorGen.py`,
    `systemGen.py`, `planetgen.db.query`, `planetgen/api/config.py`) reading its own
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
    `mysql` section (see `planetgen.util.appconfig`), then the
    hardcoded default below.
    """

    def __init__(self, host=None, port=None, user=None, password=None, database=None):
        defaults = load_config()["mysql"]
        self.host = host if host is not None else os.environ.get("PLANETGEN_MYSQL_HOST", defaults["host"])
        self.port = int(port if port is not None else os.environ.get("PLANETGEN_MYSQL_PORT", defaults["port"]))
        self.user = user if user is not None else os.environ.get("PLANETGEN_MYSQL_USER", defaults["user"])
        self.password = password if password is not None else os.environ.get("PLANETGEN_MYSQL_PASSWORD", defaults["password"])
        self.database = database if database is not None else os.environ.get("PLANETGEN_MYSQL_DATABASE", defaults["database"])

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


def add_mysql_connection_args(parser):
    """
    Adds the `--mysql-host`/`--mysql-port`/`--mysql-user`/
    `--mysql-password`/`--mysql-database` optional overrides to `parser`,
    shared by every CLI entry point in this project (`sectorGen.py`,
    `systemGen.py`, `galaxyGen.py`, `planetgen.db.query`, `planetgen.cli.migrate`) instead
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
    parser.add_argument('--mysql-port', type=int,
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
"""dict: `MySQLConfig._key() -> PooledDB`, one pool per distinct set of
connection parameters seen so far in this process. A WSGI worker or a
single CLI invocation only ever needs one (the process-wide default), but
keying by config rather than keeping a single module-level pool lets
tests (and any future multi-database use) point at a second database
within the same process without the two stomping on each other's pooled
connections."""


STATEMENT_TIMEOUT_ERRORS = (1969, 3024)
"""tuple: The errors a statement stopped by `statement_timeout_s` raises:
MariaDB's `max_statement_time` and MySQL's `MAX_EXECUTION_TIME`."""

_mariadb_servers = {}


def _is_mariadb(config):
    """Whether `config`'s server is MariaDB (cached per host and port) --
    the two name their statement time limit differently."""
    server = (config.host, config.port)
    if server not in _mariadb_servers:
        probe = pymysql.connect(host=config.host, port=config.port, user=config.user,
                                password=config.password, connect_timeout=10)
        try:
            _mariadb_servers[server] = "mariadb" in probe.get_server_info().lower()
        finally:
            probe.close()
    return _mariadb_servers[server]


def _statement_timeout_sql(config, seconds):
    if _is_mariadb(config):
        return f"SESSION max_statement_time = {float(seconds):g}"
    return f"SESSION MAX_EXECUTION_TIME = {int(float(seconds) * 1000)}"


SQL_MODE_ENV_VAR = "PLANETGEN_MYSQL_SQL_MODE"
"""str: Env var that, when set, pins every pooled connection's session
`sql_mode` to its value (TEST.7). The test suite sets it to MySQL 8's
default (`STRICT_SQL_MODE`), so a GROUP BY or a truncation MariaDB's
looser default forgives fails locally too, not only in CI. Unset (every
deployment): the server's own `sql_mode` applies."""

STRICT_SQL_MODE = ("ONLY_FULL_GROUP_BY,STRICT_TRANS_TABLES,NO_ZERO_IN_DATE,NO_ZERO_DATE,"
                   "ERROR_FOR_DIVISION_BY_ZERO,NO_ENGINE_SUBSTITUTION")
"""str: MySQL 8.0's default `sql_mode`, valid on MariaDB too."""


def _session_sql_mode():
    mode = os.environ.get(SQL_MODE_ENV_VAR, "").strip()
    if mode and not re.fullmatch(r"[A-Za-z_,]+", mode):
        raise ValueError(f"{SQL_MODE_ENV_VAR} must be a comma-separated list of sql_mode names, not {mode!r}")
    return mode


def _get_pool(config, statement_timeout_s=None):
    key = config._key() if not statement_timeout_s else (config._key(), float(statement_timeout_s))
    if key not in _pools:
        init = "SET time_zone = '+00:00'"
        if statement_timeout_s:
            init += ", " + _statement_timeout_sql(config, statement_timeout_s)
        sql_mode = _session_sql_mode()
        if sql_mode:
            init += f", SESSION sql_mode = '{sql_mode}'"
        _pools[key] = PooledDB(
            creator=pymysql,
            mincached=1,
            maxcached=5,
            maxconnections=10,
            blocking=True,
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
            # Never fail over to a fresh connection mid-use. DBUtils
            # otherwise re-runs a statement that hit any OperationalError
            # -- a deadlock included -- on a new cursor or connection, and
            # the caller carries on in a transaction MySQL has already
            # rolled back (the error 1452s of PERF.14's parallel test). A
            # dead pooled connection is still replaced when it's checked
            # out (`ping`); one lost mid-use raises.
            isfatal=_never_fail_over,
        )
    return _pools[key]


def _never_fail_over(_error):
    return False


def close_pool(config):
    """
    Closes and discards the cached `PooledDB` for `config`, if one exists.

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
        _pools.pop(pool_key).close()
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


ID_BLOCK_TABLES = frozenset((
    "system_configs", "sectors", "star_systems", "stars", "planets", "moons",
    "asteroid_belts", "comets", "black_holes", "neutron_stars", "nebulae",
    "supernova_remnants", "rogue_planets", "interstellar_comets",
    "asteroid_fields", "quasars",
))
"""frozenset: Tables whose new rows take their `id` from `id_blocks`
(schema v45, PERF.13) rather than from AUTO_INCREMENT. Every plain
`INSERT INTO <one of these> (...) VALUES (...)` that names no `id`
column gets one added by `Connection.execute`, so a sector's children
know their parents' ids before anything is written, and `cur.lastrowid`
still answers as before."""

BATCH_CHILD_TABLES = frozenset((
    "system_config_slots", "planet_evolutionary_paragraphs", "moon_evolutionary_paragraphs",
    "planet_reflection_spectrum", "moon_reflection_spectrum", "asteroid_belt_composition",
    "comet_composition", "interstellar_comet_composition", "asteroid_field_composition",
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
    raw = _get_pool(config).connection()
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
        self._txn_locks = []
        self._short_row_waits = False

    def _run(self, sql, params):
        cur = self._conn.cursor()
        if not log.debug_log_active():
            cur.execute(sql.replace("?", "%s"), params)
            return _Cursor(cur)
        start = time.perf_counter()
        try:
            cur.execute(sql.replace("?", "%s"), params)
        except Exception as exc:
            log.trace(f"SQL failed after {(time.perf_counter() - start) * 1000:.2f}ms: {_sql_for_log(sql, params)} "
                      f"-> {type(exc).__name__}: {exc}", stacklevel=4)
            raise
        log.trace(f"SQL {(time.perf_counter() - start) * 1000:.2f}ms, {cur.rowcount} row(s): "
                  f"{_sql_for_log(sql, params)}", stacklevel=4)
        return _Cursor(cur)

    def execute(self, sql, params=()):
        shape = _insert_shape(sql)
        if shape is not None:
            new_id = None
            if shape.table in ID_BLOCK_TABLES and not shape.has_id and self._config is not None:
                new_id = _allocate_id(self._config, shape.table)
            if new_id is not None:
                shape = shape._replace(columns=["id", *shape.columns], values="(?, " + shape.values[1:], has_id=True)
                params = (new_id, *params)
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
        start = time.perf_counter()
        statements = 0
        for _position, ((table, columns, values, _level), rows) in order:
            head = f"INSERT INTO {table} ({', '.join(columns)}) VALUES "
            values = values.replace("?", "%s")
            chunk, size = [], len(head)
            for row in rows:
                # A cheap upper bound on the escaped literals -- 32 bytes a
                # number or NULL, 8 a character of text (4 UTF-8 bytes,
                # doubled by escaping) -- measured exactly only for a row
                # that might not fit at all.
                row_size = len(values) + 2 + 32 * len(row) + 8 * sum(map(len, filter(str.__instancecheck__, row)))
                if len(head) + row_size > budget:
                    row_size = len(cur.mogrify(values, row).encode()) + 2
                    if len(head) + row_size > budget:
                        raise pymysql.err.OperationalError(
                            1153, f"A {table} row of {row_size} bytes doesn't fit the server's "
                                  f"max_allowed_packet ({budget + _PACKET_MARGIN} bytes); nothing was sent for it")
                if chunk and (len(chunk) >= _BATCH_ROWS or size + row_size > budget):
                    cur.execute(head + ", ".join([values] * len(chunk)), [v for r in chunk for v in r])
                    statements += 1
                    chunk, size = [], len(head)
                chunk.append(row)
                size += row_size
            if chunk:
                cur.execute(head + ", ".join([values] * len(chunk)), [v for r in chunk for v in r])
                statements += 1
        if log.debug_log_active():
            log.trace(f"SQL batch {(time.perf_counter() - start) * 1000:.2f}ms: "
                      f"{sum(len(rows) for rows in batch.values())} row(s) in {statements} INSERT(s)", stacklevel=3)

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
                      f"row(s): {_sql_for_log(sql, rows[:3])}{' ...' if len(rows) > 3 else ''}", stacklevel=3)
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

    def commit(self):
        self.flush()
        self._conn.commit()
        self._release_txn_locks()

    def rollback(self):
        if self._batch:
            self._batch = {}
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
    conn = Connection(_get_pool(config, statement_timeout_s).connection(), config)
    if ensure_schema and config._key() not in _schema_ensured:
        try:
            _ensure_schema(conn)
        except BaseException:
            conn.close()
            raise
        _schema_ensured.add(config._key())
    return conn


DB_PREFIX_ENV_VAR = "PLANETGEN_MYSQL_DATABASE_PREFIX"
"""str: Env var overriding the default schema-name prefix `list_databases`/
`resolve_database` filter by -- see `MySQLConfig`'s own `database` default.
Shared by every entry point that offers a choice among several MySQL
schemas on one server (the `html/` CGI browser's `?db=` picker, and the
Flask API's own `?db=`/`/api/databases`, both via this one implementation).
Falls back to `config.json`'s `mysql.database_prefix` (see
`planetgen.util.appconfig`), then to `DEFAULT_DB_PREFIX`, when unset."""

DEFAULT_DB_PREFIX = "planetgen"
"""str: Matches `MySQLConfig`'s own default database name -- a deployment
with just one schema names it `planetgen` and never needs to set
`DB_PREFIX_ENV_VAR` (or `config.json`'s `mysql.database_prefix`) at all;
one with several names them `planetgen_<something>` to share the prefix."""


def escape_like(value):
    """
    Escapes `\\`, `%` and `_` in `value` so it matches only itself inside
    a SQL `LIKE` pattern -- paired with `ESCAPE '\\\\'` in the query (the
    same convention as `queryDb._search_like_pattern`). Without it a
    prefix like `planetgen_` would also match `planetgenX...`, and a `%`
    would match anything.
    """
    return value.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


def configured_control_database():
    """The control schema's name: `CONTROL_DB_ENV_VAR`, else
    `config.json`'s `control_database`, else `DEFAULT_CONTROL_DATABASE`
    (what `control_mysql_config` connects to)."""
    return os.environ.get(CONTROL_DB_ENV_VAR) or load_config()["control_database"] or DEFAULT_CONTROL_DATABASE


def list_databases(base_config=None, prefix=None):
    """
    Lists every MySQL schema on `base_config`'s server whose name starts
    with `prefix` (default: `DB_PREFIX_ENV_VAR`, or `DEFAULT_DB_PREFIX`) --
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
    prefix = (
        prefix
        or os.environ.get(DB_PREFIX_ENV_VAR)
        or load_config()["mysql"]["database_prefix"]
        or DEFAULT_DB_PREFIX
    )
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
            entries.append({
                "name": name,
                "size_bytes": int(stats["size_bytes"] or 0),
                "modified_at": modified_at.strftime("%Y-%m-%d %H:%M") if modified_at else "unknown",
            })
        return entries
    finally:
        conn.close()


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
        ValueError: If `name` is empty or doesn't match a listed schema.
    """
    if not name:
        raise ValueError("No database specified.")
    if name not in {entry["name"] for entry in list_databases(base_config, prefix=prefix)}:
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


def _ensure_schema(conn):
    """
    Applies `schema.sql` to `conn`, then bootstraps `schema_migrations`
    (inserting the baseline row) if it's empty -- idempotent, safe to call
    on a database that already has some, all, or none of the schema, and
    from several connections at once (`_schema_lock`).

    The baseline is `SCHEMA_VERSION` for a new database. An existing one
    whose `schema_migrations` was emptied or lost gets the version its
    tables show instead (`detect_schema_version`, DB.4), so
    `migrate_database` still runs the steps it is missing.

    Args:
        conn (Connection): The connection to apply the schema to.

    Raises:
        SchemaTooNewError: The database is already past `SCHEMA_VERSION`
            (checked before any DDL runs).
    """
    with _schema_lock(conn):
        _apply_schema(conn)


def _apply_schema(conn):
    """`_ensure_schema`'s work, for a caller already holding the schema
    lock."""
    stored = _stored_version(conn, "schema_migrations")
    _refuse_newer(conn, stored, SCHEMA_VERSION)
    baseline = SCHEMA_VERSION if stored is not None else detect_schema_version(conn)
    with open(SCHEMA_PATH, "r", encoding="utf-8") as f:
        conn.executescript(f.read())
    if conn._config is not None:
        forget_id_blocks(conn._config._key())

    row = conn.execute("SELECT COUNT(*) AS n FROM schema_migrations").fetchone()
    if row["n"] == 0:
        conn.execute("INSERT INTO schema_migrations (version) VALUES (?)", (baseline,))
        if baseline < SCHEMA_VERSION:
            log.normal(f"Database has no schema version recorded; its tables match v{baseline}, "
                       f"so that is recorded and planetgen.cli.migrate will bring it up to v{SCHEMA_VERSION}.")
            activity_log.event("DB", "detect_version", db=conn._config.database if conn._config else "?",
                              version=baseline)
    conn.commit()


def _column_marker(table, column):
    return lambda shape: (table, column) in shape["columns"]


def _index_marker(table, index):
    return lambda shape: (table, index) in shape["indexes"]


def _table_marker(table):
    return lambda shape: table in shape["tables"]


_VERSION_MARKERS = (
    (53, _table_marker("sector_stats")),
    (52, _table_marker("generation_runs")),
    (51, _column_marker("galaxy_shape", "galaxy_seed")),
    (50, _column_marker("system_configs", "comets")),
    (49, _table_marker("bright_star_blocks")),
    (48, _column_marker("rogue_planets", "surface_regime")),
    (47, _column_marker("rogue_planets", "planet_class")),
    (46, _index_marker("sectors", "ft_sectors_name")),
    (45, _table_marker("id_blocks")),
    (44, _table_marker("species")),
    (43, _column_marker("galaxy_shape", "bright_star_min_luminosity_sol")),
    (42, _table_marker("facilities")),
    (41, _column_marker("asteroid_fields", "quadrant")),
    (40, _column_marker("system_name_registry", "first_object_table")),
    (39, _column_marker("asteroid_fields", "inside_nebula_id")),
    (38, _column_marker("nebulae", "nebula_class")),
    (37, _column_marker("star_systems", "runaway_class")),
    (36, _column_marker("black_holes", "mass_class")),
    # v35 changed no table, only which galaxy sectors exist; v34 and v35
    # look the same. Read as v35: re-running its step would delete sectors
    # placed under v35's slot rule, while skipping it at most leaves a v34
    # database's old-rule sectors in place.
    (35, lambda shape: "sector_name_registry" in shape["tables"] and "body_name_registry" not in shape["tables"]),
    (33, _table_marker("galaxy_layer")),
    (32, _column_marker("sectors", "ring_index")),
    (31, _table_marker("quasars")),
    # v30 only cleaned data, so v29 and v30 look the same; v30's step is
    # safe to repeat, so read them as v29.
    (29, lambda shape: ("rogue_planets", "center_x_pc") in shape["columns"]
         and ("star_systems", "wikitext_content") not in shape["columns"]),
    (28, _column_marker("rogue_planets", "center_x_pc")),
    (27, _column_marker("sectors", "created_at")),
    (26, _index_marker("nebulae", "idx_nebulae_center")),
    (25, _index_marker("sectors", "idx_sectors_center")),
    (24, _table_marker("sector_name_registry")),
    (23, _column_marker("sectors", "wiki_url")),
    (22, _index_marker("planets", "idx_planets_name")),
    (21, _column_marker("black_holes", "sector_id")),
    (20, _column_marker("planets", "reflex_offset_x_km")),
    (19, _table_marker("comets")),
    (18, _column_marker("nebulae", "center_x_pc")),
    (17, _column_marker("black_holes", "galactic_orbital_speed_kms")),
    (16, _table_marker("black_holes")),
    (15, _column_marker("stars", "wide_binary_a_crit_km")),
    (14, _column_marker("star_systems", "binary_mutual_position_x_km")),
    (13, _column_marker("stars", "galactic_orbital_phase_deg")),
    (12, _column_marker("planets", "min_update_interval_years")),
    (11, _column_marker("planets", "position_x_km")),
    (10, _column_marker("stars", "galactic_orbital_speed_kms")),
    (9, _column_marker("planets", "orbital_inclination_deg")),
)
"""tuple: `(version, test)` pairs, newest first, for
`detect_schema_version`: each test is true of a database's shape from that
version on (something the version added that no later one removed). Every
new migration step adds its own row at the top."""

_OLDEST_DETECTED_VERSION = 8
"""int: What `detect_schema_version` reports for a database that has
tables but none of `_VERSION_MARKERS` -- the oldest version
`migrate_database` has steps from."""


def detect_schema_version(conn):
    """
    The galaxy schema version the connected database's tables show, for a
    database whose `schema_migrations` is missing or empty (DB.4): the
    newest version in `_VERSION_MARKERS` whose mark is present, from
    `information_schema` alone (reads only, before any DDL runs).
    `SCHEMA_VERSION` for a database with no galaxy tables yet (a new one).

    Returns:
        int: The version.
    """
    tables = {row["t"].lower() for row in conn.execute(
        "SELECT table_name AS t FROM information_schema.tables WHERE table_schema = DATABASE()").fetchall()}
    if "star_systems" not in tables:
        return SCHEMA_VERSION
    shape = {
        "tables": tables,
        "columns": {(row["t"].lower(), row["c"].lower()) for row in conn.execute(
            "SELECT table_name AS t, column_name AS c FROM information_schema.columns"
            " WHERE table_schema = DATABASE()").fetchall()},
        "indexes": {(row["t"].lower(), row["i"].lower()) for row in conn.execute(
            "SELECT DISTINCT table_name AS t, index_name AS i FROM information_schema.statistics"
            " WHERE table_schema = DATABASE()").fetchall()},
    }
    for version, marked in _VERSION_MARKERS:
        if marked(shape):
            return version
    return _OLDEST_DETECTED_VERSION


def open_write(config=None):
    """
    Opens a connection for a write-capable caller against an already-
    existing content database (`ensure_schema=False` -- same reasoning as
    `queryDb.open_readonly`: if the account passed in (see
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
            `CONTROL_DB_ENV_VAR` (or `DEFAULT_CONTROL_DATABASE`).
    """
    base_config = base_config or DEFAULT_MYSQL_CONFIG
    database = configured_control_database()
    return MySQLConfig(
        host=base_config.host, port=base_config.port,
        user=base_config.user, password=base_config.password, database=database,
    )


def get_control_connection(config=None, ensure_schema=False):
    """
    Opens a pooled connection to the control schema.

    Args:
        config (MySQLConfig, optional): Connection parameters. Defaults
            to `control_mysql_config()`.
        ensure_schema (bool): Whether to run `_ensure_control_schema`
            (DDL) on this connection -- `False` by default (the normal
            runtime case: the Flask API's account may have no `CREATE`
            grant, same reasoning as `open_write` above).
            `adminAuth.bootstrap_control_schema` passes `True`, using a
            full-access account (the same one `planetgen.cli.migrate` already
            uses), to create/update the control schema once per deploy.

    Returns:
        Connection: An open connection.
    """
    config = config or control_mysql_config()
    conn = Connection(_get_pool(config).connection(), config)
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
    with open(CONTROL_SCHEMA_PATH, "r", encoding="utf-8") as f:
        conn.executescript(f.read())
    _add_control_columns(conn)

    row = conn.execute("SELECT MAX(version) AS v FROM control_schema_migrations").fetchone()
    if row["v"] is None or row["v"] < CONTROL_SCHEMA_VERSION:
        conn.execute("INSERT INTO control_schema_migrations (version) VALUES (?)", (CONTROL_SCHEMA_VERSION,))
        conn.commit()
        if row["v"] is not None:
            activity_log.event("DB", "migrate", db=configured_control_database(), from_version=row["v"],
                              to_version=CONTROL_SCHEMA_VERSION)


def _add_control_columns(conn):
    """
    Adds `_CONTROL_COLUMNS` (and the job tree's index and parent key) to
    control tables created before them -- each only when missing, so it
    is safe to run on every deploy, on MySQL as well as MariaDB (MySQL
    has no `ADD COLUMN IF NOT EXISTS`).
    """
    def existing(sql, table):
        return {row["name"] for row in conn.execute(sql, (table,)).fetchall()}

    for table, columns in _CONTROL_COLUMNS.items():
        have = existing("SELECT COLUMN_NAME AS name FROM information_schema.COLUMNS"
                        " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", table)
        missing = [f"ADD COLUMN {name} {definition}" for name, definition in columns if name not in have]
        if missing:
            conn.execute(f"ALTER TABLE {table} {', '.join(missing)}")
    if "idx_work_jobs_root" not in existing("SELECT INDEX_NAME AS name FROM information_schema.STATISTICS"
                                            " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", "work_jobs"):
        conn.execute("ALTER TABLE work_jobs ADD KEY idx_work_jobs_root (root_id)")
    if "fk_work_jobs_parent" not in existing(
            "SELECT CONSTRAINT_NAME AS name FROM information_schema.TABLE_CONSTRAINTS"
            " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", "work_jobs"):
        conn.execute("ALTER TABLE work_jobs ADD CONSTRAINT fk_work_jobs_parent"
                     " FOREIGN KEY (parent_id) REFERENCES work_jobs(id) ON DELETE CASCADE")
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
            markdown, habitable_world, asteroid_belt, comets, large_star, moons,
            max_planets, planets, star_type, name, age, intelligent_life,
            binary_system, wide_binary, num_orbits
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            int(bool(config.MARKDOWN)),
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
# planets and moons are named from their system (bodyNames.py, v34), so
# they need no registry of their own. See nameUniqueness.py's own module
# docstring for the full sector > system hierarchy.
# ---------------------------------------------------------------------------

def _rename_bodies_with_prefix(conn, star_system_id, old_prefix, new_prefix):
    """
    Renames every star, planet and moon in one system whose name carries
    `old_prefix` (`bodyNames.rename_prefix`), so derived names such as
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
    come from their star's (`bodyNames.py`), so they're unique whenever
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
    sector-vs-sector collision (`nameUniqueness.resolve_greek_roman_collision`),
    then a cross-level collision against an existing system's base name
    (renaming *that* system instead of this sector's own name -- see
    `_rename_existing_system_for_diminutive`). Draws an entirely fresh
    candidate (`stellarObjects.utils.generate_sector_name`) and starts
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
        sector_id (int): The new (or, for `planetgen.cli.dedupe`, existing)
            `sectors.id` -- recorded as the base name's first holder when
            it is one.

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
            conn.execute("UPDATE sectors SET name = ? WHERE name = ? AND id <> ?", (renamed_to, old_name, sector_id))

        return new_name, base


def _regenerate_star_name():
    """A fresh star-name candidate, for `reserve_system_names`'
    exhaustion fallback -- same generator `starData.Star` itself uses."""
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
    name (`bodyNames.py`), so they're unique whenever it is. A name whose
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
        sector_hits = set()
        for first in range(0, len(keys), _NAME_BATCH):
            chunk = [spelled[key] for key in keys[first:first + _NAME_BATCH]]
            sector_hits.update(_name_key(row["base_name"]) for row in conn.execute(
                f"SELECT base_name FROM sector_name_registry WHERE base_name IN ({', '.join('?' * len(chunk))})",
                tuple(chunk),
            ).fetchall())
        for first in range(0, len(keys), _NAME_BATCH):
            chunk = keys[first:first + _NAME_BATCH]
            conn.execute(
                "INSERT INTO system_name_registry (base_name, occurrence_count) VALUES "
                + ", ".join(["(?, ?)"] * len(chunk))
                + " ON DUPLICATE KEY UPDATE occurrence_count = occurrence_count + VALUES(occurrence_count)",
                tuple(value for key in chunk for value in (spelled[key], counts[key])),
            )
        rows = _registry_rows(conn, [spelled[key] for key in keys])

        # Keys the collation still counts as one name (a spelling
        # `_name_key` folds differently) land on one registry row: they
        # are one name's holders, counted together.
        by_row = {}
        for key in keys:
            row = rows.get(key)
            if row is None:  # stored under a spelling _name_key() doesn't match
                row = conn.execute(
                    "SELECT id, base_name, occurrence_count, diminutive_index, first_star_system_id, "
                    "first_object_table, first_object_id FROM system_name_registry WHERE base_name = ?",
                    (spelled[key],),
                ).fetchone()
            by_row.setdefault(row["id"], (row, set()))[1].add(key)

        retry = []
        for row, row_keys in by_row.values():
            uses = [i for i in todo if key_of[i] in row_keys]
            existing_before = row["occurrence_count"] - len(uses)
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
    by_base = {}
    for base, system_id, table, object_id, diminutive_index in confirmations:
        key = _name_key(base)
        if key in by_base:
            by_base[key][4] = diminutive_index
        else:
            by_base[key] = [base, system_id, table, object_id, diminutive_index]
    rows = list(by_base.values())
    vacant = "first_star_system_id IS NULL AND first_object_id IS NULL"
    for first in range(0, len(rows), _NAME_BATCH):
        chunk = rows[first:first + _NAME_BATCH]
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
    "rogue-planet": "rogue_planets",
    "black-hole": "black_holes",
    "neutron-star": "neutron_stars",
    "nebula": "nebulae",
    "supernova-remnant": "supernova_remnants",
    "quasar": "quasars",
    "comet": "interstellar_comets",
    "asteroid-field": "asteroid_fields",
    "bright-star": "star_systems",
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


def _claim_object_ids(conn, items):
    """
    Names each object by its 64-bit object ID (GEN.64, `objectId`):
    `items` is `(obj, kind, center_pc)` triples, in generation order. An
    ID already held -- by an earlier item here, or by a stored row of the
    same table -- takes the next collision number (`objectId.bump`) until
    it is free, so the
    first object in generation order keeps collision number 0. Sets each
    `obj.name`, and when `insert_sector` holds a reservation map
    (`prereserved_names`), records each object there as needing no
    registry entry.
    """
    by_table = {}
    for obj, kind, center in items:
        by_table.setdefault(OBJECT_ID_TABLES[kind], []).append((obj, objectId.pack(kind, center)))
    prereserved = getattr(conn, "prereserved_names", None)
    for table, entries in by_table.items():
        checked = {objectId.format_id(object_id) for _obj, object_id in entries}
        stored = _stored_names(conn, table, sorted(checked))
        taken = set()
        for obj, object_id in entries:
            name = objectId.format_id(object_id)
            while True:
                if name not in checked:  # a bumped ID: ask about it too
                    checked.add(name)
                    stored |= _stored_names(conn, table, [name])
                if name not in taken and name not in stored:
                    break
                object_id = objectId.bump(object_id)
                name = objectId.format_id(object_id)
            taken.add(name)
            obj.name = name
            if prereserved is not None:
                prereserved[id(obj)] = (None, None)


def _stored_names(conn, table, names):
    """The `names` (upper-cased) `table` already has a row for."""
    found = set()
    for first in range(0, len(names), _NAME_BATCH):
        chunk = names[first:first + _NAME_BATCH]
        found.update(row["name"].upper() for row in conn.execute(
            f"SELECT name FROM {table} WHERE name IN ({', '.join('?' * len(chunk))})", tuple(chunk),
        ).fetchall())
    return found


def insert_star(conn, star, star_system_id, role) -> int:
    """
    Inserts a `stars` row for one individual `Star` (never a
    `BinaryStarProxy` -- see `insert_star_system`).

    Args:
        conn (Connection): An open, schema-initialized connection.
        star: The `Star` instance.
        star_system_id (int): The owning `star_systems.id`.
        role (str): `'primary'`, `'secondary'`, or `'single'`.

    Returns:
        int: The new `stars.id`.
    """
    wide_binary_a_crit_km = (
        star.a_crit_au * physical_constants.AU_TO_KM if star.a_crit_au is not None else None
    )
    cur = conn.execute(
        """
        INSERT INTO stars (
            star_system_id, role, name, star_type, yerkes_class, mass_kg, radius_km,
            temperature_k, luminosity_w, age_gy, lifespan_gy,
            habitable_zone_inner_km, habitable_zone_outer_km,
            system_perimeter_km, heliosphere_radius_km,
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years,
            wide_binary_a_crit_km,
            reflex_offset_x_km, reflex_offset_y_km, reflex_offset_z_km
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            star_system_id, role, star.name, star.type, star.yerkes_class,
            star.mass, star.radius, star.temperature, star.luminosity, star.age,
            _lifespan_gy(star.lifespan),
            star.habitable_zone[0] * physical_constants.AU_TO_KM,
            star.habitable_zone[1] * physical_constants.AU_TO_KM,
            star.system_perimeter * physical_constants.AU_TO_KM,
            star.heliosphere_radius * physical_constants.AU_TO_KM,
            star.galactic_orbital_speed_kms,
            star.galactic_orbital_period_gy,
            star.galactic_orbital_phase_deg,
            star.galactic_min_update_interval_years,
            wide_binary_a_crit_km,
            star.reflex_offset_x * physical_constants.AU_TO_KM,
            star.reflex_offset_y * physical_constants.AU_TO_KM,
            star.reflex_offset_z * physical_constants.AU_TO_KM,
        ),
    )
    return cur.lastrowid


def _insert_paragraphs(conn, table, id_column, owner_id, paragraphs):
    """Shared body for `planet_evolutionary_paragraphs`/
    `moon_evolutionary_paragraphs` -- identical shape, different owning
    table/column. `table`/`id_column` are always one of this module's own
    literal constants below, never request-derived, so the f-string is safe."""
    for position, paragraph in enumerate(paragraphs or []):
        conn.execute(
            f"INSERT INTO {table} ({id_column}, position, paragraph) VALUES (?, ?, ?)",
            (owner_id, position, paragraph),
        )


def _insert_reflection_spectrum(conn, table, id_column, owner_id, visible, non_visible):
    """Shared body for `planet_reflection_spectrum`/
    `moon_reflection_spectrum` -- see `_insert_paragraphs`."""
    for spectrum_type, values in (("visible", visible), ("non_visible", non_visible)):
        for position, value in enumerate(values or []):
            conn.execute(
                f"INSERT INTO {table} ({id_column}, spectrum_type, position, value) VALUES (?, ?, ?, ?)",
                (owner_id, spectrum_type, position, value),
            )


_BODY_COLUMNS = (
    "body_type", "name",
    "planet_class", "distance_km", "radius_km", "mass_kg", "volume_km3", "period_years", "zone",
    "description", "gravity_g", "surface_temperature_k", "density_g_cm3", "atmosphere",
    "atm_density", "atm_molar_density", "atmospheric_pressure_pa", "composition",
    "scale_height_km", "hill_radius_km", "min_orbit_distance_km",
    "habitable_zone_inner_km", "habitable_zone_outer_km",
    "life_chemical", "evolutionary_speed", "flavor_text", "flavor_text_count",
    "orbital_inclination_deg", "orbital_ascending_node_deg", "orbital_phase_deg",
    "position_x_km", "position_y_km", "position_z_km", "orbital_speed_kms",
    "min_update_interval_years",
    "rotation_period_hours",
)
"""tuple: The generated-content columns `planets` and `moons` share, in
`body_row_values` order."""

_PLANET_ONLY_COLUMNS = ("reflex_offset_x_km", "reflex_offset_y_km", "reflex_offset_z_km")
"""tuple: The `planets` columns `moons` lacks (a moon hosts no moons)."""


def body_row_values(body):
    """
    A planet's or moon's generated-content column values, in
    `_BODY_COLUMNS` order (plus `_PLANET_ONLY_COLUMNS` for a planet), with
    every unit conversion `load_star_system` inverts: what
    `insert_planet`/`insert_moon` write and `editStore` updates in place.
    """
    min_orbit_distance_km = (
        body.min_orbit_distance * physical_constants.AU_TO_KM
        if body.min_orbit_distance is not None else None
    )
    values = [
        body.body_type, body.name,
        body.planet_class,
        body.distance * physical_constants.AU_TO_KM,
        body.radius, body.mass, body.volume, body.period, body.zone,
        body.description, body.gravity, body.surface_temperature,
        body.density, body.atmosphere,
        body.atm_density, body.atm_molar_density, body.atmospheric_pressure,
        body.composition,
        body.scale_height, body.hill_radius, min_orbit_distance_km,
        body.habitable_zone[0] * physical_constants.AU_TO_KM,
        body.habitable_zone[1] * physical_constants.AU_TO_KM,
        body.life_chemical, body.evolutionary_speed,
        body.flavor_text, body.flavor_text_count,
        body.orbital_inclination_deg, body.orbital_ascending_node_deg,
        body.orbital_phase_deg,
        body.position_x * physical_constants.AU_TO_KM,
        body.position_y * physical_constants.AU_TO_KM,
        body.position_z * physical_constants.AU_TO_KM,
        body.orbital_speed_kms,
        body.min_update_interval_years,
        body.rotation_period_hours,
    ]
    if not body.is_moon:
        values += [
            body.reflex_offset_x * physical_constants.AU_TO_KM,
            body.reflex_offset_y * physical_constants.AU_TO_KM,
            body.reflex_offset_z * physical_constants.AU_TO_KM,
        ]
    return values


def body_child_rows(conn, body, body_id):
    """Writes a planet's or moon's evolutionary paragraphs and reflection
    spectrum rows (after any old ones were deleted)."""
    prefix, id_column = ("moon", "moon_id") if body.is_moon else ("planet", "planet_id")
    _insert_paragraphs(conn, f"{prefix}_evolutionary_paragraphs", id_column, body_id, body.evolutionary_data)
    _insert_reflection_spectrum(
        conn, f"{prefix}_reflection_spectrum", id_column, body_id,
        body.reflection_spectrum_visible, body.reflection_spectrum_non_visible,
    )


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
    columns = ("star_system_id", "star_id", "orbital_index") + _BODY_COLUMNS + _PLANET_ONLY_COLUMNS
    cur = conn.execute(
        f"INSERT INTO planets ({', '.join(columns)}) VALUES ({', '.join('?' * len(columns))})",
        (star_system_id, star_id, orbital_index, *body_row_values(planet)),
    )
    planet_id = cur.lastrowid

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
    columns = ("planet_id", "star_system_id", "star_id", "orbital_index") + _BODY_COLUMNS
    cur = conn.execute(
        f"INSERT INTO moons ({', '.join(columns)}) VALUES ({', '.join('?' * len(columns))})",
        (planet_id, star_system_id, star_id, orbital_index, *body_row_values(moon)),
    )
    moon_id = cur.lastrowid

    body_child_rows(conn, moon, moon_id)

    return moon_id


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
    cur = conn.execute(
        """
        INSERT INTO asteroid_belts (
            star_system_id, star_id, orbital_index, distance_km, lower_limit_km, upper_limit_km,
            density, composition_summary
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            star_system_id, star_id, orbital_index,
            belt.distance * physical_constants.AU_TO_KM,
            belt.lower_limit * physical_constants.AU_TO_KM,
            belt.upper_limit * physical_constants.AU_TO_KM,
            belt.density,
            belt.get_composition_summary(),
        ),
    )
    belt_id = cur.lastrowid

    for position, (component, concentration) in enumerate(belt.composition):
        conn.execute(
            """
            INSERT INTO asteroid_belt_composition (belt_id, position, component, concentration)
            VALUES (?, ?, ?, ?)
            """,
            (belt_id, position, component, concentration),
        )

    return belt_id


def insert_comet(conn, comet: Comet, star_system_id, star_id=None) -> int:
    """
    Inserts a `comets` row (plus its `comet_composition` child rows).

    Modeled directly on `insert_asteroid_belt` above -- a `Comet` has no
    `orbital_index` (it isn't part of the planets/belts orbital-slot
    ordering; see `systemData.StarSystem._generate_comets`'s docstring),
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
    cur = conn.execute(
        """
        INSERT INTO comets (
            star_system_id, star_id, name, orbit_type, period_class,
            nucleus_diameter_km, composition_summary, perihelion_distance_km,
            eccentricity, inclination_deg, arg_periapsis_deg, ascending_node_deg,
            orbital_period_years, mean_anomaly_deg, parabolic_mean_anomaly,
            min_update_interval_years, primary_mass_solar, is_active,
            distance_km, position_x_km, position_y_km, position_z_km, orbital_speed_kms
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            star_system_id, star_id, comet.name, comet.orbit_type, comet.period_class,
            comet.nucleus_diameter_km, comet.get_composition_summary(),
            comet.perihelion_distance_au * physical_constants.AU_TO_KM,
            comet.eccentricity, comet.inclination_deg, comet.arg_periapsis_deg, comet.ascending_node_deg,
            comet.orbital_period_years, comet.mean_anomaly_deg, comet.parabolic_mean_anomaly,
            comet.min_update_interval_years, comet.primary_mass_solar, int(comet.is_active),
            comet.distance_au * physical_constants.AU_TO_KM,
            comet.position_x_au * physical_constants.AU_TO_KM,
            comet.position_y_au * physical_constants.AU_TO_KM,
            comet.position_z_au * physical_constants.AU_TO_KM,
            comet.orbital_speed_kms,
        ),
    )
    comet_id = cur.lastrowid

    for position, component in enumerate(comet.composition):
        conn.execute(
            "INSERT INTO comet_composition (comet_id, position, component) VALUES (?, ?, ?)",
            (comet_id, position, component),
        )

    return comet_id


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
            was generated as part of (`sectorGen.generate_sector_phenomena`)
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
    galactic_fields = (
        (None, None, None, None) if star_id is not None else (
            black_hole.galactic_orbital_speed_kms, black_hole.galactic_orbital_period_gy,
            black_hole.galactic_orbital_phase_deg, black_hole.galactic_min_update_interval_years,
        )
    )
    placement = placement or {}
    cur = conn.execute(
        """
        INSERT INTO black_holes (
            star_id, sector_id, name, mass_class, mass_solar, event_horizon_radius_km, spin,
            has_accretion_disk, temperature_k, luminosity_w, age_gy,
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years,
            center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            star_id, sector_id, black_hole.name, black_hole.mass_class, black_hole.mass_solar,
            black_hole.event_horizon_radius_km, black_hole.spin,
            int(black_hole.has_accretion_disk), black_hole.temperature,
            black_hole.luminosity, black_hole.age,
            *galactic_fields,
            placement.get("center_x_pc"), placement.get("center_y_pc"),
            placement.get("center_z_pc"), placement.get("galactic_radius_pc"),
        ),
    )
    if name_base is not None:
        confirm_object_name(conn, name_base, "black_holes", cur.lastrowid, diminutive_index)
    return cur.lastrowid


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
    galactic_fields = (
        (None, None, None, None) if star_id is not None else (
            neutron_star.galactic_orbital_speed_kms, neutron_star.galactic_orbital_period_gy,
            neutron_star.galactic_orbital_phase_deg, neutron_star.galactic_min_update_interval_years,
        )
    )
    placement = placement or {}
    cur = conn.execute(
        """
        INSERT INTO neutron_stars (
            star_id, sector_id, name, mass_solar, radius_km, spin_period_ms, magnetic_field_gauss,
            pulsar_type, surface_temperature_k, luminosity_w, age_gy,
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years,
            center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            star_id, sector_id, neutron_star.name, neutron_star.mass_solar, neutron_star.radius,
            neutron_star.spin_period_ms, neutron_star.magnetic_field_gauss,
            neutron_star.pulsar_type, neutron_star.surface_temperature_k,
            neutron_star.luminosity, neutron_star.age,
            *galactic_fields,
            placement.get("center_x_pc"), placement.get("center_y_pc"),
            placement.get("center_z_pc"), placement.get("galactic_radius_pc"),
        ),
    )
    if name_base is not None:
        confirm_object_name(conn, name_base, "neutron_stars", cur.lastrowid, diminutive_index)
    return cur.lastrowid


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
            dominant_species, density_cm3, temperature_k, extinction_av,
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years,
            center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, nebula.name, nebula.nebula_class, nebula.nebula_type, nebula.radius_ly,
            nebula.composition, nebula.formation_cause,
            nebula.dominant_species, nebula.density_cm3, nebula.temperature_k, nebula.extinction_av,
            nebula.galactic_orbital_speed_kms, nebula.galactic_orbital_period_gy,
            nebula.galactic_orbital_phase_deg, nebula.galactic_min_update_interval_years,
            placement.get("center_x_pc"), placement.get("center_y_pc"),
            placement.get("center_z_pc"), placement.get("galactic_radius_pc"),
        ),
    )
    confirm_object_name(conn, name_base, "nebulae", cur.lastrowid, diminutive_index)
    _refresh_containment_around(conn, placement, nebula.radius_ly)
    return cur.lastrowid


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
            as part of (`sectorGen.generate_sector_phenomena`) or placed
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
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years,
            center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, remnant.name, remnant.remnant_class, remnant.morphology, remnant.age_years,
            remnant.radius_ly, remnant.progenitor_type, compact_remnant_kind, black_hole_id, neutron_star_id,
            remnant.dominant_species, remnant.density_cm3, remnant.temperature_k, remnant.extinction_av,
            remnant.galactic_orbital_speed_kms, remnant.galactic_orbital_period_gy,
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
    cur = conn.execute(
        f"""
        INSERT INTO rogue_planets (
            sector_id, name, planet_type, planet_class, mass_bin, mass_kg, radius_km, composition,
            has_internal_heat, has_moons, {", ".join(ROGUE_SURFACE_FIELDS)},
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years,
            center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, {", ".join("?" * len(ROGUE_SURFACE_FIELDS))},
                  ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, planet.name, planet.planet_type, getattr(planet, "planet_class", None), planet.mass_bin, planet.mass_kg, planet.radius_km,
            planet.composition, int(planet.has_internal_heat), int(planet.has_moons),
            *_rogue_surface_values(planet),
            planet.galactic_orbital_speed_kms, planet.galactic_orbital_period_gy,
            planet.galactic_orbital_phase_deg, planet.galactic_min_update_interval_years,
            *_placement_values(placement),
        ),
    )
    confirm_object_name(conn, name_base, "rogue_planets", cur.lastrowid, diminutive_index)
    return cur.lastrowid


def sector_code(conn, sector_id):
    """
    How a designation names a sector (v40): its grid designation
    (`galaxyGeometry.provisional_sector_designation`) when it sits in the
    cylindrical grid, else its name. `None` for no sector.
    """
    if sector_id is None:
        return None
    row = conn.execute(
        "SELECT name, ring_index, layer_index, ring_slot_index FROM sectors WHERE id = ?", (sector_id,),
    ).fetchone()
    if row is None:
        return None
    if None not in (row["ring_index"], row["layer_index"], row["ring_slot_index"]):
        try:
            return provisional_sector_designation(row["ring_index"], row["layer_index"], row["ring_slot_index"])
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
    cur = conn.execute(
        """
        INSERT INTO interstellar_comets (
            sector_id, name, nucleus_diameter_km, velocity_kms, is_active, composition_summary,
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years,
            center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, comet.name, comet.nucleus_diameter_km, comet.velocity_kms,
            int(comet.is_active), comet.get_composition_summary(),
            comet.galactic_orbital_speed_kms, comet.galactic_orbital_period_gy,
            comet.galactic_orbital_phase_deg, comet.galactic_min_update_interval_years,
            *_placement_values(placement),
        ),
    )
    comet_id = cur.lastrowid

    for position, component in enumerate(comet.composition):
        conn.execute(
            "INSERT INTO interstellar_comet_composition (comet_id, position, component) VALUES (?, ?, ?)",
            (comet_id, position, component),
        )

    return comet_id


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
    placement = placement or {}
    cur = conn.execute(
        """
        INSERT INTO asteroid_fields (
            sector_id, name, field_class, composition_family, density, radius_ly, composition_summary,
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years,
            center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, field.name, field.field_class, field.composition_family, field.density,
            field.radius_ly, field.get_composition_summary(),
            field.galactic_orbital_speed_kms, field.galactic_orbital_period_gy,
            field.galactic_orbital_phase_deg, field.galactic_min_update_interval_years,
            placement.get("center_x_pc"), placement.get("center_y_pc"),
            placement.get("center_z_pc"), placement.get("galactic_radius_pc"),
        ),
    )
    field_id = cur.lastrowid

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
            center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, quasar.name, quasar.black_hole_mass_solar, quasar.event_horizon_radius_km,
            quasar.eddington_ratio, quasar.luminosity_w, quasar.accretion_rate_solar_per_year,
            quasar.broad_line_region_light_days, int(quasar.is_radio_loud), quasar.jet_length_ly,
            quasar.active_age_years,
            *_placement_values(GALACTIC_CENTER_PLACEMENT if placement is not None else None),
        ),
    )
    confirm_object_name(conn, name_base, "quasars", cur.lastrowid, diminutive_index)
    return cur.lastrowid


def _check_quasar_sector(conn, sector_id):
    """
    Refuses to place a quasar "at" a sector that doesn't host the galactic
    nucleus (every ring-0, layer-0 cell touches the origin), or into a
    galaxy that already has one -- a galaxy has a single nucleus.

    Raises:
        ValueError: If either rule would be broken.
    """
    row = conn.execute("SELECT ring_index, layer_index FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    if row is None or row["ring_index"] != 0 or row["layer_index"] != 0:
        raise ValueError(
            f"a quasar can only be placed in a ring-0, layer-0 sector (the galactic core); "
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
            call for a `sectorGen.generate_sector_phenomena`-generated
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
                return insert_star_system(conn, phenomenon, system_config)
            inserter = _PHENOMENON_INSERTERS.get(phenomenon_type)
            if inserter is None:
                raise ValueError(f"Unknown phenomenon type: {phenomenon_type!r}")
            if phenomenon_type == "quasar" and sector_id is not None:
                _check_quasar_sector(conn, sector_id)
            placement = compute_phenomenon_placement(conn, sector_id) if sector_id is not None else None
            row_id = inserter(conn, phenomenon, sector_id=sector_id, placement=placement)
            if sector_id is not None:
                refresh_containment(conn, [sector_id])
                refresh_nearest_systems(conn, [sector_id])
                # The galaxy's content stamp (`queryDb.galaxy_content_state`)
                # only sees sectors and systems; without this, a phenomenon
                # added from the command line never reaches a cached tile
                # or page (TEST.42).
                touch_sector(conn, sector_id)
            return row_id
    finally:
        conn.close()


_NULL_PROXY_ONLY_BINARY_FIELDS = (None,) * 15
"""Placeholder for the 15 `star_systems.binary_*` columns that only ever
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
    Extracts the 15 `star_systems.binary_*` column values that only ever
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
        proxy.habitable_zone[0] * physical_constants.AU_TO_KM,
        proxy.habitable_zone[1] * physical_constants.AU_TO_KM,
        proxy.system_perimeter * physical_constants.AU_TO_KM,
        proxy.heliosphere_radius * physical_constants.AU_TO_KM,
        proxy.galactic_orbital_speed_kms,
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
        proxy.binary_mutual_position_x * physical_constants.AU_TO_KM,
        proxy.binary_mutual_position_y * physical_constants.AU_TO_KM,
        proxy.binary_mutual_position_z * physical_constants.AU_TO_KM,
        proxy.binary_primary_position_x * physical_constants.AU_TO_KM,
        proxy.binary_primary_position_y * physical_constants.AU_TO_KM,
        proxy.binary_primary_position_z * physical_constants.AU_TO_KM,
        proxy.binary_secondary_position_x * physical_constants.AU_TO_KM,
        proxy.binary_secondary_position_y * physical_constants.AU_TO_KM,
        proxy.binary_secondary_position_z * physical_constants.AU_TO_KM,
        proxy.binary_secondary_mass_fraction,
    )


def _mutual_orbit_fields_from_wide_binary(pair):
    """
    The same 16 `star_systems.binary_mutual_*`/`binary_primary_position_*`/
    `binary_secondary_position_*`/`binary_secondary_mass_fraction` columns
    as `_mutual_orbit_fields_from_proxy`, from a 'wide' pair's
    `wideBinary.WideBinaryPair` instead -- these columns are reused
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
        pair.position_x_au * physical_constants.AU_TO_KM,
        pair.position_y_au * physical_constants.AU_TO_KM,
        pair.position_z_au * physical_constants.AU_TO_KM,
        pair.primary_position_x_au * physical_constants.AU_TO_KM,
        pair.primary_position_y_au * physical_constants.AU_TO_KM,
        pair.primary_position_z_au * physical_constants.AU_TO_KM,
        pair.secondary_position_x_au * physical_constants.AU_TO_KM,
        pair.secondary_position_y_au * physical_constants.AU_TO_KM,
        pair.secondary_position_z_au * physical_constants.AU_TO_KM,
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


def _location_for_entry(sector: SpaceSector, entry: SectorSystemEntry) -> str:
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

    Returns:
        str: See `_format_location_string`.
    """
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


def insert_star_system(conn, star_system: StarSystem, system_config: SystemConfig,
                        sector_id=None, position=None, location=None) -> int:
    """
    Inserts a full `StarSystem` -- the `star_systems` row, its `stars` row(s),
    and every planet/moon/asteroid belt/comet it contains -- into the database.

    No page text is stored (v29): wikitext/Markdown are rendered on demand
    from these rows by `planetgen/db/render.py` -- see
    `schema.sql`'s "v29" header note.

    `star_system.name` is reserved via `reserve_system_name` (v24,
    `nameUniqueness.py`), so every other system and sector anywhere in
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
    # discriminator (see systemData.StarSystem.__init__); is_binary is kept
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
        separation_km = star_system.star.binary_separation_au * physical_constants.AU_TO_KM
        # A close pair's mutual orbit is treated as circular (tidal
        # circularization is a legitimate simplification at its 0.05-0.25
        # AU separations) -- periapsis == apoapsis == separation, eccentricity 0.
        eccentricity = 0.0
        periapsis_km = apoapsis_km = separation_km
    elif binary_configuration == "wide":
        proxy_only_fields = _NULL_PROXY_ONLY_BINARY_FIELDS
        mutual_orbit_fields = _mutual_orbit_fields_from_wide_binary(star_system.wide_binary)
        separation_km = star_system.wide_binary.separation_au * physical_constants.AU_TO_KM
        eccentricity = star_system.wide_binary.eccentricity
        periapsis_km = star_system.wide_binary.periapsis_au * physical_constants.AU_TO_KM
        apoapsis_km = star_system.wide_binary.apoapsis_au * physical_constants.AU_TO_KM
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
            star_system.binary_planetary_wobble_x * physical_constants.AU_TO_KM,
            star_system.binary_planetary_wobble_y * physical_constants.AU_TO_KM,
            star_system.binary_planetary_wobble_z * physical_constants.AU_TO_KM,
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

    # Name-uniqueness (v24, nameUniqueness.py) -- the stars, planets and
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
            position_x_mpc, position_y_mpc, position_z_mpc, quadrant, location,
            is_binary, binary_configuration,
            binary_separation_km, binary_eccentricity, binary_periapsis_km, binary_apoapsis_km,
            binary_type, binary_temperature_k, binary_radius_km,
            binary_effective_mass_kg, binary_effective_luminosity_w, binary_age_gy, binary_lifespan_gy,
            binary_habitable_zone_inner_km, binary_habitable_zone_outer_km,
            binary_system_perimeter_km, binary_heliosphere_radius_km,
            binary_galactic_orbital_speed_kms, binary_galactic_orbital_period_gy,
            binary_galactic_orbital_phase_deg, binary_galactic_min_update_interval_years,
            binary_mutual_orbital_period_years, binary_mutual_orbital_speed_kms,
            binary_mutual_orbital_inclination_deg, binary_mutual_orbital_ascending_node_deg,
            binary_mutual_orbital_phase_deg, binary_mutual_min_update_interval_years,
            binary_mutual_position_x_km, binary_mutual_position_y_km, binary_mutual_position_z_km,
            binary_primary_position_x_km, binary_primary_position_y_km, binary_primary_position_z_km,
            binary_secondary_position_x_km, binary_secondary_position_y_km, binary_secondary_position_z_km,
            binary_secondary_mass_fraction,
            binary_planetary_wobble_x_km, binary_planetary_wobble_y_km, binary_planetary_wobble_z_km,
            system_flavor_text, runaway_class, runaway_speed_kms, schema_version
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, config_id, star_system.name,
            position_x_mpc, position_y_mpc, position_z_mpc, quadrant, location,
            int(is_binary), binary_configuration,
            separation_km, eccentricity, periapsis_km, apoapsis_km,
            *proxy_only_fields,
            *mutual_orbit_fields,
            *planetary_wobble_fields,
            star_system.system_flavor_text,
            getattr(star_system, "runaway_class", None), getattr(star_system, "runaway_speed_kms", None),
            SCHEMA_VERSION,
        ),
    )
    star_system_id = cur.lastrowid
    confirm_system_name(conn, name_base, star_system_id, diminutive_index)

    if proxy_like:
        insert_star(conn, star_system.primary_star, star_system_id, "primary")
        insert_star(conn, star_system.secondary_star, star_system_id, "secondary")
        primary_star_id = secondary_star_id = None  # planets orbit the proxy, not a stored star row -- see planets.star_id
    elif binary_configuration == "wide":
        primary_star_id = insert_star(conn, star_system.primary_star, star_system_id, "primary")
        secondary_star_id = insert_star(conn, star_system.secondary_star, star_system_id, "secondary")
    else:
        primary_star_id = insert_star(conn, star_system.star, star_system_id, "single")
        secondary_star_id = None
        # v16: a StarSystem anchored by a compact remnant (compactRemnant.py,
        # phenomenonGen.py --anchor-system) -- insert_star above already
        # wrote the stars row (BlackHole/NeutronStar's Star-compatible
        # attribute surface, with yerkes_class set to the literal 'BH'/'NS'
        # marker), so only the satellite black_holes/neutron_stars row
        # (this remnant's own extra fields) remains.
        if isinstance(star_system.star, BlackHole):
            insert_black_hole(conn, star_system.star, star_id=primary_star_id)
        elif isinstance(star_system.star, NeutronStar):
            insert_neutron_star(conn, star_system.star, star_id=primary_star_id)

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

    # Comets have their own list, own table, and no orbital_index -- see
    # insert_comet's own docstring for why they don't share the
    # planets/asteroid_belts orbital-slot loop above.
    for comet in star_system.comets:
        insert_comet(conn, comet, star_system_id, star_id=primary_star_id)
    for comet in star_system.secondary_comets:
        insert_comet(conn, comet, star_system_id, star_id=secondary_star_id)

    return star_system_id


# --- Containment (v39): what sits inside a nebula or supernova remnant ---

CONTAINABLE_TABLES = (
    "star_systems", "rogue_planets", "interstellar_comets", "black_holes",
    "neutron_stars", "asteroid_fields", "nebulae",
)
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
        dict or None: `{"type": "nebula" | "supernova_remnant", "id",
            "name", "class", "density_cm3", "temperature_k"}`, or `None` in
            open space.
    """
    if row["inside_nebula_id"] is not None:
        found = conn.execute("SELECT id, name, nebula_class AS class, density_cm3, temperature_k"
                             " FROM nebulae WHERE id = ?", (row["inside_nebula_id"],)).fetchone()
        kind = "nebula"
    elif row["inside_remnant_id"] is not None:
        found = conn.execute("SELECT id, name, remnant_class AS class, density_cm3, temperature_k"
                             " FROM supernova_remnants WHERE id = ?", (row["inside_remnant_id"],)).fetchone()
        kind = "supernova_remnant"
    else:
        return None
    if found is None:
        return None
    return {"type": kind, "id": found["id"], "name": found["name"], "class": found["class"],
            "density_cm3": found["density_cm3"], "temperature_k": found["temperature_k"]}


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
            containers.append({
                "column": column, "id": row["id"],
                "center": (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]),
                "radius_pc": ly_to_pc(row["radius_ly"]),
            })
    return containers


def innermost_container(point_pc, containers, own_radius_pc=0.0, own=None):
    """
    The smallest container in `containers` whose sphere holds `point_pc`
    -- for a nebula (`own_radius_pc` > 0), only a container larger than it,
    and never itself (`own`, a `(column, id)` pair).

    Returns:
        dict or None: The container, or `None` when the point is in open space.
    """
    best = None
    for container in containers:
        if own is not None and (container["column"], container["id"]) == own:
            continue
        if container["radius_pc"] <= own_radius_pc:
            continue
        if math.dist(point_pc, container["center"]) > container["radius_pc"]:
            continue
        if best is None or container["radius_pc"] < best["radius_pc"]:
            best = container
    return best


def refresh_containment(conn, sector_ids):
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
    """
    sector_ids = sorted(set(sector_ids))
    for start in range(0, len(sector_ids), 500):
        _refresh_containment_batch(conn, sector_ids[start:start + 500])


def _refresh_containment_batch(conn, sector_ids):
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
        best = innermost_container(point, containers, own_radius_pc, own) if point is not None else None
        new = (best["id"] if best and best["column"] == "inside_nebula_id" else None,
               best["id"] if best and best["column"] == "inside_remnant_id" else None)
        if new != (row["inside_nebula_id"], row["inside_remnant_id"]):
            changes.setdefault(table, []).append((row["id"], *new))

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


PLACED_PHENOMENON_TABLES = (
    "black_holes", "neutron_stars", "nebulae", "supernova_remnants", "rogue_planets",
    "interstellar_comets", "asteroid_fields", "quasars",
)
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
                if system_id != exclude and (distance := math.dist(point, other)) <= limit_pc
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
                            distance = math.dist(point, other)
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
        if any(math.dist(center, p) <= reach_pc for p in points):
            near[row["id"]] = center
    return near


def _edge_pc(conn):
    row = conn.execute("SELECT MAX(edge_mpc) AS edge FROM sectors").fetchone()
    return (row["edge"] or 0) / 1000.0


def refresh_nearest_systems(conn, sector_ids):
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

    Returns:
        set: The `(table, id)` of every object whose list changed.
    """
    sector_ids = sorted(set(sector_ids))
    half_diagonal = _edge_pc(conn) * math.sqrt(3) / 2
    all_changed = set()
    for start in range(0, len(sector_ids), 200):
        centers = _sector_centers(conn, sector_ids[start:start + 200])
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
    own = _sector_centers(conn, [sector_id])
    if not own:
        return
    half_diagonal = _edge_pc(conn) * math.sqrt(3) / 2
    near = _sectors_near(conn, own, 2 * half_diagonal + NEAREST_SYSTEMS_SEARCH_PC)
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
                 offset_ly=None, description=None):
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
            axes (`galaxyGeometry.sector_orientation`). The center if left
            out.
        description (str, optional): Free text.

    An asteroid facility in a belt takes no distance: it gets a random
    spot in the belt (`facilities.belt_position`) and the circular orbit
    around its star from there, which `advance_facility_orbits` moves
    like an orbital facility's.

    Returns:
        int: The new `facilities.id`.

    Raises:
        FacilityError: If the rules refuse it or the host doesn't exist.
    """
    columns, mass, radius, body_type, sphere = _facility_host(conn, host_type, host_id)
    problem = facility_rules.check_facility(kind, placement, host_type, body_type)
    if problem:
        raise FacilityError(problem)

    orbit = {}
    if placement == "orbital":
        try:
            orbit = facility_rules.orbit_for(mass, radius, distance_km, facility_rules.orbit_limits(radius, sphere)[1])
        except ValueError as exc:
            raise FacilityError(str(exc)) from exc
        if phase_deg is None:
            phase_deg = random.uniform(0.0, 360.0)
        orbit["phase_deg"] = phase_deg % 360.0
    elif distance_km is not None or phase_deg is not None:
        raise FacilityError("only an orbital facility takes an orbit distance or phase")
    elif host_type == "asteroid_belt":
        orbit = _belt_orbit(conn, host_id, mass, radius)

    placement_values = (None, None, None, None)
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
    elif offset_ly is not None:
        raise FacilityError("only a stand-alone facility has a position in space")

    cur = conn.execute(
        """
        INSERT INTO facilities (
            name, kind, placement, host_type, star_system_id, star_id, planet_id, moon_id,
            asteroid_belt_id, asteroid_field_id, sector_id,
            center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc,
            orbit_distance_km, orbit_period_years, orbital_speed_kms, orbit_phase_deg, description
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            name, kind, placement, host_type, columns.get("star_system_id"), columns.get("star_id"),
            columns.get("planet_id"), columns.get("moon_id"), columns.get("asteroid_belt_id"),
            columns.get("asteroid_field_id"), columns.get("sector_id"),
            *placement_values,
            orbit.get("distance_km"), orbit.get("period_years"), orbit.get("orbital_speed_kms"),
            orbit.get("phase_deg"), description,
        ),
    )
    if columns.get("star_system_id") is not None:
        touch_star_system(conn, columns["star_system_id"])
    return cur.lastrowid


def delete_facility(conn, facility_id):
    """Deletes one facility. Returns `False` if there was none."""
    row = conn.execute("SELECT star_system_id FROM facilities WHERE id = ?", (facility_id,)).fetchone()
    if row is None:
        return False
    conn.execute("DELETE FROM facilities WHERE id = ?", (facility_id,))
    if row["star_system_id"] is not None:
        touch_star_system(conn, row["star_system_id"])
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
    """The angle (radians) an orbit of `period_gy` sweeps in `elapsed_years`."""
    if not period_gy or period_gy <= 0 or elapsed_years <= 0:
        return 0.0
    return 2 * math.pi * elapsed_years / (period_gy * 1e9)


class _SectorIndex:
    """Every placed grid sector by address, for refiling moved objects."""

    def __init__(self, conn):
        self.edge_pc = _edge_pc(conn)
        self.by_address = {}
        self.info = {}
        for row in conn.execute(
            "SELECT id, name, ring_index, layer_index, ring_slot_index, center_x_pc, center_y_pc, center_z_pc"
            " FROM sectors WHERE center_x_pc IS NOT NULL"
        ).fetchall():
            center = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
            self.info[row["id"]] = (row["name"], center)
            if row["ring_index"] is not None:
                self.by_address[(row["ring_index"], row["layer_index"], row["ring_slot_index"])] = row["id"]

    def sector_at(self, point, current):
        """The generated sector holding `point`, or `current` when the
        cell it drifted into hasn't been generated (it stays filed where
        it was until that sector exists)."""
        if self.edge_pc <= 0:
            return current
        return self.by_address.get(sector_address_at(point, self.edge_pc), current)


def advance_galactic_positions(conn, elapsed_years):
    """
    Moves every placed star system, standalone phenomenon and stand-alone
    facility along its galactic orbit by `elapsed_years`, the same turn
    its galactic phase advances by (`advance_orbital_phases`): a system by
    its close pair's or primary star's period, a phenomenon by its own,
    a facility by the rotation curve at its radius
    (`utils.calculate_galactic_orbit`). Sectors are fixed cells, so an
    object that drifts into another generated sector is refiled there
    (`sector_id`, its sector-relative position, and later its octant,
    location text, containment and nearest systems -- see
    `refresh_after_motion`). A pure move keeps `modified_at`; a change of
    sector bumps it.

    Returns:
        dict: `moved` and `refiled` counts, and `sectors` -- the ids of
            every sector something moved in or out of.
    """
    index = _SectorIndex(conn)
    moved = refiled = 0
    touched = set()

    systems = conn.execute(
        """
        SELECT ss.id, ss.sector_id, ss.position_x_mpc, ss.position_y_mpc, ss.position_z_mpc,
               COALESCE(ss.binary_galactic_orbital_period_gy, s.galactic_orbital_period_gy) AS period_gy
        FROM star_systems ss
        LEFT JOIN stars s ON s.star_system_id = ss.id AND s.role IN ('single', 'primary')
        WHERE ss.position_x_mpc IS NOT NULL AND ss.sector_id IS NOT NULL
        """
    ).fetchall()
    updates, refile_updates = [], []
    for row in systems:
        if row["sector_id"] not in index.info:
            continue
        angle = _galactic_turn(elapsed_years, row["period_gy"])
        if angle == 0.0:
            continue
        _name, center = index.info[row["sector_id"]]
        offset = (row["position_x_mpc"] / 1000.0, row["position_y_mpc"] / 1000.0, row["position_z_mpc"] / 1000.0)
        point = _rotate_about_axis(local_to_galaxy_pc(center, offset), angle)
        sector_id = index.sector_at(point, row["sector_id"])
        local_mpc = tuple(c * 1000.0 for c in galaxy_to_local_pc(index.info[sector_id][1], point))
        quadrant, _magnitudes = classify_octant(local_mpc)
        moved += 1
        if sector_id != row["sector_id"]:
            refiled += 1
            touched.update((sector_id, row["sector_id"]))
            refile_updates.append((sector_id, *local_mpc, quadrant, row["id"]))
        else:
            updates.append((*local_mpc, quadrant, row["id"]))
    if updates:
        conn.executemany(
            "UPDATE star_systems SET position_x_mpc = ?, position_y_mpc = ?, position_z_mpc = ?, quadrant = ?,"
            " modified_at = modified_at WHERE id = ?", updates)
    if refile_updates:
        conn.executemany(
            "UPDATE star_systems SET sector_id = ?, position_x_mpc = ?, position_y_mpc = ?, position_z_mpc = ?,"
            " quadrant = ? WHERE id = ?", refile_updates)
        _refile_nearest_rows(conn, "star_systems", [(update[0], update[-1]) for update in refile_updates])

    for table in PLACED_PHENOMENON_TABLES:
        if table == "quasars":
            continue  # the galaxy's nucleus sits at the center and doesn't orbit it
        updates, refile_updates = [], []
        for row in conn.execute(
            f"SELECT id, sector_id, center_x_pc, center_y_pc, center_z_pc, galactic_orbital_period_gy"
            f" FROM {table} WHERE center_x_pc IS NOT NULL"
        ).fetchall():
            angle = _galactic_turn(elapsed_years, row["galactic_orbital_period_gy"])
            if angle == 0.0:
                continue
            point = _rotate_about_axis((row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]), angle)
            radius = math.sqrt(sum(c * c for c in point))
            moved += 1
            sector_id = index.sector_at(point, row["sector_id"]) if row["sector_id"] is not None else None
            if sector_id != row["sector_id"]:
                refiled += 1
                touched.update((sector_id, row["sector_id"]))
                refile_updates.append((sector_id, *point, radius, row["id"]))
            else:
                updates.append((*point, radius, row["id"]))
        if updates:
            conn.executemany(
                f"UPDATE {table} SET center_x_pc = ?, center_y_pc = ?, center_z_pc = ?, galactic_radius_pc = ?,"
                f" modified_at = modified_at WHERE id = ?", updates)
        if refile_updates:
            conn.executemany(
                f"UPDATE {table} SET sector_id = ?, center_x_pc = ?, center_y_pc = ?, center_z_pc = ?,"
                f" galactic_radius_pc = ? WHERE id = ?", refile_updates)
            _refile_nearest_rows(conn, table, [(update[0], update[-1]) for update in refile_updates])

    updates, refile_updates = [], []
    for row in conn.execute(
        "SELECT id, sector_id, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc FROM facilities"
        " WHERE host_type = 'space' AND center_x_pc IS NOT NULL"
    ).fetchall():
        _speed, period_gy = calculate_galactic_orbit(pc_to_ly(math.hypot(row["center_x_pc"], row["center_y_pc"])))
        angle = _galactic_turn(elapsed_years, period_gy)
        if angle == 0.0:
            continue
        point = _rotate_about_axis((row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]), angle)
        moved += 1
        sector_id = index.sector_at(point, row["sector_id"])
        if sector_id != row["sector_id"]:
            refiled += 1
            touched.update((sector_id, row["sector_id"]))
            refile_updates.append((sector_id, *point, row["galactic_radius_pc"], row["id"]))
        else:
            updates.append((*point, row["id"]))
    if updates:
        conn.executemany("UPDATE facilities SET center_x_pc = ?, center_y_pc = ?, center_z_pc = ?,"
                         " modified_at = modified_at WHERE id = ?", updates)
    if refile_updates:
        conn.executemany("UPDATE facilities SET sector_id = ?, center_x_pc = ?, center_y_pc = ?, center_z_pc = ?,"
                         " galactic_radius_pc = ? WHERE id = ?", refile_updates)

    touched.discard(None)
    for sector_id in touched:
        touch_sector(conn, sector_id)
    return {"moved": moved, "refiled": refiled, "sectors": touched}


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


def advance_facility_orbits(conn, elapsed_years):
    """Advances every orbital facility's `orbit_phase_deg` by its period,
    the way `advance_orbital_phases` does for moons, and every asteroid
    facility in a belt around its star the same way. Returns the count."""
    if elapsed_years <= 0:
        return 0
    return conn.execute(
        "UPDATE facilities SET orbit_phase_deg = MOD(orbit_phase_deg + 360.0 * ? / orbit_period_years, 360.0),"
        " modified_at = modified_at WHERE placement IN ('orbital', 'asteroid') AND orbit_period_years > 0",
        (elapsed_years,),
    ).rowcount


def refresh_after_motion(conn, refiled_sectors=()):
    """
    Brings everything that depends on where things are up to date after
    `advance_galactic_positions`: containment (`refresh_containment`),
    each phenomenon's octant and every object's nearest systems
    (`refresh_nearest_systems`) for every placed sector, then the
    `star_systems.location` text of every system whose nearest systems
    changed or that sits in a sector something moved in or out of, built
    from its sector's name and the stored nearest systems.

    Returns:
        int: How many location texts were rewritten.
    """
    placed = [row["id"] for row in conn.execute("SELECT id FROM sectors WHERE center_x_pc IS NOT NULL").fetchall()]
    refresh_containment(conn, placed)
    changed = refresh_nearest_systems(conn, placed)
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
    return rewritten


# ---------------------------------------------------------------------------
# Bright-star pre-placement (schema v43) -- storage. The scatter itself
# (generate.py plan) and the fill draw stars from Physics' sampling API.
# ---------------------------------------------------------------------------

BRIGHT_STAR_COLUMNS = (
    "ring_index", "layer_index", "ring_slot_index", "position_x_mpc", "position_y_mpc", "position_z_mpc",
    "population", "star_type", "yerkes_class", "mass_kg", "radius_km", "temperature_k", "luminosity_w",
    "age_gy", "lifespan_gy", "initial_mass_sol", "phase_end_age_gy", "seed",
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
    conn.execute("UPDATE sector_stats SET bright_level_sol = -1 WHERE bright_level_sol > 0")
    conn.execute("UPDATE sector_stats SET level_before_fill_sol = -1 WHERE level_before_fill_sol > 0")
    conn.execute("UPDATE galaxy_shape SET bright_star_min_luminosity_sol = NULL, bright_star_seed = NULL")
    conn.commit()


def record_bright_star_scatter(conn, min_luminosity_sol, seed):
    """Stores the threshold and seed a finished scatter used, so a fill
    reads them rather than today's constant."""
    conn.execute("UPDATE galaxy_shape SET bright_star_min_luminosity_sol = ?, bright_star_seed = ? WHERE id = 1",
                 (min_luminosity_sol, seed))


def bright_star_scatter_settings(conn):
    """`(min_luminosity_sol, seed)` of the galaxy's scatter, or `None`
    when none has run (a fill then behaves exactly as before)."""
    row = conn.execute(
        "SELECT bright_star_min_luminosity_sol, bright_star_seed FROM galaxy_shape WHERE id = 1").fetchone()
    if row is None or row["bright_star_min_luminosity_sol"] is None:
        return None
    return row["bright_star_min_luminosity_sol"], row["bright_star_seed"]


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
        addresses (iterable): `(ring, layer, slot)` tuples.

    Returns:
        dict: `(ring, layer, slot)` -> level, for the sectors with a row.
    """
    levels = {}
    for chunk in _address_chunks(addresses):
        marks = ", ".join("(?, ?, ?)" for _ in chunk)
        for row in conn.execute(
            "SELECT ring_index, layer_index, ring_slot_index, bright_level_sol FROM sector_stats"
            f" WHERE (ring_index, layer_index, ring_slot_index) IN ({marks})",
            tuple(value for key in chunk for value in key),
        ).fetchall():
            levels[(row["ring_index"], row["layer_index"], row["ring_slot_index"])] = row["bright_level_sol"]
    return levels


def sector_bright_level_keys(conn):
    """
    Every sector a backfill took to its own level (`bright_level_sol` above
    0, GEN.44), as `(ring, layer, slot) -> level`: what a staged scatter
    (`generate.py plan --bright-stars-down-to`) leaves out of its layers
    and tops up sector by sector instead.
    """
    rows = conn.execute("SELECT ring_index, layer_index, ring_slot_index, bright_level_sol FROM sector_stats"
                        " WHERE bright_level_sol > 0").fetchall()
    return {(row["ring_index"], row["layer_index"], row["ring_slot_index"]): row["bright_level_sol"] for row in rows}


def lock_sector_stats(conn, entries):
    """
    Takes the row locks on some sectors' `sector_stats` rows (GEN.44),
    making each row first (untouched, -1) with its expected density
    (PERF.11) if there is none, so two backfills never draw the same
    sector at once. Rows are locked in address order, so two callers
    wait rather than deadlock. Holds until the caller commits.

    Args:
        conn (Connection): An open connection.
        entries (iterable): `((ring, layer, slot), relative_density,
            expected_systems)` per sector.

    Returns:
        dict: `(ring, layer, slot)` -> the level now stored.
    """
    rows = sorted((int(address[0]), int(address[1]), int(address[2]), density, expected)
                  for address, density, expected in entries)
    for start in range(0, len(rows), 500):
        chunk = rows[start:start + 500]
        # ON DUPLICATE KEY UPDATE takes each row's exclusive lock at once;
        # INSERT IGNORE would take a shared one, and two workers both
        # upgrading it to exclusive for the SELECT below would deadlock.
        # A derived table rather than VALUES(col), deprecated on MySQL
        # (tests/test_sql_portability.py).
        selects = " UNION ALL ".join(
            ["SELECT ? AS ring_index, ? AS layer_index, ? AS ring_slot_index, ? AS density, ? AS expected"]
            + ["SELECT ?, ?, ?, ?, ?"] * (len(chunk) - 1)
        )
        conn.execute(
            "INSERT INTO sector_stats (ring_index, layer_index, ring_slot_index, relative_density, expected_systems)"
            f" SELECT * FROM ({selects}) AS incoming"
            " ON DUPLICATE KEY UPDATE relative_density = COALESCE(sector_stats.relative_density, incoming.density),"
            " expected_systems = COALESCE(sector_stats.expected_systems, incoming.expected)",
            tuple(value for row in chunk for value in row),
        )
    levels = {}
    for chunk in _address_chunks([row[:3] for row in rows]):
        marks = ", ".join("(?, ?, ?)" for _ in chunk)
        for row in conn.execute(
            "SELECT ring_index, layer_index, ring_slot_index, bright_level_sol FROM sector_stats"
            f" WHERE (ring_index, layer_index, ring_slot_index) IN ({marks}) FOR UPDATE",
            tuple(value for key in chunk for value in key),
        ).fetchall():
            levels[(row["ring_index"], row["layer_index"], row["ring_slot_index"])] = row["bright_level_sol"]
    return levels


def set_sector_bright_levels(conn, levels):
    """Records how dim some sectors' stars now go (`(ring, layer, slot) ->
    L_sun`), inside the transaction `lock_sector_stats` started."""
    conn.executemany(
        "UPDATE sector_stats SET bright_level_sol = ? WHERE ring_index = ? AND layer_index = ? AND ring_slot_index = ?",
        [(level, *address) for address, level in sorted(levels.items())],
    )


def bright_star_fill_level(conn, ring_index, layer_index, ring_slot_index):
    """
    The luminosity one sector's pre-placed stars go down to: its own
    backfill level (GEN.44, `sector_stats`) when it has one, else the
    galaxy's scatter threshold, else `None` (no bright stars placed; a
    fill caps nothing).
    """
    if ring_index is not None and layer_index is not None and ring_slot_index is not None:
        level = sector_bright_levels(conn, [(ring_index, layer_index, ring_slot_index)]).get(
            (ring_index, layer_index, ring_slot_index))
        if level is not None and level > 0:
            return level
    settings = bright_star_scatter_settings(conn)
    return settings[0] if settings else None


DENSITY_RATIO_WINDOW = 1000
"""int: The decaying average of actual against expected systems (PERF.11)
is a plain mean over the first fills, then weighs each new fill 1 in
this many, so it follows a galaxy whose model drifts."""


def record_sector_stats(conn, sector_id, address, center_pc):
    """
    Writes a sector's stats once it is generated (PERF.11, GEN.44): its
    expected density from the galaxy model, the systems and stars it got,
    the mean temperature and luminosity of those stars, and its Galaxy
    Map color worked out from them (MAP.86, `sectorLook`), and
    level 0 (filled), keeping the level it had in
    `level_before_fill_sol`; then folds its actual-to-expected systems
    into `galaxy_shape`'s decaying average. Inside the sector's own save,
    last, so the galaxy row is locked only for the commit.

    Args:
        conn (Connection): An open connection, mid-transaction.
        sector_id (int): The new `sectors.id`.
        address (tuple): Its `(ring, layer, slot)`.
        center_pc (tuple): Its galaxy-frame center, parsecs.
    """
    from planetgen.galaxy.density import relative_density
    from planetgen.galaxy.sector_look import fill_share, max_sector_systems, sector_color

    skeleton = get_galaxy_shape(conn)
    density = expected = None
    if skeleton is not None:
        density = max(relative_density(center_pc, skeleton.shape), 0.0)
        expected = density * skeleton.expected_system_count_at_density_1
    systems = conn.execute("SELECT COUNT(*) AS n FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["n"]
    stars = conn.execute(
        "SELECT COUNT(*) AS n, AVG(st.temperature_k) AS temperature, AVG(st.luminosity_w) AS luminosity"
        " FROM stars st JOIN star_systems ss ON ss.id = st.star_system_id WHERE ss.sector_id = ?", (sector_id,),
    ).fetchone()
    luminosity = None if stars["luminosity"] is None else float(stars["luminosity"]) / physical_constants.SOLAR_LUMINOSITY
    temperature = None if stars["temperature"] is None else float(stars["temperature"])
    share = fill_share(systems, max_sector_systems(skeleton))
    color = sector_color(temperature, luminosity, share) or (None, None, None)
    conn.execute(
        "INSERT INTO sector_stats (ring_index, layer_index, ring_slot_index, bright_level_sol, level_before_fill_sol,"
        " relative_density, expected_systems, actual_systems, actual_stars, mean_temperature_k, mean_luminosity_sol,"
        " fill_share, color_r, color_g, color_b, filled_at)"
        " VALUES (?, ?, ?, 0, -1, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP(3))"
        " ON DUPLICATE KEY UPDATE level_before_fill_sol = IF(bright_level_sol = 0, level_before_fill_sol,"
        " bright_level_sol), bright_level_sol = 0, relative_density = ?, expected_systems = ?, actual_systems = ?,"
        " actual_stars = ?, mean_temperature_k = ?, mean_luminosity_sol = ?, fill_share = ?, color_r = ?,"
        " color_g = ?, color_b = ?, filled_at = CURRENT_TIMESTAMP(3)",
        (*address, *(2 * (density, expected, systems, stars["n"], temperature, luminosity, share, *color))),
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
        " fill_share = NULL, color_r = NULL, color_g = NULL, color_b = NULL,"
        " filled_at = NULL WHERE ring_index = ? AND layer_index = ? AND ring_slot_index = ? AND bright_level_sol = 0",
        tuple(address),
    )


def get_sector_stats(conn, ring_index, layer_index, ring_slot_index):
    """One sector's `sector_stats` row as a dict (PERF.11), or `None`."""
    return conn.execute(
        "SELECT * FROM sector_stats WHERE ring_index = ? AND layer_index = ? AND ring_slot_index = ?",
        (ring_index, layer_index, ring_slot_index),
    ).fetchone()


def galaxy_density_ratio(conn):
    """`(average, samples)`: the decaying average of the systems sector
    fills got against the systems expected (PERF.11); the average is
    `None` before any fill."""
    row = conn.execute("SELECT density_ratio_avg, density_ratio_samples FROM galaxy_shape WHERE id = 1").fetchone()
    if row is None:
        return None, 0
    return row["density_ratio_avg"], int(row["density_ratio_samples"])


def delete_unfinished_band(conn, below_luminosity_w, keep_addresses=(), batch_size=500):
    """
    Deletes the unbuilt bright stars below `below_luminosity_w` (the
    galaxy's star-fill level) outside `keep_addresses`: what a band run
    (`generate.py plan --bright-stars-down-to`) that stopped part way left
    in the layers it got through, since the level only moves once a band
    is whole (GEN.32). Every finished scatter or band is at or above the
    level, and the stars below it that belong there (a backfilled block's,
    a filled sector's) are in `keep_addresses` or built.

    Returns:
        int: Stars deleted.
    """
    keep = set(keep_addresses)
    rows = conn.execute(
        "SELECT ring_index, layer_index, ring_slot_index FROM bright_stars"
        " WHERE luminosity_w < ? AND star_system_id IS NULL"
        " GROUP BY ring_index, layer_index, ring_slot_index", (below_luminosity_w,)).fetchall()
    cells = [(row["ring_index"], row["layer_index"], row["ring_slot_index"]) for row in rows]
    cells = [cell for cell in cells if cell not in keep]
    deleted = 0
    for start in range(0, len(cells), batch_size):
        chunk = cells[start:start + batch_size]
        placeholders = ", ".join("(?, ?, ?)" for _ in chunk)
        cur = conn.execute(
            "DELETE FROM bright_stars WHERE luminosity_w < ? AND star_system_id IS NULL"
            f" AND (ring_index, layer_index, ring_slot_index) IN ({placeholders})",
            (below_luminosity_w, *[value for cell in chunk for value in cell]))
        deleted += cur.rowcount
    return deleted


def insert_bright_stars(conn, rows, batch_size=10000):
    """
    Bulk-writes scattered bright stars (`BRIGHT_STAR_COLUMNS` order) in
    batches of `batch_size`, which need no server setting (`LOAD DATA
    LOCAL INFILE` would need `local_infile` on both ends).

    Returns:
        int: Rows written.
    """
    columns = ", ".join(BRIGHT_STAR_COLUMNS)
    marks = ", ".join("?" * len(BRIGHT_STAR_COLUMNS))
    written = 0
    batch = []
    for row in rows:
        batch.append(tuple(row))
        if len(batch) >= batch_size:
            conn.executemany(f"INSERT INTO bright_stars ({columns}) VALUES ({marks})", batch)
            written += len(batch)
            batch = []
    if batch:
        conn.executemany(f"INSERT INTO bright_stars ({columns}) VALUES ({marks})", batch)
        written += len(batch)
    return written


def bright_stars_for_sector(conn, ring_index, layer_index, ring_slot_index):
    """The pre-placed stars in one sector cell not yet built into a
    system, brightest first, as dicts of every `bright_stars` column."""
    return [dict(row) for row in conn.execute(
        "SELECT * FROM bright_stars WHERE ring_index = ? AND layer_index = ? AND ring_slot_index = ?"
        " AND star_system_id IS NULL ORDER BY luminosity_w DESC, id",
        (ring_index, layer_index, ring_slot_index),
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
            " AND ring_index IS NOT NULL ORDER BY ring_index, layer_index, ring_slot_index", (since,)
        ).fetchall()
    ]


def filled_sector_addresses(conn):
    """The `(ring, layer, slot)` of every galaxy sector already filled."""
    return {
        (row["ring_index"], row["layer_index"], row["ring_slot_index"])
        for row in conn.execute(
            "SELECT ring_index, layer_index, ring_slot_index FROM sectors WHERE ring_index IS NOT NULL"
        ).fetchall()
    }


def mark_bright_star_filled(conn, bright_star_id, star_system_id):
    """Links a pre-placed star to the system its sector's fill built
    around it."""
    conn.execute("UPDATE bright_stars SET star_system_id = ? WHERE id = ?", (star_system_id, bright_star_id))


def sectors_reached_by(conn, center_pc, radius_pc):
    """The ids of every placed sector a sphere of `radius_pc` around
    `center_pc` (galaxy-frame parsecs) overlaps."""
    edge_row = conn.execute("SELECT MAX(edge_mpc) AS edge FROM sectors").fetchone()
    if edge_row["edge"] is None:
        return []
    reach = radius_pc + _sector_half_diagonal_pc(edge_row["edge"])
    rows = conn.execute(
        "SELECT id, center_x_pc, center_y_pc, center_z_pc FROM sectors"
        " WHERE center_x_pc BETWEEN ? AND ? AND center_y_pc BETWEEN ? AND ? AND center_z_pc BETWEEN ? AND ?",
        (center_pc[0] - reach, center_pc[0] + reach, center_pc[1] - reach, center_pc[1] + reach,
         center_pc[2] - reach, center_pc[2] + reach),
    ).fetchall()
    return [
        row["id"] for row in rows
        if math.dist(center_pc, (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])) <= reach
    ]


def _refresh_containment_around(conn, placement, radius_ly):
    """After placing a nebula or remnant: refresh every sector it reaches."""
    if placement is None or placement.get("center_x_pc") is None:
        return
    center = (placement["center_x_pc"], placement["center_y_pc"], placement["center_z_pc"])
    refresh_containment(conn, sectors_reached_by(conn, center, ly_to_pc(radius_ly)))


def insert_sector(conn, sector: SpaceSector, galaxy_position=None) -> int:
    """
    Inserts a full `SpaceSector` -- the `sectors` row, every system it
    contains (with its placement), and every exotic phenomenon it contains
    (`sector.phenomena`, see `sectorGen.generate_sector_phenomena`) -- into
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
            `sectorGen.py`'s own standalone CLI keeps producing these.
            When given, must have keys `center_x_pc`, `center_y_pc`,
            `center_z_pc`, `galactic_radius_pc` (all required together --
            the schema's CHECK constraint enforces that), and optionally
            `ring_index`/`layer_index`/`ring_slot_index`, the sector's
            cell in the cylindrical grid (`None`/omitted leaves them NULL
            -- a hand-placed position with no grid address).

    Rows are written in batches (PERF.13, `Connection.batched`): one
    multi-row INSERT per table and statement shape, with ids from
    `id_blocks`, instead of one INSERT per row.

    A grid-addressed sector's stats (`sector_stats`, PERF.11, GEN.44) are
    written last, from the rows just stored (`record_sector_stats`).

    Returns:
        int: The new `sectors.id`.
    """
    with conn.batched():
        sector_id = _insert_sector_rows(conn, sector, galaxy_position)
    address = None if galaxy_position is None else tuple(
        galaxy_position.get(key) for key in ("ring_index", "layer_index", "ring_slot_index"))
    if address is not None and None not in address:
        record_sector_stats(conn, sector_id, address, (galaxy_position["center_x_pc"], galaxy_position["center_y_pc"],
                                                       galaxy_position["center_z_pc"]))
    return sector_id


def _insert_sector_rows(conn, sector, galaxy_position):
    if galaxy_position is not None:
        cur = conn.execute(
            """
            INSERT INTO sectors (
                name, edge_mpc, center_x_pc, center_y_pc, center_z_pc,
                galactic_radius_pc, ring_index, layer_index, ring_slot_index
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                sector.name, ly_to_milliparsecs(sector.edge_ly),
                galaxy_position["center_x_pc"], galaxy_position["center_y_pc"],
                galaxy_position["center_z_pc"], galaxy_position["galactic_radius_pc"],
                galaxy_position.get("ring_index"), galaxy_position.get("layer_index"),
                galaxy_position.get("ring_slot_index"),
            ),
        )
    else:
        cur = conn.execute(
            "INSERT INTO sectors (name, edge_mpc) VALUES (?, ?)",
            (sector.name, ly_to_milliparsecs(sector.edge_ly)),
        )
    sector_id = cur.lastrowid

    # Name-uniqueness (v24, nameUniqueness.py) -- reserved right after the
    # INSERT (the registry row points at it), mutating sector.name in
    # place, so every print/rendering call site that reads it afterward
    # (including this same sector's own systems, added below) sees the
    # final, collision-free name. See `reserve_sector_name`.
    final_name, _name_base = reserve_sector_name(conn, sector.name, sector_id)
    if final_name != sector.name:
        conn.execute("UPDATE sectors SET name = ? WHERE id = ?", (final_name, sector_id))
    sector.name = final_name

    # Every system and uniquely named phenomenon in the sector reserves
    # its name in one go (PERF.14); the insert functions below take their
    # reservation from `prereserved_names` and leave their registry writes
    # in `deferred_name_confirmations` for one upsert at the end.
    conn.prereserved_names = {}
    conn.deferred_name_confirmations = []
    try:
        # Every placed phenomenon and every bright-sweep system is named
        # by its object ID (GEN.64), claimed here for the whole sector in
        # generation order; only the rest go through the name registry.
        by_id = _sector_object_ids(sector, galaxy_position)
        _claim_object_ids(conn, by_id)
        claimed = {id(obj) for obj, _kind, _center in by_id}
        named = [entry.star_system for entry in sector.entries if id(entry.star_system) not in claimed]
        named += [entry.phenomenon for entry in sector.phenomena
                  if _phenomenon_registers_name(entry.phenomenon) and id(entry.phenomenon) not in claimed]
        reservations = reserve_system_names(conn, [obj.name for obj in named])
        for obj, (final, base, diminutive_index) in zip(named, reservations):
            core = getattr(obj, "compact_remnant", None) if isinstance(obj, SupernovaRemnant) else None
            if core is not None and core.name == f"{obj.name} Core":
                core.name = f"{final} Core"
            obj.name = final
            conn.prereserved_names[id(obj)] = (base, diminutive_index)

        bright_star_links = []
        for entry in sector.entries:
            star_system_id = insert_star_system(
                conn, entry.star_system, entry.system_config,
                sector_id=sector_id, position=entry.position,
                location=_location_for_entry(sector, entry),
            )
            bright_star_id = getattr(entry, "bright_star_id", None)
            if bright_star_id is not None:
                bright_star_links.append((bright_star_id, star_system_id))

        for entry in sector.phenomena:
            inserter = _PHENOMENON_INSERTERS.get(entry.phenomenon_type)
            if inserter is None:
                raise ValueError(f"Unknown phenomenon type: {entry.phenomenon_type!r}")
            placement = _galaxy_placement_from_sector_offset(galaxy_position, entry.position)
            inserter(conn, entry.phenomenon, sector_id=sector_id, placement=placement)

        confirm_system_names(conn, conn.deferred_name_confirmations)
    finally:
        conn.prereserved_names = None
        conn.deferred_name_confirmations = None
    if bright_star_links:
        _update_by_id(conn, "bright_stars", ("star_system_id",), bright_star_links, touch=True)

    if galaxy_position is not None:
        # Containment and nearest-neighbor lists read the sectors around
        # this one and rewrite theirs, so two neighbors saved at once
        # (parallel generation, PERF.7) would each miss the other. One
        # writer at a time does this last step, holding the lock until
        # commit; the slower part above still runs side by side.
        conn.lock_until_commit(_neighbor_lock_name(conn))
        _insert_field_nebulae(conn, sector, sector_id)
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
        insert_nebula(conn, nebula, sector_id=sector_id, placement=placement)


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
    for entry in sector.entries:
        if getattr(entry, "bright_star_id", None) is not None:
            center = _placement_center(_galaxy_placement_from_sector_offset(galaxy_position, entry.position))
            items.append((entry.star_system, "bright-star", center))
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
    row = conn.execute("SELECT name, edge_mpc, ring_index FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    if row is None:
        raise ValueError(f"no sectors row with id {sector_id}")
    edge_ly = milliparsecs_to_ly(row["edge_mpc"])
    cell = SectorCell.for_ring(row["ring_index"], edge_ly) if row["ring_index"] is not None else None
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
            star=SimpleNamespace(system_perimeter=perimeter_km / physical_constants.AU_TO_KM))
        position = tuple(milliparsecs_to_ly(system[f"position_{axis}_mpc"]) for axis in "xyz")
        sector.entries.append(SectorSystemEntry(stand_in, position))
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
    "binary_galactic_orbital_speed_kms", "binary_galactic_orbital_period_gy",
    "binary_galactic_orbital_phase_deg", "binary_galactic_min_update_interval_years",
    "binary_mutual_orbital_period_years", "binary_mutual_orbital_speed_kms",
    "binary_mutual_orbital_inclination_deg", "binary_mutual_orbital_ascending_node_deg",
    "binary_mutual_orbital_phase_deg", "binary_mutual_min_update_interval_years",
    "binary_mutual_position_x_km", "binary_mutual_position_y_km", "binary_mutual_position_z_km",
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


def replace_star_system_content(conn, star_system_id, star_system, system_config):
    """
    Swaps a stored system's generated content -- its stars, planets,
    moons, asteroid belts, comets and the system-level star/binary fields
    -- for a newly generated system's, keeping the system's id, name,
    sector, position, location, containment and wiki links. The new
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
    # Saved standalone first (the registry gives it a temporary name),
    # then its bodies are moved under the kept row and the stand-in row
    # dropped, which also drops its name reservation.
    temp_id = insert_star_system(conn, star_system, system_config)
    temp_name = conn.execute("SELECT name FROM star_systems WHERE id = ?", (temp_id,)).fetchone()["name"]
    conn.execute("DELETE FROM facilities WHERE star_system_id = ?", (star_system_id,))
    for table in ("moons", "planets", "asteroid_belts", "comets", "stars"):
        conn.execute(f"DELETE FROM {table} WHERE star_system_id = ?", (star_system_id,))
    for table in _SYSTEM_CONTENT_TABLES:
        conn.execute(f"UPDATE {table} SET star_system_id = ? WHERE star_system_id = ?", (star_system_id, temp_id))
    assignments = ", ".join(f"kept.{column} = made.{column}" for column in _SYSTEM_CONTENT_COLUMNS)
    conn.execute(
        f"UPDATE star_systems kept JOIN star_systems made ON made.id = ? SET {assignments},"
        " kept.modified_at = CURRENT_TIMESTAMP(3) WHERE kept.id = ?",
        (temp_id, star_system_id),
    )
    conn.execute("DELETE FROM star_systems WHERE id = ?", (temp_id,))
    _rename_bodies_with_prefix(conn, star_system_id, temp_name, row["name"])
    return True


def get_sector_galaxy_position(conn, sector_id):
    """
    Reads back a sector's stored galaxy-frame placement (see `schema.sql`'s
    "v4" header note) -- used by `galaxyGen.py`'s local-neighborhood mode
    to look up an existing sector's own center before enumerating its
    neighbors (`galaxyGeometry.enumerate_sectors_within_radius`).

    Args:
        conn (Connection): An open, schema-initialized connection.
        sector_id (int): The `sectors.id` to look up.

    Returns:
        dict or None: A dict with keys `center_x_pc`, `center_y_pc`,
            `center_z_pc`, `galactic_radius_pc`, `ring_index`,
            `layer_index`, `ring_slot_index`, or `None` if
            this sector has never been placed in a galaxy (the four
            galaxy-position columns NULL).

    Raises:
        ValueError: If no such `sectors` row exists.
    """
    row = conn.execute(
        """
        SELECT center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc,
               ring_index, layer_index, ring_slot_index
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
    radially outward from the galactic axis, +Y along the ring, +Z
    galactic north -- `galaxyGeometry.sector_orientation`), the same frame
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
            (`sectorGen.py`'s own standalone CLI).
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
    cell (`galaxyGeometry.SectorCell`; a sector with no grid address falls
    back to a `+/- edge_pc / 2` cube around its center) -- so the
    phenomenon's center always falls inside that sector, regardless of the phenomenon's own `radius_ly` (which may be
    far larger than the sector itself, and is free to spill into
    neighboring sectors -- see `queryDb.phenomena_near_sector`, which finds
    those by real distance, not by this row's `sector_id`).

    Args:
        conn (Connection): An open, schema-initialized connection.
        sector_id (int): The `sectors.id` to place this phenomenon at.

    Returns:
        dict: `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc`.

    Raises:
        ValueError: If `sector_id` doesn't exist, or exists but has never
            been placed in the galaxy itself (`galaxyGen.py`) -- there is
            no galaxy-frame position to jitter around in that case.
    """
    position = get_sector_galaxy_position(conn, sector_id)
    if position is None:
        raise ValueError(
            f"sector {sector_id} has no galaxy placement of its own -- see galaxyGen.py"
        )
    edge_row = conn.execute("SELECT edge_mpc FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    edge_pc = mpc_to_pc(edge_row["edge_mpc"])

    center_pc = (position["center_x_pc"], position["center_y_pc"], position["center_z_pc"])
    if position["ring_index"] is not None:
        offset_pc = SectorCell.for_ring(position["ring_index"], edge_pc).sample(random)
        center_x, center_y, center_z = local_to_galaxy_pc(center_pc, offset_pc)
    else:
        center_x, center_y, center_z = (c + random.uniform(-edge_pc / 2, edge_pc / 2) for c in center_pc)
    return {
        "center_x_pc": center_x, "center_y_pc": center_y, "center_z_pc": center_z,
        "galactic_radius_pc": math.sqrt(center_x ** 2 + center_y ** 2 + center_z ** 2),
    }


def get_occupied_addresses(conn, ring_indices):
    """
    Returns every already-occupied `(ring_index, layer_index,
    ring_slot_index)` address among the given rings -- used by
    `generate.py galaxy`'s batch and neighborhood modes to skip addresses a
    sector already exists at, in one query rather than one per candidate.

    Args:
        conn (Connection): An open, schema-initialized connection.
        ring_indices (iterable): Ring indices to check.

    Returns:
        set: Address tuples already present in `sectors`. Empty if
            `ring_indices` is empty.
    """
    ring_indices = sorted(set(ring_indices))
    if not ring_indices:
        return set()

    placeholders = ", ".join("?" for _ in ring_indices)
    rows = conn.execute(
        f"SELECT ring_index, layer_index, ring_slot_index FROM sectors "
        f"WHERE ring_index IN ({placeholders}) AND ring_slot_index IS NOT NULL",
        tuple(ring_indices),
    ).fetchall()
    return {(row["ring_index"], row["layer_index"], row["ring_slot_index"]) for row in rows}


def get_sector_id_at(conn, ring_index, layer_index, ring_slot_index):
    """
    Looks up the `sectors.id` already generated at one grid address, if
    any -- used by `generate.ensure_sector_generated` to check (and, on an
    `INSERT` race, re-check) one address at a time.

    Returns:
        int or None: The existing `sectors.id`, or `None`.
    """
    row = conn.execute(
        "SELECT id FROM sectors WHERE ring_index = ? AND layer_index = ? AND ring_slot_index = ?",
        (ring_index, layer_index, ring_slot_index),
    ).fetchone()
    return row["id"] if row is not None else None


GalaxySkeletonInfo = namedtuple(
    "GalaxySkeletonInfo",
    ["shape", "edge_pc", "outer_ring_index", "expected_system_count_at_density_1", "galaxy_seed"],
    defaults=(None,),
)
"""The galaxy's stored skeleton -- everything needed to recompute any
sector's exact position/density on demand (see schema.sql's "v8" header
note). `shape` is a `galaxyDensity.GalaxyShape`; the other fields are
`galaxy_shape`'s own remaining columns, `galaxy_seed` (v51, GEN.39) the
16-byte seed or `None` for a galaxy planned before it."""


def save_galaxy_shape(shape: GalaxyShape, edge_pc, outer_ring_index,
                       expected_system_count_at_density_1, config=None, galaxy_seed=None):
    """
    Replaces the galaxy's singleton `galaxy_shape` row -- there is exactly
    one galaxy, so this always overwrites whatever was there before rather
    than inserting a second row (`galaxyPlan.py` calls this once per full
    skeleton (re)build).

    Args:
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.
        edge_pc (float): The sector edge length this skeleton was built
                         at, parsecs.
        outer_ring_index (int): The last ring with any qualifying
            content (`generate.py plan`'s own discovered galaxy edge).
        expected_system_count_at_density_1 (float): See
            `galaxySkeleton.expected_system_count_at_density_1`.
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
            conn.execute(
                """
                INSERT INTO galaxy_shape (
                    id, disk_scale_length_pc, disk_scale_height_pc,
                    bulge_scale_radius_pc, bulge_amplitude, arm_count,
                    pitch_angle_rad, arm_amplitude, spiral_reference_radius_pc,
                    spiral_reference_angle_rad, k_norm, edge_pc,
                    expected_system_count_at_density_1, outer_ring_index, galaxy_seed
                ) VALUES (1, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON DUPLICATE KEY UPDATE
                    disk_scale_length_pc = VALUES(disk_scale_length_pc),
                    disk_scale_height_pc = VALUES(disk_scale_height_pc),
                    bulge_scale_radius_pc = VALUES(bulge_scale_radius_pc),
                    bulge_amplitude = VALUES(bulge_amplitude),
                    arm_count = VALUES(arm_count),
                    pitch_angle_rad = VALUES(pitch_angle_rad),
                    arm_amplitude = VALUES(arm_amplitude),
                    spiral_reference_radius_pc = VALUES(spiral_reference_radius_pc),
                    spiral_reference_angle_rad = VALUES(spiral_reference_angle_rad),
                    k_norm = VALUES(k_norm),
                    edge_pc = VALUES(edge_pc),
                    expected_system_count_at_density_1 = VALUES(expected_system_count_at_density_1),
                    outer_ring_index = VALUES(outer_ring_index),
                    bright_star_min_luminosity_sol = NULL,
                    bright_star_seed = NULL,
                    galaxy_seed = ?
                """,
                (
                    shape.disk_scale_length_pc, shape.disk_scale_height_pc,
                    shape.bulge_scale_radius_pc, shape.bulge_amplitude, shape.arm_count,
                    shape.pitch_angle_rad, shape.arm_amplitude, shape.spiral_reference_radius_pc,
                    shape.spiral_reference_angle_rad, shape.k_norm, edge_pc,
                    expected_system_count_at_density_1, outer_ring_index, bytes(galaxy_seed), bytes(galaxy_seed),
                ),
            )
            if stored_seed != bytes(galaxy_seed):
                # DB.6: a new seed is a new galaxy -- record the code making it.
                made_by = versionKey.current()
                conn.execute(
                    "UPDATE galaxy_shape SET version_key = ?, planetgen_version = ?, python_version = ?, platform = ?"
                    " WHERE id = 1",
                    (made_by["version_key"], made_by["planetgen_version"], made_by["python_version"],
                     made_by["platform"]),
                )
    finally:
        conn.close()
    return bytes(galaxy_seed)


def has_galaxy_sectors(conn):
    """Whether any sector has been placed in the galaxy grid."""
    return conn.execute("SELECT 1 FROM sectors WHERE ring_index IS NOT NULL LIMIT 1").fetchone() is not None


GalaxyMaker = namedtuple("GalaxyMaker", ["version_key", "planetgen_version", "python_version", "platform"])
"""What made the galaxy (v52, DB.6): `galaxy_shape`'s version columns."""


def get_galaxy_maker(conn):
    """The version key and versions that made the galaxy, or `None` when it
    has never been planned (or was planned before v52)."""
    row = conn.execute("SELECT * FROM galaxy_shape WHERE id = 1").fetchone()
    if row is None or row.get("version_key") is None:
        return None
    return GalaxyMaker(row["version_key"], row["planetgen_version"], row["python_version"], row["platform"])


def start_generation_run(conn, command, arguments, run_seed=None):
    """
    Records a run that changes the galaxy (v52, DB.6) as started, with the
    running code's version key and versions and the galaxy seed it runs
    against, and commits.

    Args:
        conn (Connection): An open connection to the galaxy database.
        command (str): The `generate.py` subcommand.
        arguments (list): Its command line, without the --mysql-* and
            --debug options (stored as JSON).
        run_seed (int, optional): The run's own 128-bit seed.

    Returns:
        int: The run's `generation_runs.id`.
    """
    made_by = versionKey.current()
    seed_bytes = None if run_seed is None else int(run_seed).to_bytes(galaxySeed.SEED_BYTES, "big")
    run_id = conn.execute(
        "INSERT INTO generation_runs (command, arguments, run_seed, galaxy_seed, version_key, planetgen_version,"
        " python_version, platform) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        (command, json.dumps(list(arguments)), seed_bytes, get_galaxy_seed(conn), made_by["version_key"],
         made_by["planetgen_version"], made_by["python_version"], made_by["platform"]),
    ).lastrowid
    conn.commit()
    return run_id


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
        spiral_reference_radius_pc=row["spiral_reference_radius_pc"],
        spiral_reference_angle_rad=row["spiral_reference_angle_rad"],
        k_norm=row["k_norm"],
    )
    return GalaxySkeletonInfo(
        shape=shape, edge_pc=row["edge_pc"], outer_ring_index=row["outer_ring_index"],
        expected_system_count_at_density_1=row["expected_system_count_at_density_1"],
        galaxy_seed=None if row.get("galaxy_seed") is None else bytes(row["galaxy_seed"]),
    )


def replace_galaxy_layers(layer_extents, config=None, conn=None):
    """
    Replaces the galaxy's stored outline wholesale: every `galaxy_layer`
    row (each layer's radial bound) and every `galaxy_column` row (each
    ring's stack bound, derived from the same extents) -- a full skeleton
    build always produces the whole outline in one pass, so there is no
    partial update (matches `save_galaxy_shape`).

    Args:
        layer_extents (iterable): `(layer_index, outer_ring_index)`
            pairs, any order (`galaxySkeleton.build_layer_extents`).
        config (MySQLConfig, optional): Connection parameters. Defaults
            to `DEFAULT_MYSQL_CONFIG`. Ignored when `conn` is given.
        conn (Connection, optional): Write through this connection, inside
            the caller's own transaction (the v33 migration does), instead
            of opening and committing a new one.
    """
    from planetgen.galaxy.skeleton import column_extents

    layer_extents = list(layer_extents)

    def _write(c):
        c.execute("DELETE FROM galaxy_layer")
        c.executemany(
            "INSERT INTO galaxy_layer (layer_index, outer_ring_index) VALUES (?, ?)", layer_extents,
        )
        c.execute("DELETE FROM galaxy_column")
        c.executemany(
            "INSERT INTO galaxy_column (ring_index, layer_index_min, layer_index_max) VALUES (?, ?, ?)",
            column_extents(layer_extents),
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


def get_galaxy_layer_outer_ring(conn, layer_index):
    """
    The last ring layer `layer_index` reaches, or `None` if the layer
    holds no content (including if the skeleton was never built at all).
    The layer holds rings 0 through this one.

    Returns:
        int or None: The layer's `outer_ring_index`.
    """
    row = conn.execute(
        "SELECT outer_ring_index FROM galaxy_layer WHERE layer_index = ?", (layer_index,),
    ).fetchone()
    return row["outer_ring_index"] if row is not None else None


def get_galaxy_column(conn, ring_index):
    """
    Ring `ring_index`'s stack bound: the lowest and highest layer its
    column of sectors reaches, or `None` if no layer reaches that ring.

    Returns:
        tuple or None: `(layer_index_min, layer_index_max)`, inclusive.
    """
    row = conn.execute(
        "SELECT layer_index_min, layer_index_max FROM galaxy_column WHERE ring_index = ?", (ring_index,),
    ).fetchone()
    return (row["layer_index_min"], row["layer_index_max"]) if row is not None else None


def get_galaxy_bounds(conn):
    """
    The galaxy's stored outline as a `galaxySkeleton.GalaxyBounds`, the
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
    return GalaxyBounds(get_galaxy_layers(conn), skeleton.edge_pc)


def get_galaxy_layers(conn):
    """
    The galaxy's whole stored outline, highest layer first.

    Returns:
        list[tuple]: `(layer_index, outer_ring_index)` pairs; empty if the
            skeleton was never built.
    """
    rows = conn.execute(
        "SELECT layer_index, outer_ring_index FROM galaxy_layer ORDER BY layer_index DESC"
    ).fetchall()
    return [(row["layer_index"], row["outer_ring_index"]) for row in rows]


def save_system(star_system: StarSystem, system_config: SystemConfig, config=None) -> int:
    """
    Opens the database and persists a single, standalone `StarSystem` (no
    sector -- `sector_id`/`position` are left `None`) in one transaction.
    The single-system counterpart to `save_sector`, for `systemGen.py`
    (which, unlike `sectorGen.py`, generates one system with no natural
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
    return _save_with_retries(config, names, lambda conn: insert_star_system(conn, star_system, system_config))


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


def save_sector(sector: SpaceSector, config=None, galaxy_position=None) -> int:
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

    Returns:
        int: The new `sectors.id`.
    """
    return _save_with_retries(config, _sector_names(sector),
                              lambda conn: insert_sector(conn, sector, galaxy_position=galaxy_position))


_RETRY_JITTER = random.Random()
"""The retry pause's own stream, so waiting never moves the generation's
`random` stream (GEN.39)."""


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
    state = random.getstate()
    for attempt in range(1, SECTOR_SAVE_ATTEMPTS + 1):
        conn = get_connection(config)
        try:
            conn.execute("SET TRANSACTION ISOLATION LEVEL READ COMMITTED")
            random.setstate(state)
            with conn:
                return insert(conn)
        except Exception as exc:
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
        "markdown": bool(row["markdown"]),
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
        "slots": slots,
    })


def _star_row_to_dict(row):
    """Maps a `stars` row to `Star.from_dict`'s expected dict shape,
    inverting every unit conversion `insert_star` applies.

    `temperature` is cast back to `int` when it is a whole number --
    `Star.generate_star` always sets it via `int(round(...))`, but a `REAL`/
    `DOUBLE` column hands every numeric value back as a Python `float`, and
    `f"{self.temperature} K"` (`get_table_properties`) renders `5800` vs.
    `5800.0` differently. A fractional value is kept as-is: an anchored
    `BlackHole`'s accretion-disk temperature is a float (`random.uniform`),
    and truncating it would render e.g. 3,676,064.7 K as "3,676,064 K"
    instead of the generated "3,676,065 K".
    """
    temperature = row["temperature_k"]
    if temperature is not None and float(temperature).is_integer():
        temperature = int(temperature)
    return {
        "name": row["name"],
        "type": row["star_type"],
        "yerkes_class": row["yerkes_class"],
        "mass": row["mass_kg"],
        "radius": row["radius_km"],
        "temperature": temperature,
        "luminosity": row["luminosity_w"],
        "age": row["age_gy"],
        "lifespan": row["lifespan_gy"],
        "habitable_zone": [
            row["habitable_zone_inner_km"] / physical_constants.AU_TO_KM,
            row["habitable_zone_outer_km"] / physical_constants.AU_TO_KM,
        ],
        "system_perimeter": row["system_perimeter_km"] / physical_constants.AU_TO_KM,
        "heliosphere_radius": row["heliosphere_radius_km"] / physical_constants.AU_TO_KM,
        "galactic_orbital_speed_kms": row["galactic_orbital_speed_kms"],
        "galactic_orbital_period_gy": row["galactic_orbital_period_gy"],
        "galactic_orbital_phase_deg": row["galactic_orbital_phase_deg"],
        "galactic_min_update_interval_years": row["galactic_min_update_interval_years"],
        "a_crit_au": (
            row["wide_binary_a_crit_km"] / physical_constants.AU_TO_KM
            if row["wide_binary_a_crit_km"] is not None else None
        ),
        # v20: NULL for a planet-less star (row column, not the Star's own
        # 0.0 class-level default) -- normalize back to 0.0 the same way
        # every other "0.0 in memory, NULL when not applicable" v20 column
        # does on its own read path.
        "reflex_offset_x": row["reflex_offset_x_km"] / physical_constants.AU_TO_KM if row["reflex_offset_x_km"] is not None else 0.0,
        "reflex_offset_y": row["reflex_offset_y_km"] / physical_constants.AU_TO_KM if row["reflex_offset_y_km"] is not None else 0.0,
        "reflex_offset_z": row["reflex_offset_z_km"] / physical_constants.AU_TO_KM if row["reflex_offset_z_km"] is not None else 0.0,
    }


def _black_hole_row_to_dict(star_row, black_hole_row):
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
    data = _star_row_to_dict(star_row)
    data["mass_solar"] = black_hole_row["mass_solar"]
    data["event_horizon_radius_km"] = black_hole_row["event_horizon_radius_km"]
    data["spin"] = black_hole_row["spin"]
    data["has_accretion_disk"] = bool(black_hole_row["has_accretion_disk"])
    data["mass_class"] = black_hole_row["mass_class"]
    return data


def _neutron_star_row_to_dict(star_row, neutron_star_row):
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
    data = _star_row_to_dict(star_row)
    data["mass_solar"] = neutron_star_row["mass_solar"]
    data["spin_period_ms"] = neutron_star_row["spin_period_ms"]
    data["magnetic_field_gauss"] = neutron_star_row["magnetic_field_gauss"]
    data["pulsar_type"] = neutron_star_row["pulsar_type"]
    data["surface_temperature_k"] = neutron_star_row["surface_temperature_k"]
    return data


def _load_single_star(conn, star_system_id, system_config):
    """
    Loads the `stars` row for a single (non-binary) system and reconstructs
    the correct Python class from it: a `BlackHole`/`NeutronStar` when a
    matching `black_holes`/`neutron_stars` satellite row exists for it
    (`phenomenonGen.py --anchor-system`, see `compactRemnant.py`'s module
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
        star = BlackHole.from_dict(_black_hole_row_to_dict(star_row, black_hole_row), system_config)
    elif neutron_star_row is not None:
        star = NeutronStar.from_dict(_neutron_star_row_to_dict(star_row, neutron_star_row), system_config)
    else:
        star = Star.from_dict(_star_row_to_dict(star_row), system_config)
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
            row["binary_habitable_zone_inner_km"] / physical_constants.AU_TO_KM,
            row["binary_habitable_zone_outer_km"] / physical_constants.AU_TO_KM,
        ],
        "system_perimeter": row["binary_system_perimeter_km"] / physical_constants.AU_TO_KM,
        "heliosphere_radius": row["binary_heliosphere_radius_km"] / physical_constants.AU_TO_KM,
        "galactic_orbital_speed_kms": row["binary_galactic_orbital_speed_kms"],
        "galactic_orbital_period_gy": row["binary_galactic_orbital_period_gy"],
        "galactic_orbital_phase_deg": row["binary_galactic_orbital_phase_deg"],
        "galactic_min_update_interval_years": row["binary_galactic_min_update_interval_years"],
        "binary_mutual_orbital_period_years": row["binary_mutual_orbital_period_years"],
        "binary_mutual_orbital_speed_kms": row["binary_mutual_orbital_speed_kms"],
        "binary_mutual_orbital_inclination_deg": row["binary_mutual_orbital_inclination_deg"],
        "binary_mutual_orbital_ascending_node_deg": row["binary_mutual_orbital_ascending_node_deg"],
        "binary_mutual_orbital_phase_deg": row["binary_mutual_orbital_phase_deg"],
        "binary_mutual_min_update_interval_years": row["binary_mutual_min_update_interval_years"],
        "binary_mutual_position_x": row["binary_mutual_position_x_km"] / physical_constants.AU_TO_KM,
        "binary_mutual_position_y": row["binary_mutual_position_y_km"] / physical_constants.AU_TO_KM,
        "binary_mutual_position_z": row["binary_mutual_position_z_km"] / physical_constants.AU_TO_KM,
        "binary_primary_position_x": row["binary_primary_position_x_km"] / physical_constants.AU_TO_KM,
        "binary_primary_position_y": row["binary_primary_position_y_km"] / physical_constants.AU_TO_KM,
        "binary_primary_position_z": row["binary_primary_position_z_km"] / physical_constants.AU_TO_KM,
        "binary_secondary_position_x": row["binary_secondary_position_x_km"] / physical_constants.AU_TO_KM,
        "binary_secondary_position_y": row["binary_secondary_position_y_km"] / physical_constants.AU_TO_KM,
        "binary_secondary_position_z": row["binary_secondary_position_z_km"] / physical_constants.AU_TO_KM,
        "binary_secondary_mass_fraction": row["binary_secondary_mass_fraction"],
        "_binary_separation_au": row["binary_separation_km"] / physical_constants.AU_TO_KM,
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
        "separation_au": row["binary_separation_km"] / physical_constants.AU_TO_KM,
        "eccentricity": row["binary_eccentricity"],
        "period_years": row["binary_mutual_orbital_period_years"],
        "speed_kms": row["binary_mutual_orbital_speed_kms"],
        "periapsis_au": row["binary_periapsis_km"] / physical_constants.AU_TO_KM,
        "apoapsis_au": row["binary_apoapsis_km"] / physical_constants.AU_TO_KM,
        "inclination_deg": row["binary_mutual_orbital_inclination_deg"],
        "ascending_node_deg": row["binary_mutual_orbital_ascending_node_deg"],
        "phase_deg": row["binary_mutual_orbital_phase_deg"],
        "min_update_interval_years": row["binary_mutual_min_update_interval_years"],
        "position_x_au": row["binary_mutual_position_x_km"] / physical_constants.AU_TO_KM,
        "position_y_au": row["binary_mutual_position_y_km"] / physical_constants.AU_TO_KM,
        "position_z_au": row["binary_mutual_position_z_km"] / physical_constants.AU_TO_KM,
        "primary_position_x_au": row["binary_primary_position_x_km"] / physical_constants.AU_TO_KM,
        "primary_position_y_au": row["binary_primary_position_y_km"] / physical_constants.AU_TO_KM,
        "primary_position_z_au": row["binary_primary_position_z_km"] / physical_constants.AU_TO_KM,
        "secondary_position_x_au": row["binary_secondary_position_x_km"] / physical_constants.AU_TO_KM,
        "secondary_position_y_au": row["binary_secondary_position_y_km"] / physical_constants.AU_TO_KM,
        "secondary_position_z_au": row["binary_secondary_position_z_km"] / physical_constants.AU_TO_KM,
        "secondary_mass_fraction": row["binary_secondary_mass_fraction"],
    }


def _planet_or_moon_row_to_dict(conn, row, is_moon):
    """
    Maps a `planets` or `moons` row (identical column shape -- see
    `schema.sql`'s "v2" header note) to `Planet.from_dict`'s expected dict
    shape, inverting every unit conversion `insert_planet`/`insert_moon`
    apply and pulling in the row's evolutionary-paragraph and
    reflection-spectrum child rows.

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

    paragraph_rows = conn.execute(
        f"SELECT paragraph FROM {table_prefix}_evolutionary_paragraphs "
        f"WHERE {id_column} = ? ORDER BY position",
        (row["id"],),
    ).fetchall()
    spectrum_rows = conn.execute(
        f"SELECT spectrum_type, value FROM {table_prefix}_reflection_spectrum "
        f"WHERE {id_column} = ? ORDER BY position",
        (row["id"],),
    ).fetchall()
    visible = [r["value"] for r in spectrum_rows if r["spectrum_type"] == "visible"]
    non_visible = [r["value"] for r in spectrum_rows if r["spectrum_type"] == "non_visible"]

    min_orbit_distance_km = row["min_orbit_distance_km"]

    return {
        "is_moon": is_moon,
        "zone": row["zone"],
        "description": row["description"],
        "atm_molar_density": row["atm_molar_density"],
        "gravity": row["gravity_g"],
        "atm_density": row["atm_density"],
        "surface_temperature": row["surface_temperature_k"],
        "density": row["density_g_cm3"],
        "atmospheric_pressure": row["atmospheric_pressure_pa"],
        "mass": row["mass_kg"],
        "atmosphere": row["atmosphere"],
        "composition": row["composition"],
        "radius": row["radius_km"],
        "planet_class": row["planet_class"],
        "distance": row["distance_km"] / physical_constants.AU_TO_KM,
        "body_type": row["body_type"],
        "scale_height": row["scale_height_km"],
        "hill_radius": row["hill_radius_km"],
        "min_orbit_distance": (
            min_orbit_distance_km / physical_constants.AU_TO_KM
            if min_orbit_distance_km is not None else None
        ),
        "name": row["name"],
        "life_chemical": row["life_chemical"],
        "evolutionary_speed": row["evolutionary_speed"],
        "reflection_spectrum_visible": visible or None,
        "reflection_spectrum_non_visible": non_visible or None,
        "evolutionary_data": [r["paragraph"] for r in paragraph_rows],
        "flavor_text": row["flavor_text"],
        "flavor_text_count": row["flavor_text_count"],
        "habitable_zone": [
            row["habitable_zone_inner_km"] / physical_constants.AU_TO_KM,
            row["habitable_zone_outer_km"] / physical_constants.AU_TO_KM,
        ],
        "volume": row["volume_km3"],
        "period": row["period_years"],
        "orbital_inclination_deg": row["orbital_inclination_deg"],
        "orbital_ascending_node_deg": row["orbital_ascending_node_deg"],
        "orbital_phase_deg": row["orbital_phase_deg"],
        "position_x": row["position_x_km"] / physical_constants.AU_TO_KM,
        "position_y": row["position_y_km"] / physical_constants.AU_TO_KM,
        "position_z": row["position_z_km"] / physical_constants.AU_TO_KM,
        "orbital_speed_kms": row["orbital_speed_kms"],
        "min_update_interval_years": row["min_update_interval_years"],
        "rotation_period_hours": row["rotation_period_hours"],
        # v20: only the `planets` table has these columns (a planet's own
        # wobble from its moons) -- `moons` has no such column at all
        # (moons never host their own moons), so a moon always gets the
        # same 0.0 `Planet.reflex_offset_x` class-level default instead.
        **(
            {} if is_moon else {
                "reflex_offset_x": row["reflex_offset_x_km"] / physical_constants.AU_TO_KM if row["reflex_offset_x_km"] is not None else 0.0,
                "reflex_offset_y": row["reflex_offset_y_km"] / physical_constants.AU_TO_KM if row["reflex_offset_y_km"] is not None else 0.0,
                "reflex_offset_z": row["reflex_offset_z_km"] / physical_constants.AU_TO_KM if row["reflex_offset_z_km"] is not None else 0.0,
            }
        ),
        "moons": [],
    }


def _belt_row_to_dict(row, composition_pairs):
    """Maps an `asteroid_belts` row (plus its `asteroid_belt_composition`
    child rows) to `AsteroidBelt.from_dict`'s expected dict shape."""
    return {
        "distance": row["distance_km"] / physical_constants.AU_TO_KM,
        "lower_limit": row["lower_limit_km"] / physical_constants.AU_TO_KM,
        "upper_limit": row["upper_limit_km"] / physical_constants.AU_TO_KM,
        "body_type": "a",
        "density": row["density"],
        "composition": [[pair["component"], pair["concentration"]] for pair in composition_pairs],
    }


def _comet_row_to_dict(row, composition_rows):
    """Maps a `comets` row (plus its `comet_composition` child rows) to
    `Comet.from_dict`'s expected dict shape."""
    return {
        "name": row["name"],
        "orbit_type": row["orbit_type"],
        "period_class": row["period_class"],
        "nucleus_diameter_km": row["nucleus_diameter_km"],
        "perihelion_distance_au": row["perihelion_distance_km"] / physical_constants.AU_TO_KM,
        "eccentricity": row["eccentricity"],
        "inclination_deg": row["inclination_deg"],
        "arg_periapsis_deg": row["arg_periapsis_deg"],
        "ascending_node_deg": row["ascending_node_deg"],
        "orbital_period_years": row["orbital_period_years"],
        "mean_anomaly_deg": row["mean_anomaly_deg"],
        "parabolic_mean_anomaly": row["parabolic_mean_anomaly"],
        "min_update_interval_years": row["min_update_interval_years"],
        "primary_mass_solar": row["primary_mass_solar"],
        "is_active": bool(row["is_active"]),
        "distance_au": row["distance_km"] / physical_constants.AU_TO_KM,
        "position_x_au": row["position_x_km"] / physical_constants.AU_TO_KM,
        "position_y_au": row["position_y_km"] / physical_constants.AU_TO_KM,
        "position_z_au": row["position_z_km"] / physical_constants.AU_TO_KM,
        "orbital_speed_kms": row["orbital_speed_kms"],
        "composition": [r["component"] for r in composition_rows],
    }


def load_star_system(conn, star_system_id) -> StarSystem:
    """
    Reconstructs a full `StarSystem` -- config, star(s), every planet/
    moon/asteroid belt it contains (in original orbital order), and every
    comet it contains -- from a `star_systems` row and its related rows.

    Each reconstructed star/planet/moon/belt/comet also carries `db_id`,
    its own row's `id`, so a caller rendering per-body text
    (`systemRender.render_system_sections`) can match it back to the
    rows `queryDb.system_detail` returns.

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
            row, _star_row_to_dict(primary_row), _star_row_to_dict(secondary_row)
        )
        star = BinaryStarProxy.from_dict(proxy_data, system_config)
    elif binary_configuration == "wide":
        primary_row = conn.execute(
            "SELECT * FROM stars WHERE star_system_id = ? AND role = 'primary'", (star_system_id,)
        ).fetchone()
        secondary_row = conn.execute(
            "SELECT * FROM stars WHERE star_system_id = ? AND role = 'secondary'", (star_system_id,)
        ).fetchone()
        star = Star.from_dict(_star_row_to_dict(primary_row), system_config)
        secondary_star = Star.from_dict(_star_row_to_dict(secondary_row), system_config)
        wide_binary = WideBinaryPair.from_dict(
            _wide_binary_row_to_dict(row), system_config, star, secondary_star
        )
    else:
        star = _load_single_star(conn, star_system_id, system_config)

    # The page is titled with the system's name (`star_systems.name`). A
    # single star, or a close pair's proxy, shares it -- a close pair's
    # proxy already takes it (`_binary_proxy_row_to_dict`). A wide pair's
    # `star` is its primary, which keeps its own `stars.name` (see
    # bodyNames.py).
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

    # v20: only meaningful (non-NULL) for a 'close' pair -- see
    # schema.sql's "v20" header note and StarSystem.__init__'s own comment.
    system.binary_planetary_wobble_x = (
        row["binary_planetary_wobble_x_km"] / physical_constants.AU_TO_KM
        if row["binary_planetary_wobble_x_km"] is not None else 0.0
    )
    system.binary_planetary_wobble_y = (
        row["binary_planetary_wobble_y_km"] / physical_constants.AU_TO_KM
        if row["binary_planetary_wobble_y_km"] is not None else 0.0
    )
    system.binary_planetary_wobble_z = (
        row["binary_planetary_wobble_z_km"] / physical_constants.AU_TO_KM
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
        system.planets = _build_object_list(planet_rows, belt_rows, star)
        system.secondary_planets = []
        system.comets = _build_comet_list(comet_rows)
        system.secondary_comets = []

    system.system_flavor_text = row["system_flavor_text"]
    system.runaway_class = row["runaway_class"]
    system.runaway_speed_kms = row["runaway_speed_kms"]
    system.planet_count, system.belt_count, system.moon_count = system.count_objects()
    system.hab_count, system.m_count = system.count_habitable()
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
    # determined by its ring and edge (as generate.py builds it).
    cell = SectorCell.for_ring(row["ring_index"], edge_ly) if row["ring_index"] is not None else None
    sector = SpaceSector(row["name"], edge_ly=edge_ly, cell=cell)

    system_rows = conn.execute(
        "SELECT id, position_x_mpc, position_y_mpc, position_z_mpc FROM star_systems "
        "WHERE sector_id = ? ORDER BY id",
        (sector_id,),
    ).fetchall()

    for r in system_rows:
        star_system = load_star_system(conn, r["id"])
        position = (
            milliparsecs_to_ly(r["position_x_mpc"]),
            milliparsecs_to_ly(r["position_y_mpc"]),
            milliparsecs_to_ly(r["position_z_mpc"]),
        )
        sector.entries.append(SectorSystemEntry(star_system, position, system_config=star_system.system_config))

    return sector


def _migrate_v8_to_v9(conn):
    """
    Adds v9's orbital-motion columns (`orbital_inclination_deg`,
    `orbital_ascending_node_deg`, `orbital_phase_deg`,
    `rotation_period_hours`) to `planets`/`moons` on an existing v8
    database -- see `schema.sql`'s header comment's "v9" note.

    A fresh database never reaches this function: `_ensure_schema`'s
    `CREATE TABLE IF NOT EXISTS` already creates `planets`/`moons` with
    these columns from `schema.sql` directly. This is only for a database
    whose `planets`/`moons` tables already existed at the older, v8 shape.

    Every existing row gets `0` for the three fixed-at-generation-time
    orbital-orientation columns and `24` (an arbitrary but harmless
    Earth-like placeholder) for `rotation_period_hours` -- not physically
    meaningful for those pre-existing bodies (this generator never ran its
    actual random orbital-motion generation for them), but a simple,
    always-valid default that keeps the columns `NOT NULL` and lets
    `planetgen.cli.orbits` run against them without special-casing "does this
    row predate v9". Every body generated from this point on gets real,
    randomly-generated values instead (see
    `planetPhysics.generate_orbital_motion_properties`).

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    for table in ("planets", "moons"):
        conn.execute(
            f"""
            ALTER TABLE {table}
                ADD COLUMN orbital_inclination_deg    DOUBLE NOT NULL DEFAULT 0,
                ADD COLUMN orbital_ascending_node_deg  DOUBLE NOT NULL DEFAULT 0,
                ADD COLUMN orbital_phase_deg           DOUBLE NOT NULL DEFAULT 0,
                ADD COLUMN rotation_period_hours       DOUBLE NOT NULL DEFAULT 24
            """
        )
    conn.execute(
        "CREATE TABLE IF NOT EXISTS orbit_simulation_state ("
        "    id BIGINT UNSIGNED PRIMARY KEY CHECK (id = 1),"
        "    last_updated_at TIMESTAMP NOT NULL"
        ") ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"
    )
    conn.execute("INSERT INTO schema_migrations (version) VALUES (9)")


def _migrate_v9_to_v10(conn):
    """
    Adds v10's galactic-orbit columns (`stars.galactic_orbital_speed_kms`/
    `galactic_orbital_period_gy`, `star_systems.binary_galactic_orbital_speed_kms`/
    `binary_galactic_orbital_period_gy`) to an existing v9 database -- see
    `schema.sql`'s header comment's "v10" note.

    A fresh database never reaches this function: `_ensure_schema`'s
    `CREATE TABLE IF NOT EXISTS` already creates `stars`/`star_systems` with
    these columns from `schema.sql` directly. This is only for a database
    whose tables already existed at the older, v9 shape.

    `star_systems`'s pair stays nullable, same as every other `binary_*`
    column (no default needed -- MySQL's own column-add default is `NULL`
    for a nullable column). `stars`' pair is `NOT NULL`, matching
    `system_perimeter_km`/`heliosphere_radius_km` on the same table, so
    every pre-existing row is backfilled with the value
    `utils.calculate_galactic_orbit(physical_constants.GALACTIC_CENTER_DISTANCE_LY)`
    itself produces -- the same fixed-fallback distance every pre-v10 row
    was already implicitly generated at (`Star.galactic_center_dist_ly`
    defaults to this same constant whenever a system isn't placed in a
    sector), rather than an arbitrary placeholder like
    `_migrate_v8_to_v9`'s `rotation_period_hours` default.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    conn.execute(
        "ALTER TABLE stars "
        "ADD COLUMN galactic_orbital_speed_kms DOUBLE NOT NULL DEFAULT 205.7034718241521, "
        "ADD COLUMN galactic_orbital_period_gy DOUBLE NOT NULL DEFAULT 0.2362604506547502"
    )
    conn.execute(
        "ALTER TABLE star_systems "
        "ADD COLUMN binary_galactic_orbital_speed_kms DOUBLE, "
        "ADD COLUMN binary_galactic_orbital_period_gy DOUBLE"
    )
    conn.execute("INSERT INTO schema_migrations (version) VALUES (10)")


def _migrate_v10_to_v11(conn):
    """
    Adds v11's position/speed columns (`planets`/`moons`.
    `position_x_km`/`_y_km`/`_z_km`/`orbital_speed_kms`) to an existing v10
    database -- see `schema.sql`'s header comment's "v11" note.

    A fresh database never reaches this function: `_ensure_schema`'s
    `CREATE TABLE IF NOT EXISTS` already creates `planets`/`moons` with
    these columns from `schema.sql` directly. This is only for a database
    whose tables already existed at the older, v10 shape.

    Unlike `_migrate_v9_to_v10`'s `stars` columns (which needed a fixed
    fallback distance since a pre-v10 row carries no record of which
    sector, if any, it was placed in), every pre-existing `planets`/`moons`
    row already has everything position/speed are derived from --
    `distance_km`, `period_years`, and the v9 orbital-motion columns
    (`orbital_inclination_deg`/`orbital_ascending_node_deg`/
    `orbital_phase_deg`) -- so this backfills real values via the same
    formula `advance_orbital_phases`/`utils.orbital_position_au` use,
    computed directly in SQL, rather than an arbitrary placeholder. The
    `ADD COLUMN ... DEFAULT 0` only has to satisfy `NOT NULL` for the
    instant between the `ALTER TABLE` and the `UPDATE` that immediately
    follows it -- no row is ever left at that placeholder.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    for table in ("planets", "moons"):
        conn.execute(
            f"ALTER TABLE {table} "
            f"ADD COLUMN position_x_km DOUBLE NOT NULL DEFAULT 0, "
            f"ADD COLUMN position_y_km DOUBLE NOT NULL DEFAULT 0, "
            f"ADD COLUMN position_z_km DOUBLE NOT NULL DEFAULT 0, "
            f"ADD COLUMN orbital_speed_kms DOUBLE NOT NULL DEFAULT 0"
        )
        conn.execute(
            f"""
            UPDATE {table}
            SET position_x_km = distance_km * (
                    COS(RADIANS(orbital_ascending_node_deg)) * COS(RADIANS(orbital_phase_deg))
                    - SIN(RADIANS(orbital_ascending_node_deg)) * SIN(RADIANS(orbital_phase_deg))
                      * COS(RADIANS(orbital_inclination_deg))
                ),
                position_y_km = distance_km * (
                    SIN(RADIANS(orbital_ascending_node_deg)) * COS(RADIANS(orbital_phase_deg))
                    + COS(RADIANS(orbital_ascending_node_deg)) * SIN(RADIANS(orbital_phase_deg))
                      * COS(RADIANS(orbital_inclination_deg))
                ),
                position_z_km = distance_km * SIN(RADIANS(orbital_phase_deg)) * SIN(RADIANS(orbital_inclination_deg)),
                orbital_speed_kms = (2 * PI() * distance_km) / (period_years * {physical_constants.SECONDS_PER_YEAR})
            WHERE period_years > 0
            """
        )
    conn.execute("INSERT INTO schema_migrations (version) VALUES (11)")


def _migrate_v11_to_v12(conn):
    """
    Adds v12's `planets`/`moons.min_update_interval_years` column to an
    existing v11 database -- see `schema.sql`'s header comment's "v12"
    note. Scoped to `planets`/`moons` only: `stars`' galactic-orbit values
    are fixed forever at generation time (no periodic update mechanism
    exists for them the way `advance_orbital_phases` exists for
    `orbital_phase_deg`), so there is nothing for a floating-point update
    guard to protect there.

    A fresh database never reaches this function: `_ensure_schema`'s
    `CREATE TABLE IF NOT EXISTS` already creates `planets`/`moons` with
    this column from `schema.sql` directly. This is only for a database
    whose tables already existed at the older, v11 shape.

    Every pre-existing row already has everything this is derived from --
    `period_years` -- so, like `_migrate_v10_to_v11`, this backfills real
    values via the same formula `utils.minimum_update_interval_years`
    uses, computed directly in SQL, rather than an arbitrary placeholder.
    `math.ulp(360.0)` is evaluated once in Python and spliced in as a
    literal -- SQL has no equivalent builtin, and this value is a fixed
    property of IEEE 754 double precision, not something that could ever
    legitimately differ between rows or need recomputing.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    ulp_360_deg = math.ulp(360.0)
    for table in ("planets", "moons"):
        conn.execute(f"ALTER TABLE {table} ADD COLUMN min_update_interval_years DOUBLE NOT NULL DEFAULT 0")
        conn.execute(
            f"UPDATE {table} SET min_update_interval_years = period_years * {ulp_360_deg} / 360"
        )

    conn.execute("INSERT INTO schema_migrations (version) VALUES (12)")


def _migrate_v12_to_v13(conn):
    """
    Adds v13's star-motion columns to an existing v12 database -- see
    `schema.sql`'s header comment's "v13" note: `stars.
    galactic_orbital_phase_deg`/`galactic_min_update_interval_years`, and
    `star_systems`' matching `binary_galactic_orbital_phase_deg`/
    `binary_galactic_min_update_interval_years` plus the six
    `binary_mutual_orbital_*`/`binary_mutual_min_update_interval_years`
    columns for a binary pair's own mutual orbit.

    A fresh database never reaches this function: `_ensure_schema`'s
    `CREATE TABLE IF NOT EXISTS` already creates `stars`/`star_systems`
    with these columns from `schema.sql` directly. This is only for a
    database whose tables already existed at the older, v12 shape.

    Every pre-existing row has everything the *_min_update_interval_years
    guards and the mutual orbit's period/speed are derived from
    (`galactic_orbital_period_gy`, `binary_galactic_orbital_period_gy`,
    `binary_separation_km`, `binary_effective_mass_kg`) already stored --
    so, like `_migrate_v10_to_v11`/`_migrate_v11_to_v12`, those get
    backfilled with real derived values via the same formulas
    `utils.minimum_update_interval_years`/`planetPhysics.
    calculate_orbital_period_years`/`utils.circular_orbital_speed_kms` use,
    computed directly in SQL. There is no pre-existing record of what any
    body's random *phase*/orientation roll would have been, so
    `galactic_orbital_phase_deg`, `binary_galactic_orbital_phase_deg`, and
    the three `binary_mutual_orbital_{inclination,ascending_node,phase}_deg`
    columns get an arbitrary but harmless `0` placeholder instead -- the
    same treatment `_migrate_v8_to_v9` already gives `orbital_inclination_deg`/
    `orbital_ascending_node_deg`/`orbital_phase_deg` for pre-v9 planets/
    moons.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    ulp_360_deg = math.ulp(360.0)
    au_to_km = physical_constants.AU_TO_KM
    solar_mass_to_kg = physical_constants.SOLAR_MASS_TO_KG
    seconds_per_year = physical_constants.SECONDS_PER_YEAR

    conn.execute(
        "ALTER TABLE stars "
        "ADD COLUMN galactic_orbital_phase_deg DOUBLE NOT NULL DEFAULT 0, "
        "ADD COLUMN galactic_min_update_interval_years DOUBLE NOT NULL DEFAULT 0"
    )
    conn.execute(
        f"UPDATE stars SET galactic_min_update_interval_years = "
        f"galactic_orbital_period_gy * 1e9 * {ulp_360_deg} / 360"
    )

    conn.execute(
        "ALTER TABLE star_systems "
        "ADD COLUMN binary_galactic_orbital_phase_deg DOUBLE, "
        "ADD COLUMN binary_galactic_min_update_interval_years DOUBLE, "
        "ADD COLUMN binary_mutual_orbital_period_years DOUBLE, "
        "ADD COLUMN binary_mutual_orbital_speed_kms DOUBLE, "
        "ADD COLUMN binary_mutual_orbital_inclination_deg DOUBLE, "
        "ADD COLUMN binary_mutual_orbital_ascending_node_deg DOUBLE, "
        "ADD COLUMN binary_mutual_orbital_phase_deg DOUBLE, "
        "ADD COLUMN binary_mutual_min_update_interval_years DOUBLE"
    )
    conn.execute(
        """
        UPDATE star_systems
        SET binary_galactic_orbital_phase_deg = 0,
            binary_galactic_min_update_interval_years =
                binary_galactic_orbital_period_gy * 1e9 * ? / 360,
            binary_mutual_orbital_period_years =
                SQRT(POW(binary_separation_km / ?, 3) / (binary_effective_mass_kg / ?)),
            binary_mutual_orbital_inclination_deg = 0,
            binary_mutual_orbital_ascending_node_deg = 0,
            binary_mutual_orbital_phase_deg = 0
        WHERE is_binary = 1
        """,
        (ulp_360_deg, au_to_km, solar_mass_to_kg),
    )
    # Separate UPDATE: binary_mutual_orbital_speed_kms/min_update_interval_years
    # both depend on binary_mutual_orbital_period_years, which the UPDATE
    # above just persisted -- a later statement can safely read it back,
    # unlike relying on same-statement SET left-to-right evaluation (see
    # advance_orbital_phases's docstring for why that trick works within
    # one UPDATE but isn't needed -- or reached for -- across two).
    conn.execute(
        """
        UPDATE star_systems
        SET binary_mutual_orbital_speed_kms =
                (2 * PI() * binary_separation_km) / (binary_mutual_orbital_period_years * ?),
            binary_mutual_min_update_interval_years =
                binary_mutual_orbital_period_years * ? / 360
        WHERE is_binary = 1
        """,
        (seconds_per_year, ulp_360_deg),
    )

    conn.execute("INSERT INTO schema_migrations (version) VALUES (13)")


def _migrate_v13_to_v14(conn):
    """
    Adds v14's `star_systems.binary_mutual_position_x_km`/`_y_km`/`_z_km`
    columns to an existing v13 database -- see `schema.sql`'s header
    comment's "v14" note: the secondary's Cartesian position relative to
    the primary, kept in lockstep with `binary_mutual_orbital_phase_deg`.

    A fresh database never reaches this function: `_ensure_schema`'s
    `CREATE TABLE IF NOT EXISTS` already creates `star_systems` with these
    columns from `schema.sql` directly. This is only for a database whose
    table already existed at the older, v13 shape.

    Every pre-existing binary row already has everything this is derived
    from -- `binary_separation_km` and the v13 `binary_mutual_orbital_
    {inclination,ascending_node,phase}_deg` columns -- so, like
    `_migrate_v10_to_v11`, this backfills real values via the same formula
    `utils.orbital_position_au` uses, computed directly in SQL (the same
    trig `advance_orbital_phases` already relies on for planets'/moons'
    position columns).

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    conn.execute(
        "ALTER TABLE star_systems "
        "ADD COLUMN binary_mutual_position_x_km DOUBLE, "
        "ADD COLUMN binary_mutual_position_y_km DOUBLE, "
        "ADD COLUMN binary_mutual_position_z_km DOUBLE"
    )
    conn.execute(
        """
        UPDATE star_systems
        SET binary_mutual_position_x_km = binary_separation_km * (
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
                * SIN(RADIANS(binary_mutual_orbital_phase_deg)) * SIN(RADIANS(binary_mutual_orbital_inclination_deg))
        WHERE is_binary = 1
        """
    )

    conn.execute("INSERT INTO schema_migrations (version) VALUES (14)")


def _migrate_v14_to_v15(conn):
    """
    Adds v15's S-type (wide) binary columns to an existing v14 database --
    see `schema.sql`'s header comment's "v15" note: `star_systems` gains
    `binary_configuration`/`binary_eccentricity`/`binary_periapsis_km`/
    `binary_apoapsis_km`; `stars` gains `wide_binary_a_crit_km`;
    `asteroid_belts` gains `star_id` (+ its FK/index).

    A fresh database never reaches this function: `_ensure_schema`'s
    `CREATE TABLE IF NOT EXISTS` already creates every table at this shape
    directly from `schema.sql`. This is only for a database whose tables
    already existed at the older, v14 shape.

    Every pre-existing binary row predates S-type support entirely -- it is
    necessarily a 'close' (P-type) pair, whose mutual orbit has always been
    (documented as) circular, so it's backfilled with `binary_configuration
    = 'close'`, `binary_eccentricity = 0`, and
    `binary_periapsis_km = binary_apoapsis_km = binary_separation_km`
    (a circular orbit's periapsis/apoapsis both equal its semi-major axis).
    `wide_binary_a_crit_km` and `asteroid_belts.star_id` have no equivalent
    pre-existing data to backfill from (no 'wide' pair could have existed
    yet) -- both stay `NULL` for every pre-existing row, the correct value
    regardless (a single star's or 'close' pair's asteroid belt has never
    had a specific owning star; see `planets.star_id`'s own comment for the
    identical, pre-existing convention).

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    conn.execute(
        "ALTER TABLE star_systems "
        "ADD COLUMN binary_configuration VARCHAR(8) CHECK (binary_configuration IN ('close', 'wide')), "
        "ADD COLUMN binary_eccentricity DOUBLE, "
        "ADD COLUMN binary_periapsis_km DOUBLE, "
        "ADD COLUMN binary_apoapsis_km DOUBLE"
    )
    conn.execute(
        "UPDATE star_systems SET "
        "binary_configuration = 'close', binary_eccentricity = 0, "
        "binary_periapsis_km = binary_separation_km, binary_apoapsis_km = binary_separation_km "
        "WHERE is_binary = 1"
    )
    conn.execute("ALTER TABLE stars ADD COLUMN wide_binary_a_crit_km DOUBLE")
    conn.execute("ALTER TABLE asteroid_belts ADD COLUMN star_id BIGINT UNSIGNED")
    conn.execute(
        "ALTER TABLE asteroid_belts ADD CONSTRAINT fk_asteroid_belts_star "
        "FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE SET NULL"
    )
    conn.execute("ALTER TABLE asteroid_belts ADD KEY idx_asteroid_belts_star_id (star_id)")

    conn.execute("INSERT INTO schema_migrations (version) VALUES (15)")


def _migrate_v15_to_v16(conn):
    """
    Records schema v16 -- six new exotic-phenomenon tables (`black_holes`,
    `neutron_stars`, `nebulae`, `supernova_remnants`, `rogue_planets`,
    `interstellar_comets`, plus `interstellar_comet_composition`; see
    `schema.sql`'s "v16" header note) for `phenomenonGen.py`'s separate,
    rarer generation mode.

    Unlike every earlier migration step, this one needs no `ALTER TABLE`:
    all six are brand-new tables, and `_ensure_schema`'s
    `CREATE TABLE IF NOT EXISTS` (run on every new connection) already
    creates them directly from the current `schema.sql`, even against a
    database whose `schema_migrations` bookkeeping still says v15 -- there
    is no pre-existing table whose shape needs changing the way earlier
    migrations' `ALTER TABLE` calls did. This step exists purely to keep
    that bookkeeping counter itself accurate.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    conn.execute("INSERT INTO schema_migrations (version) VALUES (16)")


def _migrate_v16_to_v17(conn):
    """
    Adds v17's galactic-orbital-motion columns to the six pre-existing
    exotic-phenomenon tables -- see `schema.sql`'s "v17" header note.
    Unlike `_migrate_v15_to_v16`, these ARE real `ALTER TABLE` steps: v16's
    tables already existed with a fixed shape by the time this step runs,
    so (unlike a brand-new table) `_ensure_schema`'s `CREATE TABLE IF NOT
    EXISTS` would never retroactively add these columns on its own.

    `asteroid_fields`/`asteroid_field_composition` (the new seventh
    phenomenon type, also added in v17) need no `ALTER TABLE` here --
    they're brand-new tables, so `_ensure_schema`'s `CREATE TABLE IF NOT
    EXISTS` already creates them directly from the current `schema.sql`,
    the same reasoning `_migrate_v15_to_v16`'s own docstring gives for
    v16's tables.

    No backfill is needed for any pre-existing row: a v16-era `black_holes`/
    `neutron_stars` row never had this data computed at all (it was
    silently dropped at insert time, the bug this version fixes going
    forward), and a v16-era `nebulae`/`supernova_remnants`/`rogue_planets`/
    `interstellar_comets` row's object is long gone by migration time, so
    there is no live value to backfill from either way -- these columns
    simply start `NULL`/default-less for pre-existing rows (nullable on
    `black_holes`/`neutron_stars`; the four standalone-only tables get a
    real one-time default of `0` on the `ADD COLUMN` itself, immediately
    followed by dropping that default, since MySQL/MariaDB require some
    value for a `NOT NULL` column added to a non-empty table).

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    for table in ("black_holes", "neutron_stars"):
        conn.execute(
            f"ALTER TABLE {table} "
            "ADD COLUMN galactic_orbital_speed_kms DOUBLE, "
            "ADD COLUMN galactic_orbital_period_gy DOUBLE, "
            "ADD COLUMN galactic_orbital_phase_deg DOUBLE, "
            "ADD COLUMN galactic_min_update_interval_years DOUBLE"
        )

    for table in ("nebulae", "supernova_remnants", "rogue_planets", "interstellar_comets"):
        conn.execute(
            f"ALTER TABLE {table} "
            "ADD COLUMN galactic_orbital_speed_kms DOUBLE NOT NULL DEFAULT 0, "
            "ADD COLUMN galactic_orbital_period_gy DOUBLE NOT NULL DEFAULT 0, "
            "ADD COLUMN galactic_orbital_phase_deg DOUBLE NOT NULL DEFAULT 0, "
            "ADD COLUMN galactic_min_update_interval_years DOUBLE NOT NULL DEFAULT 0"
        )
        # The DEFAULT above exists only to satisfy NOT NULL for any
        # pre-existing row (per this function's own docstring, there is no
        # real value to backfill); dropped immediately after so a future
        # INSERT can't silently rely on it instead of always supplying a
        # real generated value the way insert_nebula/etc. already do.
        for column in (
            "galactic_orbital_speed_kms", "galactic_orbital_period_gy",
            "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
        ):
            conn.execute(f"ALTER TABLE {table} ALTER COLUMN {column} DROP DEFAULT")

    conn.execute("INSERT INTO schema_migrations (version) VALUES (17)")


def _migrate_v17_to_v18(conn):
    """
    Adds v18's galaxy-frame placement columns (`center_x/y/z_pc`,
    `galactic_radius_pc`) to `nebulae`/`asteroid_fields` -- see
    `schema.sql`'s "v18" header note. Real `ALTER TABLE` steps, like
    `_migrate_v16_to_v17`: both tables already existed with a fixed shape.

    No backfill: a pre-v18 row was always generated fully standalone (no
    `--sector-id` option existed yet), so there is no real placement to
    recover -- the new columns simply start NULL, the same "never placed"
    state a v18-era standalone phenomenon has too.

    Also adds the same null-together CHECK constraint (`chk_nebulae_placement`/
    `chk_asteroid_fields_placement`) a freshly created v18 database already
    gets from `schema.sql`'s `CREATE TABLE` bodies directly -- a migrated
    database would otherwise silently lack it. Named explicitly in
    `schema.sql` (rather than left anonymous, as MySQL allows) specifically
    so this step can add the identical constraint by name; MySQL fully
    supports `ADD CONSTRAINT ... CHECK` via `ALTER TABLE`, unlike the
    SQLite-era multi-column `CHECK` gap `sectors`' own v4 columns still
    have (see `docs/design/galaxy-coordinate-system.md` section 4).

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    for table in ("nebulae", "asteroid_fields"):
        conn.execute(
            f"ALTER TABLE {table} "
            "ADD COLUMN center_x_pc DOUBLE, "
            "ADD COLUMN center_y_pc DOUBLE, "
            "ADD COLUMN center_z_pc DOUBLE, "
            "ADD COLUMN galactic_radius_pc DOUBLE, "
            f"ADD KEY idx_{table}_galactic_radius_pc (galactic_radius_pc), "
            f"ADD CONSTRAINT chk_{table}_placement CHECK ("
            "(center_x_pc IS NULL) = (center_y_pc IS NULL) AND "
            "(center_y_pc IS NULL) = (center_z_pc IS NULL) AND "
            "(center_z_pc IS NULL) = (galactic_radius_pc IS NULL))"
        )

    conn.execute("INSERT INTO schema_migrations (version) VALUES (18)")


def _migrate_v18_to_v19(conn):
    """
    Records schema v19 -- new tables `comets` and `comet_composition`
    (`cometData.Comet` -- see `schema.sql`'s "v19" header note) for
    star-bound comets.

    Like `_migrate_v15_to_v16`, this needs no `ALTER TABLE`: both are
    brand-new tables, and `_ensure_schema`'s `CREATE TABLE IF NOT EXISTS`
    (run on every new connection) already creates them directly from the
    current `schema.sql`, even against a database whose
    `schema_migrations` bookkeeping still says v18 -- there is no
    pre-existing table whose shape needs changing. This step exists
    purely to keep that bookkeeping counter itself accurate.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    conn.execute("INSERT INTO schema_migrations (version) VALUES (19)")


def _migrate_v19_to_v20(conn):
    """
    Adds v20's proper two-body (barycentric) trajectory columns to
    `star_systems`, `stars`, and `planets` -- see `schema.sql`'s "v20"
    header note. Real `ALTER TABLE` steps, the same as v17's/v18's.

    Unlike `_migrate_v16_to_v17` (whose objects were "long gone by
    migration time"), every value these new columns need is fully
    derivable from data already stored on existing rows -- each star's/
    planet's own `mass_kg` and already-stored `position_x/y/z_km`, and
    (for the binary columns) each pair's already-stored
    `binary_mutual_position_x/y/z_km` plus both stars' own `mass_kg` --
    so this backfills every pre-existing row with real, correct values
    using the exact same formulas `utils.calculate_reflex_offset`/this
    file's own `advance_orbital_phases` use going forward, rather than
    leaving them `NULL` until the next `planetgen.cli.orbits` run.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    conn.execute(
        "ALTER TABLE stars "
        "ADD COLUMN reflex_offset_x_km DOUBLE, "
        "ADD COLUMN reflex_offset_y_km DOUBLE, "
        "ADD COLUMN reflex_offset_z_km DOUBLE"
    )
    conn.execute(
        "ALTER TABLE planets "
        "ADD COLUMN reflex_offset_x_km DOUBLE, "
        "ADD COLUMN reflex_offset_y_km DOUBLE, "
        "ADD COLUMN reflex_offset_z_km DOUBLE"
    )
    conn.execute(
        "ALTER TABLE star_systems "
        "ADD COLUMN binary_primary_position_x_km DOUBLE, "
        "ADD COLUMN binary_primary_position_y_km DOUBLE, "
        "ADD COLUMN binary_primary_position_z_km DOUBLE, "
        "ADD COLUMN binary_secondary_position_x_km DOUBLE, "
        "ADD COLUMN binary_secondary_position_y_km DOUBLE, "
        "ADD COLUMN binary_secondary_position_z_km DOUBLE, "
        "ADD COLUMN binary_secondary_mass_fraction DOUBLE, "
        "ADD COLUMN binary_planetary_wobble_x_km DOUBLE, "
        "ADD COLUMN binary_planetary_wobble_y_km DOUBLE, "
        "ADD COLUMN binary_planetary_wobble_z_km DOUBLE"
    )

    # Backfill each star's own reflex offset from planets it already hosts
    # -- mirrors advance_orbital_phases' own "stars" UPDATE exactly.
    conn.execute(
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

    # Backfill each planet's own reflex offset from moons it already hosts.
    conn.execute(
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

    # Backfill binary_secondary_mass_fraction for every existing binary --
    # both configurations' total mass is exactly primary.mass_kg +
    # secondary.mass_kg (for a 'close' pair this equals the already-stored
    # binary_effective_mass_kg by construction), so joining the two stars
    # directly works uniformly for both configurations without needing
    # that column at all.
    conn.execute(
        """
        UPDATE star_systems ss
        JOIN stars sp ON sp.star_system_id = ss.id AND sp.role = 'primary'
        JOIN stars ssec ON ssec.star_system_id = ss.id AND ssec.role = 'secondary'
        SET ss.binary_secondary_mass_fraction = ssec.mass_kg / (sp.mass_kg + ssec.mass_kg)
        WHERE ss.is_binary = 1
        """
    )

    # Backfill binary_primary/secondary_position from the already-stored
    # binary_mutual_position (the secondary's position relative to the
    # primary) and the mass fraction just backfilled above -- the exact
    # formula doubleStar.BinaryStarProxy.__init__/wideBinary.WideBinaryPair.
    # __init__ compute at generation time going forward.
    conn.execute(
        """
        UPDATE star_systems
        SET binary_primary_position_x_km = -binary_secondary_mass_fraction * binary_mutual_position_x_km,
            binary_primary_position_y_km = -binary_secondary_mass_fraction * binary_mutual_position_y_km,
            binary_primary_position_z_km = -binary_secondary_mass_fraction * binary_mutual_position_z_km,
            binary_secondary_position_x_km = (1 - binary_secondary_mass_fraction) * binary_mutual_position_x_km,
            binary_secondary_position_y_km = (1 - binary_secondary_mass_fraction) * binary_mutual_position_y_km,
            binary_secondary_position_z_km = (1 - binary_secondary_mass_fraction) * binary_mutual_position_z_km
        WHERE is_binary = 1 AND binary_secondary_mass_fraction IS NOT NULL
        """
    )

    # Backfill binary_planetary_wobble for existing 'close' pairs from
    # already-stored circumbinary planets (star_id IS NULL), against the
    # proxy's already-stored binary_effective_mass_kg.
    conn.execute(
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
            )
        WHERE ss.binary_configuration = 'close'
        """
    )

    conn.execute("INSERT INTO schema_migrations (version) VALUES (20)")


def _migrate_v20_to_v21(conn):
    """
    Adds v21's sector-placement columns (`sector_id`, `center_x/y/z_pc`,
    `galactic_radius_pc`) to `black_holes`/`neutron_stars` -- see
    `schema.sql`'s "v21" header note. Real `ALTER TABLE` steps, the same
    shape `_migrate_v17_to_v18` already used for `nebulae`/
    `asteroid_fields`' identical addition.

    No backfill: a pre-v21 row was always generated fully standalone
    (`sectorGen.py` didn't generate phenomena yet, and `phenomenonGen.py
    --sector-id` didn't apply to these two types yet either), so there is
    no real placement to recover -- the new columns simply start NULL,
    the same "never placed" state a v21-era standalone black hole/neutron
    star has too.

    Also adds the FK to `sectors` and the same null-together CHECK
    constraint (`chk_black_holes_placement`/`chk_neutron_stars_placement`)
    a freshly created v21 database already gets from `schema.sql`'s
    `CREATE TABLE` bodies directly -- see `_migrate_v17_to_v18`'s identical
    reasoning for why this needs to be named explicitly.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    for table in ("black_holes", "neutron_stars"):
        conn.execute(
            f"ALTER TABLE {table} "
            "ADD COLUMN sector_id BIGINT UNSIGNED, "
            "ADD COLUMN center_x_pc DOUBLE, "
            "ADD COLUMN center_y_pc DOUBLE, "
            "ADD COLUMN center_z_pc DOUBLE, "
            "ADD COLUMN galactic_radius_pc DOUBLE, "
            f"ADD KEY idx_{table}_sector_id (sector_id), "
            f"ADD KEY idx_{table}_galactic_radius_pc (galactic_radius_pc), "
            f"ADD CONSTRAINT fk_{table}_sector "
            "FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE SET NULL, "
            f"ADD CONSTRAINT chk_{table}_placement CHECK ("
            "(center_x_pc IS NULL) = (center_y_pc IS NULL) AND "
            "(center_y_pc IS NULL) = (center_z_pc IS NULL) AND "
            "(center_z_pc IS NULL) = (galactic_radius_pc IS NULL))"
        )

    conn.execute("INSERT INTO schema_migrations (version) VALUES (21)")


def _migrate_v21_to_v22(conn):
    """
    Adds v22's search-facing indexes -- see `schema.sql`'s "v22" header
    note for the full reasoning (every one of these was a genuine
    full-table scan/sort on every single `GET /api/search` visit,
    confirmed in production as the API timing out once the database grew
    past a trivial size). Real `ALTER TABLE ... ADD KEY` steps, the same
    shape every earlier migration here already uses for a plain index
    addition.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    conn.execute("ALTER TABLE sectors ADD KEY idx_sectors_name (name)")
    conn.execute("ALTER TABLE star_systems ADD KEY idx_star_systems_name (name)")
    conn.execute(
        "ALTER TABLE stars "
        "ADD KEY idx_stars_name (name), "
        "ADD KEY idx_stars_yerkes_class (yerkes_class)"
    )
    conn.execute(
        "ALTER TABLE planets "
        "ADD KEY idx_planets_name (name), "
        "ADD KEY idx_planets_planet_class (planet_class), "
        "ADD KEY idx_planets_body_type (body_type), "
        "ADD KEY idx_planets_life_chemical (life_chemical)"
    )
    conn.execute(
        "ALTER TABLE moons "
        "ADD KEY idx_moons_name (name), "
        "ADD KEY idx_moons_planet_class (planet_class), "
        "ADD KEY idx_moons_body_type (body_type), "
        "ADD KEY idx_moons_life_chemical (life_chemical)"
    )
    conn.execute("ALTER TABLE asteroid_belts ADD KEY idx_asteroid_belts_density (density)")

    conn.execute("INSERT INTO schema_migrations (version) VALUES (22)")


def _has_column(conn, table, column):
    """
    Whether `table` already has a column named `column` in the connected
    database -- used only by `_migrate_v22_to_v23` below, to make each of
    its `ADD COLUMN` steps a genuine no-op (rather than a "Duplicate
    column name" error) against a database that already has one of v23's
    new columns despite `schema_migrations` still reporting an older
    version -- notably `star_systems.wikijs_url`/`mediawiki_url`, which
    shipped in `schema.sql`'s `CREATE TABLE` slightly ahead of the
    migration step that back-fills them onto an existing database (see
    that migration's own docstring), so a database created fresh in that
    window already has them.
    """
    row = conn.execute(
        "SELECT 1 FROM information_schema.columns"
        " WHERE table_schema = DATABASE() AND table_name = ? AND column_name = ?",
        (table, column),
    ).fetchone()
    return row is not None


def _has_index(conn, table, index_name):
    """
    Whether `table` already has an index named `index_name` -- `_has_column`'s
    own reasoning, applied to an index instead of a column, used by
    `_migrate_v24_to_v25` so its `ADD KEY` is a genuine no-op (rather than
    a "Duplicate key name" error) against a database that already has v25's
    index despite `schema_migrations` still reporting an older version --
    notably a database `_ensure_schema` created fresh straight from
    `schema.sql` (which has carried this index in its `CREATE TABLE` since
    before this migration step existed).
    """
    row = conn.execute(
        "SELECT 1 FROM information_schema.statistics"
        " WHERE table_schema = DATABASE() AND table_name = ? AND index_name = ?",
        (table, index_name),
    ).fetchone()
    return row is not None


def _has_constraint(conn, table, constraint_name):
    """
    Whether `table` already has a constraint named `constraint_name` --
    `_has_index`'s own reasoning, for a named CHECK (used by
    `_migrate_v27_to_v28`).
    """
    row = conn.execute(
        "SELECT 1 FROM information_schema.table_constraints"
        " WHERE table_schema = DATABASE() AND table_name = ? AND constraint_name = ?",
        (table, constraint_name),
    ).fetchone()
    return row is not None


def _split_top_level(text):
    """`text` split at the commas outside parentheses and quotes -- the
    clauses of one `ALTER TABLE`."""
    parts, depth, quote, start = [], 0, None, 0
    for index, char in enumerate(text):
        if quote:
            if char == quote:
                quote = None
        elif char in "'\"`":
            quote = char
        elif char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == "," and depth == 0:
            parts.append(text[start:index])
            start = index + 1
    parts.append(text[start:])
    return [part.strip() for part in parts if part.strip()]


_ALTER_RE = re.compile(r"^\s*ALTER\s+TABLE\s+`?(\w+)`?\s+(.*?)\s*;?\s*$", re.IGNORECASE | re.DOTALL)
_NAME = r"`?(\w+)`?"
_CLAUSE_GUARDS = (
    # (clause pattern, whether the clause still has work to do)
    (re.compile(rf"ADD\s+COLUMN\s+{_NAME}", re.I), lambda conn, table, m: not _has_column(conn, table, m[1])),
    (re.compile(rf"ADD\s+(?:UNIQUE\s+|FULLTEXT\s+|SPATIAL\s+)?(?:KEY|INDEX)\s+{_NAME}", re.I),
     lambda conn, table, m: not _has_index(conn, table, m[1])),
    (re.compile(r"ADD\s+PRIMARY\s+KEY", re.I), lambda conn, table, m: not _has_index(conn, table, "PRIMARY")),
    (re.compile(rf"ADD\s+CONSTRAINT\s+{_NAME}", re.I), lambda conn, table, m: not _has_constraint(conn, table, m[1])),
    (re.compile(rf"DROP\s+COLUMN\s+{_NAME}", re.I), lambda conn, table, m: _has_column(conn, table, m[1])),
    (re.compile(rf"DROP\s+(?:INDEX|KEY)\s+{_NAME}", re.I), lambda conn, table, m: _has_index(conn, table, m[1])),
    (re.compile(rf"DROP\s+(?:FOREIGN\s+KEY|CHECK|CONSTRAINT)\s+{_NAME}", re.I),
     lambda conn, table, m: _has_constraint(conn, table, m[1])),
    (re.compile(rf"(?:CHANGE\s+COLUMN|RENAME\s+COLUMN)\s+{_NAME}\s+(?:TO\s+)?{_NAME}", re.I),
     lambda conn, table, m: m[1].lower() == m[2].lower() or _has_column(conn, table, m[1])),
)
_TABLE_OPTION_RE = re.compile(r"^(ALGORITHM|LOCK)\s*=", re.IGNORECASE)


def _rerunnable_sql(conn, sql):
    """
    `sql` as a migration step may safely run it again (TEST.8, TEST.9):
    an `ALTER TABLE` loses each clause whose work is already there (a
    column, index or constraint it adds that exists, or one it drops that
    is gone), and a plain `CREATE TABLE`/`DROP TABLE` gains `IF NOT
    EXISTS`/`IF EXISTS`. `None` when nothing is left to run.

    MySQL commits DDL as it goes, so a step that stopped halfway leaves its
    first changes behind, and a database migrated from an old version
    already has every brand-new table in its newest shape (`_ensure_schema`
    makes them before the steps run). Either way a step meets some of its
    own work done, which a bare `ADD COLUMN` would refuse.
    """
    match = _ALTER_RE.match(sql)
    if match:
        table, kept, actions = match[1], [], 0
        for clause in _split_top_level(match[2]):
            if _TABLE_OPTION_RE.match(clause):
                kept.append(clause)
                continue
            for pattern, needed in _CLAUSE_GUARDS:
                found = pattern.match(clause)
                if found:
                    if needed(conn, table, found):
                        kept.append(clause)
                        actions += 1
                    break
            else:
                kept.append(clause)
                actions += 1
        return f"ALTER TABLE {table} {', '.join(kept)}" if actions else None
    sql = re.sub(r"^(\s*CREATE\s+TABLE\s+)(?!IF\s+NOT\s+EXISTS)", r"\1IF NOT EXISTS ", sql, flags=re.I)
    sql = re.sub(r"^(\s*DROP\s+TABLE\s+)(?!IF\s+EXISTS)", r"\1IF EXISTS ", sql, flags=re.I)
    return re.sub(r"^(\s*)INSERT\s+INTO\s+schema_migrations\b", r"\1INSERT IGNORE INTO schema_migrations", sql,
                  flags=re.I)


class _MigrationConnection:
    """
    The connection each migration step runs on: `execute` passes every
    statement through `_rerunnable_sql` first, so a step can run against
    a database that already has part of its work -- after a crash halfway
    through it, or twice -- and finish the job rather than fail. Anything
    else goes to the real `Connection`.
    """

    def __init__(self, conn):
        self._conn = conn

    def execute(self, sql, params=()):
        rerunnable = _rerunnable_sql(self._conn, sql)
        if rerunnable is None:
            return self._conn.execute("DO 0")
        return self._conn.execute(rerunnable, params)

    def __getattr__(self, name):
        return getattr(self._conn, name)


def _migrate_v22_to_v23(conn):
    """
    Adds v23's wiki-publishing link columns -- see `schema.sql`'s header
    comment's "v23" note. `star_systems.wikijs_url`/`mediawiki_url` were
    already present in `schema.sql`'s `CREATE TABLE` (so a brand-new
    database already has them), but were never added to an existing
    database by any earlier migration step -- this is the one that
    actually does that. `sectors.wiki_url` is new outright. Every column
    is added through `_has_column`'s guard (see that function's own
    docstring) rather than a bare `ADD COLUMN`, since a database that
    already has any of them (despite `schema_migrations` still reporting
    a pre-v23 version) is exactly the case this whole migration exists to
    handle correctly rather than erroring on.

    No backfill for either table: nothing before v23 ever uploaded a page
    or recorded a link, so every row's new column(s) simply start NULL.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    if not _has_column(conn, "sectors", "wiki_url"):
        conn.execute("ALTER TABLE sectors ADD COLUMN wiki_url VARCHAR(2048)")

    for column in ("mediawiki_url", "wikijs_url"):
        if not _has_column(conn, "star_systems", column):
            conn.execute(f"ALTER TABLE star_systems ADD COLUMN {column} VARCHAR(2048)")

    conn.execute("INSERT INTO schema_migrations (version) VALUES (23)")


def _migrate_v23_to_v24(conn):
    """
    Adds v24's three name-uniqueness registry tables (`sector_name_registry`/
    `system_name_registry`/`body_name_registry`) -- see `schema.sql`'s
    "v24" header note and `planetgen/names/uniqueness.py`. Real `CREATE
    TABLE` steps (not `ALTER TABLE` -- these are new tables, not new
    columns on an existing one), copied verbatim from `schema.sql` so a
    migrated database ends up with exactly the same shape a fresh one
    gets from `_ensure_schema`.

    No backfill: an existing database may already hold duplicate names
    from before this feature existed, and this migration doesn't scan for
    or fix them -- run `planetgen.cli.dedupe` once, separately, for that
    (safe to run on a database this migration has already brought
    current, and idempotent on repeat runs).

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS sector_name_registry (
            id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
            base_name         VARCHAR(255) NOT NULL,
            occurrence_count  INT NOT NULL,
            first_sector_id   BIGINT UNSIGNED NOT NULL,

            UNIQUE (base_name),
            CONSTRAINT fk_sector_name_registry_first_sector
                FOREIGN KEY (first_sector_id) REFERENCES sectors(id) ON DELETE CASCADE
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
        """
    )
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS system_name_registry (
            id                     BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
            base_name              VARCHAR(255) NOT NULL,
            occurrence_count       INT NOT NULL,
            first_star_system_id   BIGINT UNSIGNED NOT NULL,
            diminutive_index       INT,

            UNIQUE (base_name),
            CONSTRAINT fk_system_name_registry_first_star_system
                FOREIGN KEY (first_star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
        """
    )
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS body_name_registry (
            id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
            base_name         VARCHAR(255) NOT NULL,
            occurrence_count  INT NOT NULL,
            first_body_kind   VARCHAR(8) NOT NULL CHECK (first_body_kind IN ('planet', 'moon')),
            first_body_id     BIGINT UNSIGNED NOT NULL,
            suffix_index      INT,

            UNIQUE (base_name)
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
        """
    )

    conn.execute("INSERT INTO schema_migrations (version) VALUES (24)")


def _migrate_v24_to_v25(conn):
    """
    Adds v25's spatial index on `sectors.center_x_pc`/`center_y_pc`/
    `center_z_pc` -- see `schema.sql`'s "v25" header note (the interactive
    3D Galaxy Map's live-viewport query, `queryDb.galaxy_sectors_in_view`,
    was doing a genuine full-table scan on every call without it --
    confirmed in production, the same failure mode v22's note documents
    for the pre-v22 `/api/search`). Guarded through `_has_index` (`_has_
    column`'s own reasoning, for an index) rather than a bare `ADD KEY`,
    since `schema.sql`'s `CREATE TABLE` has carried this index from the
    start of this migration's own existence -- a database `_ensure_schema`
    creates fresh already has it despite `schema_migrations` reporting an
    older version.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    if not _has_index(conn, "sectors", "idx_sectors_center"):
        conn.execute(
            "ALTER TABLE sectors ADD KEY idx_sectors_center (center_x_pc, center_y_pc, center_z_pc)"
        )

    conn.execute("INSERT INTO schema_migrations (version) VALUES (25)")


_V26_CENTER_INDEX_TABLES = ("nebulae", "asteroid_fields", "black_holes", "neutron_stars")
"""tuple: The four exotic-phenomenon tables `_migrate_v25_to_v26` adds a
spatial index to -- see that function's own docstring."""


def _migrate_v25_to_v26(conn):
    """
    Adds v26's spatial indexes on `nebulae`/`asteroid_fields`/
    `black_holes`/`neutron_stars`' own `center_x_pc`/`center_y_pc`/
    `center_z_pc` -- see `schema.sql`'s "v26" header note
    (`queryDb.phenomena_near_sector`, called on every `GET
    /api/sectors/<id>`, was doing a genuine full-table-scan-times-four on
    every call without these -- the same failure mode v25's own note
    documents for `sectors`, confirmed in production). Guarded through
    `_has_index` (`_has_column`'s own reasoning, for an index) rather than
    a bare `ADD KEY`, since `schema.sql`'s `CREATE TABLE` has carried these
    indexes from the start of this migration's own existence -- a database
    `_ensure_schema` creates fresh already has them despite
    `schema_migrations` reporting an older version.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    for table in _V26_CENTER_INDEX_TABLES:
        index_name = f"idx_{table}_center"
        if not _has_index(conn, table, index_name):
            conn.execute(f"ALTER TABLE {table} ADD KEY {index_name} (center_x_pc, center_y_pc, center_z_pc)")

    conn.execute("INSERT INTO schema_migrations (version) VALUES (26)")


TIMESTAMPED_TABLES = (
    "sectors", "star_systems",
    "black_holes", "neutron_stars", "nebulae", "supernova_remnants",
    "rogue_planets", "interstellar_comets", "asteroid_fields", "quasars",
)
"""tuple: The top-level tables carrying v27's `created_at`/`modified_at`
row timestamps -- see `schema.sql`'s "v27" header note. Child rows
(stars, planets, moons, belts, comets, ...) have no timestamps of their
own; a change to one bumps its parent system's `modified_at` instead
(`touch_star_system`)."""

_V27_CREATED_AT_DDL = "created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP"
_V27_MODIFIED_AT_DDL = (
    "modified_at TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3)"
)


def _alter_table_online(conn, table, alteration, algorithms):
    """
    Runs `ALTER TABLE {table} {alteration}` with the cheapest of
    `algorithms` the server accepts, falling back to a plain `ALTER TABLE`
    (the server's own choice) only if none of them is supported.

    For v27's migration on a large production table: `ALGORITHM=INSTANT`
    (MySQL 8.0.12+/MariaDB 10.3+) adds a column as a metadata-only change
    with no table rebuild at all, and `ALGORITHM=INPLACE, LOCK=NONE` builds
    an index (or, where INSTANT isn't available, a column) while still
    letting reads and writes through. Asking for them explicitly makes the
    server refuse rather than silently fall back to a blocking table copy,
    which is what lets this function try the next option instead. An
    older server that doesn't know a clause at all rejects it the same
    way (a syntax error), so it's handled by the same fallback.

    Args:
        conn (Connection): An open connection, mid-migration.
        table (str): Table name (a module-level constant, never input).
        alteration (str): The `ALTER TABLE` body, e.g. `ADD COLUMN ...`.
        algorithms (tuple[str]): Clauses to try in order, e.g.
            `("ALGORITHM=INSTANT", "ALGORITHM=INPLACE, LOCK=NONE")`.
    """
    for algorithm in algorithms:
        try:
            conn.execute(f"ALTER TABLE {table} {alteration}, {algorithm}")
            return
        except pymysql.MySQLError:
            continue
    conn.execute(f"ALTER TABLE {table} {alteration}")


def _migrate_v26_to_v27(conn):
    """
    Adds v27's `created_at`/`modified_at` row timestamps (and an index on
    `modified_at`) to every table in `TIMESTAMPED_TABLES` -- see
    `schema.sql`'s "v27" header note. `star_systems` already had
    `created_at`, so it only gains `modified_at`.

    Written to be safe on a large, live table: each column is added with
    `ALGORITHM=INSTANT` where the server supports it (no rebuild, no
    lock), falling back to an online `INPLACE` rebuild, and each index is
    built online (`ALGORITHM=INPLACE, LOCK=NONE`) -- see
    `_alter_table_online`. Every step is guarded through `_has_column`/
    `_has_index`, so a database `_ensure_schema` created fresh (which
    already has all of this) makes the whole step a no-op.

    Then backfills what can be recovered, via `_backfill_v27_timestamps`:
    a system's `modified_at` starts at its own (pre-existing, real)
    `created_at`, and a sector's `created_at`/`modified_at` both start at
    its oldest system's `created_at`. Nothing recorded when a phenomenon
    was made, so those rows (and any sector with no systems) keep the
    time this migration ran.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    for table in TIMESTAMPED_TABLES:
        columns = []
        if not _has_column(conn, table, "created_at"):
            columns.append(_V27_CREATED_AT_DDL)
        if not _has_column(conn, table, "modified_at"):
            columns.append(_V27_MODIFIED_AT_DDL)
        if columns:
            _alter_table_online(
                conn, table, ", ".join(f"ADD COLUMN {column}" for column in columns),
                ("ALGORITHM=INSTANT", "ALGORITHM=INPLACE, LOCK=NONE"),
            )

        index_name = f"idx_{table}_modified_at"
        if not _has_index(conn, table, index_name):
            _alter_table_online(
                conn, table, f"ADD KEY {index_name} (modified_at)", ("ALGORITHM=INPLACE, LOCK=NONE",),
            )

    _backfill_v27_timestamps(conn)

    conn.execute("INSERT INTO schema_migrations (version) VALUES (27)")


_V27_BACKFILL_BATCH_SIZE = 5000
"""int: How many `id`s `_backfill_v27_timestamps` covers per `UPDATE`,
committing after each one so no single statement holds row locks on a
large table for long."""


def _id_batches(conn, table):
    """
    Yields `(first_id, last_id)` ranges covering every `id` in `table`,
    `_V27_BACKFILL_BATCH_SIZE` ids at a time (ranges over the id space, so
    gaps from deleted rows just make a batch smaller).
    """
    row = conn.execute(f"SELECT MIN(id) AS lo, MAX(id) AS hi FROM {table}").fetchone()
    if row["lo"] is None:
        return
    for first_id in range(row["lo"], row["hi"] + 1, _V27_BACKFILL_BATCH_SIZE):
        yield first_id, first_id + _V27_BACKFILL_BATCH_SIZE - 1


def _backfill_v27_timestamps(conn):
    """
    `_migrate_v26_to_v27`'s backfill -- see that function's docstring for
    what's recovered. Works through each table in primary-key batches,
    committing after each, so a large live table is never locked by one
    long `UPDATE`. Every value is derived from `star_systems.created_at`,
    which this never changes, so a backfill interrupted partway just
    redoes the same work when `migrate_database` is run again (the
    `schema_migrations` row isn't written until it finishes).

    Setting `modified_at` explicitly here also keeps `ON UPDATE` from
    overwriting it with the current time.
    """
    for first_id, last_id in _id_batches(conn, "star_systems"):
        conn.execute(
            "UPDATE star_systems SET modified_at = created_at WHERE id BETWEEN ? AND ?",
            (first_id, last_id),
        )
        conn.commit()

    for first_id, last_id in _id_batches(conn, "sectors"):
        conn.execute(
            """
            UPDATE sectors s
            JOIN (
                SELECT sector_id, MIN(created_at) AS first_created
                FROM star_systems
                WHERE sector_id BETWEEN ? AND ?
                GROUP BY sector_id
            ) f ON f.sector_id = s.id
            SET s.created_at = f.first_created, s.modified_at = f.first_created
            """,
            (first_id, last_id),
        )
        conn.commit()


V28_PLACED_TABLES = ("supernova_remnants", "rogue_planets", "interstellar_comets")
"""tuple: The phenomenon tables v28 gave galaxy-frame placement columns --
see `schema.sql`'s "v28" header note."""

_V28_PLACEMENT_COLUMNS_DDL = (
    "ADD COLUMN center_x_pc DOUBLE, ADD COLUMN center_y_pc DOUBLE, "
    "ADD COLUMN center_z_pc DOUBLE, ADD COLUMN galactic_radius_pc DOUBLE"
)


def _migrate_v27_to_v28(conn):
    """
    Adds v28's galaxy-frame placement columns (plus their null-together
    CHECK and spatial/radius indexes) to every table in
    `V28_PLACED_TABLES`, then backfills a position for every existing row
    it can -- see `schema.sql`'s "v28" header note and
    `_backfill_v28_placements`.

    Safe on a large, live table the same way `_migrate_v26_to_v27` is:
    columns go in with `ALGORITHM=INSTANT` where supported and indexes are
    built online (`_alter_table_online`). Every step is guarded through
    `_has_column`/`_has_index`/`_has_constraint`, so a database
    `_ensure_schema` created fresh (which already has all of this) makes
    the whole step a no-op.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    for table in V28_PLACED_TABLES:
        if not _has_column(conn, table, "center_x_pc"):
            _alter_table_online(
                conn, table, _V28_PLACEMENT_COLUMNS_DDL,
                ("ALGORITHM=INSTANT", "ALGORITHM=INPLACE, LOCK=NONE"),
            )
        for index_name, columns in (
            (f"idx_{table}_galactic_radius_pc", "galactic_radius_pc"),
            (f"idx_{table}_center", "center_x_pc, center_y_pc, center_z_pc"),
        ):
            if not _has_index(conn, table, index_name):
                _alter_table_online(
                    conn, table, f"ADD KEY {index_name} ({columns})", ("ALGORITHM=INPLACE, LOCK=NONE",),
                )

    _backfill_v28_placements(conn)

    # The CHECKs go on last: MySQL validates every existing row when one is
    # added, and the backfill above only ever writes all four columns
    # together, so they pass.
    for table in V28_PLACED_TABLES:
        if not _has_constraint(conn, table, f"chk_{table}_placement"):
            conn.execute(
                f"ALTER TABLE {table} ADD CONSTRAINT chk_{table}_placement CHECK ("
                "(center_x_pc IS NULL) = (center_y_pc IS NULL) AND "
                "(center_y_pc IS NULL) = (center_z_pc IS NULL) AND "
                "(center_z_pc IS NULL) = (galactic_radius_pc IS NULL))"
            )

    conn.execute("INSERT INTO schema_migrations (version) VALUES (28)")


def _random_placement_in_sector(sector_row, seed):
    """
    A galaxy-frame placement at a uniformly random point inside one
    galaxy-placed sector's own (rotated) cube -- what `sectorGen` would
    have picked, for a pre-v28 row whose real in-sector position was never
    saved. Seeded, so the same row always lands on the same point.

    Args:
        sector_row (dict-like): `center_x_pc`/`center_y_pc`/`center_z_pc`/
            `edge_mpc` of the owning sector.
        seed (str): Seeds the random offset (the row's table and id).

    Returns:
        dict: `_galaxy_placement_from_sector_offset`'s return shape.
    """
    rng = random.Random(seed)
    half_edge_ly = milliparsecs_to_ly(sector_row["edge_mpc"]) / 2
    offset_ly = tuple(rng.uniform(-half_edge_ly, half_edge_ly) for _ in range(3))
    return _galaxy_placement_from_sector_offset(sector_row, offset_ly)


def _backfill_v28_placements(conn):
    """
    `_migrate_v27_to_v28`'s backfill. Every `V28_PLACED_TABLES` row that
    is linked to a galaxy-placed sector but has no center yet gets
    `_random_placement_in_sector`'s point; then every supernova remnant's
    embedded black hole/neutron star that still has no center takes its
    remnant's sector and center (it sits at the remnant's middle). Works
    in primary-key batches, committing after each, like
    `_backfill_v27_timestamps`; only still-NULL rows are touched, so an
    interrupted run just carries on where it stopped when re-run.
    """
    for table in V28_PLACED_TABLES:
        for first_id, last_id in _id_batches(conn, table):
            rows = conn.execute(
                f"""
                SELECT t.id, s.center_x_pc, s.center_y_pc, s.center_z_pc, s.edge_mpc
                FROM {table} t
                JOIN sectors s ON s.id = t.sector_id
                WHERE t.id BETWEEN ? AND ?
                  AND t.center_x_pc IS NULL AND s.center_x_pc IS NOT NULL
                """,
                (first_id, last_id),
            ).fetchall()
            for row in rows:
                placement = _random_placement_in_sector(row, f"{table}:{row['id']}")
                conn.execute(
                    f"UPDATE {table} SET center_x_pc = ?, center_y_pc = ?, center_z_pc = ?,"
                    " galactic_radius_pc = ? WHERE id = ?",
                    (*_placement_values(placement), row["id"]),
                )
            conn.commit()

    for remnant_table, id_column in (
        ("black_holes", "compact_remnant_black_hole_id"),
        ("neutron_stars", "compact_remnant_neutron_star_id"),
    ):
        conn.execute(
            f"""
            UPDATE {remnant_table} r
            JOIN supernova_remnants snr ON snr.{id_column} = r.id
            SET r.sector_id = snr.sector_id,
                r.center_x_pc = snr.center_x_pc, r.center_y_pc = snr.center_y_pc,
                r.center_z_pc = snr.center_z_pc, r.galactic_radius_pc = snr.galactic_radius_pc
            WHERE r.center_x_pc IS NULL AND r.star_id IS NULL AND snr.center_x_pc IS NOT NULL
            """
        )
        conn.commit()


V29_DROPPED_COLUMNS = ("wikitext_content", "markdown_content")
"""tuple: The `star_systems` page-text columns v29 dropped -- see
`schema.sql`'s "v29" header note."""


def _migrate_v28_to_v29(conn):
    """
    Drops `star_systems.wikitext_content`/`markdown_content` -- v29 renders
    both on demand from the rest of the system's rows instead
    (`planetgen/db/render.py`; see `schema.sql`'s "v29" header
    note). This deletes data, so take a backup first if you want the old
    stored copies (the PR that added this step describes how).

    Tries `ALGORITHM=INSTANT` (a metadata-only drop on MySQL 8.0.29+/
    MariaDB 10.4+), then an online `INPLACE` rebuild, via
    `_alter_table_online`. Guarded through `_has_column`, so a database
    `_ensure_schema` created fresh (which never has these columns) makes
    this a no-op.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    present = [column for column in V29_DROPPED_COLUMNS if _has_column(conn, "star_systems", column)]
    if present:
        _alter_table_online(
            conn, "star_systems", ", ".join(f"DROP COLUMN {column}" for column in present),
            ("ALGORITHM=INSTANT", "ALGORITHM=INPLACE, LOCK=NONE"),
        )

    conn.execute("INSERT INTO schema_migrations (version) VALUES (29)")


def _migrate_v29_to_v30(conn):
    """
    Cleans up two kinds of bad surface-condition values older generator
    code saved on `planets`/`moons` -- see `schema.sql`'s "v30" header
    note. No columns change, so a database `_ensure_schema` created fresh
    has nothing to fix and every `UPDATE` here matches no rows.

    - Airless bodies (`atmosphere = 'None'`) that kept a previous class's
      `atm_density`/`atm_molar_density`/`scale_height_km` after
      `planetPhysics.reconcile_zone_and_class` moved them into an airless
      class: set back to NULL, like every other airless body.
    - `surface_temperature_k` below the cosmic microwave background
      (`physical_constants.COSMIC_BACKGROUND_TEMPERATURE_K`): raised to it.

    Deliberately leaves `star_systems.modified_at` alone, like the orbit
    ticks do (see `schema.sql`'s "v27" note): nothing the Galaxy Map tile
    cache draws depends on these values.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    floor_k = physical_constants.COSMIC_BACKGROUND_TEMPERATURE_K
    for table in ("planets", "moons"):
        conn.execute(
            f"UPDATE {table} SET atm_density = NULL, atm_molar_density = NULL, scale_height_km = NULL "
            "WHERE atmosphere = 'None' AND (atm_density IS NOT NULL OR atm_molar_density IS NOT NULL "
            "OR scale_height_km IS NOT NULL)"
        )
        conn.execute(
            f"UPDATE {table} SET surface_temperature_k = ? WHERE surface_temperature_k < ?",
            (floor_k, floor_k),
        )

    conn.execute("INSERT INTO schema_migrations (version) VALUES (30)")


def _migrate_v30_to_v31(conn):
    """
    Records schema v31 -- the new `quasars` table (see `schema.sql`'s
    "v31" header note). Like `_migrate_v15_to_v16`, a brand-new table
    needs no `ALTER TABLE`: `_ensure_schema` already created it, so this
    step only keeps the bookkeeping counter accurate.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    conn.execute("INSERT INTO schema_migrations (version) VALUES (31)")


V32_PLACED_CONTENT_TABLES = (
    "quasars", "supernova_remnants", "rogue_planets", "interstellar_comets", "nebulae", "asteroid_fields",
    "black_holes", "neutron_stars", "star_systems",
)
"""tuple: The tables whose rows in a galaxy-placed sector v32 deletes, in
an order that respects their foreign keys (supernova remnants point at
black holes and neutron stars, so they go first)."""


def _migrate_v31_to_v32(conn):
    """
    Moves galaxy placement from spherical shells to the cylindrical
    ring/layer/slot grid -- see `schema.sql`'s "v32" header note. A
    shell-addressed sector has no matching cell, so every galaxy-placed
    sector is **deleted**, together with its star systems (their planets,
    moons and stars go with them through `ON DELETE CASCADE`) and every
    phenomenon filed under it; visiting the galaxy regenerates them.
    Sectors that were never placed in the galaxy, and their contents, are
    untouched. Then:

    - `sector_vertices` is dropped (a cell's corners are closed-form).
    - `sectors` swaps `shell_index`/`shell_slot_index` for `ring_index`/
      `layer_index`/`ring_slot_index` plus their UNIQUE key.
    - `galaxy_shell_band` is dropped and `galaxy_shape.outer_shell_index`
      becomes `outer_ring_index`; `_migrate_v32_to_v33` rebuilds the
      skeleton.

    Guarded on `sectors.shell_index`, so a database `_ensure_schema`
    created fresh (already the new shape) makes this a no-op apart from
    its bookkeeping row.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    if _has_column(conn, "sectors", "shell_index"):
        placed = "SELECT id FROM sectors WHERE center_x_pc IS NOT NULL"
        for table in V32_PLACED_CONTENT_TABLES:
            conn.execute(f"DELETE FROM {table} WHERE sector_id IN ({placed})")
        conn.execute("DROP TABLE IF EXISTS sector_vertices")
        conn.execute("DELETE FROM sectors WHERE center_x_pc IS NOT NULL")

        if _has_index(conn, "sectors", "idx_sectors_shell_index"):
            conn.execute("ALTER TABLE sectors DROP INDEX idx_sectors_shell_index")
        for index_row in conn.execute(
            "SHOW INDEX FROM sectors WHERE Column_name = 'shell_index' AND Key_name <> 'PRIMARY'"
        ).fetchall():
            if _has_index(conn, "sectors", index_row["Key_name"]):
                conn.execute(f"ALTER TABLE sectors DROP INDEX `{index_row['Key_name']}`")
        conn.execute(
            "ALTER TABLE sectors DROP COLUMN shell_index, DROP COLUMN shell_slot_index, "
            "ADD COLUMN ring_index INT, ADD COLUMN layer_index INT, ADD COLUMN ring_slot_index INT, "
            "ADD UNIQUE KEY uq_sectors_address (ring_index, layer_index, ring_slot_index)"
        )

        conn.execute("DROP TABLE IF EXISTS galaxy_shell_band")
        if _has_column(conn, "galaxy_shape", "outer_shell_index"):
            conn.execute("ALTER TABLE galaxy_shape CHANGE COLUMN outer_shell_index outer_ring_index INT NOT NULL")
        # The skeleton itself is rebuilt by `_migrate_v32_to_v33`, which
        # always runs right after this step.

    conn.execute("INSERT INTO schema_migrations (version) VALUES (32)")


def _migrate_v32_to_v33(conn):
    """
    Moves the galaxy to the one sector standard -- see `schema.sql`'s
    "v33" header note: a whole-parsec edge
    (`tuning.DEFAULT_SECTOR_EDGE_PC`), `round(2*pi*(i + 1/2))`
    slots per ring, and a per-layer skeleton (`galaxy_layer`). Nearly
    every cell's address and position changes, so, as in v32, every
    galaxy-placed sector is **deleted** together with its star systems and
    every phenomenon filed under it; visiting the galaxy regenerates them.
    Sectors that were never placed in the galaxy are untouched.

    Guarded on `galaxy_ring_band`, which only a real v32 database has: a
    database `_ensure_schema` created fresh (already the new shape), or
    one v32 just cleared of shell-addressed sectors, keeps its sectors.
    `galaxy_ring_band` is dropped. If a galaxy shape is stored, the
    skeleton is rebuilt from it at the standard edge (`galaxy_shape.
    edge_pc`, `expected_system_count_at_density_1` and `outer_ring_index`
    are rewritten and `galaxy_layer` is filled), so the galaxy keeps its
    shape without re-running `generate.py plan`.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    from planetgen.galaxy.skeleton import build_layer_extents, expected_system_count_at_density_1

    has_ring_bands = conn.execute(
        "SELECT 1 FROM information_schema.tables WHERE table_schema = DATABASE() AND table_name = 'galaxy_ring_band'"
    ).fetchone() is not None
    if has_ring_bands:
        placed = "SELECT id FROM sectors WHERE center_x_pc IS NOT NULL"
        for table in V32_PLACED_CONTENT_TABLES:
            conn.execute(f"DELETE FROM {table} WHERE sector_id IN ({placed})")
        conn.execute("DELETE FROM sectors WHERE center_x_pc IS NOT NULL")
        conn.execute("DROP TABLE galaxy_ring_band")

    skeleton = get_galaxy_shape(conn)
    if skeleton is not None:
        edge_pc = float(tuning.DEFAULT_SECTOR_EDGE_PC)
        e_value = expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY)
        extents, outer_ring_index, _confirmed = build_layer_extents(skeleton.shape, edge_pc, 1.0 / e_value)
        replace_galaxy_layers(extents, conn=conn)
        conn.execute(
            "UPDATE galaxy_shape SET edge_pc = ?, expected_system_count_at_density_1 = ?, "
            "outer_ring_index = ? WHERE id = 1",
            (edge_pc, e_value, outer_ring_index),
        )

    conn.execute("INSERT INTO schema_migrations (version) VALUES (33)")


def _migrate_v33_to_v34(conn):
    """
    Drops `body_name_registry` -- see `schema.sql`'s "v34" header note.
    Planet and moon names now derive from their system's (`bodyNames.py`),
    so they need no registry. Existing rows keep their names; regenerating
    the galaxy gives them the new ones.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    conn.execute("DROP TABLE IF EXISTS body_name_registry")
    conn.execute("INSERT INTO schema_migrations (version) VALUES (34)")


def _migrate_v34_to_v35(conn):
    """
    Moves ring slot counts to the hybrid master-wedge rule
    (`galaxyGeometry.ring_sector_count`) -- see `schema.sql`'s "v35"
    header note. Slot counts change in all but 15 of the default galaxy's
    3,856 rings, so a stored `ring_slot_index` no longer names the same
    wedge: as in v32 and v33, every grid-addressed sector in a ring whose
    count changed is **deleted** with its star systems and every
    phenomenon filed under it, and regenerating the galaxy refills them.
    Sectors in the 15 unchanged rings (0, 1, 9, 10, 11, 30, ...), and sectors never
    placed on the grid, are untouched. `galaxy_layer` and `galaxy_column` are keyed by
    ring and layer, which keep their meaning, so the skeleton stays; its
    candidate counts are computed on the fly from `ring_sector_count`.

    A database `_ensure_schema` creates fresh starts at `SCHEMA_VERSION`,
    so this only ever runs on a database that holds old-rule addresses.
    The Galaxy Map's tile caches drop everything on their own: the
    deleted sectors (and the new release) make `queryDb.galaxy_changes`
    answer `full`.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    from planetgen.galaxy.geometry import ring_sector_count

    max_ring = conn.execute("SELECT MAX(ring_index) AS ring FROM sectors").fetchone()["ring"]
    changed = [
        ring for ring in range((max_ring if max_ring is not None else -1) + 1)
        if ring_sector_count(ring) != max(1, round(2 * math.pi * (ring + 0.5)))
    ]
    for start in range(0, len(changed), 1000):
        rings = changed[start:start + 1000]
        marks = ", ".join("?" * len(rings))
        placed = f"SELECT id FROM sectors WHERE ring_index IN ({marks})"
        for table in V32_PLACED_CONTENT_TABLES:
            conn.execute(f"DELETE FROM {table} WHERE sector_id IN ({placed})", tuple(rings))
        conn.execute(f"DELETE FROM sectors WHERE ring_index IN ({marks})", tuple(rings))
    conn.execute("INSERT INTO schema_migrations (version) VALUES (35)")


def _migrate_v35_to_v36(conn):
    """
    Adds `black_holes.mass_class` -- see `schema.sql`'s "v36" header note
    -- and fills it from each row's mass (`compactRemnant.
    infer_black_hole_mass_class`'s thresholds). Guarded on the column, so
    a database already at the new shape only gets its bookkeeping row.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    if not _has_column(conn, "black_holes", "mass_class"):
        conn.execute(
            "ALTER TABLE black_holes ADD COLUMN mass_class VARCHAR(16) NOT NULL DEFAULT 'stellar' AFTER name, "
            "ADD CONSTRAINT chk_black_holes_mass_class "
            "CHECK (mass_class IN ('stellar', 'intermediate', 'supermassive'))"
        )
        conn.execute(
            "UPDATE black_holes SET mass_class = CASE "
            "WHEN mass_solar >= ? THEN 'supermassive' WHEN mass_solar >= ? THEN 'intermediate' "
            "ELSE 'stellar' END, modified_at = modified_at",
            (tuning.BLACK_HOLE_SUPERMASSIVE_MASS_RANGE_SOLAR[0],
             tuning.BLACK_HOLE_INTERMEDIATE_MASS_RANGE_SOLAR[0]),
        )
    conn.execute("INSERT INTO schema_migrations (version) VALUES (36)")


def _migrate_v36_to_v37(conn):
    """
    Adds `rogue_planets.mass_bin` and `star_systems.runaway_class`/
    `runaway_speed_kms` -- see `schema.sql`'s "v37" header note -- and
    fills `mass_bin` from each existing rogue's mass
    (`roguePlanetData.infer_rogue_mass_bin`'s thresholds). Guarded per
    column.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    if not _has_column(conn, "rogue_planets", "mass_bin"):
        conn.execute(
            "ALTER TABLE rogue_planets ADD COLUMN mass_bin VARCHAR(16) NOT NULL DEFAULT 'terrestrial' AFTER planet_type"
        )
        earth = physical_constants.EARTH_MASS_TO_KG
        bins = tuning.ROGUE_PLANET_MASS_BINS
        brown_dwarf_kg = tuning.ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER[0] * physical_constants.JUPITER_MASS_TO_KG
        conn.execute(
            "UPDATE rogue_planets SET mass_bin = CASE "
            "WHEN mass_kg >= ? THEN 'brown-dwarf' WHEN mass_kg < ? THEN 'terrestrial' "
            "WHEN mass_kg < ? THEN 'sub-neptune' WHEN mass_kg < ? THEN 'saturn' ELSE 'jupiter' END, "
            "modified_at = modified_at",
            (brown_dwarf_kg, bins["terrestrial"][1] * earth, bins["sub-neptune"][1] * earth,
             bins["saturn"][1] * earth),
        )
    if not _has_column(conn, "star_systems", "runaway_class"):
        conn.execute(
            "ALTER TABLE star_systems ADD COLUMN runaway_class VARCHAR(16) AFTER system_flavor_text, "
            "ADD COLUMN runaway_speed_kms DOUBLE AFTER runaway_class"
        )
    conn.execute("INSERT INTO schema_migrations (version) VALUES (37)")


def _drop_checks_mentioning(conn, table, column):
    """
    Drops every CHECK constraint on `table` whose clause mentions
    `column`, whatever the server named it (an inline column CHECK gets an
    automatic name: `<table>_chk_<n>` on MySQL, the column's own name on
    MariaDB).
    """
    # MariaDB's unnamed table CHECKs are CONSTRAINT_<n> per table, so its
    # rows must match by table too (only its check_constraints has one).
    same_table = " AND cc.table_name = tc.table_name" if _is_mariadb_connection(conn) else ""
    rows = conn.execute(
        "SELECT tc.constraint_name AS name FROM information_schema.table_constraints tc"
        " JOIN information_schema.check_constraints cc"
        "   ON cc.constraint_schema = tc.constraint_schema AND cc.constraint_name = tc.constraint_name"
        + same_table +
        " WHERE tc.table_schema = DATABASE() AND tc.table_name = ? AND tc.constraint_type = 'CHECK'"
        "   AND cc.check_clause LIKE ?",
        (table, f"%{column}%"),
    ).fetchall()
    for name in {row["name"] for row in rows}:
        definition = _mariadb_column_check_definition(conn, table, name)
        if definition is not None:
            # MariaDB keeps an inline CHECK with its column, and DROP
            # CONSTRAINT can't reach it (error 1091); redefining the
            # column without it drops it (TEST.8).
            conn.execute(f"ALTER TABLE {table} MODIFY {definition}")
        else:
            conn.execute(f"ALTER TABLE {table} DROP CONSTRAINT `{name}`")


def _mariadb_column_check_definition(conn, table, name):
    """The definition of column `name` without its inline CHECK, when
    `name` is a MariaDB column-level CHECK (whose name is its column's);
    `None` for any other constraint."""
    row = conn.execute(
        "SELECT 1 FROM information_schema.check_constraints"
        " WHERE constraint_schema = DATABASE() AND table_name = ? AND constraint_name = ? AND level = 'Column'",
        (table, name),
    ).fetchone() if _is_mariadb_connection(conn) else None
    if row is None:
        return None
    create = conn.execute(f"SHOW CREATE TABLE {table}").fetchone()["Create Table"]
    for line in create.splitlines():
        line = line.strip().rstrip(",")
        if line.startswith(f"`{name}` "):
            start = line.index(" CHECK (")
            depth = 0
            for end in range(start + len(" CHECK "), len(line)):
                depth += {"(": 1, ")": -1}.get(line[end], 0)
                if depth == 0:
                    return (line[:start] + line[end + 1:]).rstrip()
    raise ValueError(f"no inline CHECK on {table}.{name} in SHOW CREATE TABLE")


def _is_mariadb_connection(conn):
    """Whether `conn` is to a MariaDB server (whose information_schema
    has columns MySQL's lacks, such as `check_constraints.level`)."""
    return "mariadb" in conn.execute("SELECT VERSION() AS v").fetchone()["v"].lower()


def _migrate_v37_to_v38(conn):
    """
    Adds letter classes and contents to nebulae, supernova remnants and
    asteroid fields -- see `schema.sql`'s "v38" header note. Existing rows
    get the class their stored family, shape and size most likely mean
    (`nebulaData.infer_nebula_class`, `supernovaRemnantData.
    infer_remnant_class`; asteroid fields become the "mixed" family, since
    their random mineral lists were exactly that) and that class's typical
    contents. `nebulae.nebula_type`'s CHECK widens to the five families
    (adding `diffuse`). Guarded per column.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    from planetgen.generation.phenomena.asteroid_field import asteroid_field_class
    from planetgen.generation.phenomena.nebula import infer_nebula_class, typical_class_contents
    from planetgen.generation.phenomena.supernova_remnant import infer_remnant_class

    contents_columns = (
        "ADD COLUMN dominant_species VARCHAR(255) NOT NULL DEFAULT '', "
        "ADD COLUMN density_cm3 DOUBLE NOT NULL DEFAULT 0, "
        "ADD COLUMN temperature_k DOUBLE NOT NULL DEFAULT 0, "
        "ADD COLUMN extinction_av DOUBLE NOT NULL DEFAULT 0"
    )
    contents_update = (
        "dominant_species = ?, density_cm3 = ?, temperature_k = ?, extinction_av = ?, modified_at = modified_at"
    )

    if not _has_column(conn, "nebulae", "nebula_class"):
        _drop_checks_mentioning(conn, "nebulae", "nebula_type")
        conn.execute(
            "ALTER TABLE nebulae ADD COLUMN nebula_class CHAR(1) NOT NULL DEFAULT 'D' AFTER name, "
            + contents_columns + ", "
            "ADD CONSTRAINT chk_nebulae_type CHECK "
            "(nebula_type IN ('diffuse', 'emission', 'reflection', 'planetary', 'dark'))"
        )
        for row in conn.execute("SELECT id, nebula_type, radius_ly FROM nebulae").fetchall():
            letter = infer_nebula_class(row["nebula_type"], row["radius_ly"])
            conn.execute(
                f"UPDATE nebulae SET nebula_class = ?, {contents_update} WHERE id = ?",
                (letter, *typical_class_contents(letter), row["id"]),
            )

    if not _has_column(conn, "supernova_remnants", "remnant_class"):
        conn.execute(
            "ALTER TABLE supernova_remnants ADD COLUMN remnant_class CHAR(1) NOT NULL DEFAULT 'S' AFTER name, "
            + contents_columns
        )
        rows = conn.execute("SELECT id, morphology, progenitor_type, age_years FROM supernova_remnants").fetchall()
        for row in rows:
            letter = infer_remnant_class(row["morphology"], row["progenitor_type"], row["age_years"])
            conn.execute(
                f"UPDATE supernova_remnants SET remnant_class = ?, {contents_update} WHERE id = ?",
                (letter, *typical_class_contents(letter), row["id"]),
            )

    if not _has_column(conn, "asteroid_fields", "field_class"):
        conn.execute(
            "ALTER TABLE asteroid_fields ADD COLUMN field_class VARCHAR(4) NOT NULL DEFAULT 'R1' AFTER name, "
            "ADD COLUMN composition_family VARCHAR(16) NOT NULL DEFAULT 'mixed' AFTER field_class"
        )
        for row in conn.execute("SELECT id, density, radius_ly FROM asteroid_fields").fetchall():
            conn.execute(
                "UPDATE asteroid_fields SET field_class = ?, modified_at = modified_at WHERE id = ?",
                (asteroid_field_class("mixed", row["density"], row["radius_ly"]), row["id"]),
            )

    conn.execute("INSERT INTO schema_migrations (version) VALUES (38)")


def _migrate_v38_to_v39(conn):
    """
    Adds `inside_nebula_id`/`inside_remnant_id` to every containable table
    (`CONTAINABLE_TABLES`) with their foreign keys -- see `schema.sql`'s
    "v39" header note -- then fills them with `refresh_containment` for
    every placed sector a nebula or supernova remnant reaches. Guarded per
    table.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    for table in CONTAINABLE_TABLES:
        if _has_column(conn, table, "inside_nebula_id"):
            continue
        conn.execute(
            f"ALTER TABLE {table} ADD COLUMN inside_nebula_id BIGINT UNSIGNED, "
            "ADD COLUMN inside_remnant_id BIGINT UNSIGNED, "
            f"ADD CONSTRAINT fk_{table}_inside_nebula FOREIGN KEY (inside_nebula_id) "
            "REFERENCES nebulae(id) ON DELETE SET NULL, "
            f"ADD CONSTRAINT fk_{table}_inside_remnant FOREIGN KEY (inside_remnant_id) "
            "REFERENCES supernova_remnants(id) ON DELETE SET NULL"
        )
    reached = set()
    for table, _column in CONTAINER_TABLES:
        for row in conn.execute(
            f"SELECT center_x_pc, center_y_pc, center_z_pc, radius_ly FROM {table} WHERE center_x_pc IS NOT NULL"
        ).fetchall():
            center = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
            reached.update(sectors_reached_by(conn, center, ly_to_pc(row["radius_ly"])))
    refresh_containment(conn, reached)
    conn.execute("INSERT INTO schema_migrations (version) VALUES (39)")


def _migrate_v39_to_v40(conn):
    """
    Names under one standard -- see `schema.sql`'s "v40" header note.
    Lets `system_name_registry` point at a phenomenon, indexes each
    uniquely named phenomenon table's `name`, then reserves every
    existing uniquely named phenomenon's name (renaming the ones that
    clash) and gives every comet and asteroid field its designation.
    Guarded per step.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    if not _has_column(conn, "system_name_registry", "first_object_table"):
        conn.execute(
            "ALTER TABLE system_name_registry MODIFY first_star_system_id BIGINT UNSIGNED NULL, "
            "ADD COLUMN first_object_table VARCHAR(32), ADD COLUMN first_object_id BIGINT UNSIGNED"
        )
    for table in NAMED_PHENOMENON_TABLES:
        if not _has_index(conn, table, f"idx_{table}_name"):
            conn.execute(f"ALTER TABLE {table} ADD KEY idx_{table}_name (name)")

    registered = conn.execute(
        "SELECT first_object_table, first_object_id FROM system_name_registry WHERE first_object_table IS NOT NULL"
    ).fetchall()
    registered = {(row["first_object_table"], row["first_object_id"]) for row in registered}
    cores = {("black_holes", row["compact_remnant_black_hole_id"]) for row in conn.execute(
        "SELECT compact_remnant_black_hole_id FROM supernova_remnants WHERE compact_remnant_black_hole_id IS NOT NULL"
    ).fetchall()}
    cores |= {("neutron_stars", row["compact_remnant_neutron_star_id"]) for row in conn.execute(
        "SELECT compact_remnant_neutron_star_id FROM supernova_remnants WHERE compact_remnant_neutron_star_id IS NOT NULL"
    ).fetchall()}
    for table in NAMED_PHENOMENON_TABLES:
        anchored = " WHERE star_id IS NULL" if table in ("black_holes", "neutron_stars") else ""
        for row in conn.execute(f"SELECT id, name FROM {table}{anchored} ORDER BY id").fetchall():
            if (table, row["id"]) in cores or (table, row["id"]) in registered:
                continue
            final_name, base, diminutive_index = reserve_system_name(conn, row["name"])
            if final_name != row["name"]:
                rename_phenomenon(conn, table, row["id"], final_name)
            confirm_object_name(conn, base, table, row["id"], diminutive_index)

    _designate_existing_comets(conn)
    for table in ("interstellar_comets", "asteroid_fields"):
        extra = ", field_class" if table == "asteroid_fields" else ""
        counts = {}
        codes = {}
        for row in conn.execute(f"SELECT id, sector_id{extra} FROM {table} ORDER BY id").fetchall():
            sector_id = row["sector_id"]
            counts[sector_id] = counts.get(sector_id, 0) + 1
            if sector_id not in codes:
                codes[sector_id] = sector_code(conn, sector_id)
            if table == "interstellar_comets":
                name = interstellar_comet_designation(codes[sector_id], counts[sector_id])
            else:
                name = asteroid_field_designation(row["field_class"], codes[sector_id], counts[sector_id])
            conn.execute(f"UPDATE {table} SET name = ? WHERE id = ?", (name, row["id"]))
    conn.execute("INSERT INTO schema_migrations (version) VALUES (40)")


def _migrate_v40_to_v41(conn):
    """
    Adds each placeable phenomenon's `quadrant` and the `nearest_systems`
    table -- see `schema.sql`'s "v41" header note -- then fills both for
    every placed sector (`refresh_nearest_systems`). Guarded per step.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    for table in PLACED_PHENOMENON_TABLES:
        if not _has_column(conn, table, "quadrant"):
            conn.execute(
                f"ALTER TABLE {table} ADD COLUMN quadrant VARCHAR(4) "
                "CHECK (quadrant IN ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII'))"
            )
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS nearest_systems (
            id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
            sector_id            BIGINT UNSIGNED NOT NULL,
            object_table         VARCHAR(32) NOT NULL,
            object_id            BIGINT UNSIGNED NOT NULL,
            star_system_id       BIGINT UNSIGNED,
            neighbor_rank        TINYINT NOT NULL,
            neighbor_system_id   BIGINT UNSIGNED NOT NULL,
            distance_pc          DOUBLE NOT NULL,

            UNIQUE KEY uq_nearest_systems_object_rank (object_table, object_id, neighbor_rank),
            KEY idx_nearest_systems_neighbor (neighbor_system_id),
            CONSTRAINT fk_nearest_systems_sector
                FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
            CONSTRAINT fk_nearest_systems_system
                FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
            CONSTRAINT fk_nearest_systems_neighbor
                FOREIGN KEY (neighbor_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
        """
    )
    placed = [row["id"] for row in conn.execute("SELECT id FROM sectors WHERE center_x_pc IS NOT NULL").fetchall()]
    refresh_nearest_systems(conn, placed)
    conn.execute("INSERT INTO schema_migrations (version) VALUES (41)")


def _migrate_v41_to_v42(conn):
    """
    Creates `facilities` -- see `schema.sql`'s "v42" header note. Reads
    the table's definition from `schema.sql` itself so the two can't
    drift.

    Args:
        conn (Connection): An open connection, mid-migration.
    """
    conn.execute(_schema_statement("facilities"))
    conn.execute("INSERT INTO schema_migrations (version) VALUES (42)")


def _migrate_v42_to_v43(conn):
    """
    Adds bright-star pre-placement's storage -- see `schema.sql`'s "v43"
    header note: the `bright_stars` table (empty; the next
    `generate.py plan` scatters) and `galaxy_shape`'s threshold and seed.

    Args:
        conn (Connection): An open connection, mid-migration.
    """
    if not _has_column(conn, "galaxy_shape", "bright_star_min_luminosity_sol"):
        conn.execute("ALTER TABLE galaxy_shape ADD COLUMN bright_star_min_luminosity_sol DOUBLE, "
                     "ADD COLUMN bright_star_seed BIGINT UNSIGNED")
    conn.execute(_schema_statement("bright_stars"))
    conn.execute("INSERT INTO schema_migrations (version) VALUES (43)")


def _migrate_v43_to_v44(conn):
    """
    Adds population and politics' storage (POP.1 to POP.4) -- see
    `schema.sql`'s "v44" header note: the `species`, `polities`,
    `system_owners` and `population_state` tables, empty until the next
    `generate.py population` pass.

    Args:
        conn (Connection): An open connection, mid-migration.
    """
    for table in ("species", "polities", "system_owners", "population_state"):
        conn.execute(_schema_statement(table))
    conn.execute("INSERT INTO schema_migrations (version) VALUES (44)")


def _migrate_v44_to_v45(conn):
    """
    Adds `id_blocks` (PERF.13) -- see `schema.sql`'s "v45" header note.
    Empty: each table's row is made on its first block, starting above
    the table's current `MAX(id)`.

    Args:
        conn (Connection): An open connection, mid-migration.
    """
    conn.execute(_schema_statement("id_blocks"))
    conn.execute("INSERT INTO schema_migrations (version) VALUES (45)")
    if conn._config is not None:
        forget_id_blocks(conn._config._key())


FULLTEXT_NAME_TABLES = ("sectors", "star_systems", "stars", "planets", "moons")
"""tuple: Tables with a FULLTEXT index on `name` (v46, PERF.16)."""


def _migrate_v45_to_v46(conn):
    """
    Adds a FULLTEXT index on `name` to each of `FULLTEXT_NAME_TABLES`
    (PERF.16) -- see `schema.sql`'s "v46" header note. InnoDB rebuilds a
    table for its first FULLTEXT index and blocks writes (not reads)
    meanwhile, so on a large galaxy this step takes a while. Guarded on
    `_has_index`, so a database `_ensure_schema` created fresh just
    records the version.

    Args:
        conn (Connection): An open connection, mid-migration.
    """
    for table in FULLTEXT_NAME_TABLES:
        if not _has_index(conn, table, f"ft_{table}_name"):
            _alter_table_online(conn, table, f"ADD FULLTEXT KEY ft_{table}_name (name)",
                                ("ALGORITHM=INPLACE, LOCK=SHARED",))
    conn.execute("INSERT INTO schema_migrations (version) VALUES (46)")


ROGUE_CLASS_BACKFILL_BATCH = 20000
"""int: Rows per read when `_migrate_v46_to_v47` gives stored rogue
planets their class."""


def _migrate_v46_to_v47(conn):
    """
    Adds `rogue_planets.planet_class` (GEN.8) -- see `schema.sql`'s "v47"
    header note -- and fills it for every stored rogue planet with
    `roguePlanetData.default_rogue_planet_class` (deterministic, so a
    rerun gives the same classes). Brown dwarfs stay NULL. Read in id
    order a batch at a time, one UPDATE per class per batch.

    Args:
        conn (Connection): An open connection, mid-migration.
    """
    if not _has_column(conn, "rogue_planets", "planet_class"):
        conn.execute("ALTER TABLE rogue_planets ADD COLUMN planet_class VARCHAR(4) AFTER planet_type")
    last_id = 0
    while True:
        rows = conn.execute(
            "SELECT id, planet_type, mass_bin, mass_kg, radius_km FROM rogue_planets "
            "WHERE id > ? AND planet_class IS NULL ORDER BY id LIMIT ?",
            (last_id, ROGUE_CLASS_BACKFILL_BATCH),
        ).fetchall()
        if not rows:
            break
        by_class = {}
        for row in rows:
            code = default_rogue_planet_class(row["planet_type"], row["radius_km"], row["mass_kg"], row["mass_bin"])
            if code is not None:
                by_class.setdefault(code, []).append(row["id"])
        for code, ids in by_class.items():
            conn.execute(
                f"UPDATE rogue_planets SET planet_class = ? WHERE id IN ({', '.join('?' * len(ids))})",
                (code, *ids),
            )
        last_id = rows[-1]["id"]
    conn.execute("INSERT INTO schema_migrations (version) VALUES (47)")


ROGUE_SURFACE_COLUMNS = (
    ("age_gy", "DOUBLE"), ("internal_heat_flux_w_m2", "DOUBLE"), ("effective_temperature_k", "DOUBLE"),
    ("surface_regime", "VARCHAR(24)"), ("surface_temperature_k", "DOUBLE"), ("surface_pressure_pa", "DOUBLE"),
    ("ice_shell_thickness_km", "DOUBLE"), ("ocean_depth_km", "DOUBLE"), ("has_liquid_water", "TINYINT(1)"),
)
"""tuple: `(column, type)` of each v48 `rogue_planets` column, in
`ROGUE_SURFACE_FIELDS` order."""


def _migrate_v47_to_v48(conn):
    """
    Adds rogue planet surface conditions -- see `schema.sql`'s "v48"
    header note -- and fills them for every stored rogue planet with
    `rogueSurface.rogue_surface_conditions`, its draws seeded by the
    rogue's name (as `RoguePlanet.from_dict` does), so a rerun gives the
    same answers. `has_internal_heat` is reset to the computed answer,
    so it agrees with the new heat flow. Read in id order a batch at a
    time.

    Args:
        conn (Connection): An open connection, mid-migration.
    """
    missing = [f"ADD COLUMN {column} {kind}" for column, kind in ROGUE_SURFACE_COLUMNS
               if not _has_column(conn, "rogue_planets", column)]
    if missing:
        conn.execute(f"ALTER TABLE rogue_planets {', '.join(missing)}")
    sets = ", ".join(f"{column} = ?" for column, _kind in ROGUE_SURFACE_COLUMNS)
    last_id = 0
    while True:
        rows = conn.execute(
            "SELECT id, name, planet_type, mass_bin, mass_kg, radius_km, has_moons FROM rogue_planets "
            "WHERE id > ? AND surface_regime IS NULL ORDER BY id LIMIT ?",
            (last_id, ROGUE_CLASS_BACKFILL_BATCH),
        ).fetchall()
        if not rows:
            break
        params = []
        for row in rows:
            conditions = rogue_surface_conditions(
                row["mass_kg"], row["radius_km"], row["planet_type"], row["mass_bin"], bool(row["has_moons"]),
                random.Random(row["name"]))
            values = [conditions[field] for field in ROGUE_SURFACE_FIELDS]
            values[ROGUE_SURFACE_FIELDS.index("has_liquid_water")] = int(conditions["has_liquid_water"])
            params.append((*values, int(conditions["has_internal_heat"]), row["id"]))
        conn.executemany(f"UPDATE rogue_planets SET {sets}, has_internal_heat = ? WHERE id = ?", params)
        last_id = rows[-1]["id"]
    conn.execute("INSERT INTO schema_migrations (version) VALUES (48)")


def _migrate_v48_to_v49(conn):
    """
    Adds `bright_star_blocks` (GEN.23) -- see `schema.sql`'s "v49" header
    note. Empty: every block starts at the galaxy's scatter level.

    Args:
        conn (Connection): An open connection, mid-migration.
    """
    # Inline, not `_schema_statement`: v53 drops the table from schema.sql.
    conn.execute(
        "CREATE TABLE IF NOT EXISTS bright_star_blocks (block_ring INT NOT NULL, block_wedge INT NOT NULL,"
        " block_slab SMALLINT NOT NULL, min_luminosity_sol DOUBLE,"
        " updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,"
        " PRIMARY KEY (block_ring, block_wedge, block_slab))"
        " ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"
    )
    conn.execute("INSERT INTO schema_migrations (version) VALUES (49)")


_V50_PLACEHOLDER_DEFAULTS = {
    "planets": ("orbital_inclination_deg", "orbital_ascending_node_deg", "orbital_phase_deg",
                "rotation_period_hours", "position_x_km", "position_y_km", "position_z_km",
                "orbital_speed_kms", "min_update_interval_years"),
    "moons": ("orbital_inclination_deg", "orbital_ascending_node_deg", "orbital_phase_deg",
              "rotation_period_hours", "position_x_km", "position_y_km", "position_z_km",
              "orbital_speed_kms", "min_update_interval_years"),
    "stars": ("galactic_orbital_speed_kms", "galactic_orbital_period_gy", "galactic_orbital_phase_deg",
              "galactic_min_update_interval_years"),
    "nebulae": ("nebula_class", "dominant_species", "density_cm3", "temperature_k", "extinction_av"),
    "supernova_remnants": ("remnant_class", "dominant_species", "density_cm3", "temperature_k", "extinction_av"),
    "asteroid_fields": ("field_class", "composition_family"),
}
"""dict: `{table: columns}` the v9-v13 and v38 steps added `NOT NULL`
with a placeholder DEFAULT (to fill existing rows) and never dropped it;
`schema.sql` gives them none."""

_V50_SET_NULL_FOREIGN_KEYS = (
    ("nebulae", "fk_nebulae_sector"),
    ("asteroid_fields", "fk_asteroid_fields_sector"),
)
"""tuple: `(table, constraint)` of the `sector_id` foreign keys v16/v17
created `ON DELETE CASCADE`; `schema.sql` has had them `ON DELETE SET
NULL` since v18 with no step changing an existing database."""


def _migrate_v49_to_v50(conn):
    """
    Brings a database migrated from an old version to exactly the shape
    `schema.sql` gives a new one (TEST.8, which migrates every released
    schema and compares) -- see `schema.sql`'s "v50" header note: drops
    the placeholder DEFAULTs in `_V50_PLACEHOLDER_DEFAULTS`, and remakes
    the two `_V50_SET_NULL_FOREIGN_KEYS` as `ON DELETE SET NULL`, so
    deleting a sector keeps its nebulae and asteroid fields (unplaced)
    there too. Checks each first; a database that never had them is
    untouched. Also adds `system_configs.comets` and `.wide_binary`, so a
    stored recipe keeps `--comets`/`--wide-binary` (TEST.12 found both
    dropped on save: regenerating from the stored config lost them).

    Args:
        conn (Connection): An open connection, mid-migration.
    """
    for table, columns in _V50_PLACEHOLDER_DEFAULTS.items():
        rows = conn.execute(
            "SELECT column_name AS name FROM information_schema.columns WHERE table_schema = DATABASE()"
            f" AND table_name = ? AND column_default IS NOT NULL AND column_name IN ({', '.join('?' * len(columns))})",
            (table, *columns),
        ).fetchall()
        if rows:
            conn.execute(f"ALTER TABLE {table} "
                         + ", ".join(f"ALTER COLUMN {row['name']} DROP DEFAULT" for row in rows))
    for table, constraint in _V50_SET_NULL_FOREIGN_KEYS:
        row = conn.execute(
            "SELECT delete_rule AS delete_rule FROM information_schema.referential_constraints"
            " WHERE constraint_schema = DATABASE() AND table_name = ? AND constraint_name = ?",
            (table, constraint),
        ).fetchone()
        if row is not None and row["delete_rule"] != "SET NULL":
            conn.execute(f"ALTER TABLE {table} DROP FOREIGN KEY {constraint}")
            conn.execute(f"ALTER TABLE {table} ADD CONSTRAINT {constraint} "
                         "FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE SET NULL")
    conn.execute("ALTER TABLE system_configs"
                 " ADD COLUMN comets TINYINT(1) CHECK (comets IN (0, 1)) AFTER asteroid_belt,"
                 " ADD COLUMN wide_binary TINYINT(1) CHECK (wide_binary IN (0, 1)) AFTER binary_system")
    conn.execute("INSERT INTO schema_migrations (version) VALUES (50)")


def _migrate_v50_to_v51(conn):
    """
    Adds `galaxy_shape.galaxy_seed` (GEN.39) -- see `schema.sql`'s "v51"
    header note. Left NULL: a galaxy planned before has no seed, and its
    sectors can't be given one after the fact.

    Args:
        conn (Connection): An open connection, mid-migration.
    """
    if not _has_column(conn, "galaxy_shape", "galaxy_seed"):
        conn.execute("ALTER TABLE galaxy_shape ADD COLUMN galaxy_seed BINARY(16)")
    conn.execute("INSERT INTO schema_migrations (version) VALUES (51)")


def _migrate_v51_to_v52(conn):
    """
    Adds DB.6's record of what made the galaxy -- see `schema.sql`'s "v52"
    header note: `galaxy_shape`'s version columns (left NULL: a galaxy
    planned before doesn't say what made it) and the empty
    `generation_runs` table.

    Args:
        conn (Connection): An open connection, mid-migration.
    """
    if not _has_column(conn, "galaxy_shape", "version_key"):
        conn.execute("ALTER TABLE galaxy_shape ADD COLUMN version_key CHAR(22), ADD COLUMN planetgen_version VARCHAR(32),"
                     " ADD COLUMN python_version VARCHAR(16), ADD COLUMN platform VARCHAR(64)")
    conn.execute(_schema_statement("generation_runs"))
    conn.execute("INSERT INTO schema_migrations (version) VALUES (52)")


def _migrate_v52_to_v53(conn):
    """
    Adds the per-sector stats (GEN.44, PERF.11) -- see `schema.sql`'s
    "v53" header note: the `sector_stats` table and `galaxy_shape`'s
    density-ratio columns. Every filled grid sector gets a row at level 0;
    each `bright_star_blocks` level moves onto the block's sectors (the
    unfilled ones' level, the filled ones' level before their fill), and
    the table is dropped. The density stats start empty: no sector is
    measured after the fact (GEN.39 starts a fresh galaxy).

    Args:
        conn (Connection): An open connection, mid-migration.
    """
    conn.execute(_schema_statement("sector_stats"))
    if not _has_column(conn, "galaxy_shape", "density_ratio_avg"):
        conn.execute("ALTER TABLE galaxy_shape ADD COLUMN density_ratio_avg DOUBLE,"
                     " ADD COLUMN density_ratio_samples BIGINT UNSIGNED NOT NULL DEFAULT 0")
    conn.execute(
        "INSERT INTO sector_stats (ring_index, layer_index, ring_slot_index, bright_level_sol)"
        " SELECT ring_index, layer_index, ring_slot_index, 0 FROM sectors WHERE ring_index IS NOT NULL"
        " ON DUPLICATE KEY UPDATE bright_level_sol = 0"
    )
    has_blocks = conn.execute(
        "SELECT 1 FROM information_schema.tables WHERE table_schema = DATABASE() AND table_name = 'bright_star_blocks'"
    ).fetchone() is not None
    if has_blocks:
        from planetgen.galaxy.drill import drill_block_sectors, drill_slabs

        bounds = get_galaxy_bounds(conn)
        rows = []
        for row in conn.execute("SELECT block_ring, block_wedge, block_slab, min_luminosity_sol FROM bright_star_blocks"
                                " WHERE min_luminosity_sol IS NOT NULL").fetchall():
            block = DrillBlock(3, row["block_ring"], row["block_wedge"], row["block_slab"])
            for slab in drill_slabs(block):
                for sector in drill_block_sectors(block, slab):
                    if bounds is None or bounds.contains(sector.ring, sector.slab):
                        rows.append((sector.ring, sector.slab, sector.wedge, row["min_luminosity_sol"]))
        for start in range(0, len(rows), 1000):
            chunk = rows[start:start + 1000]
            conn.executemany(
                "INSERT IGNORE INTO sector_stats (ring_index, layer_index, ring_slot_index, bright_level_sol)"
                " VALUES (?, ?, ?, ?)",
                chunk,
            )
            # A filled sector's row (level 0) keeps the block's level as
            # the one a delete puts back.
            conn.executemany(
                "UPDATE sector_stats SET level_before_fill_sol = ? WHERE ring_index = ? AND layer_index = ?"
                " AND ring_slot_index = ? AND bright_level_sol = 0",
                [(level, ring, layer, slot) for ring, layer, slot, level in chunk],
            )
        conn.execute("DROP TABLE bright_star_blocks")
    conn.execute("INSERT INTO schema_migrations (version) VALUES (53)")


def _schema_statement(table):
    """`schema.sql`'s own `CREATE TABLE IF NOT EXISTS <table>` statement."""
    with open(SCHEMA_PATH, "r", encoding="utf-8") as handle:
        text = handle.read()
    start = text.index(f"CREATE TABLE IF NOT EXISTS {table} (")
    end = text.index(";", text.index(") ENGINE=InnoDB", start))
    return text[start:end]


def _designate_existing_comets(conn):
    """Gives every stored star-bound comet its designation
    (`cometData.comet_designation`) -- `_designate_comets`' counterpart
    for rows, for `_migrate_v39_to_v40`."""
    rows = conn.execute(
        "SELECT c.id, c.star_system_id, c.star_id, c.orbit_type, c.orbital_period_years, "
        "ss.name AS system_name, s.name AS star_name, s.role AS star_role "
        "FROM comets c JOIN star_systems ss ON ss.id = c.star_system_id "
        "LEFT JOIN stars s ON s.id = c.star_id ORDER BY c.id"
    ).fetchall()
    counts = {}
    for row in rows:
        host = row["star_name"] if row["star_id"] is not None and row["star_role"] != "single" else row["system_name"]
        key = (row["star_system_id"], host)
        counts[key] = counts.get(key, 0) + 1
        comet = SimpleNamespace(orbit_type=row["orbit_type"], orbital_period_years=row["orbital_period_years"])
        conn.execute("UPDATE comets SET name = ? WHERE id = ?",
                     (comet_designation(host, counts[key], comet), row["id"]))


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
    conn.execute("UPDATE sectors SET modified_at = CURRENT_TIMESTAMP(3) WHERE id = ?", (sector_id,))


def _migration_steps():
    """Every migration step `migrate_database` knows, oldest first, as
    `(version it brings the database to, step function)`."""
    return [
        (9, _migrate_v8_to_v9),
        (10, _migrate_v9_to_v10),
        (11, _migrate_v10_to_v11),
        (12, _migrate_v11_to_v12),
        (13, _migrate_v12_to_v13),
        (14, _migrate_v13_to_v14),
        (15, _migrate_v14_to_v15),
        (16, _migrate_v15_to_v16),
        (17, _migrate_v16_to_v17),
        (18, _migrate_v17_to_v18),
        (19, _migrate_v18_to_v19),
        (20, _migrate_v19_to_v20),
        (21, _migrate_v20_to_v21),
        (22, _migrate_v21_to_v22),
        (23, _migrate_v22_to_v23),
        (24, _migrate_v23_to_v24),
        (25, _migrate_v24_to_v25),
        (26, _migrate_v25_to_v26),
        (27, _migrate_v26_to_v27),
        (28, _migrate_v27_to_v28),
        (29, _migrate_v28_to_v29),
        (30, _migrate_v29_to_v30),
        (31, _migrate_v30_to_v31),
        (32, _migrate_v31_to_v32),
        (33, _migrate_v32_to_v33),
        (34, _migrate_v33_to_v34),
        (35, _migrate_v34_to_v35),
        (36, _migrate_v35_to_v36),
        (37, _migrate_v36_to_v37),
        (38, _migrate_v37_to_v38),
        (39, _migrate_v38_to_v39),
        (40, _migrate_v39_to_v40),
        (41, _migrate_v40_to_v41),
        (42, _migrate_v41_to_v42),
        (43, _migrate_v42_to_v43),
        (44, _migrate_v43_to_v44),
        (45, _migrate_v44_to_v45),
        (46, _migrate_v45_to_v46),
        (47, _migrate_v46_to_v47),
        (48, _migrate_v47_to_v48),
        (49, _migrate_v48_to_v49),
        (50, _migrate_v49_to_v50),
        (51, _migrate_v50_to_v51),
        (52, _migrate_v51_to_v52),
        (53, _migrate_v52_to_v53),
    ]


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
    # pending -- unless it already has tables, whose shape says which
    # version it is (DB.4).
    conn = get_connection(config, ensure_schema=False)
    try:
        version = _stored_version(conn, "schema_migrations")
        if version is None:
            version = detect_schema_version(conn)
        _refuse_newer(conn, version, SCHEMA_VERSION)
        return version, sum(1 for target, _ in _migration_steps() if version < target)
    finally:
        conn.close()


def migrate_database(config=None, on_step=None):
    """
    Brings a database's `schema_migrations` bookkeeping up to
    `SCHEMA_VERSION`, applying any migration step in between.

    `get_connection`/`_ensure_schema` always create a brand-new database
    already at `SCHEMA_VERSION` (every `CREATE TABLE IF NOT EXISTS` in
    `schema.sql` reflects the current shape directly), so this function
    only has real work to do against a database created by an older
    version of this project -- `_migrate_v8_to_v9` (added for the v9
    orbital-motion columns), `_migrate_v9_to_v10` (added for the v10
    galactic-orbit columns), `_migrate_v10_to_v11` (added for the v11
    planet/moon position columns), `_migrate_v11_to_v12` (added for
    the v12 floating-point update-guard column), `_migrate_v12_to_v13`
    (added for the v13 star-motion columns), `_migrate_v13_to_v14`
    (added for the v14 binary-mutual-orbit position columns),
    `_migrate_v14_to_v15` (added for v15's S-type/wide-binary columns),
    `_migrate_v15_to_v16` (added for v16's exotic-phenomenon tables),
    `_migrate_v16_to_v17` (added for v17's galactic-orbital-motion columns
    on those tables plus the new `asteroid_fields` phenomenon),
    `_migrate_v17_to_v18` (added for v18's galaxy-frame placement columns
    on `nebulae`/`asteroid_fields`), `_migrate_v18_to_v19` (added for
    v19's new `comets`/`comet_composition` tables), `_migrate_v19_to_v20`
    (added for v20's proper two-body/barycentric trajectory columns on
    `star_systems`/`stars`/`planets`), `_migrate_v20_to_v21` (added for
    v21's sector-placement columns on `black_holes`/`neutron_stars`),
    `_migrate_v21_to_v22` (added for v22's search-facing indexes),
    `_migrate_v22_to_v23` (added for v23's wiki-publishing link columns on
    `sectors`/`star_systems`), `_migrate_v23_to_v24` (added for v24's
    name-uniqueness registry tables), `_migrate_v24_to_v25` (added for
    v25's spatial index on `sectors`), `_migrate_v25_to_v26` (added for
    v26's spatial indexes on `nebulae`/`asteroid_fields`/`black_holes`/
    `neutron_stars`), `_migrate_v26_to_v27` (added for v27's
    `created_at`/`modified_at` row timestamps), and `_migrate_v27_to_v28`
    (added for v28's placement columns on `supernova_remnants`/
    `rogue_planets`/`interstellar_comets`), `_migrate_v28_to_v29`
    (dropping the stored wikitext/Markdown page text v29 renders on
    demand instead), `_migrate_v29_to_v30` (clearing stale atmosphere
    values on airless bodies and flooring surface temperatures at the
    cosmic background), `_migrate_v30_to_v31` (recording v31's new
    `quasars` table), and `_migrate_v31_to_v32` (moving galaxy placement
    to the cylindrical sector grid, which deletes every galaxy-placed
    sector and its contents), and `_migrate_v32_to_v33` (the whole-parsec
    sector standard and per-layer skeleton, which deletes them again), and
    `_migrate_v33_to_v34` (dropping the planet/moon name registry), and
    `_migrate_v34_to_v35` (the hybrid master-wedge slot rule, which deletes
    galaxy-placed sectors once more), and `_migrate_v35_to_v36` (black
    hole mass classes) are the migration steps so far; see
    `schema.sql`'s header comment for the versioning convention, and
    `planetgen.cli.migrate` for the CLI wrapper around this.

    Args:
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.
        on_step (callable, optional): Called as `on_step(number, total,
            from_version, to_version)` just before each pending step
            runs (`number` counts from 1 to `total`), so a caller can
            show progress (`planetgen.cli.migrate`'s progress bar).

    Returns:
        int: The database's `schema_migrations` version (always
            `SCHEMA_VERSION` after this call).

    Raises:
        SchemaTooNewError: The database is past `SCHEMA_VERSION` already
            (a newer planetGen migrated it); nothing is changed.
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
            pending = [(target, step) for target, step in _migration_steps() if version < target]
            for number, (target, step) in enumerate(pending, start=1):
                if on_step is not None:
                    on_step(number, len(pending), version, target)
                step(_MigrationConnection(conn))
                activity_log.event("DB", "migrate", db=(config or DEFAULT_MYSQL_CONFIG).database,
                                  from_version=version, to_version=target)
                version = target
            conn.commit()
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
    return row["elapsed_seconds"] / physical_constants.SECONDS_PER_YEAR


# Quasars have no galactic orbit: the nucleus sits at the center. Galactic
# positions follow the phases in `advance_galactic_positions`
# (updateOrbits.main runs both).
def advance_orbital_phases(conn, elapsed_years):
    """
    Advances every planet's and moon's `orbital_phase_deg` in place by the
    fraction of a full revolution `elapsed_years` represents, given each
    body's own already-stored `period_years` -- one set-based `UPDATE` per
    table rather than a per-row Python loop, so this stays fast regardless
    of how many bodies the database holds (see `planetgen.cli.orbits`).
    `position_x/y/z_km` are recomputed in lockstep from the *new* phase --
    position is a pure function of distance/inclination/ascending-node/
    phase, so it has no independent update of its own; it just has to move
    whenever phase does. `orbital_phase_deg` is assigned first in the
    `SET` list and every position expression below reads it back
    afterward, relying on documented single-table `UPDATE` behavior (MySQL/
    MariaDB evaluate a single table's `SET` assignments left to right, so a
    later expression sees an earlier assignment's *new* value) rather than
    a CTE -- tried first, but MariaDB (unlike MySQL 8) doesn't allow a CTE
    to be joined into a multi-table `UPDATE`; see
    `utils.orbital_position_au`'s docstring for the same formula in its
    Python form.

    A row is skipped entirely (not just a no-op write, no `UPDATE` attempt
    at all) when `elapsed_years < min_update_interval_years` -- that body's
    own precomputed floor below which the phase delta added is smaller
    than `orbital_phase_deg`'s own floating-point resolution, so the write
    is guaranteed to round back to the exact value already stored (see
    `utils.minimum_update_interval_years`'s docstring). In practice this
    floor sits many orders of magnitude below any realistic `elapsed_years`
    (`planetgen.cli.orbits` runs "once a month or so"), so the guard exists for
    correctness against a caller advancing time in much smaller steps
    (e.g. a fast-forward simulation), not because today's actual usage
    pattern comes close to tripping it.

    `orbital_inclination_deg`/`orbital_ascending_node_deg`/
    `rotation_period_hours` are untouched -- fixed at generation time, per
    `planetPhysics.generate_orbital_motion_properties`. `orbital_speed_kms`/
    `min_update_interval_years` are also untouched -- both constant around
    a circular orbit, only changing if `distance_km`/`period_years`
    themselves do (never from phase advancing alone).

    Also advances `stars.galactic_orbital_phase_deg` (guarded by that row's
    own `galactic_min_update_interval_years`, `galactic_orbital_period_gy`
    converted from Gy to years the same way `Star.__init__` computes the
    guard itself) and, for a binary `star_systems` row,
    `binary_mutual_orbital_phase_deg` (both binary configurations --
    `binary_configuration` 'close' and 'wide' alike, see `schema.sql`'s
    "v15" header note on why these columns are shared) and, separately,
    `binary_galactic_orbital_phase_deg` (a 'close' pair only -- a 'wide'
    pair's two stars already each advance their own galactic phase via the
    per-row `stars` `UPDATE` above, individually, since they're real
    stored `stars` rows rather than a merged proxy). These are TWO
    separate `UPDATE`s, not one combined statement: an earlier version
    combined them, guarded by "either interval qualifies" against a single
    `WHERE` that (incorrectly) required BOTH `binary_galactic_orbital_period_gy`
    and `binary_mutual_orbital_period_years` to be positive -- for a 'wide'
    pair, `binary_galactic_orbital_period_gy` is always NULL (see
    `schema.sql`'s "v15" note), which made that combined `WHERE` silently
    false for every 'wide' row, forever, so its mutual-orbit phase (and
    the P-type-only-in-spirit galactic phase this generator never intended
    to apply to it in the first place) would never advance. Splitting
    the mutual-orbit update out with its own, independent guard fixes
    this: it now applies correctly to both configurations, matching the
    exact same "own guard interval" pattern the per-table planets/moons/
    stars `UPDATE`s above already use, rather than the removed combined
    approach's non-independent one.
    `binary_mutual_position_x/y/z_km` are recomputed in lockstep from the
    *new* `binary_mutual_orbital_phase_deg` the same way a planet's/moon's
    position is -- see `schema.sql`'s "v14" note; `binary_mutual_orbital_phase_deg`
    is assigned earlier in this same `SET` list so the position expressions
    read back its new value, the identical left-to-right trick the
    planets/moons `UPDATE`s above use. The galactic orbit has no such
    position to keep in lockstep -- it's treated as planar (no
    inclination/ascending node to resolve a 3D position from), unlike the
    mutual orbit's full orbital-element set.

    Also advances every standalone exotic phenomenon's own
    `galactic_orbital_phase_deg` the identical way `stars`' is advanced
    above -- `black_holes`/`neutron_stars` (only the rows with `star_id IS
    NULL`; an anchored remnant's motion already lives on its own `stars`
    row, updated by the `stars` `UPDATE` above instead) and `nebulae`/
    `supernova_remnants`/`rogue_planets`/`interstellar_comets`/
    `asteroid_fields` (always, every row there is standalone) -- see
    `schema.sql`'s "v17" header note. Each gets its own independent
    `UPDATE` with its own `galactic_min_update_interval_years` guard, the
    same "one table, one guard" pattern every other `UPDATE` in this
    function already follows.

    v27: every `UPDATE` here on a table with a `modified_at` column
    (`star_systems` and the phenomenon tables) sets `modified_at =
    modified_at` explicitly, which stops MySQL's `ON UPDATE
    CURRENT_TIMESTAMP` from firing -- an orbit tick is the simulation
    clock moving, not the row being edited, and would otherwise mark
    every row in the galaxy as changed on each run. Anything that needs
    to know when the simulation last moved reads
    `orbit_simulation_state.last_updated_at` instead. See `schema.sql`'s
    "v27" header note.

    Also upserts `orbit_simulation_state.last_updated_at` to `NOW()` (the
    reference point the *next* call's `elapsed_years` should be measured
    from), in the same transaction, so a caller can never advance phases
    without also recording that it did.

    v20 additionally recomputes three "reflex offset"/"wobble" values --
    a proper two-body (barycentric) treatment layered on top of the
    existing relative-position model, never changing what any existing
    column means (see `schema.sql`'s "v20" header note and
    `utils.calculate_reflex_offset`'s docstring for the underlying
    formula):
      - `stars.reflex_offset_x/y/z_km`, from each star's own hosted
        planets (`planets.star_id`).
      - `planets.reflex_offset_x/y/z_km`, from each planet's own hosted
        moons (`moons.planet_id`).
      - `star_systems.binary_primary_position_*_km`/
        `binary_secondary_position_*_km`, recomputed from the same-`SET`-
        list's freshly-advanced `binary_mutual_position_*_km` and the
        stored constant `binary_secondary_mass_fraction`, folded into the
        existing mutual-orbit `UPDATE` rather than a separate statement.
      - `star_systems.binary_planetary_wobble_*_km`, a 'close' pair's
        combined pull from its own circumbinary planets (`star_id IS
        NULL`).
    Unlike every phase-advancing `UPDATE` above, these three are cheap
    values *derived from* other rows' just-advanced positions rather than
    an independently advancing phase of their own, so they're recomputed
    unconditionally on every call -- no `min_update_interval_years`-style
    guard of their own. The two per-child-table ones use a correlated
    subquery (`SELECT SUM(...) FROM <children> WHERE <parent link>`) to
    sum a parent's pull from *multiple* children in one set-based
    statement -- a different, well-supported mechanism from the CTE-in-
    multi-table-UPDATE approach flagged as unsupported above; this one
    works because it's a plain correlated scalar subquery in a
    single-table `UPDATE`'s own `SET` clause, not a join.

    Args:
        conn (Connection): An open, schema-initialized, read-write
                           connection.
        elapsed_years (float): How much simulated time has passed since
                               the reference point `elapsed_years` was
                               computed from (typically
                               `get_last_orbit_update`'s return value).
                               Must be >= 0.

    Returns:
        dict: `{table_name: rows_updated}` for every table this function
            touches -- `"planets"`, `"moons"`, `"stars"`,
            `"star_reflex_offsets"`, `"planet_reflex_offsets"`,
            `"binary_mutual_orbits"`, `"binary_planetary_wobbles"`,
            `"binary_galactic_orbits"`, `"black_holes"`, `"neutron_stars"`,
            `"nebulae"`, `"supernova_remnants"`, `"rogue_planets"`,
            `"interstellar_comets"`, `"asteroid_fields"`. A dict rather
            than a positional tuple (this function's shape before v17)
            specifically because this list keeps growing as new phenomena
            gain their own tracked motion -- a name-keyed result stays
            self-describing and immune to callers silently unpacking the
            wrong position as the list grows further.

    Raises:
        ValueError: If `elapsed_years` is negative.
    """
    if elapsed_years < 0:
        raise ValueError(f"elapsed_years must be >= 0, got {elapsed_years}")

    counts = {}
    for table in ("planets", "moons"):
        cur = conn.execute(
            f"""
            UPDATE {table}
            SET orbital_phase_deg = MOD(orbital_phase_deg + (? / period_years) * 360, 360),
                position_x_km = distance_km * (
                    COS(RADIANS(orbital_ascending_node_deg)) * COS(RADIANS(orbital_phase_deg))
                    - SIN(RADIANS(orbital_ascending_node_deg)) * SIN(RADIANS(orbital_phase_deg))
                      * COS(RADIANS(orbital_inclination_deg))
                ),
                position_y_km = distance_km * (
                    SIN(RADIANS(orbital_ascending_node_deg)) * COS(RADIANS(orbital_phase_deg))
                    + COS(RADIANS(orbital_ascending_node_deg)) * SIN(RADIANS(orbital_phase_deg))
                      * COS(RADIANS(orbital_inclination_deg))
                ),
                position_z_km = distance_km * SIN(RADIANS(orbital_phase_deg)) * SIN(RADIANS(orbital_inclination_deg))
            WHERE period_years > 0 AND ? >= min_update_interval_years
            """,
            (elapsed_years, elapsed_years),
        )
        counts[table] = cur.rowcount

    cur = conn.execute(
        """
        UPDATE stars
        SET galactic_orbital_phase_deg =
            MOD(galactic_orbital_phase_deg + (? / (galactic_orbital_period_gy * 1e9)) * 360, 360)
        WHERE galactic_orbital_period_gy > 0 AND ? >= galactic_min_update_interval_years
        """,
        (elapsed_years, elapsed_years),
    )
    counts["stars"] = cur.rowcount

    # v20: each star's own reflex-offset "wobble" from the planets it
    # hosts (planets.star_id) -- a correlated subquery summing every
    # hosted planet's individual pairwise pull, the same
    # utils.calculate_reflex_offset formula generation time uses (see its
    # docstring). Recomputed unconditionally on every run, using each
    # planet's just-advanced position_x/y/z_km above -- this is a cheap
    # derived value, not an independently advancing phase, so (unlike
    # every other UPDATE in this function) it has no min-update-interval
    # guard of its own. NULL/0 rows (no planets) are simply left alone by
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

    # Mutual orbit: shared by both binary configurations (see this
    # function's own docstring on why this is now a separate UPDATE from
    # the galactic-phase one below, guarded independently). v20: also
    # recomputes each star's own offset from the pair's barycenter
    # (binary_primary/secondary_position_*_km) from the freshly-advanced
    # binary_mutual_position_*_km above and the constant
    # binary_secondary_mass_fraction, using the same left-to-right
    # single-table SET evaluation trick binary_mutual_position_*_km's own
    # computation already relies on (each position expression here reads
    # binary_mutual_position_*_km's *new* value, assigned earlier in this
    # same SET list).
    cur = conn.execute(
        """
        UPDATE star_systems
        SET binary_mutual_orbital_phase_deg =
                MOD(binary_mutual_orbital_phase_deg + (? / binary_mutual_orbital_period_years) * 360, 360),
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
            modified_at = modified_at
        WHERE is_binary = 1
          AND binary_mutual_orbital_period_years > 0
          AND ? >= binary_mutual_min_update_interval_years
        """,
        (elapsed_years, elapsed_years),
    )
    counts["binary_mutual_orbits"] = cur.rowcount

    # v20: circumbinary (P-type) planets' combined pull on the whole pair
    # -- same correlated-subquery shape as the stars/planets reflex-offset
    # UPDATEs above, grouped by star_system_id instead (circumbinary
    # planets have star_id IS NULL, so there's no stars row to attach this
    # to -- see schema.sql's "v20" header note on why it's modeled as one
    # shared wobble rather than split between primary/secondary). Also
    # recomputed unconditionally, no guard interval of its own.
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

    # Galactic phase: a 'close' pair only -- a 'wide' pair's two stars
    # already each advance their own galactic phase individually via the
    # per-row `stars` UPDATE above (real stored rows, not a merged proxy).
    cur = conn.execute(
        """
        UPDATE star_systems
        SET binary_galactic_orbital_phase_deg =
                MOD(binary_galactic_orbital_phase_deg
                    + (? / (binary_galactic_orbital_period_gy * 1e9)) * 360, 360),
            modified_at = modified_at
        WHERE binary_configuration = 'close'
          AND binary_galactic_orbital_period_gy > 0
          AND ? >= binary_galactic_min_update_interval_years
        """,
        (elapsed_years, elapsed_years),
    )
    counts["binary_galactic_orbits"] = cur.rowcount

    # v17: standalone exotic phenomena -- black_holes/neutron_stars only
    # for their star_id IS NULL rows (an anchored remnant's motion already
    # advanced via the stars UPDATE above); the other five tables are
    # always standalone, so every row there qualifies.
    for table in ("black_holes", "neutron_stars"):
        cur = conn.execute(
            f"""
            UPDATE {table}
            SET galactic_orbital_phase_deg =
                MOD(galactic_orbital_phase_deg + (? / (galactic_orbital_period_gy * 1e9)) * 360, 360),
                modified_at = modified_at
            WHERE star_id IS NULL AND galactic_orbital_period_gy > 0 AND ? >= galactic_min_update_interval_years
            """,
            (elapsed_years, elapsed_years),
        )
        counts[table] = cur.rowcount

    for table in ("nebulae", "supernova_remnants", "rogue_planets", "interstellar_comets", "asteroid_fields"):
        cur = conn.execute(
            f"""
            UPDATE {table}
            SET galactic_orbital_phase_deg =
                MOD(galactic_orbital_phase_deg + (? / (galactic_orbital_period_gy * 1e9)) * 360, 360),
                modified_at = modified_at
            WHERE galactic_orbital_period_gy > 0 AND ? >= galactic_min_update_interval_years
            """,
            (elapsed_years, elapsed_years),
        )
        counts[table] = cur.rowcount

    conn.execute(
        "INSERT INTO orbit_simulation_state (id, last_updated_at) VALUES (1, NOW()) "
        "ON DUPLICATE KEY UPDATE last_updated_at = NOW()"
    )
    conn.commit()
    return counts


def advance_comet_orbits(conn, elapsed_years):
    """
    Advances every comet's own orbital anomaly (`mean_anomaly_deg` for an
    elliptical comet, `parabolic_mean_anomaly` for a parabolic one) by
    `elapsed_years`, and recomputes `distance_km`/`position_x/y/z_km`/
    `orbital_speed_kms` from the new anomaly via
    `kepler.comet_orbital_state` -- the Kepler/Barker-equation
    analog of `advance_orbital_phases`'s planet/moon handling, called
    separately by `planetgen.cli.orbits` alongside it.

    Unlike `advance_orbital_phases` (a pure, set-based SQL `UPDATE` for
    every table it touches -- `orbital_phase_deg` is a LINEAR function of
    elapsed time for a circular orbit, so MySQL/MariaDB can compute the
    resulting position directly), a comet's position is NOT a linear SQL
    expression: turning an advanced anomaly into a distance/position
    requires solving Kepler's equation (Newton-Raphson, elliptical) or
    Barker's equation (a real-cube-root closed form, parabolic) -- neither
    expressible in standard SQL. So this fetches every `comets` row and
    does that computation in Python -- one Python loop instead of one
    set-based statement, the necessary tradeoff for correctness here (in
    practice a small table -- see `tuning.SYSTEM_COMET_COUNT_RANGE`
    -- so this isn't the scaling concern it would be for `planets`/`moons`).
    The resulting rows are still written back in one batched `UPDATE` via
    `executemany` (like `insert_sector`'s per-vertex rows), not one
    `execute` per row -- the per-row work has to stay in Python, but the
    round trips to the database don't.

    An elliptical comet's `mean_anomaly_deg` advances the same
    `MOD(current + (elapsed_years / period_years) * 360, 360)` way
    `orbital_phase_deg` does, guarded by its own `min_update_interval_years`
    the identical way (see `advance_orbital_phases`'s docstring) -- skipped
    entirely, not just a no-op write, when `elapsed_years` is below it. A
    parabolic comet's `parabolic_mean_anomaly` instead advances LINEARLY
    (via `kepler.parabolic_mean_anomaly`) and does NOT wrap (a
    parabolic pass is a one-shot event, not periodic -- see
    `cometData.Comet`'s own `parabolic_mean_anomaly` docstring), and has no
    `min_update_interval_years` floor to guard against (same docstring) --
    so every parabolic row is always updated, regardless of `elapsed_years`.

    Args:
        conn (Connection): An open, schema-initialized, read-write
                           connection.
        elapsed_years (float): How much simulated time has passed since
                               the reference point `elapsed_years` was
                               computed from (the same value passed to
                               `advance_orbital_phases` -- both share one
                               `orbit_simulation_state` clock, which only
                               that function updates). Must be >= 0.

    Returns:
        int: The number of `comets` rows actually updated (a skipped,
            below-guard elliptical row doesn't count).

    Raises:
        ValueError: If `elapsed_years` is negative.
    """
    if elapsed_years < 0:
        raise ValueError(f"elapsed_years must be >= 0, got {elapsed_years}")

    rows = conn.execute(
        "SELECT id, orbit_type, perihelion_distance_km, eccentricity, inclination_deg, "
        "arg_periapsis_deg, ascending_node_deg, orbital_period_years, mean_anomaly_deg, "
        "parabolic_mean_anomaly, min_update_interval_years, primary_mass_solar FROM comets"
    ).fetchall()

    update_params = []
    for row in rows:
        perihelion_distance_au = row["perihelion_distance_km"] / physical_constants.AU_TO_KM

        if row["orbit_type"] == "elliptical":
            if elapsed_years < row["min_update_interval_years"]:
                continue
            new_mean_anomaly_deg = (
                row["mean_anomaly_deg"] + (elapsed_years / row["orbital_period_years"]) * 360
            ) % 360
            mean_anomaly_rad = math.radians(new_mean_anomaly_deg)
            new_parabolic_mean_anomaly = None
            parabolic_mean_anomaly_value = None
        else:
            new_mean_anomaly_deg = None
            mean_anomaly_rad = None
            new_parabolic_mean_anomaly = row["parabolic_mean_anomaly"] + kepler.parabolic_mean_anomaly(
                elapsed_years, perihelion_distance_au, row["primary_mass_solar"]
            )
            parabolic_mean_anomaly_value = new_parabolic_mean_anomaly

        state = kepler.comet_orbital_state(
            row["orbit_type"], perihelion_distance_au, row["eccentricity"],
            row["inclination_deg"], row["arg_periapsis_deg"], row["ascending_node_deg"],
            row["primary_mass_solar"],
            mean_anomaly_rad=mean_anomaly_rad,
            parabolic_mean_anomaly_value=parabolic_mean_anomaly_value,
            orbital_period_years=row["orbital_period_years"],
        )

        update_params.append((
            new_mean_anomaly_deg, new_parabolic_mean_anomaly,
            state["distance_au"] * physical_constants.AU_TO_KM,
            state["position_x_au"] * physical_constants.AU_TO_KM,
            state["position_y_au"] * physical_constants.AU_TO_KM,
            state["position_z_au"] * physical_constants.AU_TO_KM,
            state["orbital_speed_kms"],
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
                orbital_speed_kms = ?
            WHERE id = ?
            """,
            update_params,
        )

    conn.commit()
    return len(update_params)
