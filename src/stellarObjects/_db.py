# stellarObjects/_db.py

"""
Database Persistence (private)
===============================

Writes already-generated `StarSystem`/`SpaceSector` objects into the
MySQL database described by `stellarObjects/schema.sql` and
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
allowlists/reconstruction logic from `stellarObjects/serialization.py`
(Phase 1) without a redundant object-to-dict-to-object round trip for data
that was never nested to begin with.

MySQL port (TODO.md Phase 5): this module used to talk directly to
`sqlite3`. It now goes through `pymysql` (pure-Python driver, no system
libraries to build against on the deployment host) via a small
`Connection` wrapper (below) that keeps every existing call site's
`conn.execute(sql, params)` shape working unchanged -- including this
module's and `queryDb.py`'s `?` positional placeholders, which the
wrapper rewrites to `pymysql`'s `%s` at the point of execution, and
`sqlite3.Row`-style `row["column"]` access, which `pymysql`'s
`DictCursor` already provides natively. Real concurrent access uses a
connection pool (`DBUtils.PooledDB`) rather than opening a fresh TCP
connection per call, per TODO.md's "add real connection pooling" note.
"""

import math
import os
from collections import namedtuple

import pymysql
import pymysql.cursors
from dbutils.pooled_db import PooledDB

from . import keplerMotion, physical_constants
from .appconfig import load_config
from .asteroidData import AsteroidBelt
from .asteroidFieldData import AsteroidField
from .compactRemnant import BlackHole, NeutronStar
from .cometData import Comet
from .config import SystemConfig
from .doubleStar import BinaryStarProxy
from .galaxyDensity import GalaxyShape
from .nebulaData import Nebula
from .planetData import Planet
from .roguePlanetData import InterstellarComet, RoguePlanet
from .spaceSector import SectorSystemEntry, SpaceSector, classify_octant, distance_between
from .starData import Star
from .supernovaRemnantData import SupernovaRemnant
from .systemData import StarSystem
from .utils import ly_to_milliparsecs, milliparsecs_to_ly
from .wideBinary import WideBinaryPair

SCHEMA_VERSION = 18
"""int: Matches `star_systems.schema_version` and the highest row in the
`schema_migrations` table (see `stellarObjects/schema.sql`'s header
comment). Also the target version `migrate_database` brings a database's
`schema_migrations` bookkeeping up to."""

_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

SCHEMA_PATH = os.path.join(_PACKAGE_DIR, "schema.sql")
"""str: Path to the DDL file applied by `_ensure_schema`."""

CONTROL_SCHEMA_VERSION = 1
"""int: Version counter for `control_schema.sql`, independent of
`SCHEMA_VERSION` above -- see that file's header comment for why the
control plane (admin identities/sessions/API keys/audit log) is a
separate schema with its own versioning."""

CONTROL_SCHEMA_PATH = os.path.join(_PACKAGE_DIR, "control_schema.sql")
"""str: Path to the DDL file applied by `_ensure_control_schema`."""

CONTROL_DB_ENV_VAR = "PLANETGEN_CONTROL_DATABASE"
"""str: Env var naming the one MySQL schema the control plane lives in
(admin identities are global to a deployment, not per-galaxy -- see
`control_schema.sql`'s header comment). Falls back to `config.json`'s
`control_database` (see `stellarObjects.appconfig`), then
`DEFAULT_CONTROL_DATABASE`, when unset."""

DEFAULT_CONTROL_DATABASE = "planetgen_control"
"""str: Default control-schema name when neither `CONTROL_DB_ENV_VAR` nor
`config.json`'s `control_database` is set."""


class MySQLConfig:
    """
    MySQL connection parameters, read from environment variables --
    mirrors every other entry point in this project (`sectorGen.py`,
    `systemGen.py`, `queryDb.py`, `html/api/config.py`) reading its own
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
    `mysql` section (see `stellarObjects.appconfig`), then the
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
    `systemGen.py`, `galaxyGen.py`, `queryDb.py`, `migrateDb.py`) instead
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


def _get_pool(config):
    key = config._key()
    if key not in _pools:
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
        )
    return _pools[key]


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


class Connection:
    """
    Thin wrapper around a pooled `pymysql` connection that keeps this
    module's (and `queryDb.py`'s) existing `conn.execute(sql, params)` /
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
    """

    def __init__(self, pooled_conn):
        self._conn = pooled_conn

    def execute(self, sql, params=()):
        cur = self._conn.cursor()
        cur.execute(sql.replace("?", "%s"), params)
        return _Cursor(cur)

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
        cur = self._conn.cursor()
        cur.executemany(sql.replace("?", "%s"), list(seq_of_params))
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

    def commit(self):
        self._conn.commit()

    def rollback(self):
        self._conn.rollback()

    def close(self):
        self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is None:
            self.commit()
        else:
            self.rollback()
        return False


def get_connection(config=None, ensure_schema=True):
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
        ensure_schema (bool): Whether to run `_ensure_schema` (DDL) on
            this connection. `True` (the default) suits every read-write
            caller (`save_sector`/`save_system`/the generation CLIs) --
            the whole point of `CREATE ... IF NOT EXISTS` is that a fresh
            database gets its schema the first time anything connects.
            Every read-only caller (`queryDb.py`'s `open_readonly`, the
            Flask API's `get_db`, `html/lib/dbutil.py`) should instead
            pass `False`: this project's own docs recommend pointing
            those tools at a database account with `SELECT`-only grants
            (see `queryDb.py`'s module docstring), and DDL requires
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
    conn = Connection(_get_pool(config).connection())
    if ensure_schema:
        _ensure_schema(conn)
    return conn


DB_PREFIX_ENV_VAR = "PLANETGEN_MYSQL_DATABASE_PREFIX"
"""str: Env var overriding the default schema-name prefix `list_databases`/
`resolve_database` filter by -- see `MySQLConfig`'s own `database` default.
Shared by every entry point that offers a choice among several MySQL
schemas on one server (the `html/` CGI browser's `?db=` picker, and the
Flask API's own `?db=`/`/api/databases`, both via this one implementation).
Falls back to `config.json`'s `mysql.database_prefix` (see
`stellarObjects.appconfig`), then to `DEFAULT_DB_PREFIX`, when unset."""

DEFAULT_DB_PREFIX = "planetgen"
"""str: Matches `MySQLConfig`'s own default database name -- a deployment
with just one schema names it `planetgen` and never needs to set
`DB_PREFIX_ENV_VAR` (or `config.json`'s `mysql.database_prefix`) at all;
one with several names them `planetgen_<something>` to share the prefix."""


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


def _ensure_schema(conn):
    """
    Applies `schema.sql` to `conn`, then bootstraps `schema_migrations`
    (inserting `SCHEMA_VERSION` as the baseline row) if it's empty --
    idempotent, safe to call on a database that already has some, all, or
    none of the schema.

    Args:
        conn (Connection): The connection to apply the schema to.
    """
    with open(SCHEMA_PATH, "r", encoding="utf-8") as f:
        conn.executescript(f.read())

    row = conn.execute("SELECT COUNT(*) AS n FROM schema_migrations").fetchone()
    if row["n"] == 0:
        conn.execute("INSERT INTO schema_migrations (version) VALUES (?)", (SCHEMA_VERSION,))
        conn.commit()


def open_write(config=None):
    """
    Opens a connection for a write-capable caller against an already-
    existing content database (`ensure_schema=False` -- same reasoning as
    `queryDb.open_readonly`: a write-capable account (see
    `docs/apache-deployment.md`'s `PLANETGEN_MYSQL_WRITE_*`) deliberately
    has no `CREATE`/`ALTER` grant, so attempting `_ensure_schema`'s DDL
    here would fail every connection instead of just skipping a step a
    full-access account has already done once, via `migrateDb.py`).

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
    host/port/user/password -- typically the same write-capable account
    `open_write` uses (`docs/apache-deployment.md`'s `PLANETGEN_MYSQL_WRITE_*`),
    since the control schema needs the same `SELECT`/`INSERT`/`UPDATE`/
    `DELETE` grants, just on a different schema name.

    Args:
        base_config (MySQLConfig, optional): Connection parameters to
            reuse host/port/user/password from. Defaults to
            `DEFAULT_MYSQL_CONFIG`.

    Returns:
        MySQLConfig: `base_config` with `database` replaced by
            `CONTROL_DB_ENV_VAR` (or `DEFAULT_CONTROL_DATABASE`).
    """
    base_config = base_config or DEFAULT_MYSQL_CONFIG
    database = os.environ.get(CONTROL_DB_ENV_VAR) or load_config()["control_database"] or DEFAULT_CONTROL_DATABASE
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
            runtime case: the Flask API's write-capable account has no
            `CREATE` grant, same reasoning as `open_write` above).
            `adminAuth.bootstrap_control_schema` passes `True`, using a
            full-access account (the same one `migrateDb.py` already
            uses), to create/update the control schema once per deploy.

    Returns:
        Connection: An open connection.
    """
    config = config or control_mysql_config()
    conn = Connection(_get_pool(config).connection())
    if ensure_schema:
        _ensure_control_schema(conn)
    return conn


def _ensure_control_schema(conn):
    """
    Applies `control_schema.sql` to `conn`, then bootstraps
    `control_schema_migrations` (inserting `CONTROL_SCHEMA_VERSION` as the
    baseline row) if it's empty -- idempotent, mirrors `_ensure_schema`
    above exactly, just against the control schema's own DDL file/version
    counter.

    Args:
        conn (Connection): The connection to apply the schema to.
    """
    with open(CONTROL_SCHEMA_PATH, "r", encoding="utf-8") as f:
        conn.executescript(f.read())

    row = conn.execute("SELECT COUNT(*) AS n FROM control_schema_migrations").fetchone()
    if row["n"] == 0:
        conn.execute("INSERT INTO control_schema_migrations (version) VALUES (?)", (CONTROL_SCHEMA_VERSION,))
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
            markdown, habitable_world, asteroid_belt, large_star, moons,
            max_planets, planets, star_type, name, age, intelligent_life,
            binary_system, num_orbits
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            int(bool(config.MARKDOWN)),
            _tristate(config.HABITABLE_WORLD),
            _tristate(config.ASTEROID_BELT),
            _tristate(config.LARGE_STAR),
            _tristate(config.MOONS),
            _tristate(config.MAX_PLANETS),
            _tristate(config.PLANETS),
            config.STAR_TYPE,
            config.NAME,
            config.AGE,
            _tristate(config.INTELLIGENT_LIFE),
            _tristate(config.BINARY_SYSTEM),
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
            wide_binary_a_crit_km
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
    min_orbit_distance_km = (
        planet.min_orbit_distance * physical_constants.AU_TO_KM
        if planet.min_orbit_distance is not None else None
    )
    cur = conn.execute(
        """
        INSERT INTO planets (
            star_system_id, star_id, orbital_index, body_type, name,
            planet_class, distance_km, radius_km, mass_kg, volume_km3, period_years, zone,
            description, gravity_g, surface_temperature_k, density_g_cm3, atmosphere,
            atm_density, atm_molar_density, atmospheric_pressure_pa, composition,
            scale_height_km, hill_radius_km, min_orbit_distance_km,
            habitable_zone_inner_km, habitable_zone_outer_km,
            life_chemical, evolutionary_speed, flavor_text, flavor_text_count,
            orbital_inclination_deg, orbital_ascending_node_deg, orbital_phase_deg,
            position_x_km, position_y_km, position_z_km, orbital_speed_kms,
            min_update_interval_years,
            rotation_period_hours
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            star_system_id, star_id, orbital_index, planet.body_type, planet.name,
            planet.planet_class,
            planet.distance * physical_constants.AU_TO_KM,
            planet.radius, planet.mass, planet.volume, planet.period, planet.zone,
            planet.description, planet.gravity, planet.surface_temperature,
            planet.density, planet.atmosphere,
            planet.atm_density, planet.atm_molar_density, planet.atmospheric_pressure,
            planet.composition,
            planet.scale_height, planet.hill_radius, min_orbit_distance_km,
            planet.habitable_zone[0] * physical_constants.AU_TO_KM,
            planet.habitable_zone[1] * physical_constants.AU_TO_KM,
            planet.life_chemical, planet.evolutionary_speed,
            planet.flavor_text, planet.flavor_text_count,
            planet.orbital_inclination_deg, planet.orbital_ascending_node_deg,
            planet.orbital_phase_deg,
            planet.position_x * physical_constants.AU_TO_KM,
            planet.position_y * physical_constants.AU_TO_KM,
            planet.position_z * physical_constants.AU_TO_KM,
            planet.orbital_speed_kms,
            planet.min_update_interval_years,
            planet.rotation_period_hours,
        ),
    )
    planet_id = cur.lastrowid

    _insert_paragraphs(conn, "planet_evolutionary_paragraphs", "planet_id", planet_id, planet.evolutionary_data)
    _insert_reflection_spectrum(
        conn, "planet_reflection_spectrum", "planet_id", planet_id,
        planet.reflection_spectrum_visible, planet.reflection_spectrum_non_visible,
    )

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
    min_orbit_distance_km = (
        moon.min_orbit_distance * physical_constants.AU_TO_KM
        if moon.min_orbit_distance is not None else None
    )
    cur = conn.execute(
        """
        INSERT INTO moons (
            planet_id, star_system_id, star_id, orbital_index, body_type, name,
            planet_class, distance_km, radius_km, mass_kg, volume_km3, period_years, zone,
            description, gravity_g, surface_temperature_k, density_g_cm3, atmosphere,
            atm_density, atm_molar_density, atmospheric_pressure_pa, composition,
            scale_height_km, hill_radius_km, min_orbit_distance_km,
            habitable_zone_inner_km, habitable_zone_outer_km,
            life_chemical, evolutionary_speed, flavor_text, flavor_text_count,
            orbital_inclination_deg, orbital_ascending_node_deg, orbital_phase_deg,
            position_x_km, position_y_km, position_z_km, orbital_speed_kms,
            min_update_interval_years,
            rotation_period_hours
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            planet_id, star_system_id, star_id, orbital_index, moon.body_type, moon.name,
            moon.planet_class,
            moon.distance * physical_constants.AU_TO_KM,
            moon.radius, moon.mass, moon.volume, moon.period, moon.zone,
            moon.description, moon.gravity, moon.surface_temperature,
            moon.density, moon.atmosphere,
            moon.atm_density, moon.atm_molar_density, moon.atmospheric_pressure,
            moon.composition,
            moon.scale_height, moon.hill_radius, min_orbit_distance_km,
            moon.habitable_zone[0] * physical_constants.AU_TO_KM,
            moon.habitable_zone[1] * physical_constants.AU_TO_KM,
            moon.life_chemical, moon.evolutionary_speed,
            moon.flavor_text, moon.flavor_text_count,
            moon.orbital_inclination_deg, moon.orbital_ascending_node_deg,
            moon.orbital_phase_deg,
            moon.position_x * physical_constants.AU_TO_KM,
            moon.position_y * physical_constants.AU_TO_KM,
            moon.position_z * physical_constants.AU_TO_KM,
            moon.orbital_speed_kms,
            moon.min_update_interval_years,
            moon.rotation_period_hours,
        ),
    )
    moon_id = cur.lastrowid

    _insert_paragraphs(conn, "moon_evolutionary_paragraphs", "moon_id", moon_id, moon.evolutionary_data)
    _insert_reflection_spectrum(
        conn, "moon_reflection_spectrum", "moon_id", moon_id,
        moon.reflection_spectrum_visible, moon.reflection_spectrum_non_visible,
    )

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


def insert_black_hole(conn, black_hole: BlackHole, star_id=None) -> int:
    """
    Inserts a `black_holes` row (see `schema.sql`'s "v16"/"v17" header
    notes).

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
            truth for the same object's position.

    Returns:
        int: The new `black_holes.id`.
    """
    galactic_fields = (
        (None, None, None, None) if star_id is not None else (
            black_hole.galactic_orbital_speed_kms, black_hole.galactic_orbital_period_gy,
            black_hole.galactic_orbital_phase_deg, black_hole.galactic_min_update_interval_years,
        )
    )
    cur = conn.execute(
        """
        INSERT INTO black_holes (
            star_id, name, mass_solar, event_horizon_radius_km, spin,
            has_accretion_disk, temperature_k, luminosity_w, age_gy,
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            star_id, black_hole.name, black_hole.mass_solar,
            black_hole.event_horizon_radius_km, black_hole.spin,
            int(black_hole.has_accretion_disk), black_hole.temperature,
            black_hole.luminosity, black_hole.age,
            *galactic_fields,
        ),
    )
    return cur.lastrowid


def insert_neutron_star(conn, neutron_star: NeutronStar, star_id=None) -> int:
    """
    Inserts a `neutron_stars` row (see `schema.sql`'s "v16"/"v17" header
    notes).

    Args:
        conn (Connection): An open, schema-initialized connection.
        neutron_star (NeutronStar): The neutron star to persist.
        star_id (int, optional): The owning `stars.id`, when this neutron
            star anchors a `StarSystem` (`phenomenonGen.py --anchor-system`
            -- see `insert_star_system`'s single-star branch, the only
            caller that passes this). `None` for one generated standalone.
            When set, the `galactic_orbital_*` columns are left `NULL` --
            see `insert_black_hole`'s identical reasoning.

    Returns:
        int: The new `neutron_stars.id`.
    """
    galactic_fields = (
        (None, None, None, None) if star_id is not None else (
            neutron_star.galactic_orbital_speed_kms, neutron_star.galactic_orbital_period_gy,
            neutron_star.galactic_orbital_phase_deg, neutron_star.galactic_min_update_interval_years,
        )
    )
    cur = conn.execute(
        """
        INSERT INTO neutron_stars (
            star_id, name, mass_solar, radius_km, spin_period_ms, magnetic_field_gauss,
            pulsar_type, surface_temperature_k, luminosity_w, age_gy,
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            star_id, neutron_star.name, neutron_star.mass_solar, neutron_star.radius,
            neutron_star.spin_period_ms, neutron_star.magnetic_field_gauss,
            neutron_star.pulsar_type, neutron_star.surface_temperature_k,
            neutron_star.luminosity, neutron_star.age,
            *galactic_fields,
        ),
    )
    return cur.lastrowid


def insert_nebula(conn, nebula: Nebula, sector_id=None) -> int:
    """
    Inserts a `nebulae` row (see `schema.sql`'s "v16" header note).

    Args:
        conn (Connection): An open, schema-initialized connection.
        nebula (Nebula): The nebula to persist.
        sector_id (int, optional): Reserved for a future sector-context
            encounter. `None` (always, today -- `phenomenonGen.py` never
            creates or attaches a sector).

    Returns:
        int: The new `nebulae.id`.
    """
    cur = conn.execute(
        """
        INSERT INTO nebulae (
            sector_id, name, nebula_type, radius_ly, composition, formation_cause,
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, nebula.name, nebula.nebula_type, nebula.radius_ly, nebula.composition, nebula.formation_cause,
            nebula.galactic_orbital_speed_kms, nebula.galactic_orbital_period_gy,
            nebula.galactic_orbital_phase_deg, nebula.galactic_min_update_interval_years,
        ),
    )
    return cur.lastrowid


def insert_supernova_remnant(conn, remnant: SupernovaRemnant, sector_id=None) -> int:
    """
    Inserts a `supernova_remnants` row (see `schema.sql`'s "v16" header
    note), plus (for a core-collapse progenitor whose collapsed core is
    still detectable) its embedded `black_holes`/`neutron_stars` row.

    Args:
        conn (Connection): An open, schema-initialized connection.
        remnant (SupernovaRemnant): The remnant to persist.
        sector_id (int, optional): Reserved for a future sector-context
            encounter. `None` (always, today).

    Returns:
        int: The new `supernova_remnants.id`.
    """
    compact_remnant_kind = None
    black_hole_id = None
    neutron_star_id = None
    if isinstance(remnant.compact_remnant, BlackHole):
        compact_remnant_kind = "black_hole"
        black_hole_id = insert_black_hole(conn, remnant.compact_remnant)
    elif isinstance(remnant.compact_remnant, NeutronStar):
        compact_remnant_kind = "neutron_star"
        neutron_star_id = insert_neutron_star(conn, remnant.compact_remnant)

    cur = conn.execute(
        """
        INSERT INTO supernova_remnants (
            sector_id, name, morphology, age_years, radius_ly, progenitor_type,
            compact_remnant_kind, compact_remnant_black_hole_id, compact_remnant_neutron_star_id,
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, remnant.name, remnant.morphology, remnant.age_years, remnant.radius_ly,
            remnant.progenitor_type, compact_remnant_kind, black_hole_id, neutron_star_id,
            remnant.galactic_orbital_speed_kms, remnant.galactic_orbital_period_gy,
            remnant.galactic_orbital_phase_deg, remnant.galactic_min_update_interval_years,
        ),
    )
    return cur.lastrowid


def insert_rogue_planet(conn, planet: RoguePlanet, sector_id=None) -> int:
    """
    Inserts a `rogue_planets` row (see `schema.sql`'s "v16" header note).

    Args:
        conn (Connection): An open, schema-initialized connection.
        planet (RoguePlanet): The rogue planet to persist.
        sector_id (int, optional): Reserved for a future sector-context
            encounter. `None` (always, today).

    Returns:
        int: The new `rogue_planets.id`.
    """
    cur = conn.execute(
        """
        INSERT INTO rogue_planets (
            sector_id, name, planet_type, mass_kg, radius_km, composition, has_internal_heat, has_moons,
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, planet.name, planet.planet_type, planet.mass_kg, planet.radius_km,
            planet.composition, int(planet.has_internal_heat), int(planet.has_moons),
            planet.galactic_orbital_speed_kms, planet.galactic_orbital_period_gy,
            planet.galactic_orbital_phase_deg, planet.galactic_min_update_interval_years,
        ),
    )
    return cur.lastrowid


def insert_interstellar_comet(conn, comet: InterstellarComet, sector_id=None) -> int:
    """
    Inserts an `interstellar_comets` row (plus its
    `interstellar_comet_composition` child rows; see `schema.sql`'s "v16"
    header note).

    Args:
        conn (Connection): An open, schema-initialized connection.
        comet (InterstellarComet): The comet to persist.
        sector_id (int, optional): Reserved for a future sector-context
            encounter. `None` (always, today).

    Returns:
        int: The new `interstellar_comets.id`.
    """
    cur = conn.execute(
        """
        INSERT INTO interstellar_comets (
            sector_id, name, nucleus_diameter_km, velocity_kms, is_active, composition_summary,
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, comet.name, comet.nucleus_diameter_km, comet.velocity_kms,
            int(comet.is_active), comet.get_composition_summary(),
            comet.galactic_orbital_speed_kms, comet.galactic_orbital_period_gy,
            comet.galactic_orbital_phase_deg, comet.galactic_min_update_interval_years,
        ),
    )
    comet_id = cur.lastrowid

    for position, component in enumerate(comet.composition):
        conn.execute(
            "INSERT INTO interstellar_comet_composition (comet_id, position, component) VALUES (?, ?, ?)",
            (comet_id, position, component),
        )

    return comet_id


def insert_asteroid_field(conn, field: AsteroidField, sector_id=None) -> int:
    """
    Inserts an `asteroid_fields` row (plus its `asteroid_field_composition`
    child rows; see `schema.sql`'s "v17" header note) -- the standalone
    counterpart to `insert_asteroid_belt`, following the identical
    belt-plus-child-rows shape.

    Args:
        conn (Connection): An open, schema-initialized connection.
        field (AsteroidField): The asteroid field to persist.
        sector_id (int, optional): Reserved for a future sector-context
            encounter. `None` (always, today).

    Returns:
        int: The new `asteroid_fields.id`.
    """
    cur = conn.execute(
        """
        INSERT INTO asteroid_fields (
            sector_id, name, density, radius_ly, composition_summary,
            galactic_orbital_speed_kms, galactic_orbital_period_gy,
            galactic_orbital_phase_deg, galactic_min_update_interval_years
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, field.name, field.density, field.radius_ly, field.get_composition_summary(),
            field.galactic_orbital_speed_kms, field.galactic_orbital_period_gy,
            field.galactic_orbital_phase_deg, field.galactic_min_update_interval_years,
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


def save_phenomenon(phenomenon, system_config: SystemConfig, phenomenon_type: str, config=None) -> int:
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
            `program_constants.PHENOMENON_TYPE_CHOICES`, naming which
            table `phenomenon` belongs in (ignored when `phenomenon` is a
            `StarSystem`, which always goes through `insert_star_system`
            regardless of whether it's anchored by a black hole or a
            neutron star).
        config (MySQLConfig, optional): Connection parameters. Defaults to
            `DEFAULT_MYSQL_CONFIG`.

    Returns:
        int: The new row's id, in whichever table `phenomenon`/
            `phenomenon_type` maps to.

    Raises:
        ValueError: If `phenomenon_type` isn't one of the recognized
                   choices (and `phenomenon` isn't a `StarSystem`).
    """
    conn = get_connection(config)
    try:
        with conn:
            if isinstance(phenomenon, StarSystem):
                # An anchored compact remnant: insert_star_system's own
                # single-star branch already inserts the satellite
                # black_holes/neutron_stars row -- see that function.
                return insert_star_system(conn, phenomenon, system_config)
            if phenomenon_type == "black-hole":
                return insert_black_hole(conn, phenomenon)
            if phenomenon_type == "neutron-star":
                return insert_neutron_star(conn, phenomenon)
            if phenomenon_type == "nebula":
                return insert_nebula(conn, phenomenon)
            if phenomenon_type == "supernova-remnant":
                return insert_supernova_remnant(conn, phenomenon)
            if phenomenon_type == "rogue-planet":
                return insert_rogue_planet(conn, phenomenon)
            if phenomenon_type == "comet":
                return insert_interstellar_comet(conn, phenomenon)
            if phenomenon_type == "asteroid-field":
                return insert_asteroid_field(conn, phenomenon)
            raise ValueError(f"Unknown phenomenon type: {phenomenon_type!r}")
    finally:
        conn.close()


_NULL_PROXY_ONLY_BINARY_FIELDS = (None,) * 15
"""Placeholder for the 15 `star_systems.binary_*` columns that only ever
describe a merged `BinaryStarProxy` (a 'close'/P-type pair) -- always NULL
for a 'wide'/S-type pair or a single star, since no merged effective star
exists to describe in either of those cases. See `schema.sql`'s "v15"
header note and `_proxy_only_binary_fields`."""

_NULL_MUTUAL_ORBIT_FIELDS = (None,) * 9
"""Placeholder for the 9 `star_systems.binary_mutual_*` columns (excluding
`binary_separation_km`, handled separately in `insert_star_system` since it
sits earlier in column order, alongside the eccentricity/periapsis/apoapsis
columns it's grouped with) -- NULL for a single (non-binary) star. See
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
    Extracts the 9 `star_systems.binary_mutual_*` column values (period,
    speed, inclination, ascending node, phase, update-guard interval, and
    x/y/z position -- everything except `binary_separation_km`, handled
    separately) from a 'close' pair's `BinaryStarProxy`.

    Returns:
        tuple: 9 values, ready to splice into the `INSERT` parameters.
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
    )


def _mutual_orbit_fields_from_wide_binary(pair):
    """
    The same 9 `star_systems.binary_mutual_*` columns as
    `_mutual_orbit_fields_from_proxy`, from a 'wide' pair's
    `doubleStar.WideBinaryPair` instead -- these columns are reused
    unchanged across both binary configurations (see `schema.sql`'s "v15"
    header note): a wide pair's own (circular-approximation) mutual orbit
    fits the exact same shape a close pair's already occupies.

    Args:
        pair (WideBinaryPair): The system's wide-binary orbital pair.

    Returns:
        tuple: 9 values, ready to splice into the `INSERT` parameters.
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
        str: e.g. `"Voranthis Kelmoor -- nearest: Alpha Prime (4.2 ly), ..."`,
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
        (neighbor.star_system.star.name, distance_between(entry, neighbor))
        for neighbor in neighbors
    ]
    return _format_location_string(sector.name, neighbor_info)


def insert_star_system(conn, star_system: StarSystem, system_config: SystemConfig,
                        sector_id=None, position=None, location=None) -> int:
    """
    Inserts a full `StarSystem` -- the `star_systems` row, its `stars` row(s),
    and every planet/moon/asteroid belt/comet it contains -- into the database.

    Both `wikitext_content` and `markdown_content` are rendered here, from
    this same already-generated `star_system` object, back-to-back
    (toggling `system_config.MARKDOWN` and restoring it afterward) -- see
    `schema.sql`'s header comment for why they can't be independently
    regenerated later and still match.

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

    if position is not None:
        position_x_mpc = ly_to_milliparsecs(position[0])
        position_y_mpc = ly_to_milliparsecs(position[1])
        position_z_mpc = ly_to_milliparsecs(position[2])
        quadrant, _magnitudes = classify_octant(position)
    else:
        position_x_mpc = position_y_mpc = position_z_mpc = None
        quadrant = None
        location = None

    # Render both formats from this same generated object -- rendering is
    # idempotent (Phase 0), so toggling MARKDOWN here has no other effect
    # on the object and doesn't re-roll anything.
    original_markdown = system_config.MARKDOWN
    system_config.MARKDOWN = False
    wikitext_content = str(star_system)
    system_config.MARKDOWN = True
    markdown_content = str(star_system)
    system_config.MARKDOWN = original_markdown

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
            system_flavor_text, schema_version, wikitext_content, markdown_content
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, config_id, star_system.star.name,
            position_x_mpc, position_y_mpc, position_z_mpc, quadrant, location,
            int(is_binary), binary_configuration,
            separation_km, eccentricity, periapsis_km, apoapsis_km,
            *proxy_only_fields,
            *mutual_orbit_fields,
            star_system.system_flavor_text, SCHEMA_VERSION, wikitext_content, markdown_content,
        ),
    )
    star_system_id = cur.lastrowid

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


def insert_sector(conn, sector: SpaceSector, galaxy_position=None) -> int:
    """
    Inserts a full `SpaceSector` -- the `sectors` row and every system it
    contains (with its placement) -- into the database.

    Args:
        conn (Connection): An open, schema-initialized connection.
        sector (SpaceSector): The sector to persist.
        galaxy_position (dict, optional): This sector's galaxy-frame
            placement (see `schema.sql`'s "v4"/"v6/v7" header notes), or
            `None` (the default) for a sector never placed in a galaxy --
            `sectorGen.py`'s own standalone CLI keeps producing these.
            When given, must have keys `center_x_pc`, `center_y_pc`,
            `center_z_pc`, `galactic_radius_pc`, `vertices_pc` (all
            required together -- the schema's CHECK constraint enforces
            the first four on every other path, but this function trusts
            the caller rather than re-deriving `galactic_radius_pc`
            itself; `vertices_pc` isn't part of that CHECK since it lives
            in a separate table SQLite can't cross-reference in a CHECK,
            but is expected NULL-together with the other four all the
            same), and optionally `shell_index`/`shell_slot_index` (each
            independently optional -- `None`/omitted leaves that one
            column NULL, per `schema.sql`'s note that they aren't implied
            by a center point the way the other five are). `vertices_pc`
            is a `{"inner": [...], "outer": [...]}` dict, each a list of
            `(x, y, z)` tuples -- variable length, not fixed at 8 (see
            `sectorGeometry.prism_vertices`) -- written as one
            `sector_vertices` row per vertex, ordinary columns throughout,
            never serialized.

    Returns:
        int: The new `sectors.id`.
    """
    if galaxy_position is not None:
        cur = conn.execute(
            """
            INSERT INTO sectors (
                name, edge_mpc, center_x_pc, center_y_pc, center_z_pc,
                galactic_radius_pc, shell_index, shell_slot_index
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                sector.name, ly_to_milliparsecs(sector.edge_ly),
                galaxy_position["center_x_pc"], galaxy_position["center_y_pc"],
                galaxy_position["center_z_pc"], galaxy_position["galactic_radius_pc"],
                galaxy_position.get("shell_index"), galaxy_position.get("shell_slot_index"),
            ),
        )
        sector_id = cur.lastrowid
        for ring in ("inner", "outer"):
            conn.executemany(
                "INSERT INTO sector_vertices (sector_id, ring, vertex_index, x_pc, y_pc, z_pc) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (
                    (sector_id, ring, i, x, y, z)
                    for i, (x, y, z) in enumerate(galaxy_position["vertices_pc"][ring])
                ),
            )
    else:
        cur = conn.execute(
            "INSERT INTO sectors (name, edge_mpc) VALUES (?, ?)",
            (sector.name, ly_to_milliparsecs(sector.edge_ly)),
        )
        sector_id = cur.lastrowid

    for entry in sector.entries:
        insert_star_system(
            conn, entry.star_system, entry.system_config,
            sector_id=sector_id, position=entry.position,
            location=_location_for_entry(sector, entry),
        )

    return sector_id


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
            `center_z_pc`, `galactic_radius_pc`, `shell_index`,
            `shell_slot_index`, `vertices_pc` (rebuilt from `sector_vertices`
            rows into a `{"inner": [...], "outer": [...]}` dict of
            `[x, y, z]` lists, ordered by `vertex_index`), or `None` if
            this sector has never been placed in a galaxy (the four
            galaxy-position columns NULL).

    Raises:
        ValueError: If no such `sectors` row exists.
    """
    row = conn.execute(
        """
        SELECT center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc,
               shell_index, shell_slot_index
        FROM sectors WHERE id = ?
        """,
        (sector_id,),
    ).fetchone()
    if row is None:
        raise ValueError(f"no sectors row with id {sector_id}")
    if row["center_x_pc"] is None:
        return None

    result = dict(row)
    vertex_rows = conn.execute(
        """
        SELECT ring, x_pc, y_pc, z_pc FROM sector_vertices
        WHERE sector_id = ? ORDER BY ring, vertex_index
        """,
        (sector_id,),
    ).fetchall()
    vertices_pc = {"inner": [], "outer": []}
    for vrow in vertex_rows:
        vertices_pc[vrow["ring"]].append([vrow["x_pc"], vrow["y_pc"], vrow["z_pc"]])
    result["vertices_pc"] = vertices_pc
    return result


def get_occupied_shell_slots(conn, shell_indices):
    """
    Returns every already-occupied `(shell_index, shell_slot_index)`
    address among the given shell indices -- used by `galaxyGen.py` (both
    batch and local-neighborhood modes) to skip addresses a sector already
    exists at, in one query per batch of candidate shells rather than one
    query per candidate slot.

    Args:
        conn (Connection): An open, schema-initialized connection.
        shell_indices (iterable): Shell indices to check.

    Returns:
        set: `(shell_index, shell_slot_index)` tuples already present in
            `sectors`. Empty if `shell_indices` is empty.
    """
    shell_indices = list(shell_indices)
    if not shell_indices:
        return set()

    placeholders = ", ".join("?" for _ in shell_indices)
    rows = conn.execute(
        f"SELECT shell_index, shell_slot_index FROM sectors "
        f"WHERE shell_index IN ({placeholders}) AND shell_slot_index IS NOT NULL",
        tuple(shell_indices),
    ).fetchall()
    return {(row["shell_index"], row["shell_slot_index"]) for row in rows}


def get_sector_id_at(conn, shell_index, shell_slot_index):
    """
    Looks up the `sectors.id` already generated at a specific galaxy
    address, if any -- the single-address counterpart to
    `get_occupied_shell_slots`'s batch-of-a-shell query, used by
    `galaxyGen.ensure_sector_generated` to check (and, on an `INSERT`
    race, re-check) one address at a time.

    Args:
        conn (Connection): An open, schema-initialized connection.
        shell_index (int): The shell index to look up.
        shell_slot_index (int): The slot index within that shell.

    Returns:
        int or None: The existing `sectors.id`, or `None` if no sector has
            been generated at this address yet.
    """
    row = conn.execute(
        "SELECT id FROM sectors WHERE shell_index = ? AND shell_slot_index = ?",
        (shell_index, shell_slot_index),
    ).fetchone()
    return row["id"] if row is not None else None


GalaxySkeletonInfo = namedtuple(
    "GalaxySkeletonInfo",
    ["shape", "edge_pc", "outer_shell_index", "expected_system_count_at_density_1"],
)
"""The galaxy's stored skeleton -- everything needed to recompute any
sector's exact position/density on demand (see schema.sql's "v8" header
note). `shape` is a `galaxyDensity.GalaxyShape`; the other three fields
are `galaxy_shape`'s own remaining columns."""


def save_galaxy_shape(shape: GalaxyShape, edge_pc, outer_shell_index,
                       expected_system_count_at_density_1, config=None):
    """
    Replaces the galaxy's singleton `galaxy_shape` row -- there is exactly
    one galaxy, so this always overwrites whatever was there before rather
    than inserting a second row (`galaxyPlan.py` calls this once per full
    skeleton (re)build).

    Args:
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.
        edge_pc (float): The sector edge length this skeleton was built
                         at, parsecs.
        outer_shell_index (int): The last shell index with any qualifying
            content (`galaxyPlan.py`'s own discovered galaxy edge).
        expected_system_count_at_density_1 (float): See
            `galaxySkeleton.expected_system_count_at_density_1`.
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.
    """
    conn = get_connection(config)
    try:
        with conn:
            conn.execute(
                """
                INSERT INTO galaxy_shape (
                    id, disk_scale_length_pc, disk_scale_height_pc,
                    bulge_scale_radius_pc, bulge_amplitude, arm_count,
                    pitch_angle_rad, arm_amplitude, spiral_reference_radius_pc,
                    spiral_reference_angle_rad, k_norm, edge_pc,
                    expected_system_count_at_density_1, outer_shell_index
                ) VALUES (1, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
                    outer_shell_index = VALUES(outer_shell_index)
                """,
                (
                    shape.disk_scale_length_pc, shape.disk_scale_height_pc,
                    shape.bulge_scale_radius_pc, shape.bulge_amplitude, shape.arm_count,
                    shape.pitch_angle_rad, shape.arm_amplitude, shape.spiral_reference_radius_pc,
                    shape.spiral_reference_angle_rad, shape.k_norm, edge_pc,
                    expected_system_count_at_density_1, outer_shell_index,
                ),
            )
    finally:
        conn.close()


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
    row = conn.execute(
        """
        SELECT disk_scale_length_pc, disk_scale_height_pc, bulge_scale_radius_pc,
               bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude,
               spiral_reference_radius_pc, spiral_reference_angle_rad, k_norm,
               edge_pc, expected_system_count_at_density_1, outer_shell_index
        FROM galaxy_shape WHERE id = 1
        """
    ).fetchone()
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
        shape=shape, edge_pc=row["edge_pc"], outer_shell_index=row["outer_shell_index"],
        expected_system_count_at_density_1=row["expected_system_count_at_density_1"],
    )


def replace_galaxy_shell_bands(shell_bands, config=None):
    """
    Replaces every `galaxy_shell_band` row wholesale -- `galaxyPlan.py`'s
    own full-galaxy skeleton build is the only writer, and it always
    produces a complete, coherent set for the whole galaxy in one pass, so
    there is no notion of an incremental/partial update here (matches
    `save_galaxy_shape`'s own "replace the one true answer" behavior).

    Args:
        shell_bands (iterable): `(shell_index, band_index, slot_index_min,
            slot_index_max)` tuples, any order -- `band_index` is the
            0-based position of that band within its own shell (almost
            always just `0`; see `galaxySkeleton.find_shell_bands`).
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.
    """
    conn = get_connection(config)
    try:
        with conn:
            conn.execute("DELETE FROM galaxy_shell_band")
            conn.executemany(
                "INSERT INTO galaxy_shell_band (shell_index, band_index, slot_index_min, slot_index_max) "
                "VALUES (?, ?, ?, ?)",
                shell_bands,
            )
    finally:
        conn.close()


def get_galaxy_shell_bands(conn, shell_index):
    """
    This shell's stored candidate band(s), in `band_index` order.

    Args:
        conn (Connection): An open, schema-initialized connection.
        shell_index (int): The shell index to look up.

    Returns:
        list: `(slot_index_min, slot_index_max)` tuples, in ascending
              `band_index` order -- empty if this shell has no stored
              qualifying content (including if the skeleton was never
              built at all).
    """
    rows = conn.execute(
        "SELECT slot_index_min, slot_index_max FROM galaxy_shell_band "
        "WHERE shell_index = ? ORDER BY band_index",
        (shell_index,),
    ).fetchall()
    return [(row["slot_index_min"], row["slot_index_max"]) for row in rows]


def save_system(star_system: StarSystem, system_config: SystemConfig, config=None) -> int:
    """
    Opens the database and persists a single, standalone `StarSystem` (no
    sector -- `sector_id`/`position` are left `None`) in one transaction.
    The single-system counterpart to `save_sector`, for `systemGen.py`
    (which, unlike `sectorGen.py`, generates one system with no natural
    sector placement of its own).

    Args:
        star_system (StarSystem): The generated system to persist.
        system_config (SystemConfig): The config it was generated from.
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.

    Returns:
        int: The new `star_systems.id`.
    """
    conn = get_connection(config)
    try:
        with conn:
            star_system_id = insert_star_system(conn, star_system, system_config)
        return star_system_id
    finally:
        conn.close()


def save_sector(sector: SpaceSector, config=None, galaxy_position=None) -> int:
    """
    Opens the database and persists a full `SpaceSector` to it in one
    transaction.

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
    conn = get_connection(config)
    try:
        with conn:
            sector_id = insert_sector(conn, sector, galaxy_position=galaxy_position)
        return sector_id
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
        "large_star": _tristate_from_db(row["large_star"]),
        "moons": _tristate_from_db(row["moons"]),
        "max_planets": _tristate_from_db(row["max_planets"]),
        "planets": _tristate_from_db(row["planets"]),
        "star_type": row["star_type"],
        "name": row["name"],
        "age": row["age"],
        "intelligent_life": _tristate_from_db(row["intelligent_life"]),
        "binary_system": _tristate_from_db(row["binary_system"]),
        "num_orbits": row["num_orbits"],
        "slots": slots,
    })


def _star_row_to_dict(row):
    """Maps a `stars` row to `Star.from_dict`'s expected dict shape,
    inverting every unit conversion `insert_star` applies.

    `temperature` is cast back to `int` -- `Star.generate_star` always sets
    it via `int(round(...))`, but SQLite's `REAL` column type hands every
    numeric value back as a Python `float` regardless of what was stored,
    and `f"{self.temperature} K"` (`get_table_properties`) renders `5800`
    vs. `5800.0` differently -- the one place a plain round-trip through a
    `REAL` column would otherwise silently break render fidelity.
    """
    return {
        "name": row["name"],
        "type": row["star_type"],
        "yerkes_class": row["yerkes_class"],
        "mass": row["mass_kg"],
        "radius": row["radius_km"],
        "temperature": int(row["temperature_k"]),
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
    if black_hole_row is not None:
        return BlackHole.from_dict(_black_hole_row_to_dict(star_row, black_hole_row), system_config)

    neutron_star_row = conn.execute("SELECT * FROM neutron_stars WHERE star_id = ?", (star_row["id"],)).fetchone()
    if neutron_star_row is not None:
        return NeutronStar.from_dict(_neutron_star_row_to_dict(star_row, neutron_star_row), system_config)

    return Star.from_dict(_star_row_to_dict(star_row), system_config)


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

    system = object.__new__(StarSystem)
    system.system_config = system_config
    system.star = star
    system.binary_type = binary_configuration
    system.wide_binary = wide_binary

    if binary_configuration == "close":
        system.primary_star = star._primary
        system.secondary_star = star._secondary
        system.stars = [star._primary, star._secondary]
    elif binary_configuration == "wide":
        system.primary_star = star
        system.secondary_star = secondary_star
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
            comets.append(Comet.from_dict(_comet_row_to_dict(r, comp_rows), system_config))
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
                objects.append(AsteroidBelt.from_dict(_belt_row_to_dict(r, comp_rows), system_config))
            else:
                planet_data = _planet_or_moon_row_to_dict(conn, r, is_moon=False)
                moon_rows = conn.execute(
                    "SELECT * FROM moons WHERE planet_id = ? ORDER BY orbital_index", (r["id"],)
                ).fetchall()
                planet_data["moons"] = [_planet_or_moon_row_to_dict(conn, mr, is_moon=True) for mr in moon_rows]
                objects.append(Planet.from_dict(planet_data, owning_star, system_config))
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

    sector = SpaceSector(row["name"], edge_ly=milliparsecs_to_ly(row["edge_mpc"]))

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
    `updateOrbits.py` run against them without special-casing "does this
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
    Records schema v18 -- new tables `comets` and `comet_composition`
    (`cometData.Comet` -- see `schema.sql`'s "v18" header note) for
    star-bound comets.

    Like `_migrate_v15_to_v16`, this needs no `ALTER TABLE`: both are
    brand-new tables, and `_ensure_schema`'s `CREATE TABLE IF NOT EXISTS`
    (run on every new connection) already creates them directly from the
    current `schema.sql`, even against a database whose
    `schema_migrations` bookkeeping still says v17 -- there is no
    pre-existing table whose shape needs changing. This step exists
    purely to keep that bookkeeping counter itself accurate.

    Args:
        conn (Connection): An open connection, mid-migration (not yet
                           committed -- the caller commits once every step
                           up to `SCHEMA_VERSION` has run).
    """
    conn.execute("INSERT INTO schema_migrations (version) VALUES (18)")


def migrate_database(config=None):
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
    on those tables plus the new `asteroid_fields` phenomenon), and
    `_migrate_v17_to_v18` (added for v18's new `comets`/`comet_composition`
    tables) are the migration steps so far; see `schema.sql`'s header
    comment for the versioning convention, and `migrateDb.py` for the CLI
    wrapper around this.

    Args:
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.

    Returns:
        int: The database's `schema_migrations` version (always
            `SCHEMA_VERSION` after this call).
    """
    conn = get_connection(config)
    try:
        row = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()
        version = row["version"]

        if version < 9:
            _migrate_v8_to_v9(conn)
            version = 9

        if version < 10:
            _migrate_v9_to_v10(conn)
            version = 10

        if version < 11:
            _migrate_v10_to_v11(conn)
            version = 11

        if version < 12:
            _migrate_v11_to_v12(conn)
            version = 12

        if version < 13:
            _migrate_v12_to_v13(conn)
            version = 13

        if version < 14:
            _migrate_v13_to_v14(conn)
            version = 14

        if version < 15:
            _migrate_v14_to_v15(conn)
            version = 15

        if version < 16:
            _migrate_v15_to_v16(conn)
            version = 16

        if version < 17:
            _migrate_v16_to_v17(conn)
            version = 17

        if version < 18:
            _migrate_v17_to_v18(conn)
            version = 18

        conn.commit()
        return version
    finally:
        conn.close()


def get_orbit_update_elapsed_years(conn):
    """
    Returns how many years have elapsed since `updateOrbits.py` last
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


def advance_orbital_phases(conn, elapsed_years):
    """
    Advances every planet's and moon's `orbital_phase_deg` in place by the
    fraction of a full revolution `elapsed_years` represents, given each
    body's own already-stored `period_years` -- one set-based `UPDATE` per
    table rather than a per-row Python loop, so this stays fast regardless
    of how many bodies the database holds (see `updateOrbits.py`).
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
    (`updateOrbits.py` runs "once a month or so"), so the guard exists for
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

    Also upserts `orbit_simulation_state.last_updated_at` to `NOW()` (the
    reference point the *next* call's `elapsed_years` should be measured
    from), in the same transaction, so a caller can never advance phases
    without also recording that it did.

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
            `"binary_mutual_orbits"`, `"binary_galactic_orbits"`,
            `"black_holes"`, `"neutron_stars"`, `"nebulae"`,
            `"supernova_remnants"`, `"rogue_planets"`,
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

    # Mutual orbit: shared by both binary configurations (see this
    # function's own docstring on why this is now a separate UPDATE from
    # the galactic-phase one below, guarded independently).
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
                * SIN(RADIANS(binary_mutual_orbital_phase_deg)) * SIN(RADIANS(binary_mutual_orbital_inclination_deg))
        WHERE is_binary = 1
          AND binary_mutual_orbital_period_years > 0
          AND ? >= binary_mutual_min_update_interval_years
        """,
        (elapsed_years, elapsed_years),
    )
    counts["binary_mutual_orbits"] = cur.rowcount

    # Galactic phase: a 'close' pair only -- a 'wide' pair's two stars
    # already each advance their own galactic phase individually via the
    # per-row `stars` UPDATE above (real stored rows, not a merged proxy).
    cur = conn.execute(
        """
        UPDATE star_systems
        SET binary_galactic_orbital_phase_deg =
                MOD(binary_galactic_orbital_phase_deg
                    + (? / (binary_galactic_orbital_period_gy * 1e9)) * 360, 360)
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
                MOD(galactic_orbital_phase_deg + (? / (galactic_orbital_period_gy * 1e9)) * 360, 360)
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
                MOD(galactic_orbital_phase_deg + (? / (galactic_orbital_period_gy * 1e9)) * 360, 360)
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
    `keplerMotion.comet_orbital_state` -- the Kepler/Barker-equation
    analog of `advance_orbital_phases`'s planet/moon handling, called
    separately by `updateOrbits.py` alongside it.

    Unlike `advance_orbital_phases` (a pure, set-based SQL `UPDATE` for
    every table it touches -- `orbital_phase_deg` is a LINEAR function of
    elapsed time for a circular orbit, so MySQL/MariaDB can compute the
    resulting position directly), a comet's position is NOT a linear SQL
    expression: turning an advanced anomaly into a distance/position
    requires solving Kepler's equation (Newton-Raphson, elliptical) or
    Barker's equation (a real-cube-root closed form, parabolic) -- neither
    expressible in standard SQL. So this fetches every `comets` row, does
    that computation in Python, and writes each result back with its own
    `UPDATE` -- one Python loop instead of one set-based statement, the
    necessary tradeoff for correctness here (in practice a small table --
    see `program_constants.SYSTEM_COMET_COUNT_RANGE` -- so this isn't the
    scaling concern it would be for `planets`/`moons`).

    An elliptical comet's `mean_anomaly_deg` advances the same
    `MOD(current + (elapsed_years / period_years) * 360, 360)` way
    `orbital_phase_deg` does, guarded by its own `min_update_interval_years`
    the identical way (see `advance_orbital_phases`'s docstring) -- skipped
    entirely, not just a no-op write, when `elapsed_years` is below it. A
    parabolic comet's `parabolic_mean_anomaly` instead advances LINEARLY
    (via `keplerMotion.parabolic_mean_anomaly`) and does NOT wrap (a
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

    updated = 0
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
            new_parabolic_mean_anomaly = row["parabolic_mean_anomaly"] + keplerMotion.parabolic_mean_anomaly(
                elapsed_years, perihelion_distance_au, row["primary_mass_solar"]
            )
            parabolic_mean_anomaly_value = new_parabolic_mean_anomaly

        state = keplerMotion.comet_orbital_state(
            row["orbit_type"], perihelion_distance_au, row["eccentricity"],
            row["inclination_deg"], row["arg_periapsis_deg"], row["ascending_node_deg"],
            row["primary_mass_solar"],
            mean_anomaly_rad=mean_anomaly_rad,
            parabolic_mean_anomaly_value=parabolic_mean_anomaly_value,
            orbital_period_years=row["orbital_period_years"],
        )

        conn.execute(
            """
            UPDATE comets
            SET mean_anomaly_deg = ?, parabolic_mean_anomaly = ?,
                distance_km = ?, position_x_km = ?, position_y_km = ?, position_z_km = ?,
                orbital_speed_kms = ?
            WHERE id = ?
            """,
            (
                new_mean_anomaly_deg, new_parabolic_mean_anomaly,
                state["distance_au"] * physical_constants.AU_TO_KM,
                state["position_x_au"] * physical_constants.AU_TO_KM,
                state["position_y_au"] * physical_constants.AU_TO_KM,
                state["position_z_au"] * physical_constants.AU_TO_KM,
                state["orbital_speed_kms"],
                row["id"],
            ),
        )
        updated += 1

    conn.commit()
    return updated
