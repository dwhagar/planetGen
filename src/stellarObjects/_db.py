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

import os
from collections import namedtuple

import pymysql
import pymysql.cursors
from dbutils.pooled_db import PooledDB

from . import physical_constants
from .asteroidData import AsteroidBelt
from .config import SystemConfig
from .doubleStar import BinaryStarProxy
from .galaxyDensity import GalaxyShape
from .planetData import Planet
from .spaceSector import SectorSystemEntry, SpaceSector, classify_octant, distance_between
from .starData import Star
from .systemData import StarSystem
from .utils import ly_to_milliparsecs, milliparsecs_to_ly

SCHEMA_VERSION = 9
"""int: Matches `star_systems.schema_version` and the highest row in the
`schema_migrations` table (see `stellarObjects/schema.sql`'s header
comment). Also the target version `migrate_database` brings a database's
`schema_migrations` bookkeeping up to."""

_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

SCHEMA_PATH = os.path.join(_PACKAGE_DIR, "schema.sql")
"""str: Path to the DDL file applied by `_ensure_schema`."""


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
    """

    def __init__(self, host=None, port=None, user=None, password=None, database=None):
        self.host = host if host is not None else os.environ.get("PLANETGEN_MYSQL_HOST", "127.0.0.1")
        self.port = int(port if port is not None else os.environ.get("PLANETGEN_MYSQL_PORT", 3306))
        self.user = user if user is not None else os.environ.get("PLANETGEN_MYSQL_USER", "planetgen")
        self.password = password if password is not None else os.environ.get("PLANETGEN_MYSQL_PASSWORD", "")
        self.database = database if database is not None else os.environ.get("PLANETGEN_MYSQL_DATABASE", "planetgen")

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
    non-standard JSON `Infinity` token.

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
    cur = conn.execute(
        """
        INSERT INTO stars (
            star_system_id, role, name, star_type, yerkes_class, mass_kg, radius_km,
            temperature_k, luminosity_w, age_gy, lifespan_gy,
            habitable_zone_inner_km, habitable_zone_outer_km,
            system_perimeter_km, heliosphere_radius_km
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            star_system_id, role, star.name, star.type, star.yerkes_class,
            star.mass, star.radius, star.temperature, star.luminosity, star.age,
            _lifespan_gy(star.lifespan),
            star.habitable_zone[0] * physical_constants.AU_TO_KM,
            star.habitable_zone[1] * physical_constants.AU_TO_KM,
            star.system_perimeter * physical_constants.AU_TO_KM,
            star.heliosphere_radius * physical_constants.AU_TO_KM,
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
            rotation_period_hours
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
            planet.orbital_phase_deg, planet.rotation_period_hours,
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
            rotation_period_hours
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
            moon.orbital_phase_deg, moon.rotation_period_hours,
        ),
    )
    moon_id = cur.lastrowid

    _insert_paragraphs(conn, "moon_evolutionary_paragraphs", "moon_id", moon_id, moon.evolutionary_data)
    _insert_reflection_spectrum(
        conn, "moon_reflection_spectrum", "moon_id", moon_id,
        moon.reflection_spectrum_visible, moon.reflection_spectrum_non_visible,
    )

    return moon_id


def insert_asteroid_belt(conn, belt: AsteroidBelt, star_system_id, orbital_index) -> int:
    """
    Inserts an `asteroid_belts` row (plus its `asteroid_belt_composition`
    child rows).

    Args:
        conn (Connection): An open, schema-initialized connection.
        belt (AsteroidBelt): The belt to persist.
        star_system_id (int): The owning `star_systems.id`.
        orbital_index (int): This belt's position in the star's `planets`
                             list (shared index space with `Planet`
                             entries, so orbital order across both types is
                             preserved).

    Returns:
        int: The new `asteroid_belts.id`.
    """
    cur = conn.execute(
        """
        INSERT INTO asteroid_belts (
            star_system_id, orbital_index, distance_km, lower_limit_km, upper_limit_km,
            density, composition_summary
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        (
            star_system_id, orbital_index,
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


_NULL_BINARY_FIELDS = (None,) * 12
"""Placeholder for every `binary_*` column when `star_system.star` isn't a
`BinaryStarProxy` -- see `_binary_fields`."""


def _binary_fields(proxy: BinaryStarProxy):
    """
    Extracts the 12 `star_systems.binary_*` column values from a
    `BinaryStarProxy`, in the exact order `insert_star_system`'s `INSERT`
    lists them.

    Args:
        proxy (BinaryStarProxy): The system's combined-pair proxy.

    Returns:
        tuple: 12 values, ready to splice into the `INSERT` parameters.
    """
    return (
        proxy.binary_separation_au * physical_constants.AU_TO_KM,
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
    and every planet/moon/asteroid belt it contains -- into the database.

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

    is_binary = isinstance(star_system.star, BinaryStarProxy)
    binary_fields = _binary_fields(star_system.star) if is_binary else _NULL_BINARY_FIELDS

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
            is_binary,
            binary_separation_km, binary_type, binary_temperature_k, binary_radius_km,
            binary_effective_mass_kg, binary_effective_luminosity_w, binary_age_gy, binary_lifespan_gy,
            binary_habitable_zone_inner_km, binary_habitable_zone_outer_km,
            binary_system_perimeter_km, binary_heliosphere_radius_km,
            system_flavor_text, schema_version, wikitext_content, markdown_content
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, config_id, star_system.star.name,
            position_x_mpc, position_y_mpc, position_z_mpc, quadrant, location,
            int(is_binary),
            *binary_fields,
            star_system.system_flavor_text, SCHEMA_VERSION, wikitext_content, markdown_content,
        ),
    )
    star_system_id = cur.lastrowid

    if is_binary:
        insert_star(conn, star_system.primary_star, star_system_id, "primary")
        insert_star(conn, star_system.secondary_star, star_system_id, "secondary")
        planet_star_id = None  # planets orbit the proxy, not a stored star row -- see planets.star_id
    else:
        planet_star_id = insert_star(conn, star_system.star, star_system_id, "single")

    for orbital_index, obj in enumerate(star_system.planets):
        if obj.body_type == "a":
            insert_asteroid_belt(conn, obj, star_system_id, orbital_index)
        else:
            insert_planet(conn, obj, star_system_id, planet_star_id, orbital_index)

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
    }


def _binary_proxy_row_to_dict(star_system_row, primary_dict, secondary_dict):
    """Maps a `star_systems` row's `binary_*` columns to
    `BinaryStarProxy.from_dict`'s expected dict shape, inverting every unit
    conversion `_binary_fields` applies."""
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
        "_binary_separation_au": row["binary_separation_km"] / physical_constants.AU_TO_KM,
        "_effective_mass": row["binary_effective_mass_kg"],
        "_effective_luminosity": row["binary_effective_luminosity_w"],
        "primary": primary_dict,
        "secondary": secondary_dict,
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


def load_star_system(conn, star_system_id) -> StarSystem:
    """
    Reconstructs a full `StarSystem` -- config, star(s), and every planet/
    moon/asteroid belt it contains, in original orbital order -- from a
    `star_systems` row and its related rows.

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

    if row["is_binary"]:
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
    else:
        single_row = conn.execute(
            "SELECT * FROM stars WHERE star_system_id = ? AND role = 'single'", (star_system_id,)
        ).fetchone()
        star = Star.from_dict(_star_row_to_dict(single_row), system_config)

    system = object.__new__(StarSystem)
    system.system_config = system_config
    system.star = star
    if isinstance(star, BinaryStarProxy):
        system.primary_star = star._primary
        system.secondary_star = star._secondary
        system.stars = [star._primary, star._secondary]
    else:
        system.primary_star = star
        system.stars = [star]

    planet_rows = conn.execute(
        "SELECT * FROM planets WHERE star_system_id = ? ORDER BY orbital_index", (star_system_id,)
    ).fetchall()
    belt_rows = conn.execute(
        "SELECT * FROM asteroid_belts WHERE star_system_id = ? ORDER BY orbital_index", (star_system_id,)
    ).fetchall()

    # planets/belts share one orbital_index space (see insert_star_system's
    # enumerate over star_system.planets) -- merge and re-sort by it to
    # restore that original interleaved order.
    combined = [("p", r) for r in planet_rows] + [("b", r) for r in belt_rows]
    combined.sort(key=lambda item: item[1]["orbital_index"])

    system.planets = []
    for kind, r in combined:
        if kind == "b":
            comp_rows = conn.execute(
                "SELECT component, concentration FROM asteroid_belt_composition "
                "WHERE belt_id = ? ORDER BY position",
                (r["id"],),
            ).fetchall()
            system.planets.append(AsteroidBelt.from_dict(_belt_row_to_dict(r, comp_rows), system_config))
        else:
            planet_data = _planet_or_moon_row_to_dict(conn, r, is_moon=False)
            moon_rows = conn.execute(
                "SELECT * FROM moons WHERE planet_id = ? ORDER BY orbital_index", (r["id"],)
            ).fetchall()
            planet_data["moons"] = [_planet_or_moon_row_to_dict(conn, mr, is_moon=True) for mr in moon_rows]
            system.planets.append(Planet.from_dict(planet_data, star, system_config))

    system.system_flavor_text = row["system_flavor_text"]
    system.planet_count, system.belt_count, system.moon_count = system.count_objects()
    system.hab_count, system.m_count = system.count_habitable()

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


def migrate_database(config=None):
    """
    Brings a database's `schema_migrations` bookkeeping up to
    `SCHEMA_VERSION`, applying any migration step in between.

    `get_connection`/`_ensure_schema` always create a brand-new database
    already at `SCHEMA_VERSION` (every `CREATE TABLE IF NOT EXISTS` in
    `schema.sql` reflects the current shape directly), so this function
    only has real work to do against a database created by an older
    version of this project -- `_migrate_v8_to_v9` (added for the v9
    orbital-motion columns) is the first such step; see `schema.sql`'s
    header comment for the versioning convention, and `migrateDb.py` for
    the CLI wrapper around this.

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

    `orbital_inclination_deg`/`orbital_ascending_node_deg`/
    `rotation_period_hours` are untouched -- fixed at generation time, per
    `planetPhysics.generate_orbital_motion_properties`.

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
        tuple: (planets_updated, moons_updated) -- row counts, straight
              from each `UPDATE`'s own affected-row count.

    Raises:
        ValueError: If `elapsed_years` is negative.
    """
    if elapsed_years < 0:
        raise ValueError(f"elapsed_years must be >= 0, got {elapsed_years}")

    counts = []
    for table in ("planets", "moons"):
        cur = conn.execute(
            f"""
            UPDATE {table}
            SET orbital_phase_deg = MOD(orbital_phase_deg + (? / period_years) * 360, 360)
            WHERE period_years > 0
            """,
            (elapsed_years,),
        )
        counts.append(cur.rowcount)

    conn.execute(
        "INSERT INTO orbit_simulation_state (id, last_updated_at) VALUES (1, NOW()) "
        "ON DUPLICATE KEY UPDATE last_updated_at = NOW()"
    )
    conn.commit()
    return tuple(counts)
