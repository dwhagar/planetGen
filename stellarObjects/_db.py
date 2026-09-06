# stellarObjects/_db.py

"""
Database Persistence (private)
===============================

Writes already-generated `StarSystem`/`SpaceSector` objects into the
SQLite database described by `stellarObjects/schema.sql` and
`db/README.md`. Leading underscore -- this module is an internal
implementation detail of the persistence boundary, not part of the
package's public generation API (`StarSystem`, `SpaceSector`, `Planet`,
etc. remain the public surface).

This module owns every unit conversion at the point of writing a value
into the database -- generation/physics code elsewhere in the package
keeps its own native units (km/AU/ly) throughout and never imports this
module. See `db/README.md` for the two-tier distance-unit convention
(milliparsecs for sector-scale placement, kilometers for everything else)
and the full column reference.

Nothing here reads the database back into live objects yet -- that's
Phase 3 territory (a query/listing tool) once there's real data to query.
"""

import os
import sqlite3

from . import physical_constants
from .asteroidData import AsteroidBelt
from .config import SystemConfig
from .doubleStar import BinaryStarProxy
from .spaceSector import SpaceSector, classify_octant
from .systemData import StarSystem
from .utils import ly_to_milliparsecs

SCHEMA_VERSION = 1
"""int: Matches `star_systems.schema_version` and `PRAGMA user_version` in
`stellarObjects/schema.sql` -- see that file's header comment."""

_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_PACKAGE_DIR)

SCHEMA_PATH = os.path.join(_PACKAGE_DIR, "schema.sql")
"""str: Path to the DDL file applied by `_ensure_schema`."""

DEFAULT_DB_PATH = os.path.join(_PROJECT_ROOT, "db", "planetgen.db")
"""str: Where the database lives by default -- the `db/` directory
scaffolded at the project root specifically for this file (gitignored)."""


def get_connection(db_path=None):
    """
    Opens a SQLite connection to `db_path` (or `DEFAULT_DB_PATH`),
    creating the containing directory and applying the schema if needed.

    Safe to call repeatedly against the same file -- `_ensure_schema` uses
    `CREATE ... IF NOT EXISTS` throughout, so an existing database is left
    untouched beyond having any missing tables/indexes/views added.

    Args:
        db_path (str, optional): Path to the `.db` file. Defaults to
                                 `DEFAULT_DB_PATH`.

    Returns:
        sqlite3.Connection: An open connection with foreign keys enabled.
    """
    path = db_path or DEFAULT_DB_PATH
    os.makedirs(os.path.dirname(path), exist_ok=True)
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA foreign_keys = ON")
    _ensure_schema(conn)
    return conn


def _ensure_schema(conn):
    """
    Applies `schema.sql` to `conn`. Idempotent -- every `CREATE TABLE`/
    `CREATE INDEX`/`CREATE VIEW` statement in the file uses
    `IF NOT EXISTS`, so this is safe to call on a database that already
    has some or all of the schema.

    Args:
        conn (sqlite3.Connection): The connection to apply the schema to.
    """
    with open(SCHEMA_PATH, "r", encoding="utf-8") as f:
        conn.executescript(f.read())


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
        conn (sqlite3.Connection): An open, schema-initialized connection.
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
        conn (sqlite3.Connection): An open, schema-initialized connection.
        star: The `Star` instance.
        star_system_id (int): The owning `star_systems.id`.
        role (str): `'primary'`, `'secondary'`, or `'single'`.

    Returns:
        int: The new `stars.id`.
    """
    props = star.get_table_properties()
    cur = conn.execute(
        """
        INSERT INTO stars (
            star_system_id, role, name, star_type, yerkes_class, mass_kg, radius_km,
            temperature_k, luminosity_w, age_gy, lifespan_gy,
            habitable_zone_inner_km, habitable_zone_outer_km,
            system_perimeter_km, heliosphere_radius_km,
            table_type, table_radius, table_mass, table_temp, table_lum, table_hab, table_loc
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            star_system_id, role, star.name, star.type, star.yerkes_class,
            star.mass, star.radius, star.temperature, star.luminosity, star.age,
            _lifespan_gy(star.lifespan),
            star.habitable_zone[0] * physical_constants.AU_TO_KM,
            star.habitable_zone[1] * physical_constants.AU_TO_KM,
            star.system_perimeter * physical_constants.AU_TO_KM,
            star.heliosphere_radius * physical_constants.AU_TO_KM,
            props["type"], props["radius"], props["mass"], props["temp"],
            props["lum"], props["hab"], props["loc"],
        ),
    )
    return cur.lastrowid


def insert_planet(conn, planet, star_system_id, star_id, orbital_index, parent_planet_id=None) -> int:
    """
    Inserts a `planets` row for one `Planet` (a top-level planet or a
    moon), then recurses over its `moons` list.

    Args:
        conn (sqlite3.Connection): An open, schema-initialized connection.
        planet (Planet): The planet or moon to persist.
        star_system_id (int): The owning `star_systems.id`.
        star_id (int or None): The specific `stars.id` this body orbits,
                               or `None` for a binary system (see the
                               `planets.star_id` column comment in
                               `schema.sql`). Threaded unchanged into every
                               recursive moon call.
        orbital_index (int): This body's position in its parent's ordered
                             list (the star's `planets` list for a
                             top-level body, or a planet's `moons` list).
        parent_planet_id (int, optional): The owning `planets.id`, for a
                                          moon. `None` for a top-level body.

    Returns:
        int: The new `planets.id`.
    """
    props = planet.get_table_properties()
    min_orbit_distance_km = (
        planet.min_orbit_distance * physical_constants.AU_TO_KM
        if planet.min_orbit_distance is not None else None
    )
    cur = conn.execute(
        """
        INSERT INTO planets (
            star_system_id, star_id, parent_planet_id, orbital_index, is_moon, body_type, name,
            planet_class, distance_km, radius_km, mass_kg, volume_km3, period_years, zone,
            description, gravity_g, surface_temperature_k, density_g_cm3, atmosphere,
            atm_density, atm_molar_density, atmospheric_pressure_pa, composition,
            scale_height_km, hill_radius_km, min_orbit_distance_km,
            habitable_zone_inner_km, habitable_zone_outer_km,
            life_chemical, evolutionary_speed, flavor_text, flavor_text_count,
            table_class, table_distance, table_period, table_radius, table_gravity
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            star_system_id, star_id, parent_planet_id, orbital_index,
            int(planet.is_moon), planet.body_type, planet.name,
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
            props.get("class"), props["distance"], props["period"],
            props["radius"], props.get("gravity"),
        ),
    )
    planet_id = cur.lastrowid

    for position, paragraph in enumerate(planet.evolutionary_data or []):
        conn.execute(
            "INSERT INTO planet_evolutionary_paragraphs (planet_id, position, paragraph) VALUES (?, ?, ?)",
            (planet_id, position, paragraph),
        )

    for spectrum_type, values in (
        ("visible", planet.reflection_spectrum_visible),
        ("non_visible", planet.reflection_spectrum_non_visible),
    ):
        for position, value in enumerate(values or []):
            conn.execute(
                """
                INSERT INTO planet_reflection_spectrum (planet_id, spectrum_type, position, value)
                VALUES (?, ?, ?, ?)
                """,
                (planet_id, spectrum_type, position, value),
            )

    for position, moon in enumerate(planet.moons or []):
        insert_planet(conn, moon, star_system_id, star_id, position, parent_planet_id=planet_id)

    return planet_id


def insert_asteroid_belt(conn, belt: AsteroidBelt, star_system_id, orbital_index) -> int:
    """
    Inserts an `asteroid_belts` row (plus its `asteroid_belt_composition`
    child rows).

    Args:
        conn (sqlite3.Connection): An open, schema-initialized connection.
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


_NULL_BINARY_FIELDS = (None,) * 18
"""Placeholder for every `binary_*` column when `star_system.star` isn't a
`BinaryStarProxy` -- see `_binary_fields`."""


def _binary_fields(proxy: BinaryStarProxy):
    """
    Extracts the 18 `star_systems.binary_*` column values from a
    `BinaryStarProxy`, in the exact order `insert_star_system`'s `INSERT`
    lists them.

    Args:
        proxy (BinaryStarProxy): The system's combined-pair proxy.

    Returns:
        tuple: 18 values, ready to splice into the `INSERT` parameters.
    """
    table_props = proxy.get_table_properties()
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
        table_props["type"], table_props["mass"], table_props["lum"],
        table_props["hab"], table_props["separation"], table_props["loc"],
    )


def insert_star_system(conn, star_system: StarSystem, system_config: SystemConfig,
                        sector_id=None, position=None) -> int:
    """
    Inserts a full `StarSystem` -- the `star_systems` row, its `stars` row(s),
    and every planet/moon/asteroid belt it contains -- into the database.

    Both `wikitext_content` and `markdown_content` are rendered here, from
    this same already-generated `star_system` object, back-to-back
    (toggling `system_config.MARKDOWN` and restoring it afterward) -- see
    `schema.sql`'s header comment for why they can't be independently
    regenerated later and still match.

    Args:
        conn (sqlite3.Connection): An open, schema-initialized connection.
        star_system (StarSystem): The generated system to persist.
        system_config (SystemConfig): The config it was generated from.
        sector_id (int, optional): The owning `sectors.id`, if this system
                                   belongs to a sector. `None` for a
                                   standalone system.
        position (tuple, optional): `(x, y, z)` in light-years, relative to
                                    the sector's center (as stored on
                                    `SectorSystemEntry.position`). `None`
                                    if the system isn't placed in a sector
                                    -- `position_x/y/z_mpc` and `quadrant`
                                    are then left `NULL`.

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
            position_x_mpc, position_y_mpc, position_z_mpc, quadrant,
            is_binary,
            binary_separation_km, binary_type, binary_temperature_k, binary_radius_km,
            binary_effective_mass_kg, binary_effective_luminosity_w, binary_age_gy, binary_lifespan_gy,
            binary_habitable_zone_inner_km, binary_habitable_zone_outer_km,
            binary_system_perimeter_km, binary_heliosphere_radius_km,
            binary_table_type, binary_table_mass, binary_table_lum, binary_table_hab,
            binary_table_separation, binary_table_loc,
            system_flavor_text, schema_version, wikitext_content, markdown_content
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sector_id, config_id, star_system.star.name,
            position_x_mpc, position_y_mpc, position_z_mpc, quadrant,
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


def insert_sector(conn, sector: SpaceSector) -> int:
    """
    Inserts a full `SpaceSector` -- the `sectors` row and every system it
    contains (with its placement) -- into the database.

    Args:
        conn (sqlite3.Connection): An open, schema-initialized connection.
        sector (SpaceSector): The sector to persist.

    Returns:
        int: The new `sectors.id`.
    """
    cur = conn.execute(
        "INSERT INTO sectors (name, edge_mpc) VALUES (?, ?)",
        (sector.name, ly_to_milliparsecs(sector.edge_ly)),
    )
    sector_id = cur.lastrowid

    for entry in sector.entries:
        insert_star_system(
            conn, entry.star_system, entry.system_config,
            sector_id=sector_id, position=entry.position,
        )

    return sector_id


def save_sector(sector: SpaceSector, db_path=None) -> int:
    """
    Opens (or creates) the database and persists a full `SpaceSector` to
    it in one transaction.

    Args:
        sector (SpaceSector): The sector to persist.
        db_path (str, optional): Path to the `.db` file. Defaults to
                                 `DEFAULT_DB_PATH`.

    Returns:
        int: The new `sectors.id`.
    """
    conn = get_connection(db_path)
    try:
        with conn:
            sector_id = insert_sector(conn, sector)
        return sector_id
    finally:
        conn.close()
