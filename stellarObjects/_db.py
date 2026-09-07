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
"""

import os
import shutil
import sqlite3
import time

from . import physical_constants
from .asteroidData import AsteroidBelt
from .config import SystemConfig
from .doubleStar import BinaryStarProxy
from .planetData import Planet
from .spaceSector import SectorSystemEntry, SpaceSector, classify_octant, distance_between
from .starData import Star
from .systemData import StarSystem
from .utils import ly_to_milliparsecs, milliparsecs_to_ly

SCHEMA_VERSION = 3
"""int: Matches `star_systems.schema_version` and `PRAGMA user_version` in
`stellarObjects/schema.sql` -- see that file's header comment. Also the
target version `migrate_database` converts an older database up to."""

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
        sqlite3.Connection: An open connection with foreign keys enabled
                            and `row_factory` set to `sqlite3.Row` (supports
                            both `row[0]` and `row["column"]` access, so
                            this is backward compatible with existing
                            positional-index callers).
    """
    path = db_path or DEFAULT_DB_PATH
    os.makedirs(os.path.dirname(path), exist_ok=True)
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.row_factory = sqlite3.Row
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
        conn (sqlite3.Connection): An open, schema-initialized connection.
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
    props = planet.get_table_properties()
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
            table_class, table_distance, table_period, table_radius, table_gravity
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
            props.get("class"), props["distance"], props["period"],
            props["radius"], props.get("gravity"),
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
        conn (sqlite3.Connection): An open, schema-initialized connection.
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
    props = moon.get_table_properties()
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
            table_class, table_distance, table_period, table_radius, table_gravity
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
            props.get("class"), props["distance"], props["period"],
            props["radius"], props.get("gravity"),
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


def _format_location_string(sector_name, neighbors):
    """
    Builds the `star_systems.location` display string: the owning sector's
    name, followed by up to 3 nearest neighbors and their distances in
    light-years, nearest first -- see `schema.sql`'s "v3" header note.

    Shared by both the live-object write path (`_location_for_entry`, used
    from `insert_sector`) and `migrate_database`'s v3 backfill
    (`_backfill_v3_locations`), so the two never drift into different
    formats.

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
            binary_table_type, binary_table_mass, binary_table_lum, binary_table_hab,
            binary_table_separation, binary_table_loc,
            system_flavor_text, schema_version, wikitext_content, markdown_content
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
            location=_location_for_entry(sector, entry),
        )

    return sector_id


def save_system(star_system: StarSystem, system_config: SystemConfig, db_path=None) -> int:
    """
    Opens (or creates) the database and persists a single, standalone
    `StarSystem` (no sector -- `sector_id`/`position` are left `None`) in
    one transaction. The single-system counterpart to `save_sector`, for
    `systemGen.py` (which, unlike `sectorGen.py`, generates one system with
    no natural sector placement of its own).

    Args:
        star_system (StarSystem): The generated system to persist.
        system_config (SystemConfig): The config it was generated from.
        db_path (str, optional): Path to the `.db` file. Defaults to
                                 `DEFAULT_DB_PATH`.

    Returns:
        int: The new `star_systems.id`.
    """
    conn = get_connection(db_path)
    try:
        with conn:
            star_system_id = insert_star_system(conn, star_system, system_config)
        return star_system_id
    finally:
        conn.close()


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
        conn (sqlite3.Connection): An open, schema-initialized connection.
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
        conn (sqlite3.Connection): An open, schema-initialized connection.
        row (sqlite3.Row): The `planets` or `moons` row.
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
        conn (sqlite3.Connection): An open, schema-initialized connection.
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
        conn (sqlite3.Connection): An open, schema-initialized connection.
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


class UnsupportedSchemaVersionError(Exception):
    """Raised by `migrate_database` for a `PRAGMA user_version` this
    module has neither a migration path from nor can use as-is."""


# The v1 `planets` table's columns, minus `is_moon`/`parent_planet_id`
# (the discriminator `_migrate_v1_to_v2` splits on) -- identical to both
# v2 `planets` and v2 `moons`' own column sets (see schema.sql), which is
# what makes a single shared column list usable for copying into either
# one below. `id` is included deliberately: a migrated moon keeps its
# original `planets.id` value as its new `moons.id`, so
# `planet_evolutionary_paragraphs`/`planet_reflection_spectrum` rows can
# move to their moon-owned counterparts by that same, unchanged id --
# no id remapping table needed.
_PLANET_MOON_COLUMNS = (
    "id", "star_system_id", "star_id", "orbital_index", "body_type", "name",
    "planet_class", "distance_km", "radius_km", "mass_kg", "volume_km3",
    "period_years", "zone", "description", "gravity_g", "surface_temperature_k",
    "density_g_cm3", "atmosphere", "atm_density", "atm_molar_density",
    "atmospheric_pressure_pa", "composition", "scale_height_km", "hill_radius_km",
    "min_orbit_distance_km", "habitable_zone_inner_km", "habitable_zone_outer_km",
    "life_chemical", "evolutionary_speed", "flavor_text", "flavor_text_count",
    "table_class", "table_distance", "table_period", "table_radius", "table_gravity",
)

# The full `star_systems` column set as of schema v1/v2 -- i.e. every
# column except v3's new `location` (see schema.sql's "v3" header note).
# v1 -> v2 made no structural change to `star_systems` at all (only
# `planets`/`moons` changed -- see the "v2" header note), so this one list
# covers both source versions' actual column layout, and is used by both
# `_migrate_v1_to_v2` and `_migrate_v2_to_v3` to copy `star_systems`
# explicitly column-by-column instead of `INSERT ... SELECT *` -- required
# now that the current (v3) schema has one more column than either v1 or
# v2 did. `location` is left unset by this copy and backfilled afterward
# by `_backfill_v3_locations`.
_STAR_SYSTEMS_PRE_V3_COLUMNS = (
    "id", "sector_id", "system_config_id", "name",
    "position_x_mpc", "position_y_mpc", "position_z_mpc", "quadrant",
    "is_binary",
    "binary_separation_km", "binary_type", "binary_temperature_k", "binary_radius_km",
    "binary_effective_mass_kg", "binary_effective_luminosity_w", "binary_age_gy", "binary_lifespan_gy",
    "binary_habitable_zone_inner_km", "binary_habitable_zone_outer_km",
    "binary_system_perimeter_km", "binary_heliosphere_radius_km",
    "binary_table_type", "binary_table_mass", "binary_table_lum", "binary_table_hab",
    "binary_table_separation", "binary_table_loc",
    "system_flavor_text", "schema_version", "wikitext_content", "markdown_content",
    "mediawiki_url", "wikijs_url", "created_at",
)


def _backfill_v3_locations(conn):
    """
    Populates `star_systems.location` (see `schema.sql`'s "v3" header note)
    for every already-migrated row that has a sector position, after
    `_migrate_v1_to_v2`/`_migrate_v2_to_v3` copy `star_systems` across
    without it -- a pre-v3 database predates the column, so it was never
    computed for these rows the way `insert_sector`/`_location_for_entry`
    compute it for a freshly generated sector.

    Works directly against the rows already written to `main` (there's no
    live `SpaceSector`/`StarSystem` object graph to call
    `SpaceSector.nearest_neighbors` on during a migration) -- same
    nearest-3, nearest-first logic, using plain SQL-row Euclidean distance
    on the stored milliparsec positions (converted to light-years for
    display via `milliparsecs_to_ly`, since that conversion is linear and
    applies equally to a distance magnitude as to a coordinate) and
    `_format_location_string` for the shared text format.

    Args:
        conn (sqlite3.Connection): Connection to the new database, with
                                   the v3 schema and `sectors`/`star_systems`
                                   rows already populated (`location` still
                                   `NULL` on every row).
    """
    sectors = conn.execute("SELECT id, name FROM main.sectors").fetchall()
    for sector_id, sector_name in sectors:
        systems = conn.execute(
            """
            SELECT id, name, position_x_mpc, position_y_mpc, position_z_mpc
            FROM main.star_systems
            WHERE sector_id = ? AND position_x_mpc IS NOT NULL
            """,
            (sector_id,),
        ).fetchall()

        for system_id, system_name, x, y, z in systems:
            others = [s for s in systems if s[0] != system_id]
            others.sort(key=lambda s: distance_between((x, y, z), (s[2], s[3], s[4])))
            neighbor_info = [
                (other[1], milliparsecs_to_ly(distance_between((x, y, z), (other[2], other[3], other[4]))))
                for other in others[:3]
            ]
            location = _format_location_string(sector_name, neighbor_info)
            conn.execute("UPDATE main.star_systems SET location = ? WHERE id = ?", (location, system_id))


# Tables schema v1 -> v2 leaves structurally untouched -- copied verbatim,
# column-for-column, from the attached old database. Split into "before"/
# "after" `star_systems` groups (rather than one flat list) purely for
# insert ORDER: `stars`/`asteroid_belts` carry a `star_system_id` FK, so
# `star_systems` itself must already exist in `main` before they're copied
# -- and `star_systems` needs the explicit `_STAR_SYSTEMS_PRE_V3_COLUMNS`
# column list instead of `SELECT *`, since the current schema's `location`
# column means `star_systems` can't just be one more entry in this list.
_V1_VERBATIM_TABLES_BEFORE_STAR_SYSTEMS = ("sectors", "system_configs", "system_config_slots")
_V1_VERBATIM_TABLES_AFTER_STAR_SYSTEMS = ("stars", "asteroid_belts", "asteroid_belt_composition")


def _migrate_v1_to_v2(conn):
    """
    Populates a fresh, empty database (`conn`'s main schema, always the
    *current* `schema.sql` -- there's no separate "as of v2" snapshot of
    it) from a schema-v1 database attached as `old` (see `migrate_database`
    for the attach/backup/swap this runs inside of, and for why this one
    function now migrates a v1 database all the way to the current schema
    in a single hop rather than stopping at an intermediate v2 shape).

    The structural changes this bridges are moons moving out of the shared
    `planets` table (`is_moon`/`parent_planet_id`) into their own `moons`
    table (`planet_id`) -- see `schema.sql`'s "v2" header note -- and the
    v3 `star_systems.location` column not existing yet, backfilled
    afterward by `_backfill_v3_locations`. Every other table is unaffected
    and copied straight across.

    Args:
        conn (sqlite3.Connection): Connection to the new database, with
                                   the old one already `ATTACH`ed as `old`
                                   and the current schema already applied.
    """
    for table in _V1_VERBATIM_TABLES_BEFORE_STAR_SYSTEMS:
        conn.execute(f"INSERT INTO main.{table} SELECT * FROM old.{table}")

    star_systems_columns = ", ".join(_STAR_SYSTEMS_PRE_V3_COLUMNS)
    conn.execute(
        f"INSERT INTO main.star_systems ({star_systems_columns}) "
        f"SELECT {star_systems_columns} FROM old.star_systems"
    )

    for table in _V1_VERBATIM_TABLES_AFTER_STAR_SYSTEMS:
        conn.execute(f"INSERT INTO main.{table} SELECT * FROM old.{table}")

    columns = ", ".join(_PLANET_MOON_COLUMNS)
    conn.execute(f"INSERT INTO main.planets ({columns}) SELECT {columns} FROM old.planets WHERE is_moon = 0")
    conn.execute(
        f"INSERT INTO main.moons (planet_id, {columns}) "
        f"SELECT parent_planet_id, {columns} FROM old.planets WHERE is_moon = 1"
    )

    conn.execute("""
        INSERT INTO main.planet_evolutionary_paragraphs (id, planet_id, position, paragraph)
        SELECT pep.id, pep.planet_id, pep.position, pep.paragraph
        FROM old.planet_evolutionary_paragraphs pep
        JOIN old.planets p ON p.id = pep.planet_id
        WHERE p.is_moon = 0
    """)
    conn.execute("""
        INSERT INTO main.moon_evolutionary_paragraphs (id, moon_id, position, paragraph)
        SELECT pep.id, pep.planet_id, pep.position, pep.paragraph
        FROM old.planet_evolutionary_paragraphs pep
        JOIN old.planets p ON p.id = pep.planet_id
        WHERE p.is_moon = 1
    """)
    conn.execute("""
        INSERT INTO main.planet_reflection_spectrum (id, planet_id, spectrum_type, position, value)
        SELECT prs.id, prs.planet_id, prs.spectrum_type, prs.position, prs.value
        FROM old.planet_reflection_spectrum prs
        JOIN old.planets p ON p.id = prs.planet_id
        WHERE p.is_moon = 0
    """)
    conn.execute("""
        INSERT INTO main.moon_reflection_spectrum (id, moon_id, spectrum_type, position, value)
        SELECT prs.id, prs.planet_id, prs.spectrum_type, prs.position, prs.value
        FROM old.planet_reflection_spectrum prs
        JOIN old.planets p ON p.id = prs.planet_id
        WHERE p.is_moon = 1
    """)

    _backfill_v3_locations(conn)


# Tables schema v2 -> v3 leaves structurally untouched -- copied verbatim,
# column-for-column, from the attached old database. Split into "before"/
# "after" `star_systems` groups for the same insert-order reason as
# `_V1_VERBATIM_TABLES_BEFORE_STAR_SYSTEMS`/`_AFTER_STAR_SYSTEMS` above
# (several of these tables carry a `star_system_id` FK); `star_systems`
# itself is handled separately (see `_STAR_SYSTEMS_PRE_V3_COLUMNS`).
_V2_VERBATIM_TABLES_BEFORE_STAR_SYSTEMS = ("sectors", "system_configs", "system_config_slots")
_V2_VERBATIM_TABLES_AFTER_STAR_SYSTEMS = (
    "stars", "planets", "planet_evolutionary_paragraphs", "planet_reflection_spectrum",
    "moons", "moon_evolutionary_paragraphs", "moon_reflection_spectrum",
    "asteroid_belts", "asteroid_belt_composition",
)


def _migrate_v2_to_v3(conn):
    """
    Populates a fresh, empty database (`conn`'s main schema, the current
    `schema.sql`) from a schema-v2 database attached as `old` (see
    `migrate_database` for the attach/backup/swap this runs inside of).

    The only structural change between v2 and v3 is the new
    `star_systems.location` column (see `schema.sql`'s "v3" header note).
    Every other table is copied straight across; `star_systems` needs the
    explicit `_STAR_SYSTEMS_PRE_V3_COLUMNS` column list instead of
    `INSERT ... SELECT *`, since the current schema has one more column
    than a v2 database does. `location` is then backfilled for every
    migrated row that has a sector position (`_backfill_v3_locations`) --
    a v2 database predates the column, so it was never computed for these
    rows the way `insert_sector` computes it for a freshly generated
    sector.

    Args:
        conn (sqlite3.Connection): Connection to the new database, with
                                   the old one already `ATTACH`ed as `old`
                                   and the current schema already applied.
    """
    for table in _V2_VERBATIM_TABLES_BEFORE_STAR_SYSTEMS:
        conn.execute(f"INSERT INTO main.{table} SELECT * FROM old.{table}")

    star_systems_columns = ", ".join(_STAR_SYSTEMS_PRE_V3_COLUMNS)
    conn.execute(
        f"INSERT INTO main.star_systems ({star_systems_columns}) "
        f"SELECT {star_systems_columns} FROM old.star_systems"
    )

    for table in _V2_VERBATIM_TABLES_AFTER_STAR_SYSTEMS:
        conn.execute(f"INSERT INTO main.{table} SELECT * FROM old.{table}")

    _backfill_v3_locations(conn)


def migrate_database(db_path):
    """
    Migrates one database file to `SCHEMA_VERSION`, in place, backing up
    the original first. A no-op (returns `None`) if the database is
    already current.

    Safe by construction: the migration reads only from a backup copy
    (never the live file) and builds the new database at a temporary
    path, so `db_path` itself is never touched until the final atomic
    `os.replace` -- a crash or error partway through leaves the original
    file exactly as it was, plus the backup copy already made.

    Per TODO.md's "any future structural change gets a small sequential
    `_migrate_vN_to_vN+1()` function" decision: this knows two source
    versions now -- v1 (`_migrate_v1_to_v2`) and v2 (`_migrate_v2_to_v3`),
    picked by one `if`/`elif` below based on the file's detected version.
    Still deliberately not a generic chaining dispatcher/loop: `migrate_step`
    is invoked exactly once either way, since `_ensure_schema` always
    applies the one current `schema.sql` regardless of source version (there
    is no separate "as of v2" schema snapshot to stop at partway) -- so each
    `_migrate_vN_to_vN+1` function's job is really "map data shaped like
    version N straight into the current schema", not "map N to N+1 and let
    the next hop take it from there". `_migrate_v1_to_v2` in particular
    already goes all the way from v1 to the current (v3) schema in this one
    call, including the v3 `location` backfill -- it doesn't hand off to
    `_migrate_v2_to_v3` partway. A future v4 would need every existing
    `_migrate_vN_to_vN+1` updated to also produce whatever v4 adds (as this
    v3 change did to `_migrate_v1_to_v2`), or a genuine chaining rewrite if
    that keeps growing.

    Args:
        db_path (str): Path to the `.db` file.

    Returns:
        str or None: Path to the backup copy made before migrating, or
                     `None` if no migration was needed.

    Raises:
        UnsupportedSchemaVersionError: If `PRAGMA user_version` is neither
                                       `SCHEMA_VERSION` nor a version this
                                       module has a migration path from.
    """
    probe = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        version = probe.execute("PRAGMA user_version").fetchone()[0]
    finally:
        probe.close()

    if version == SCHEMA_VERSION:
        return None
    if version == 1:
        migrate_step = _migrate_v1_to_v2
    elif version == 2:
        migrate_step = _migrate_v2_to_v3
    else:
        raise UnsupportedSchemaVersionError(
            f"{db_path}: unsupported schema version {version} (expected 1, 2, or {SCHEMA_VERSION})"
        )

    backup_path = f"{db_path}.v{version}-backup-{time.strftime('%Y%m%dT%H%M%S')}.db"
    shutil.copy2(db_path, backup_path)

    tmp_path = f"{db_path}.migrating.tmp"
    if os.path.exists(tmp_path):
        os.remove(tmp_path)
    new_conn = sqlite3.connect(tmp_path)
    try:
        _ensure_schema(new_conn)
        new_conn.execute("ATTACH DATABASE ? AS old", (backup_path,))
        with new_conn:
            migrate_step(new_conn)
        new_conn.execute("DETACH DATABASE old")
    finally:
        new_conn.close()

    os.replace(tmp_path, db_path)
    return backup_path
