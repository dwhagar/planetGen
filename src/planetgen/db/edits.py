# planetgen/db/edits.py

"""
Writes an admin's edits back to the content database (TODO ADM.1): a
star system changed in memory (a body regenerated, removed, reclassified,
or a new star), a phenomenon regenerated in place, and a sector deleted
with everything in it.

A system edit keeps every row it can: a body loaded from the database
(`store.load_star_system` gives each one its row `db_id`) is updated in
place, so its id, its facilities and anything else pointing at it stay.
Only a body that is new gets a new row, and only one that is gone loses
its row (with what points at it).
"""

import uuid

from planetgen.db import store
from planetgen.physics import constants

PHENOMENON_TABLES = {
    "black_hole": "black_holes",
    "neutron_star": "neutron_stars",
    "nebula": "nebulae",
    "supernova_remnant": "supernova_remnants",
    "rogue_planet": "rogue_planets",
    "interstellar_comet": "interstellar_comets",
    "asteroid_field": "asteroid_fields",
    "quasar": "quasars",
}
"""dict: The web's phenomenon type (`/phenomenon/<type>/<id>`) -> its table."""

GENERATOR_TYPES = {
    "black_hole": "black-hole",
    "neutron_star": "neutron-star",
    "nebula": "nebula",
    "supernova_remnant": "supernova-remnant",
    "rogue_planet": "rogue-planet",
    "interstellar_comet": "comet",
    "asteroid_field": "asteroid-field",
    "quasar": "quasar",
}
"""dict: The web's phenomenon type -> `generate.generate_phenomenon`'s."""

_KEPT_PHENOMENON_COLUMNS = {
    "id", "uid", "name", "designation", "sector_id", "star_id",
    "center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc", "quadrant",
    "inside_nebula_id", "inside_remnant_id", "created_at", "modified_at",
}
"""set: What a regenerated phenomenon keeps of its old row: its identity,
name, uid, where it is and what it sits in. Everything else is the new roll."""


class EditError(Exception):
    """An edit that can't be made (the message says why, in plain words)."""


def _star_ids(conn, system_id, system):
    """`[(star_id, planets), ...]`: the `planets.star_id` each of the
    system's planet lists is stored under (`None` for a close pair's
    circumbinary planets)."""
    if system.binary_type == "wide":
        return [(system.primary_star.db_id, system.planets),
                (system.secondary_star.db_id, system.secondary_planets)]
    if system.binary_type == "close":
        return [(None, system.planets)]
    row = conn.execute("SELECT id FROM stars WHERE star_system_id = ? AND role = 'single'", (system_id,)).fetchone()
    return [(row["id"] if row else None, system.planets)]


def _update_body(conn, table, body, assignments, key_values):
    """UPDATEs one `planets`/`moons` row from `body` and rewrites its
    paragraph and spectrum rows."""
    columns = list(assignments) + list(store._BODY_COLUMNS)
    if table == "planets":
        columns += list(store._PLANET_ONLY_COLUMNS)
    values = list(key_values) + store.body_row_values(body)
    # An edited orbit is worked out again by the next orbit update (GEN.106).
    conn.execute(
        f"UPDATE {table} SET {', '.join(f'{c} = ?' for c in columns)}, next_update_due = NULL"
        " WHERE id = ?",
        (*values, body.db_id),
    )
    prefix, id_column = ("moon", "moon_id") if table == "moons" else ("planet", "planet_id")
    conn.execute(f"DELETE FROM {prefix}_evolutionary_paragraphs WHERE {id_column} = ?", (body.db_id,))
    conn.execute(f"DELETE FROM {prefix}_reflection_spectrum WHERE {id_column} = ?", (body.db_id,))
    store.body_child_rows(conn, body, body.db_id)


def _save_moons(conn, system_id, star_id, planet):
    """Updates, inserts and deletes one planet's moon rows."""
    kept = []
    for index, moon in enumerate(planet.moons):
        if getattr(moon, "db_id", None) is not None:
            _update_body(conn, "moons", moon, ("planet_id", "star_id", "orbital_index"),
                         (planet.db_id, star_id, index))
        else:
            moon.db_id = store.insert_moon(conn, moon, system_id, star_id, planet.db_id, index)
        kept.append(moon.db_id)
    _delete_missing(conn, "moons", "planet_id", planet.db_id, kept)


def _delete_missing(conn, table, column, owner_id, kept_ids):
    """Deletes `table`'s rows under `owner_id` whose ids aren't kept."""
    if kept_ids:
        marks = ", ".join("?" * len(kept_ids))
        conn.execute(f"DELETE FROM {table} WHERE {column} = ? AND id NOT IN ({marks})", (owner_id, *kept_ids))
    else:
        conn.execute(f"DELETE FROM {table} WHERE {column} = ?", (owner_id,))


def _save_belt(conn, system_id, star_id, index, belt):
    lower_km = belt.lower_limit * constants.AU_TO_KM
    upper_km = belt.upper_limit * constants.AU_TO_KM
    distance_km = belt.distance * constants.AU_TO_KM
    if getattr(belt, "db_id", None) is None:
        belt.db_id = store.insert_asteroid_belt(conn, belt, system_id, index, star_id=star_id)
        return
    conn.execute(
        "UPDATE asteroid_belts SET star_id = ?, orbital_index = ?, distance_km = ?, lower_limit_km = ?,"
        " upper_limit_km = ?, density = ?, composition_summary = ?"
        " WHERE id = ?",
        (star_id, index, distance_km, lower_km, upper_km, belt.density, belt.get_composition_summary(),
         belt.db_id),
    )
    conn.execute("DELETE FROM asteroid_belt_composition WHERE belt_id = ?", (belt.db_id,))
    for position, (component, concentration) in enumerate(belt.composition):
        conn.execute(
            "INSERT INTO asteroid_belt_composition (belt_id, position, component, concentration)"
            " VALUES (?, ?, ?, ?)",
            (belt.db_id, position, component, concentration),
        )


_STAR_COLUMNS = (
    "name", "star_type", "yerkes_class", "mass_kg", "radius_km",
    "temperature_k", "luminosity_w", "age_gy", "lifespan_gy",
    "habitable_zone_inner_km", "habitable_zone_outer_km",
    "system_perimeter_km", "heliosphere_radius_km",
    "galactic_orbital_speed_kms", "galactic_orbital_period_gy",
    "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
    "wide_binary_a_crit_km",
    "reflex_offset_x_km", "reflex_offset_y_km", "reflex_offset_z_km",
)


def _star_values(star):
    au = constants.AU_TO_KM
    return (
        star.name, star.type, star.yerkes_class, star.mass, star.radius,
        star.temperature, star.luminosity, star.age, store._lifespan_gy(star.lifespan),
        star.habitable_zone[0] * au, star.habitable_zone[1] * au,
        star.system_perimeter * au, star.heliosphere_radius * au,
        star.galactic_orbital_speed_kms, star.galactic_orbital_period_gy,
        star.galactic_orbital_phase_deg, star.galactic_min_update_interval_years,
        star.a_crit_au * au if star.a_crit_au is not None else None,
        star.reflex_offset_x * au, star.reflex_offset_y * au, star.reflex_offset_z * au,
    )


def save_star(conn, star):
    """UPDATEs a loaded star's row (`star.db_id`) from the object; its
    system's galactic speed may have changed, so the next orbit update works
    out the system's due time again (GEN.106)."""
    conn.execute(
        f"UPDATE stars SET {', '.join(f'{c} = ?' for c in _STAR_COLUMNS)}"
        " WHERE id = ?",
        (*_star_values(star), star.db_id),
    )
    conn.execute("UPDATE star_systems ss JOIN stars s ON s.star_system_id = ss.id SET ss.next_update_due = NULL,"
                 " ss.modified_at = ss.modified_at WHERE s.id = ?", (star.db_id,))


def save_system_edits(conn, system_id, system, stars=()):
    """
    Writes an edited `StarSystem` (loaded with `store.load_star_system`) back
    to its rows: every planet, moon and belt it still holds is updated in
    place or inserted, in its current orbital order, and every one it no
    longer holds is deleted. `stars` lists the stars whose own rows
    changed (ADM.7). Comets are written back too (their orbits follow a
    changed star).

    Args:
        conn (Connection): Inside a transaction.
        system_id (int): The `star_systems.id`.
        system (StarSystem): The edited system.
        stars (iterable): Loaded `Star`s (with `db_id`) to update.
    """
    kept_planets, kept_belts = [], []
    for star_id, planets in _star_ids(conn, system_id, system):
        for index, body in enumerate(planets):
            if body.body_type == 'a':
                _save_belt(conn, system_id, star_id, index, body)
                kept_belts.append(body.db_id)
                continue
            if getattr(body, "db_id", None) is None:
                body.db_id = store.insert_planet(conn, body, system_id, star_id, index)
                for moon in body.moons:
                    moon.db_id = None
                rows = conn.execute("SELECT id FROM moons WHERE planet_id = ? ORDER BY orbital_index",
                                    (body.db_id,)).fetchall()
                for moon, row in zip(body.moons, rows):
                    moon.db_id = row["id"]
            else:
                _update_body(conn, "planets", body, ("star_id", "orbital_index"), (star_id, index))
                _save_moons(conn, system_id, star_id, body)
            kept_planets.append(body.db_id)
    _delete_missing(conn, "planets", "star_system_id", system_id, kept_planets)
    _delete_missing(conn, "asteroid_belts", "star_system_id", system_id, kept_belts)
    for star in stars:
        save_star(conn, star)
    for comet in list(system.comets) + list(getattr(system, "secondary_comets", [])):
        if getattr(comet, "db_id", None) is not None:
            _save_comet_orbit(conn, comet)
    store.touch_star_system(conn, system_id)


def _save_comet_orbit(conn, comet):
    au = constants.AU_TO_KM
    conn.execute(
        "UPDATE comets SET perihelion_distance_km = ?, orbital_period_years = ?, primary_mass_solar = ?,"
        " distance_km = ?, position_x_km = ?, position_y_km = ?, position_z_km = ?, orbital_speed_kms = ?,"
        " velocity_x_kms = ?, velocity_y_kms = ?, velocity_z_kms = ?,"
        " min_update_interval_years = ?, next_update_due = NULL WHERE id = ?",
        (comet.perihelion_distance_au * au, comet.orbital_period_years, comet.primary_mass_solar,
         comet.distance_au * au, comet.position_x_au * au, comet.position_y_au * au, comet.position_z_au * au,
         comet.orbital_speed_kms, comet.velocity_x_kms, comet.velocity_y_kms, comet.velocity_z_kms,
         comet.min_update_interval_years, comet.db_id),
    )


# ---------------------------------------------------------------------
# Phenomena
# ---------------------------------------------------------------------

def _table_columns(conn, table):
    rows = conn.execute(
        "SELECT COLUMN_NAME AS name FROM information_schema.COLUMNS"
        " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ? ORDER BY ORDINAL_POSITION",
        (table,),
    ).fetchall()
    return [row["name"] for row in rows]


def _child_tables(conn, table):
    """`[(child_table, column)]`: the tables whose rows belong to a row of
    `table` (a cascading foreign key), other than facilities, which stay
    with the kept row."""
    rows = conn.execute(
        "SELECT k.TABLE_NAME AS child, k.COLUMN_NAME AS col FROM information_schema.KEY_COLUMN_USAGE k"
        " JOIN information_schema.REFERENTIAL_CONSTRAINTS r"
        "   ON r.CONSTRAINT_SCHEMA = k.CONSTRAINT_SCHEMA AND r.CONSTRAINT_NAME = k.CONSTRAINT_NAME"
        " WHERE k.TABLE_SCHEMA = DATABASE() AND k.REFERENCED_TABLE_NAME = ? AND r.DELETE_RULE = 'CASCADE'",
        (table,),
    ).fetchall()
    return [(row["child"], row["col"]) for row in rows if row["child"] != "facilities"]


def phenomenon_row(conn, phenomenon_type, phenomenon_id):
    """The phenomenon's row, or `None`. Raises `EditError` for an unknown
    type."""
    table = PHENOMENON_TABLES.get(phenomenon_type)
    if table is None:
        raise EditError(f"unknown phenomenon type: {phenomenon_type}")
    return conn.execute(f"SELECT * FROM {table} WHERE id = ?", (phenomenon_id,)).fetchone()


def replace_phenomenon_content(conn, phenomenon_type, phenomenon_id, phenomenon):
    """
    Swaps a stored phenomenon's generated content for `phenomenon`'s (a
    freshly generated object of the same type), keeping its id, name,
    sector, galaxy position and containment. The new object is inserted
    as a stand-in row first, its content copied onto the kept row, its
    child rows (composition) moved over, then the stand-in deleted.

    Raises:
        EditError: For a black hole or neutron star that is a star
            system's star (regenerate the system instead).
    """
    table = PHENOMENON_TABLES[phenomenon_type]
    row = phenomenon_row(conn, phenomenon_type, phenomenon_id)
    if row is None:
        raise EditError("not found")
    if row.get("star_id") is not None:
        raise EditError("this is the star of a star system; regenerate the system instead")
    inserter = store._PHENOMENON_INSERTERS[GENERATOR_TYPES[phenomenon_type]]
    placement = {key: row.get(key) for key in ("center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc")}
    # A throwaway name for the stand-in row, so reserving it never
    # collides with (and renames) the kept row's own name.
    phenomenon.name = f"Regenerating {uuid.uuid4().hex}"
    temp_id = inserter(conn, phenomenon, sector_id=row.get("sector_id"),
                       placement=placement if placement["center_x_pc"] is not None else None)
    columns = [c for c in _table_columns(conn, table) if c not in _KEPT_PHENOMENON_COLUMNS]
    if columns:
        assignments = ", ".join(f"kept.{c} = made.{c}" for c in columns)
        conn.execute(
            f"UPDATE {table} kept JOIN {table} made ON made.id = ? SET {assignments} WHERE kept.id = ?",
            (temp_id, phenomenon_id),
        )
    for child, column in _child_tables(conn, table):
        conn.execute(f"DELETE FROM {child} WHERE {column} = ?", (phenomenon_id,))
        conn.execute(f"UPDATE {child} SET {column} = ? WHERE {column} = ?", (phenomenon_id, temp_id))
    conn.execute(f"DELETE FROM {table} WHERE id = ?", (temp_id,))
    conn.execute("DELETE FROM system_name_registry WHERE first_object_table = ? AND first_object_id = ?",
                 (table, temp_id))
    conn.execute(f"UPDATE {table} SET modified_at = CURRENT_TIMESTAMP(3) WHERE id = ?", (phenomenon_id,))


def delete_phenomenon(conn, phenomenon_type, phenomenon_id):
    """Deletes one standalone phenomenon. Returns `False` when there is no
    such row; raises `EditError` for a system's own star."""
    row = phenomenon_row(conn, phenomenon_type, phenomenon_id)
    if row is None:
        return False
    if row.get("star_id") is not None:
        raise EditError("this is the star of a star system; delete the system instead")
    table = PHENOMENON_TABLES[phenomenon_type]
    conn.execute("DELETE FROM nearest_systems WHERE object_table = ? AND object_id = ?", (table, phenomenon_id))
    conn.execute(f"DELETE FROM {table} WHERE id = ?", (phenomenon_id,))
    store.touch_sector(conn, row.get("sector_id"))
    return True


# ---------------------------------------------------------------------
# Sectors
# ---------------------------------------------------------------------

def sector_address(conn, sector_id):
    """`(ring, layer, slot)` of a sector, `None` for one off the grid, or
    raises `EditError` when there's no such sector."""
    row = conn.execute("SELECT ring_index, layer_index, ring_slot_index FROM sectors WHERE id = ?",
                       (sector_id,)).fetchone()
    if row is None:
        raise EditError("not found")
    if row["ring_index"] is None:
        return None
    return row["ring_index"], row["layer_index"], row["ring_slot_index"]


def delete_sector_with_contents(conn, sector_id):
    """
    Deletes a sector and everything filed under it: its star systems
    (with their bodies and facilities), the standalone phenomena whose
    home is this sector, and facilities parked in it. Its grid slot is
    left unfilled, so it can be generated again. The neighboring
    sectors' nearest-system lists are refreshed.

    Returns:
        dict: `{"systems": n, "phenomena": n}` deleted, or `None` when
        there was no such sector.
    """
    row = conn.execute("SELECT center_x_pc, center_y_pc, center_z_pc, ring_index, layer_index, ring_slot_index"
                       " FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    if row is None:
        return None
    systems = conn.execute("DELETE FROM star_systems WHERE sector_id = ?", (sector_id,)).rowcount
    phenomena = 0
    for table in PHENOMENON_TABLES.values():
        where = "sector_id = ?" + (" AND star_id IS NULL" if table in ("black_holes", "neutron_stars") else "")
        phenomena += conn.execute(f"DELETE FROM {table} WHERE {where}", (sector_id,)).rowcount
    conn.execute("DELETE FROM facilities WHERE sector_id = ?", (sector_id,))
    conn.execute("DELETE FROM nearest_systems WHERE sector_id = ?", (sector_id,))
    conn.execute("DELETE FROM sectors WHERE id = ?", (sector_id,))
    # GEN.44: the slot goes back to the bright-star level it had unfilled.
    store.forget_sector_fill(conn, (row["ring_index"], row["layer_index"], row["ring_slot_index"]))
    if row["center_x_pc"] is not None:
        refresh_nearest_around(conn, (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]))
    return {"systems": systems, "phenomena": phenomena}


def refresh_nearest_around(conn, center_pc):
    """Refreshes the nearest-system lists of every sector close enough to
    `center_pc` to have listed something there."""
    half_diagonal = store._edge_pc(conn) * 3 ** 0.5 / 2
    near = store._sectors_near(conn, {0: center_pc}, 2 * half_diagonal + store.NEAREST_SYSTEMS_SEARCH_PC)
    if near:
        store.refresh_nearest_systems(conn, near.keys())


def system_sector_center(conn, system_id):
    """The galaxy-frame center (pc) of a system's sector, or `None`."""
    row = conn.execute(
        "SELECT s.center_x_pc AS x, s.center_y_pc AS y, s.center_z_pc AS z FROM star_systems ss"
        " JOIN sectors s ON s.id = ss.sector_id WHERE ss.id = ?",
        (system_id,),
    ).fetchone()
    if row is None or row["x"] is None:
        return None
    return row["x"], row["y"], row["z"]
