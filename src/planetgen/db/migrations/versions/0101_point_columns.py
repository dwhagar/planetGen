"""The point of every kind of object, in one additive revision (DB.55, for GEN.294, GEN.289, GEN.142, GEN.110, GEN.109).

Schema v101. Every object in the galaxy is one `SpatialPosition3D`: a position, a point mass and a movement vector at an
epoch (`planetgen/db/point_columns.py`). This revision adds, nullable and backfilled, what the tables lack of that, so
the code items that fill them (GEN.336, GEN.337, GEN.289) do not each need a revision of their own:

* A movement vector (`velocity_x/y/z_kms`, galactic axes) on black holes, neutron stars, rogue planets, interstellar
  comets, asteroid fields, nebulae, remnants, belts and the quasar. The scalar `galactic_orbital_speed_kms` stays
  until the code stops reading it (a later small revision drops it); the vector is backfilled from it as the tangent of
  the galaxy's rotation at the object's place, counterclockwise from galactic north, at the same speed
  (`galaxy/system_position.galactic_velocity_ms`). The quasar and a belt, which do not move, take a stated zero. A star
  has no columns of its own: it moves with its system, whose point, vector and epoch it reads.
* A mass (`mass_kg`) on interstellar comets, asteroid fields, nebulae, remnants, belts, comets and facilities. A comet
  nucleus (interstellar or bound) takes its sphere's mass at the density of 67P, and a nebula or remnant the mass of its bounding sphere at its
  particle density (an upper bound for a cloud that fills part of its sphere); a body whose mass is not known (a field's
  or belt's bulk) takes a stated zero, which the gravity sum skips by rule; a facility takes the mass of
  what it sits on (planet, moon or star), else zero.
* An epoch and a due time (`epoch_unix`, `next_update_due`) on belts and the quasar, and an epoch on scatter members;
  NULL until the first orbit update, as on the tables that have them.
* A system-frame position (`position_x/y/z_km`, from the system centre) on belts, which sit at the centre (zero).
* `sectors.max_mass_kg` and `sectors.max_mass_id` (the heaviest object of the sector and its uid, for the influence
  search of GEN.109; NULL until it is worked out), `asteroid_fields.formed_at` and `asteroid_fields.origin_event_id`
  (when a field formed and the event that made it, GEN.110; the event table is DB.56, so no foreign key), and for
  planets and moons the state GEN.289 needs: `parent_type` (`star`, `planet` or `system`), `parent_id` and
  `dynamical_state` (`bound` or `unbound`), backfilled from `star_id` and `planet_id`.

Each step is skipped when it is already done, so a stopped run can be started again; the `sectors` columns go last and
mark the revision done. Columns go on the end of their table, which MySQL adds without copying the rows.
"""

import sqlalchemy as sa
from alembic import op

revision = "0101"
down_revision = "0100"
branch_labels = None
depends_on = None

VECTOR = tuple((f"velocity_{axis}_kms", "DOUBLE NULL") for axis in "xyz")
POSITION = tuple((f"position_{axis}_km", "DOUBLE NULL") for axis in "xyz")
MASS = (("mass_kg", "DOUBLE NULL"),)
CLOCK = (("epoch_unix", "DOUBLE NULL"), ("next_update_due", "DOUBLE NULL"))

PLANET_STATE = (
    ("parent_type", "VARCHAR(8) NOT NULL DEFAULT 'star'"),
    ("parent_id", "BIGINT UNSIGNED NULL"),
    ("dynamical_state", "VARCHAR(8) NOT NULL DEFAULT 'bound'"),
)
MOON_STATE = (
    ("parent_type", "VARCHAR(8) NOT NULL DEFAULT 'planet'"),
    ("parent_id", "BIGINT UNSIGNED NULL"),
    ("dynamical_state", "VARCHAR(8) NOT NULL DEFAULT 'bound'"),
)

ADDITIONS = {
    "black_holes": VECTOR,
    "neutron_stars": VECTOR,
    "rogue_planets": VECTOR,
    "interstellar_comets": VECTOR + MASS,
    "asteroid_fields": VECTOR + MASS + (("formed_at", "DOUBLE NULL"), ("origin_event_id", "BIGINT UNSIGNED NULL")),
    "nebulae": VECTOR + MASS,
    "supernova_remnants": VECTOR + MASS,
    "quasars": VECTOR + CLOCK,
    "asteroid_belts": POSITION + VECTOR + MASS + CLOCK,
    "comets": MASS,
    "facilities": MASS,
    "scatter_members": (("epoch_unix", "DOUBLE NULL"),),
    "planets": PLANET_STATE,
    "moons": MOON_STATE,
}
"""dict: Table -> `(column, definition)` added, in the order they are added."""

SECTOR_COLUMNS = (("max_mass_kg", "DOUBLE NULL"), ("max_mass_id", "BINARY(12) NULL"))
"""tuple: The `sectors` columns, added after every backfill: they mark the revision done."""

CHECKS = {
    "planets": (("chk_planets_parent_type", "parent_type IN ('star', 'planet', 'system')"),
                ("chk_planets_dynamical_state", "dynamical_state IN ('bound', 'unbound')")),
    "moons": (("chk_moons_parent_type", "parent_type IN ('star', 'planet', 'system')"),
              ("chk_moons_dynamical_state", "dynamical_state IN ('bound', 'unbound')")),
}
"""dict: Table -> `(constraint name, clause)` added with the state columns."""

# Free objects that carry the scalar speed and a centre in parsecs: the vector is the tangent at the same speed.
TANGENT_TABLES = ("black_holes", "neutron_stars", "rogue_planets", "interstellar_comets", "asteroid_fields", "nebulae",
                  "supernova_remnants")
# Objects that do not move on their own: a stated zero.
ZERO_VECTOR_TABLES = ("quasars", "asteroid_belts")
# Bodies whose mass is not known: a stated zero, which the gravity sum skips.
ZERO_MASS_TABLES = ("asteroid_fields", "asteroid_belts")
# A cloud's mass is its bounding sphere at its density, with the interstellar gas's helium counted; a comet nucleus's is
# its sphere at 67P's density (`physics/bulk_mass.py`, `tuning.COMET_NUCLEUS_DENSITY_KG_M3`). The constants are written
# out so the revision does not change with a later retuning.
LIGHT_YEAR_M = 9.4607304725808e15
HYDROGEN_ATOM_KG = 1.67353286206015e-27
GAS_MASS_PER_HYDROGEN = 1.4
COMET_NUCLEUS_DENSITY = 533.0


def _has_column(connection, table, column):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar())


def _has_check(connection, name):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.TABLE_CONSTRAINTS WHERE TABLE_SCHEMA = DATABASE()"
        " AND CONSTRAINT_NAME = :name AND CONSTRAINT_TYPE = 'CHECK'"), {"name": name}).scalar())


def _add_columns(connection, table, columns):
    missing = [(name, definition) for name, definition in columns if not _has_column(connection, table, name)]
    if missing:
        adds = ", ".join(f"ADD COLUMN {name} {definition}" for name, definition in missing)
        op.execute(f"ALTER TABLE {table} {adds}")


def _tangent(table):
    """UPDATE giving the vector of a free object: `speed` along the tangent at its centre, counterclockwise seen from
    galactic north, at rest on the axis or without a speed."""
    r = f"SQRT({table}.center_x_pc * {table}.center_x_pc + {table}.center_y_pc * {table}.center_y_pc)"
    moving = f"{table}.galactic_orbital_speed_kms > 0 AND {r} > 0"
    return (
        f"UPDATE {table} SET"
        f" velocity_x_kms = CASE WHEN {moving} THEN -{table}.galactic_orbital_speed_kms * {table}.center_y_pc / {r}"
        f" ELSE 0 END,"
        f" velocity_y_kms = CASE WHEN {moving} THEN {table}.galactic_orbital_speed_kms * {table}.center_x_pc / {r}"
        f" ELSE 0 END,"
        f" velocity_z_kms = 0 WHERE velocity_x_kms IS NULL")


def upgrade():
    connection = op.get_bind()
    for table, columns in ADDITIONS.items():
        _add_columns(connection, table, columns)

    for table, checks in CHECKS.items():
        for name, clause in checks:
            if not _has_check(connection, name):
                op.execute(f"ALTER TABLE {table} ADD CONSTRAINT {name} CHECK ({clause})")

    # Each backfill writes only the rows still NULL, so a run that is started again picks up where it stopped.
    for table in TANGENT_TABLES:
        op.execute(_tangent(table))
    for table in ZERO_VECTOR_TABLES:
        op.execute(f"UPDATE {table} SET velocity_x_kms = 0, velocity_y_kms = 0, velocity_z_kms = 0"
                   f" WHERE velocity_x_kms IS NULL")
    for table in ZERO_MASS_TABLES:
        op.execute(f"UPDATE {table} SET mass_kg = 0 WHERE mass_kg IS NULL")
    for table in ("nebulae", "supernova_remnants"):
        op.execute(
            f"UPDATE {table} SET mass_kg = 4.0 / 3.0 * PI() * POW(radius_ly * {LIGHT_YEAR_M!r}, 3)"
            f" * density_cm3 * 1000000.0 * {HYDROGEN_ATOM_KG!r} * {GAS_MASS_PER_HYDROGEN!r} WHERE mass_kg IS NULL")
    for table in ("interstellar_comets", "comets"):
        op.execute(
            f"UPDATE {table} SET mass_kg = 4.0 / 3.0 * PI() * POW(nucleus_diameter_km * 500.0, 3)"
            f" * {COMET_NUCLEUS_DENSITY!r} WHERE mass_kg IS NULL")

    op.execute("UPDATE asteroid_belts SET position_x_km = 0, position_y_km = 0, position_z_km = 0"
               " WHERE position_x_km IS NULL")
    op.execute(
        "UPDATE facilities f SET f.mass_kg = CASE"
        " WHEN f.planet_id IS NOT NULL THEN COALESCE((SELECT mass_kg FROM planets WHERE id = f.planet_id), 0)"
        " WHEN f.moon_id IS NOT NULL THEN COALESCE((SELECT mass_kg FROM moons WHERE id = f.moon_id), 0)"
        " WHEN f.star_id IS NOT NULL THEN COALESCE((SELECT mass_kg FROM stars WHERE id = f.star_id), 0)"
        " ELSE 0 END WHERE f.mass_kg IS NULL")
    op.execute("UPDATE planets SET parent_id = star_id WHERE parent_id IS NULL AND star_id IS NOT NULL")
    op.execute("UPDATE moons SET parent_id = planet_id WHERE parent_id IS NULL AND planet_id IS NOT NULL")
    _add_columns(connection, "sectors", SECTOR_COLUMNS)
