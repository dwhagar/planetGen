"""96-bit object IDs, the grid's reach, and a point in space for the scatter (GEN.230, GEN.233, GEN.223).

Schema v93. Three changes that renumber or reshape the same rows, so they share one migration:

* GEN.230. The object ID goes from 80 to 96 bits. The birth sector's three fields were 12, 12 and 16 bits (shell, biased
  band, wedge); they are 19, 18 and 19 bits now, so the 12-bit band field no longer stops the grid at shell 1,303 and
  the grid can be named out to the 262 kpc reach (shell 65,535) and a little past it (`spherical.check_reach`). Each
  `uid` (the fifteen object tables and `nebula_field`) is `BINARY(12)`; the sector fields move to their new places and
  the serial and body number are untouched, so a body keeps its number in its system. `sectors.uid` (the designation)
  moves the same way (a sector with no grid address counts up from `1 << 56` instead of `1 << 40`), and the scopes of
  `id_counters` follow (a sector's is 7 bytes, a system's 12).
* GEN.233. `galaxy_shape.reach_shell_index`, 65,535 for a galaxy planned before it.
* The companions of a scattered system (docs/design/systems/multi-star-systems.md section 5): `member_count` and
  `system_mass_solar` on both scatters and the tables `scatter_members` and `scatter_nodes`, empty until the
  member-tree generator writes them.
* GEN.223. `bright_stars` and `phenomenon_scatter` keep their point in space the way every other object does: position
  (milliparsecs, now DOUBLE so a slow drift is not rounded away), velocity vector, `epoch_unix` and an indexed
  `next_update_due`. A row without a velocity is on the rotation curve until the orbit update gives it one.

* Region queries (docs/design/galaxy/spherical-region-queries.md). Every table that holds a point in the galaxy frame gets
  `cell_key`, a stored generated column with the Morton key of its 0.25 pc cube, and the covering index
  `(cell_key, x, y, z)`. A star system keeps its position relative to its sector, so it also gets `galaxy_x/y/z_pc`
  (its point in the galaxy frame), written here for the systems already stored and by every save and move after.
* The stored form of the member trees also carries what the multi-star catalogue needs (docs/design/systems/multi-star-systems.md
  sections 13 and 14): a member's `state`, `flags`, and, for a young group, its offset, velocity and `decay_age_yr`; a
  node's `interaction`; and a SMALLINT `slot` for a group of up to 12.

A stopped run can be started again: each step skips what it has done, and the last one (the reach column) marks the whole as done.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner
from planetgen.galaxy import cellkey

revision = "0093"
down_revision = "0092"
branch_labels = None
depends_on = None

CHUNK = 20000

UID_TABLES = (
    # table, the name of its index on uid, whether that index is UNIQUE
    ("star_systems", "uq_star_systems_uid", True),
    ("stars", "uq_stars_uid", True),
    ("planets", "uq_planets_uid", True),
    ("moons", "uq_moons_uid", True),
    ("asteroid_belts", "uq_asteroid_belts_uid", True),
    ("comets", "uq_comets_uid", True),
    ("black_holes", "uq_black_holes_uid", True),
    ("neutron_stars", "uq_neutron_stars_uid", True),
    ("nebulae", "uq_nebulae_uid", True),
    ("supernova_remnants", "uq_supernova_remnants_uid", True),
    ("quasars", "uq_quasars_uid", True),
    ("rogue_planets", "uq_rogue_planets_uid", True),
    ("interstellar_comets", "uq_interstellar_comets_uid", True),
    ("asteroid_fields", "uq_asteroid_fields_uid", True),
    ("facilities", "uq_facilities_uid", True),
    ("nebula_field", "idx_nebula_field_uid", False),
)

SHELL_SHIFT = 37
"""int: The shell sits above the 18-bit band and 19-bit wedge fields."""
BAND_SHIFT = 19
BAND_REBIAS = (1 << 17) - (1 << 11)
"""int: The new band bias (2 ** 17) less the old one (2 ** 11), added to an old biased band."""
UNPLACED_OLD = 1 << 40
UNPLACED_NEW = 1 << 56


def _sector_sql(hex_digits):
    """SQL for the new 56-bit sector field of an old 10-hex-digit one (`hex_digits` is an expression for it)."""
    shell = f"CAST(CONV(SUBSTR({hex_digits}, 1, 3), 16, 10) AS UNSIGNED)"
    band = f"CAST(CONV(SUBSTR({hex_digits}, 4, 3), 16, 10) AS UNSIGNED)"
    wedge = f"CAST(CONV(SUBSTR({hex_digits}, 7, 4), 16, 10) AS UNSIGNED)"
    return f"(({shell} << {SHELL_SHIFT}) | (({band} + {BAND_REBIAS}) << {BAND_SHIFT}) | {wedge})"


def _hex14(sector_sql):
    return f"LPAD(HEX({sector_sql}), 14, '0')"


def _new_uid_sql(column):
    """SQL for the 12-byte ID of an old 10-byte one in `column`."""
    old = f"HEX({column})"
    return f"UNHEX(CONCAT({_hex14(_sector_sql(f'SUBSTR({old}, 1, 10)'))}, SUBSTR({old}, 11, 10)))"


def _has_column(connection, table, column):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :t AND COLUMN_NAME = :c"), {"t": table, "c": column}).scalar()


def _column_type(connection, table, column):
    return connection.execute(sa.text(
        "SELECT COLUMN_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :t AND COLUMN_NAME = :c"), {"t": table, "c": column}).scalar()


def _has_index(connection, table, index):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.STATISTICS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :t AND INDEX_NAME = :i"), {"t": table, "i": index}).scalar()


def _touch(connection, table):
    """`, modified_at = modified_at` for a table that has the column, so an UPDATE leaves the time alone."""
    return ", modified_at = modified_at" if _has_column(connection, table, "modified_at") else ""


def _renumber_column(connection, table, index, unique, column_type, expression, label, where="uid IS NOT NULL"):
    """
    Replaces `table.uid` by the value of `expression` (SQL over the old `uid`), as `column_type`: the new value goes
    into a `uid_new` column in pieces, then the old column and its index are swapped for it in one rebuild of the
    table. Nothing is lost when this stops half way: `uid_new` is filled only where it is empty, and the swap is
    one statement.
    """
    if _column_type(connection, table, "uid") == column_type.lower():
        return
    if not _has_column(connection, table, "uid_new"):
        connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN uid_new {column_type} NULL"))
    touch = _touch(connection, table)
    if not _has_column(connection, table, "id"):
        # nebula_field: a cache keyed by cell, small beside the object tables.
        connection.execute(sa.text(
            f"UPDATE {table} SET uid_new = {expression}{touch} WHERE {where} AND uid_new IS NULL"))
    else:
        low, high = connection.execute(sa.text(f"SELECT MIN(id), MAX(id) FROM {table}")).fetchone()
        if low is not None:
            total = high - low + 1
            for first in range(low, high + 1, CHUNK):
                connection.execute(sa.text(
                    f"UPDATE {table} SET uid_new = {expression}{touch}"
                    f" WHERE id BETWEEN :a AND :b AND {where} AND uid_new IS NULL"),
                    {"a": first, "b": min(first + CHUNK - 1, high)})
                alembic_runner.report_progress(label, min(first + CHUNK - low, total), total)
    drop_index = f"DROP INDEX {index}, " if _has_index(connection, table, index) else ""
    add_index = f"ADD UNIQUE KEY {index} (uid)" if unique else f"ADD KEY {index} (uid)"
    connection.execute(sa.text(
        f"ALTER TABLE {table} {drop_index}DROP COLUMN uid, CHANGE COLUMN uid_new uid {column_type} NULL, {add_index}"))


def _renumber_sector_uids(connection):
    """`sectors.uid`, the designation, as the 56-bit field; a sector with no grid address keeps its offset."""
    designation = _sector_sql("LPAD(HEX(uid), 10, '0')")
    expression = (f"CASE WHEN uid >= {UNPLACED_OLD} THEN uid - {UNPLACED_OLD} + {UNPLACED_NEW}"
                  f" ELSE {designation} END")
    if not _has_column(connection, "sectors", "uid_new"):
        connection.execute(sa.text("ALTER TABLE sectors ADD COLUMN uid_new BIGINT UNSIGNED NULL"))
    low, high = connection.execute(sa.text("SELECT MIN(id), MAX(id) FROM sectors")).fetchone()
    if low is not None:
        total = high - low + 1
        for first in range(low, high + 1, CHUNK):
            connection.execute(sa.text(
                f"UPDATE sectors SET uid_new = {expression}, modified_at = modified_at"
                " WHERE id BETWEEN :a AND :b AND uid IS NOT NULL AND uid_new IS NULL"),
                {"a": first, "b": min(first + CHUNK - 1, high)})
            alembic_runner.report_progress("Sector IDs", min(first + CHUNK - low, total), total)
    drop_index = "DROP INDEX uq_sectors_uid, " if _has_index(connection, "sectors", "uq_sectors_uid") else ""
    connection.execute(sa.text(
        f"ALTER TABLE sectors {drop_index}DROP COLUMN uid, CHANGE COLUMN uid_new uid BIGINT UNSIGNED NULL,"
        " ADD UNIQUE KEY uq_sectors_uid (uid)"))


def _renumber_counters(connection):
    """A sector's counter is keyed by its 5-byte address, a system's by its 10-byte ID; both widen."""
    sector = _hex14(_sector_sql("HEX(scope)"))
    connection.execute(sa.text(
        f"INSERT IGNORE INTO id_counters (kind, scope, next_value) SELECT kind, UNHEX({sector}), next_value"
        " FROM id_counters WHERE kind = 'sector' AND LENGTH(scope) = 5"))
    connection.execute(sa.text("DELETE FROM id_counters WHERE kind = 'sector' AND LENGTH(scope) = 5"))
    system = f"UNHEX(CONCAT({_hex14(_sector_sql('SUBSTR(HEX(scope), 1, 10)'))}, SUBSTR(HEX(scope), 11, 10)))"
    connection.execute(sa.text(
        f"INSERT IGNORE INTO id_counters (kind, scope, next_value) SELECT kind, {system}, next_value"
        " FROM id_counters WHERE kind = 'system' AND LENGTH(scope) = 10"))
    connection.execute(sa.text("DELETE FROM id_counters WHERE kind = 'system' AND LENGTH(scope) = 10"))


def _reach(connection):
    if not _has_column(connection, "galaxy_shape", "reach_shell_index"):
        connection.execute(sa.text(
            "ALTER TABLE galaxy_shape ADD COLUMN reach_shell_index INT NOT NULL DEFAULT 65535 AFTER outer_shell_index"))


MEMBER_TABLES = (
    """CREATE TABLE IF NOT EXISTS scatter_members (
    source               TINYINT NOT NULL,
    source_id            BIGINT UNSIGNED NOT NULL,
    member_no            SMALLINT NOT NULL,
    parent_node_no       SMALLINT NOT NULL,
    slot                 SMALLINT NOT NULL,
    kind                 VARCHAR(20) NOT NULL,
    state                VARCHAR(16),
    flags                SMALLINT UNSIGNED,
    offset_x_au          DOUBLE,
    offset_y_au          DOUBLE,
    offset_z_au          DOUBLE,
    velocity_x_kms       DOUBLE,
    velocity_y_kms       DOUBLE,
    velocity_z_kms       DOUBLE,
    decay_age_yr         DOUBLE,
    mass_solar           DOUBLE NOT NULL,
    initial_mass_solar   DOUBLE NOT NULL,
    seed                 BIGINT UNSIGNED NOT NULL,
    PRIMARY KEY (source, source_id, member_no),
    CONSTRAINT chk_scatter_members_source CHECK (source IN (1, 2)),
    CONSTRAINT chk_scatter_members_kind CHECK (kind IN (
        'main-sequence', 'giant', 'white-dwarf', 'neutron-star', 'black-hole', 'brown-dwarf',
        'pre-main-sequence', 'protostar', 'stripped-helium', 'thorne-zytkow', 'wolf-rayet'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;""",
    """CREATE TABLE IF NOT EXISTS scatter_nodes (
    source               TINYINT NOT NULL,
    source_id            BIGINT UNSIGNED NOT NULL,
    node_no              SMALLINT NOT NULL,
    parent_node_no       SMALLINT,
    slot                 SMALLINT,
    mass_solar           DOUBLE NOT NULL,
    separation_km        DOUBLE NOT NULL,
    eccentricity         DOUBLE NOT NULL,
    inclination_deg      DOUBLE NOT NULL,
    ascending_node_deg   DOUBLE NOT NULL,
    arg_periapsis_deg    DOUBLE NOT NULL,
    mean_anomaly_deg     DOUBLE NOT NULL,
    epoch_unix           DOUBLE NOT NULL,
    period_years         DOUBLE NOT NULL,
    kind                 VARCHAR(16) NOT NULL,
    interaction          VARCHAR(20),
    flags                TINYINT UNSIGNED NOT NULL DEFAULT 0,
    PRIMARY KEY (source, source_id, node_no),
    CONSTRAINT chk_scatter_nodes_source CHECK (source IN (1, 2)),
    CONSTRAINT chk_scatter_nodes_kind CHECK (kind IN ('close', 'wide', 'compact', 'group')),
    CONSTRAINT chk_scatter_nodes_interaction CHECK (interaction IS NULL OR interaction IN (
        'detached', 'semi-detached', 'contact', 'accreting', 'post-common-envelope', 'merger-product'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;""",
)
"""tuple: DDL of the tables that hold the companions of a scattered system (multi-star note, section 5)."""


def _scatter(connection):
    """The point in space of the two scatters: DOUBLE positions, a velocity vector, an epoch and a due time."""
    for ddl in MEMBER_TABLES:
        connection.execute(sa.text(ddl))
    if not _has_column(connection, "bright_stars", "next_update_due"):
        connection.execute(sa.text(
            "ALTER TABLE bright_stars MODIFY COLUMN position_x_mpc DOUBLE NOT NULL,"
            " MODIFY COLUMN position_y_mpc DOUBLE NOT NULL, MODIFY COLUMN position_z_mpc DOUBLE NOT NULL,"
            " ADD COLUMN velocity_x_kms DOUBLE NULL, ADD COLUMN velocity_y_kms DOUBLE NULL,"
            " ADD COLUMN velocity_z_kms DOUBLE NULL, ADD COLUMN epoch_unix DOUBLE NULL,"
            " ADD COLUMN next_update_due DOUBLE NULL, ADD COLUMN member_count SMALLINT UNSIGNED NULL,"
            " ADD COLUMN system_mass_solar DOUBLE NULL,"
            " ADD KEY idx_bright_stars_next_update_due (next_update_due)"))
    if not _has_column(connection, "phenomenon_scatter", "next_update_due"):
        connection.execute(sa.text(
            "ALTER TABLE phenomenon_scatter MODIFY COLUMN position_x_mpc DOUBLE NOT NULL,"
            " MODIFY COLUMN position_y_mpc DOUBLE NOT NULL, MODIFY COLUMN position_z_mpc DOUBLE NOT NULL,"
            " ADD COLUMN next_update_due DOUBLE NULL, ADD COLUMN member_count SMALLINT UNSIGNED NULL,"
            " ADD COLUMN system_mass_solar DOUBLE NULL,"
            " ADD KEY idx_phenomenon_scatter_next_update_due (next_update_due)"))

PC_COLUMNS = ("center_x_pc", "center_y_pc", "center_z_pc")
MPC_COLUMNS = ("position_x_mpc", "position_y_mpc", "position_z_mpc")
GALAXY_COLUMNS = ("galaxy_x_pc", "galaxy_y_pc", "galaxy_z_pc")
CELL_TABLES = (
    # table, the columns of its point and the SQL for them in parsecs (the cell key is computed from those)
    ("black_holes", PC_COLUMNS, PC_COLUMNS),
    ("neutron_stars", PC_COLUMNS, PC_COLUMNS),
    ("nebulae", PC_COLUMNS, PC_COLUMNS),
    ("supernova_remnants", PC_COLUMNS, PC_COLUMNS),
    ("quasars", PC_COLUMNS, PC_COLUMNS),
    ("rogue_planets", PC_COLUMNS, PC_COLUMNS),
    ("interstellar_comets", PC_COLUMNS, PC_COLUMNS),
    ("asteroid_fields", PC_COLUMNS, PC_COLUMNS),
    ("facilities", PC_COLUMNS, PC_COLUMNS),
    ("nebula_field", PC_COLUMNS, PC_COLUMNS),
    ("bright_stars", MPC_COLUMNS, tuple(f"{c} / 1000.0" for c in MPC_COLUMNS)),
    ("phenomenon_scatter", MPC_COLUMNS, tuple(f"{c} / 1000.0" for c in MPC_COLUMNS)),
    ("star_systems", GALAXY_COLUMNS, GALAXY_COLUMNS),
)


def _galaxy_points(connection):
    """`star_systems.galaxy_*_pc`: the sector's center plus the system's offset, turned to the galaxy's axes."""
    for column in GALAXY_COLUMNS:
        if not _has_column(connection, "star_systems", column):
            connection.execute(sa.text(f"ALTER TABLE star_systems ADD COLUMN {column} DOUBLE NULL"))
    radius = "SQRT(sec.center_x_pc * sec.center_x_pc + sec.center_y_pc * sec.center_y_pc)"
    cos = f"IF({radius} < 1e-12, 1, sec.center_x_pc / {radius})"
    sin = f"IF({radius} < 1e-12, 0, sec.center_y_pc / {radius})"
    px, py, pz = ("ss.position_x_mpc / 1000.0", "ss.position_y_mpc / 1000.0", "ss.position_z_mpc / 1000.0")
    low, high = connection.execute(sa.text("SELECT MIN(id), MAX(id) FROM star_systems")).fetchone()
    if low is None:
        return
    total = high - low + 1
    for first in range(low, high + 1, CHUNK):
        connection.execute(sa.text(
            "UPDATE star_systems ss JOIN sectors sec ON sec.id = ss.sector_id SET"
            f" ss.galaxy_x_pc = sec.center_x_pc + {px} * {cos} - {py} * {sin},"
            f" ss.galaxy_y_pc = sec.center_y_pc + {px} * {sin} + {py} * {cos},"
            f" ss.galaxy_z_pc = sec.center_z_pc + {pz}, ss.modified_at = ss.modified_at"
            " WHERE ss.id BETWEEN :a AND :b AND ss.position_x_mpc IS NOT NULL AND sec.center_x_pc IS NOT NULL"),
            {"a": first, "b": min(first + CHUNK - 1, high)})
        alembic_runner.report_progress("System points", min(first + CHUNK - low, total), total)


def _cell_keys(connection):
    """The cell key (region queries): a generated column and a covering index on every table that holds a point."""
    for table, columns, expressions in CELL_TABLES:
        if _has_column(connection, table, "cell_key"):
            continue
        if table == "star_systems":
            _galaxy_points(connection)
        key = cellkey.sql_cell_key(*expressions)
        connection.execute(sa.text(
            f"ALTER TABLE {table} ADD COLUMN cell_key BIGINT UNSIGNED GENERATED ALWAYS AS ({key}) STORED,"
            f" ADD KEY idx_{table}_cell (cell_key, {', '.join(columns)})"))


def upgrade():
    connection = op.get_bind()
    if _has_column(connection, "galaxy_shape", "reach_shell_index"):
        return  # the last step below: everything else is done
    for table, index, unique in UID_TABLES:
        _renumber_column(connection, table, index, unique, "BINARY(12)", _new_uid_sql("uid"), f"Object IDs ({table})")
    _renumber_counters(connection)
    _scatter(connection)
    _renumber_sector_uids(connection)
    _cell_keys(connection)
    _reach(connection)
