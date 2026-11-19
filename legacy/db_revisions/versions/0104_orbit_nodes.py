"""Built-system star storage: the orbit tree, the member columns of a star and the disk descriptor (DB.53).

Schema v104. A built multi-star system had nowhere to keep its tree of orbits (`scatter_nodes` holds the scatter side
only), and the items of the multi-star train (GEN.242, GEN.264, GEN.269, GEN.273, GEN.274, GEN.275) each assumed columns
on tables that did not exist. This is the one revision for all of them:

* `orbit_nodes`, with the names and units of `scatter_nodes` (the key is the system, `star_system_id`, instead of the
  scatter `source` and `source_id`) so a fill copies a node column for column, and the disk descriptor
  (`disk_kind`, `disk_mass_solar`, `disk_inner_au`, `disk_outer_au`, nullable, no bodies in it).
* On `stars`: `parent_node_id` (the node a star sits in), the member columns of `scatter_members` (`state`, `flags`,
  `offset_x/y/z_au`, `velocity_x/y/z_kms`, `decay_age_yr`, `epoch_unix`) and the same disk descriptor.
* On `planets`: `host_node_id` (the node a circumbinary planet orbits; NULL for its star).

The existing binaries convert: a `close` or `wide` pair becomes one node (number 0, no parent) from its stored mutual
orbit, its two stars take `parent_node_id` 0, and every star's offset from its system's centre is the stored binary
position (zero for a single star). A star's velocity and epoch stay NULL until the code that moves members fills them.
Each step is skipped when it is done, so a stopped run can be started again; `planets.host_node_id` goes last and marks
the revision done.
"""

import sqlalchemy as sa
from alembic import op

revision = "0104"
down_revision = "0103"
branch_labels = None
depends_on = None

AU_KM = 149597870.7
SOLAR_MASS_KG = 1.988409870698051e30

CREATE_ORBIT_NODES = """
CREATE TABLE IF NOT EXISTS orbit_nodes (
    star_system_id       BIGINT UNSIGNED NOT NULL,
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
    disk_kind            VARCHAR(16),
    disk_mass_solar      DOUBLE,
    disk_inner_au        DOUBLE,
    disk_outer_au        DOUBLE,
    PRIMARY KEY (star_system_id, node_no),
    CONSTRAINT fk_orbit_nodes_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT chk_orbit_nodes_kind CHECK (kind IN ('close', 'wide', 'compact', 'group')),
    CONSTRAINT chk_orbit_nodes_interaction CHECK (interaction IS NULL OR interaction IN (
        'detached', 'semi-detached', 'contact', 'accreting', 'post-common-envelope', 'merger-product'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
"""

DISK = (("disk_kind", "VARCHAR(16) NULL"), ("disk_mass_solar", "DOUBLE NULL"), ("disk_inner_au", "DOUBLE NULL"),
        ("disk_outer_au", "DOUBLE NULL"))
STAR_COLUMNS = (
    ("parent_node_id", "SMALLINT NULL"), ("state", "VARCHAR(16) NULL"), ("flags", "SMALLINT UNSIGNED NULL"),
    ("offset_x_au", "DOUBLE NULL"), ("offset_y_au", "DOUBLE NULL"), ("offset_z_au", "DOUBLE NULL"),
    ("velocity_x_kms", "DOUBLE NULL"), ("velocity_y_kms", "DOUBLE NULL"), ("velocity_z_kms", "DOUBLE NULL"),
    ("decay_age_yr", "DOUBLE NULL"), ("epoch_unix", "DOUBLE NULL"),
) + DISK
PLANET_COLUMNS = (("host_node_id", "SMALLINT NULL"),)
"""Added last, on `planets`: the marker of the revision."""


def _has_column(connection, table, column):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar())


def _add_columns(connection, table, columns):
    missing = [(name, definition) for name, definition in columns if not _has_column(connection, table, name)]
    if missing:
        op.execute(f"ALTER TABLE {table} " + ", ".join(f"ADD COLUMN {name} {definition}" for name, definition in missing))


def upgrade():
    connection = op.get_bind()
    op.execute(CREATE_ORBIT_NODES)
    _add_columns(connection, "stars", STAR_COLUMNS)

    # A close or wide pair is one node of the tree: its mutual orbit as stored, the mass of its two stars.
    op.execute(
        "INSERT IGNORE INTO orbit_nodes (star_system_id, node_no, parent_node_no, slot, mass_solar, separation_km,"
        " eccentricity, inclination_deg, ascending_node_deg, arg_periapsis_deg, mean_anomaly_deg, epoch_unix,"
        " period_years, kind, interaction, flags)"
        " SELECT ss.id, 0, NULL, 0,"
        f" COALESCE((SELECT SUM(s.mass_kg) FROM stars s WHERE s.star_system_id = ss.id), 0) / {SOLAR_MASS_KG!r},"
        " COALESCE(ss.binary_separation_km, 0), COALESCE(ss.binary_eccentricity, 0),"
        " COALESCE(ss.binary_mutual_orbital_inclination_deg, 0), COALESCE(ss.binary_mutual_orbital_ascending_node_deg, 0),"
        " 0, COALESCE(ss.binary_mutual_orbital_phase_deg, 0), COALESCE(ss.binary_epoch_unix, 0),"
        " COALESCE(ss.binary_mutual_orbital_period_years, 0), ss.binary_configuration, NULL, 0"
        " FROM star_systems ss WHERE ss.binary_configuration IN ('close', 'wide')")
    op.execute(
        "UPDATE stars s JOIN star_systems ss ON ss.id = s.star_system_id SET s.parent_node_id = 0"
        " WHERE ss.binary_configuration IN ('close', 'wide') AND s.role IN ('primary', 'secondary')"
        " AND s.parent_node_id IS NULL")
    # Each star's offset from its system's centre: the stored binary position, nothing for a single star.
    for role, column in (("primary", "primary"), ("secondary", "secondary")):
        axes = ", ".join(
            f"s.offset_{axis}_au = COALESCE(ss.binary_{column}_position_{axis}_km, 0) / {AU_KM!r}" for axis in "xyz")
        op.execute(
            f"UPDATE stars s JOIN star_systems ss ON ss.id = s.star_system_id SET {axes}"
            f" WHERE s.role = '{role}' AND s.offset_x_au IS NULL")
    op.execute("UPDATE stars SET offset_x_au = 0, offset_y_au = 0, offset_z_au = 0 WHERE offset_x_au IS NULL")

    _add_columns(connection, "planets", PLANET_COLUMNS)
