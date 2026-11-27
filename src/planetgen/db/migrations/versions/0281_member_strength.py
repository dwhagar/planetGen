"""The roundness ratio, yield strength and binding energy of a belt or field member (GEN.1159), and the tilt and the fallen moons a planet keeps (GEN.1874).

Schema v281 (a placeholder: the train that merges this renumbers it to the next free revision). `belt_members` and `field_members` gain `roundness_ratio` (DOUBLE, the mean radius over the round
radius, 1 or more is round), `yield_strength_pa` (DOUBLE, Pa) and `binding_energy_j` (DOUBLE, J), worked out by `physics/body_strength.py` from the member's own mass, size and makeup at the cold
strength, the same code a comet goes through. All three are NULL here: a belt or a field made before gets them when its sector is regenerated, and a migration does not draw.

GEN.1874: `body_history` gains `tilt_chaotic` (1 when the planet's tilt wanders chaotically, the quantity oddities W6 and W61 read) and a table `fallen_moons` keeps the moons a planet's tide brought inside its
Roche limit with the years after the planet formed that each fell, which its rings are made from. A galaxy made earlier has NULL and no rows, which reads as 'not kept'. A column or table already there is
skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0281"
down_revision = "0280"
branch_labels = None
depends_on = None

TABLES = ("belt_members", "field_members")
"""tuple: The tables that gain the columns."""

COLUMNS = (("roundness_ratio", "DOUBLE NULL"), ("yield_strength_pa", "DOUBLE NULL"), ("binding_energy_j", "DOUBLE NULL"))
"""tuple: The columns added to each table, each with its definition."""


HISTORY_COLUMN = ("tilt_chaotic", "TINYINT UNSIGNED NULL CHECK (tilt_chaotic IS NULL OR tilt_chaotic IN (0, 1))", "forced_eccentricity")
"""tuple: The column `body_history` gains (GEN.1874): its name, definition and the column it goes after."""

FALLEN = """
    CREATE TABLE IF NOT EXISTS fallen_moons (
        star_system_id      BIGINT UNSIGNED NOT NULL,
        planet_id           BIGINT UNSIGNED NOT NULL,
        fall_index          SMALLINT UNSIGNED NOT NULL,
        mass_kg             DOUBLE NOT NULL CHECK (mass_kg > 0),
        time_yr             DOUBLE NOT NULL CHECK (time_yr >= 0),
        PRIMARY KEY (star_system_id, planet_id, fall_index),
        CONSTRAINT fk_fallen_moons_star_system
            FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        for column, definition in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} {definition}"))
    name, definition, after = HISTORY_COLUMN
    if not helpers.exists_column(connection, "body_history", name):
        connection.execute(sa.text(f"ALTER TABLE body_history ADD COLUMN {name} {definition} AFTER {after}"))
    connection.execute(sa.text(FALLEN))


def downgrade():
    connection = op.get_bind()
    connection.execute(sa.text("DROP TABLE IF EXISTS fallen_moons"))
    if helpers.exists_column(connection, "body_history", HISTORY_COLUMN[0]):
        connection.execute(sa.text(f"ALTER TABLE body_history DROP COLUMN {HISTORY_COLUMN[0]}"))
    for table in reversed(TABLES):
        for column, _definition in reversed(COLUMNS):
            if helpers.exists_column(connection, table, column):
                connection.execute(sa.text(f"ALTER TABLE {table} DROP COLUMN {column}"))
