"""The adaptations the life of a world needed (GEN.1809).

Schema v274 (a placeholder: the train that merges this renumbers it to the next free revision). `planet_adaptations` and `moon_adaptations` are new: one row for each thing the life of a planet or moon had to
become able to do to live where it does (`physics/adaptation_record.py`, GEN.1808): the order of the rows from 1, the axis of the table of adaptations it is on (`t_cold`, `t_hot`, `light_floor`,
`water_hours`, `ground_water_depth`, `water_activity`; none for a habitat row such as a pigment that works in no light), the habitat it comes from, what the places demand on the axis (in the unit of the axis),
the lowest grade (2 or 3) of the table that reaches it (none when even the extreme grade falls short, or for a habitat row) and the ids of the rows of the table that give it, comma separated. A body stored before
has no rows, and nothing is backfilled (a migration does not draw): the next regeneration makes them. A table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0274"
down_revision = "0273"
branch_labels = None
depends_on = None

TABLE = """
CREATE TABLE {prefix}_adaptations (
    id           BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    {prefix}_id   BIGINT UNSIGNED NOT NULL,
    entry_number SMALLINT UNSIGNED NOT NULL CHECK (entry_number >= 1),
    axis         VARCHAR(24) NULL,
    habitat      VARCHAR(24) NOT NULL,
    demand       DOUBLE NULL,
    grade        TINYINT UNSIGNED NULL CHECK (grade IS NULL OR grade BETWEEN 2 AND 3),
    row_ids      VARCHAR(64) NOT NULL,
    CONSTRAINT fk_{prefix}_adaptations_{prefix}
        FOREIGN KEY ({prefix}_id) REFERENCES {prefix}s(id) ON DELETE CASCADE,
    UNIQUE KEY uq_{prefix}_adaptations ({prefix}_id, entry_number),
    KEY idx_{prefix}_adaptations_axis (axis, grade)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
"""


def upgrade():
    connection = op.get_bind()
    for prefix in ("planet", "moon"):
        if not helpers.exists_table(connection, f"{prefix}_adaptations"):
            connection.execute(sa.text(TABLE.format(prefix=prefix)))


def downgrade():
    connection = op.get_bind()
    for prefix in ("moon", "planet"):
        if helpers.exists_table(connection, f"{prefix}_adaptations"):
            connection.execute(sa.text(f"DROP TABLE {prefix}_adaptations"))
