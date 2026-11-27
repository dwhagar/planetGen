"""The exosphere of an airless body (GEN.455, part B).

Schema v280 (a placeholder: the train that merges this renumbers it to the next free revision). `planet_exosphere` and `moon_exosphere` are new: one row for each species a source of the surface supplies to a body
with no air at all (`physics/exosphere.py`): the order of the rows from 1 (the largest column first), the species (the gas table's formula), the column in atoms per square centimetre, the source that gives most of it
and the scale height, km. A body stored before has no rows, and nothing is backfilled (a migration does not draw): the next regeneration makes them. A table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0280"
down_revision = "0279"
branch_labels = None
depends_on = None

TABLE = """
CREATE TABLE {prefix}_exosphere (
    id              BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    {prefix}_id     BIGINT UNSIGNED NOT NULL,
    entry_number    SMALLINT UNSIGNED NOT NULL CHECK (entry_number >= 1),
    species         VARCHAR(16) NOT NULL,
    column_cm2      DOUBLE NOT NULL CHECK (column_cm2 > 0),
    source          VARCHAR(40) NOT NULL,
    scale_height_km DOUBLE NOT NULL,
    CONSTRAINT fk_{prefix}_exosphere_{prefix}
        FOREIGN KEY ({prefix}_id) REFERENCES {prefix}s(id) ON DELETE CASCADE,
    UNIQUE KEY uq_{prefix}_exosphere ({prefix}_id, entry_number),
    KEY idx_{prefix}_exosphere_species (species, column_cm2)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
"""


def upgrade():
    connection = op.get_bind()
    for prefix in ("planet", "moon"):
        if not helpers.exists_table(connection, f"{prefix}_exosphere"):
            connection.execute(sa.text(TABLE.format(prefix=prefix)))


def downgrade():
    connection = op.get_bind()
    for prefix in ("moon", "planet"):
        if helpers.exists_table(connection, f"{prefix}_exosphere"):
            connection.execute(sa.text(f"DROP TABLE {prefix}_exosphere"))
