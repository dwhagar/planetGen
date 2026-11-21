"""The life ages of a world (GEN.1746 slice 3).

Schema v273 (a placeholder: the train that merges this renumbers it to the next free revision). `planet_life_ages` and `moon_life_ages` are new: one row for each life age of a planet or moon, the stretch of its life
between two mass extinctions (`physics/life_ages.py`): its number from 1, its start and end age, the cause of the extinction that ended it and the share of species the extinction killed. The age the world is in
has no end, cause or kill. A body stored before has no rows, and nothing is backfilled (a migration does not draw): the next regeneration makes them. A table already there is skipped, so a stopped run can be
started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0273"
down_revision = "0272"
branch_labels = None
depends_on = None

TABLE = """
CREATE TABLE {prefix}_life_ages (
    id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    {prefix}_id  BIGINT UNSIGNED NOT NULL,
    age_number  SMALLINT UNSIGNED NOT NULL CHECK (age_number >= 1),
    start_gyr   DOUBLE NOT NULL CHECK (start_gyr >= 0),
    end_gyr     DOUBLE NULL,
    cause       VARCHAR(16) NULL,
    kill_share  DOUBLE NULL CHECK (kill_share IS NULL OR (kill_share > 0 AND kill_share <= 1)),
    CONSTRAINT fk_{prefix}_life_ages_{prefix}
        FOREIGN KEY ({prefix}_id) REFERENCES {prefix}s(id) ON DELETE CASCADE,
    UNIQUE KEY uq_{prefix}_life_ages ({prefix}_id, age_number),
    KEY idx_{prefix}_life_ages_number (age_number, end_gyr)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
"""


def upgrade():
    connection = op.get_bind()
    for prefix in ("planet", "moon"):
        if not helpers.exists_table(connection, f"{prefix}_life_ages"):
            connection.execute(sa.text(TABLE.format(prefix=prefix)))


def downgrade():
    connection = op.get_bind()
    for prefix in ("moon", "planet"):
        if helpers.exists_table(connection, f"{prefix}_life_ages"):
            connection.execute(sa.text(f"DROP TABLE {prefix}_life_ages"))
