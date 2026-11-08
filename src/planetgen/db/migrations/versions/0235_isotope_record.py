"""The isotope record of a body (GEN.1377, GEN.1420).

Schema v235. Three child tables, `planet_isotopes`, `moon_isotopes` and `rogue_planet_isotopes`, hold one row per isotope of each element a body holds (`element`, `isotope`, `fraction`, the
atom fraction of the element, and `origin`: natural, physical, decay or radiogenic). A body saved before v235 has no rows, which reads as 'no record'; nothing is backfilled. A table already there is
skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0235"
down_revision = "0234"
branch_labels = None
depends_on = None

TABLES = {"planets": "planet", "moons": "moon", "rogue_planets": "rogue_planet"}

CHILD = """
    CREATE TABLE IF NOT EXISTS {prefix}_isotopes (
        id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
        {prefix}_id  BIGINT UNSIGNED NOT NULL,
        element     VARCHAR(3) NOT NULL,
        isotope     VARCHAR(8) NOT NULL,
        fraction    DOUBLE NOT NULL CHECK (fraction >= 0 AND fraction <= 1),
        origin      VARCHAR(12) NOT NULL,
        UNIQUE KEY uq_{prefix}_isotopes ({prefix}_id, isotope),
        CONSTRAINT fk_{prefix}_isotopes_{prefix}
            FOREIGN KEY ({prefix}_id) REFERENCES {table}(id) ON DELETE CASCADE
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def upgrade():
    connection = op.get_bind()
    for table, prefix in TABLES.items():
        connection.execute(sa.text(CHILD.format(prefix=prefix, table=table)))
