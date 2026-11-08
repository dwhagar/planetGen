"""The crust of a rocky body and the magmatic ore types it was assessed for (GEN.1380).

Schema v231. `planets`, `moons` and `rogue_planets` gain `crust_thickness_km`, `craton_fraction` and `crust_sulfur_index` (DOUBLE, nullable), and three child tables,
`planet_deposit_types`, `moon_deposit_types` and `rogue_planet_deposit_types`, hold one row per magmatic ore type a body was assessed for (`deposit_type`, `status` present,
absent_by_state or absent_by_chance, `endowment_kg`, `grade`, `key_element`, the age window and whether the host capped the grade). A body saved before v228 has NULL columns and no rows,
which reads as 'not assessed'; nothing is backfilled. A column or table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0231"
down_revision = "0230"
branch_labels = None
depends_on = None

TABLES = {"planets": "planet", "moons": "moon", "rogue_planets": "rogue_planet"}
COLUMNS = ("crust_thickness_km", "craton_fraction", "crust_sulfur_index")

CHILD = """
    CREATE TABLE IF NOT EXISTS {prefix}_deposit_types (
        id                    BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
        {prefix}_id           BIGINT UNSIGNED NOT NULL,
        deposit_type          VARCHAR(24) NOT NULL,
        status                VARCHAR(20) NOT NULL,
        endowment_kg          DOUBLE,
        grade                 DOUBLE,
        key_element           VARCHAR(8) NOT NULL,
        age_window_gyr_start  DOUBLE,
        age_window_gyr_end    DOUBLE,
        capped_by_host        TINYINT NOT NULL DEFAULT 0,
        UNIQUE KEY uq_{prefix}_deposit_types ({prefix}_id, deposit_type),
        CONSTRAINT fk_{prefix}_deposit_types_{prefix}
            FOREIGN KEY ({prefix}_id) REFERENCES {table}(id) ON DELETE CASCADE
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def upgrade():
    connection = op.get_bind()
    for table, prefix in TABLES.items():
        missing = [name for name in COLUMNS if not helpers.exists_column(connection, table, name)]
        if missing:
            # No ALGORITHM clause: `planets` and `moons` carry a FULLTEXT key, which refuses ALGORITHM=INSTANT; one rebuild adds all of the columns.
            connection.execute(sa.text(f"ALTER TABLE {table} " + ", ".join(f"ADD COLUMN {name} DOUBLE NULL" for name in missing)))
        connection.execute(sa.text(CHILD.format(prefix=prefix, table=table)))
