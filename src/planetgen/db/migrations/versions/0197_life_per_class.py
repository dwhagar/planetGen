"""Several life classes on one body (GEN.1209).

Schema v197. `planet_life` and `moon_life` are new: one row for each life class that lives on a planet or moon (the class id, its share of the food
production, the stage its own schedule reached, its animal grade and complexity, the five ages of its schedule and the four kind scores it gives). A body
stored before has no rows: its one class is the dominant one in the body's own columns, and nothing is backfilled (the next regeneration makes the rows).
A table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0197"
down_revision = "0196"
branch_labels = None
depends_on = None

TABLE = """
CREATE TABLE {prefix}_life (
    id                                      BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    {prefix}_id                               BIGINT UNSIGNED NOT NULL,
    life_class                              VARCHAR(2) NOT NULL,
    share                                   DOUBLE NOT NULL CHECK (share > 0 AND share <= 1),
    life_stage                              VARCHAR(32) NULL,
    life_milestone_age_gyr                  DOUBLE NULL,
    life_system_age_gyr                     DOUBLE NULL,
    life_animal_grade                       TINYINT UNSIGNED NULL,
    life_complexity                         DOUBLE NULL,
    life_age_abiogenesis_gyr                DOUBLE NULL,
    life_age_photosynthesis_gyr             DOUBLE NULL,
    life_age_complex_cells_gyr              DOUBLE NULL,
    life_age_multicellularity_gyr           DOUBLE NULL,
    life_age_technological_civilization_gyr DOUBLE NULL,
    life_microbial                          TINYINT UNSIGNED NULL,
    life_flora                              TINYINT UNSIGNED NULL,
    life_fauna                              TINYINT UNSIGNED NULL,
    life_sentient                           TINYINT UNSIGNED NULL,
    CONSTRAINT fk_{prefix}_life_{prefix}
        FOREIGN KEY ({prefix}_id) REFERENCES {prefix}s(id) ON DELETE CASCADE,
    UNIQUE KEY uq_{prefix}_life ({prefix}_id, life_class),
    KEY idx_{prefix}_life_class (life_class, life_stage)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
"""


def upgrade():
    connection = op.get_bind()
    for prefix in ("planet", "moon"):
        if not helpers.exists_table(connection, f"{prefix}_life"):
            connection.execute(sa.text(TABLE.format(prefix=prefix)))
