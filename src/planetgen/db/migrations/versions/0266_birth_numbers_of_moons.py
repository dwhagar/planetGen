"""Where a body formed and the moons the history took away (GEN.1753).

Schema v266. `body_history` gains `formation_distance_au` (where the body formed, which keys its end-stage draw streams), `birth_distance_au` (the orbit a moon was born on) and `born_index` (its place in its
planet's born moons), and a table `lost_moons` keeps the moons a planet was born with that its history took away (mass, birth orbit, density, direction), so the history of a restored planet runs again from
what it started with (`physics/body_history.py`, `physics/moon_history.py`). A galaxy made earlier has NULL in the new columns and no lost moons, which reads as 'not kept'; nothing is backfilled (the numbers
of a body made earlier are stored nowhere). A column or table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0266"
down_revision = "0265"
branch_labels = None
depends_on = None

COLUMNS = (
    ("formation_distance_au", "DOUBLE NULL CHECK (formation_distance_au IS NULL OR formation_distance_au >= 0)", "capture_draw_2"),
    ("birth_distance_au", "DOUBLE NULL CHECK (birth_distance_au IS NULL OR birth_distance_au >= 0)", "formation_distance_au"),
    ("born_index", "SMALLINT UNSIGNED NULL", "birth_distance_au"),
)

LOST = """
    CREATE TABLE IF NOT EXISTS lost_moons (
        star_system_id      BIGINT UNSIGNED NOT NULL,
        planet_id           BIGINT UNSIGNED NOT NULL,
        born_index          SMALLINT UNSIGNED NOT NULL,
        mass_kg             DOUBLE NOT NULL CHECK (mass_kg > 0),
        birth_distance_au   DOUBLE NOT NULL CHECK (birth_distance_au > 0),
        density_kg_m3       DOUBLE NOT NULL CHECK (density_kg_m3 > 0),
        retrograde          TINYINT UNSIGNED NOT NULL CHECK (retrograde IN (0, 1)),
        PRIMARY KEY (star_system_id, planet_id, born_index),
        CONSTRAINT fk_lost_moons_star_system
            FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def _has_column(name):
    return bool(op.get_bind().execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'body_history' AND COLUMN_NAME = :name"), {"name": name}).scalar())


def upgrade():
    for name, definition, after in COLUMNS:
        if not _has_column(name):
            op.execute(f"ALTER TABLE body_history ADD COLUMN {name} {definition} AFTER {after}")
    op.get_bind().execute(sa.text(LOST))


def downgrade():
    op.execute("DROP TABLE IF EXISTS lost_moons")
    for name, _definition, _after in reversed(COLUMNS):
        if _has_column(name):
            op.execute(f"ALTER TABLE body_history DROP COLUMN {name}")
