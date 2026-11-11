"""The numbers a body's history leaves (GEN.1497).

Schema v243. One table, `body_history`: for a planet or moon, the eccentricity and the spin it was born with, the two uniforms that decided its spin resonance captures and
the version of the formulas that wrote the row (see `schema.sql` and `physics/body_history.py`). A galaxy made before v243 has no rows, which reads as 'not kept'; nothing is
backfilled (the draws of a body made earlier are not stored anywhere). A table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0243"
down_revision = "0242"
branch_labels = None
depends_on = None

HISTORY = """
    CREATE TABLE IF NOT EXISTS body_history (
        star_system_id      BIGINT UNSIGNED NOT NULL,
        body_type           TINYINT UNSIGNED NOT NULL,
        body_id             BIGINT UNSIGNED NOT NULL,
        birth_eccentricity  DOUBLE NULL CHECK (birth_eccentricity IS NULL OR birth_eccentricity >= 0),
        birth_spin_rate     DOUBLE NULL,
        birth_obliquity_deg DOUBLE NULL CHECK (birth_obliquity_deg IS NULL OR birth_obliquity_deg BETWEEN 0 AND 180),
        capture_draw_1      FLOAT NULL CHECK (capture_draw_1 IS NULL OR capture_draw_1 BETWEEN 0 AND 1),
        capture_draw_2      FLOAT NULL CHECK (capture_draw_2 IS NULL OR capture_draw_2 BETWEEN 0 AND 1),
        formula_version     TINYINT UNSIGNED NOT NULL,
        PRIMARY KEY (star_system_id, body_type, body_id),
        CONSTRAINT fk_body_history_star_system
            FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def upgrade():
    op.get_bind().execute(sa.text(HISTORY))
