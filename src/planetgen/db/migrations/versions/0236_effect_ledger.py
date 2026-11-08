"""The effect ledger (GEN.1465).

Schema v236. One table, `effect_ledger`: for every body of a system and every effect that acts on it, the share of each source's contribution in the largest source's (see `schema.sql`
and `physics/effect_ledger.py`). A galaxy made before v236 has no rows, which reads as 'no ledger'; nothing is backfilled. A table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0236"
down_revision = "0235"
branch_labels = None
depends_on = None

LEDGER = """
    CREATE TABLE IF NOT EXISTS effect_ledger (
        star_system_id  BIGINT UNSIGNED NOT NULL,
        body_type       TINYINT UNSIGNED NOT NULL,
        body_id         BIGINT UNSIGNED NOT NULL,
        effect          TINYINT UNSIGNED NOT NULL,
        source_type     TINYINT UNSIGNED NOT NULL,
        source_id       BIGINT UNSIGNED NOT NULL,
        share           FLOAT NOT NULL CHECK (share >= 0),
        n_below         INT UNSIGNED NULL,
        largest_value   DOUBLE NULL CHECK (largest_value IS NULL OR largest_value > 0),
        PRIMARY KEY (star_system_id, body_type, body_id, effect, source_type, source_id),
        CONSTRAINT fk_effect_ledger_star_system
            FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def upgrade():
    op.get_bind().execute(sa.text(LEDGER))
