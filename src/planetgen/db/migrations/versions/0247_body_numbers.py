"""The body number each slot of a system has had (DB.75).

Schema v247. One table, `body_numbers`: for a system, the slot path of a body (`store._slot_key`) and the number it was given, so a draw of the system's content again gives a body at a
remembered slot the same number (and so the same Object ID) it had, even when an earlier draw left the slot empty (see `schema.sql`). A galaxy made before v247 has no rows; nothing is
backfilled, because the first draw of a system writes the numbers its bodies hold then. A table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0247"
down_revision = "0246"
branch_labels = None
depends_on = None

NUMBERS = """
    CREATE TABLE IF NOT EXISTS body_numbers (
        star_system_id      BIGINT UNSIGNED NOT NULL,
        slot_key            VARCHAR(96) NOT NULL,
        body_number         SMALLINT UNSIGNED NOT NULL CHECK (body_number BETWEEN 1 AND 4095),
        PRIMARY KEY (star_system_id, slot_key),
        UNIQUE KEY uq_body_numbers_number (star_system_id, body_number),
        CONSTRAINT fk_body_numbers_star_system
            FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def upgrade():
    op.get_bind().execute(sa.text(NUMBERS))
