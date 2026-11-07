"""How an asteroid field spreads and where it came from (GEN.1103).

Schema v221. `asteroid_fields` gains `rms_speed_kms` (DOUBLE), `birth_age_yr` (DOUBLE), `shape` (VARCHAR(8)), `origin` (VARCHAR(32)), `history` (TEXT) and `largest_body_kg` (DOUBLE), all
nullable: the rms speed of the field's pieces, its age at the sector's epoch (its radius follows from the two, `field_dispersal.rms_radius_m`), what shape it was made as, how it came to be, the
sentence of where it came from and the mass of its largest piece. A row made before v221 has them NULL and keeps the radius it has; nothing is backfilled. Each column is skipped when it is
already there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0221"
down_revision = "0220"
branch_labels = None
depends_on = None

COLUMNS = (
    ("rms_speed_kms", "DOUBLE NULL"),
    ("birth_age_yr", "DOUBLE NULL"),
    ("shape", "VARCHAR(8) NULL"),
    ("origin", "VARCHAR(32) NULL"),
    ("history", "TEXT NULL"),
    ("largest_body_kg", "DOUBLE NULL"),
)


def upgrade():
    connection = op.get_bind()
    for name, definition in COLUMNS:
        if not helpers.exists_column(connection, "asteroid_fields", name):
            connection.execute(sa.text(f"ALTER TABLE asteroid_fields ADD COLUMN {name} {definition}, ALGORITHM=INSTANT"))
