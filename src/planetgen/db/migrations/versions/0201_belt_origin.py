"""Where a belt came from and its largest body, on every asteroid belt (GEN.1069).

Schema v201. `asteroid_belts` gain what a belt knows about its origin: `origin_cause` (what formed it), `co_causes` (JSON list), `formation_snow_line_au`, `primordial_mass_earth`,
`cleared_mass_earth`, `scattered_mass_earth`, `stirred`, `origin_history` (JSON list of sentences), `origin_sentence` and `scatter_fates` (JSON object), and `largest_body_km`, the largest body
of the belt (Dc) settled with its mass (`physics/belt_ledger.py`, GEN.1000).

NULL here: a belt made before v201 has no origin to read back (the generator epoch makes it), so the regenerated galaxy replaces it. A column already there is skipped, so a stopped run
can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0201"
down_revision = "0200"
branch_labels = None
depends_on = None

_COLUMNS = (
    ("origin_cause", "VARCHAR(24) NULL"),
    ("co_causes", "TEXT NULL"),
    ("formation_snow_line_au", "DOUBLE NULL"),
    ("primordial_mass_earth", "DOUBLE NULL"),
    ("cleared_mass_earth", "DOUBLE NULL"),
    ("scattered_mass_earth", "DOUBLE NULL"),
    ("stirred", "TINYINT(1) NULL"),
    ("origin_history", "TEXT NULL"),
    ("origin_sentence", "TEXT NULL"),
    ("scatter_fates", "TEXT NULL"),
    ("largest_body_km", "DOUBLE NULL"),
)
"""tuple: The `(column, definition)` added to `asteroid_belts`."""


def upgrade():
    connection = op.get_bind()
    for column, definition in _COLUMNS:
        if not helpers.exists_column(connection, "asteroid_belts", column):
            connection.execute(sa.text(f"ALTER TABLE asteroid_belts ADD COLUMN {column} {definition}"))
