"""The coma of an interstellar comet (GEN.1227).

Schema v275 (a placeholder: the train that merges this renumbers it to the next free revision). `interstellar_comets` gains `coma_class` (VARCHAR(16): `CO2-dominated` or `CO-dominated`),
`coma_volatiles` (VARCHAR(48): the dominant volatiles after water, the most abundant first, as a comma-separated list of formulas) and `carbon_chain_depleted` (TINYINT(1): drawn with little CN and C2), read
off the volatiles the nucleus was drawn with (`physics/comet_coma.py`). They are NULL here: a row made before v275 gets them when the galaxy is regenerated (a migration does not draw, and the volatiles of an
old row were never stored). A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0275"
down_revision = "0274"
branch_labels = None
depends_on = None

COLUMNS = (
    ("coma_class", "VARCHAR(16) NULL"),
    ("coma_volatiles", "VARCHAR(48) NULL"),
    ("carbon_chain_depleted", "TINYINT(1) NULL CHECK (carbon_chain_depleted IS NULL OR carbon_chain_depleted IN (0, 1))"),
)
"""tuple: The columns added to `interstellar_comets`, each with its definition."""


def upgrade():
    connection = op.get_bind()
    for name, definition in COLUMNS:
        if not helpers.exists_column(connection, "interstellar_comets", name):
            connection.execute(sa.text(f"ALTER TABLE interstellar_comets ADD COLUMN {name} {definition}"))


def downgrade():
    connection = op.get_bind()
    for name, _definition in reversed(COLUMNS):
        if helpers.exists_column(connection, "interstellar_comets", name):
            connection.execute(sa.text(f"ALTER TABLE interstellar_comets DROP COLUMN {name}"))
