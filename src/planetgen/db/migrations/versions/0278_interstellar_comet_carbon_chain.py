"""Whether an interstellar comet was drawn depleted in carbon chains (GEN.1227).

Schema v278 (the number of the revision the train that merges this takes if one has gone in before it). `interstellar_comets` gains `carbon_chain_depleted` (TINYINT(1): the nucleus was drawn with little CN
and C2, so its coma has no green), the flag of the volatiles the nucleus was drawn with, which the coma rows of v275 do not hold (they hold the rates of the gases, from which the coma class and the dominant
gases are read, `physics/comet_coma.py`). It is NULL for a row made before v278: a migration does not draw and the volatiles of an old row were never stored. A column already there is skipped, so a
stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0278"
down_revision = "0277"
branch_labels = None
depends_on = None

COLUMN = "carbon_chain_depleted"
"""str: The column `interstellar_comets` gains."""

DEFINITION = "TINYINT(1) NULL CHECK (carbon_chain_depleted IS NULL OR carbon_chain_depleted IN (0, 1))"
"""str: Its definition."""


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_column(connection, "interstellar_comets", COLUMN):
        connection.execute(sa.text(f"ALTER TABLE interstellar_comets ADD COLUMN {COLUMN} {DEFINITION}"))


def downgrade():
    connection = op.get_bind()
    if helpers.exists_column(connection, "interstellar_comets", COLUMN):
        connection.execute(sa.text(f"ALTER TABLE interstellar_comets DROP COLUMN {COLUMN}"))
