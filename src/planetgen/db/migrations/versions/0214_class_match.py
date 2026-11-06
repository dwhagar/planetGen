"""How well the class fits a body (GEN.1207).

Schema v214. `planets`, `moons` and `rogue_planets` gain `class_match` (DOUBLE, the percent match of the stored class as a fraction, 0 when no class fits the body and
the stored one is the nearest), `runner_up_class` (CHAR(1), the next best class) and `runner_up_match` (DOUBLE). All three are NULL here: a row made before v214 gets them
when the galaxy is regenerated, because the best-fit classifier reads the finished body and a stored letter is not rescored by a migration. A column already there is
skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0214"
down_revision = "0213"
branch_labels = None
depends_on = None

COLUMNS = (("class_match", "DOUBLE NULL"), ("runner_up_class", "CHAR(1) NULL"), ("runner_up_match", "DOUBLE NULL"))
"""tuple: The columns added to `planets`, `moons` and `rogue_planets`, each with its definition."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons", "rogue_planets"):
        for column, definition in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} {definition}"))
