"""How much a class was fudged to fit a body (GEN.1314).

Schema v244. `planets`, `moons` and `rogue_planets` gain `class_fudge_values` (VARCHAR(120), the values fudged to fit the stored class, comma separated quantity names, the
temperatures as `temperature_k`; NULL for a full fit) and `class_fudge_amount` (DOUBLE, the total move of those values in widths of the fudge table, `physics/data/fudge_table.json`;
0 for a full fit). Both are NULL here: a row made before v244 gets them when the galaxy is regenerated, because the class is read from the finished body and a stored letter is not
rescored by a migration. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0244"
down_revision = "0243"
branch_labels = None
depends_on = None

COLUMNS = (("class_fudge_values", "VARCHAR(120) NULL"), ("class_fudge_amount", "DOUBLE NULL"))
"""tuple: The columns added to `planets`, `moons` and `rogue_planets`, each with its definition."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons", "rogue_planets"):
        for column, definition in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} {definition}"))
