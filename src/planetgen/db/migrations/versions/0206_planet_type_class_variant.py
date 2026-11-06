"""The body type and the digit of the class code (GEN.1219, GEN.1220), and the envelope share of every body (GEN.1226).

Schema v206. `planets` and `moons` gain `planet_type` (VARCHAR(12): terrestrial, transitional or gas, from the envelope share of the mass,
`physics/planet_type.py`) and `class_variant` (TINYINT, 0 to 9, the digit of the class code; 0 is not written, `physics/class_code.py`). Both are
NULL here: a row made before v206 gets them when the galaxy is regenerated. `envelope_mass_fraction` was there already and now holds the true
share of a hydrogen-rich air (`physics/envelope_share.py`) for every rocky body. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0206"
down_revision = "0205"
branch_labels = None
depends_on = None

COLUMNS = (("planet_type", "VARCHAR(12) NULL"), ("class_variant", "TINYINT NULL"))
"""tuple: The columns added to `planets` and `moons`, each with its definition."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        for column, definition in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} {definition}"))
