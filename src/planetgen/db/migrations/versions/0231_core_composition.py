"""The make-up of a metal core (GEN.1364).

Schema v231. `planets`, `moons` and `rogue_planets` gain `core_ni_fraction`, `core_s_fraction`, `core_si_fraction`, `core_o_fraction`, `core_c_fraction`, `core_h_fraction` (DOUBLE,
NULL: each element's mass fraction of the metal core, the iron is the rest) and `core_light_fraction` (DOUBLE, NULL: the five light elements together). A row made before v231 has NULL
in them and its core takes the Earth's reference mix until the galaxy is regenerated (a migration does not draw; `physics/core_composition.py`). A column already there is skipped,
so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0231"
down_revision = "0230"
branch_labels = None
depends_on = None

FIELDS = ('core_ni_fraction', 'core_light_fraction', 'core_s_fraction', 'core_si_fraction', 'core_o_fraction', 'core_c_fraction', 'core_h_fraction')
TABLES = ("planets", "moons", "rogue_planets")
COLUMNS = tuple((table, field, "DOUBLE NULL") for table in TABLES for field in FIELDS)
"""tuple: `(table, column, definition)` of every column the migration adds."""


def upgrade():
    connection = op.get_bind()
    for table, column, definition in COLUMNS:
        if not helpers.exists_column(connection, table, column):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} {definition}"))
