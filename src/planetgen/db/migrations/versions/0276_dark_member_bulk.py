"""The bulk density and binding energy of a dark member of an orbit tree (GEN.1204).

Schema v276 (a placeholder: the train that merges this renumbers it to the next free revision). `dark_members` (a brown dwarf, a neutron star, a black hole or a stripped helium star that hangs in a
system's tree) gain `density_g_cm3` (DOUBLE, g/cm3) and `binding_energy_j` (DOUBLE, J; -1 for the infinite binding energy of a black hole), worked out by `physics/body_strength.py` from the member's
mass and the radius the tree drew it with, the same code the single brown dwarf, neutron star and black hole go through (GEN.1162). Both are NULL here: a system made before gets them when its sector is
regenerated, and a migration does not draw. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0276"
down_revision = "0275"
branch_labels = None
depends_on = None

COLUMNS = (("density_g_cm3", "DOUBLE NULL"), ("binding_energy_j", "DOUBLE NULL"))
"""tuple: The columns added to `dark_members`, each with its definition."""


def upgrade():
    connection = op.get_bind()
    for column, definition in COLUMNS:
        if not helpers.exists_column(connection, "dark_members", column):
            connection.execute(sa.text(f"ALTER TABLE dark_members ADD COLUMN {column} {definition}"))
