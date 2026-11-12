"""A free planet's moons are `moons` rows, so the one-bit flag is gone (GEN.558).

Schema v244. `rogue_planets.has_moons` was one boolean drawn at 20 percent for a brown dwarf and NULL for a row of the physics rogue (migration 0207), whose moons are `moons` rows the ordinary
moon routine makes inside the outer edge the ejection left (`generation/phenomena/starless_host.py`). The brown dwarf no longer draws it, so nothing writes it and nothing reads it, and the column is
dropped. A column already gone is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0244"
down_revision = "0243"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_column(connection, "rogue_planets", "has_moons"):
        return
    try:
        helpers.drop_column_instant(connection, "rogue_planets", "has_moons")
    except Exception:  # a column check on a fresh schema can refuse the instant form; the plain drop takes it with the column
        connection.execute(sa.text("ALTER TABLE rogue_planets DROP COLUMN has_moons"))
