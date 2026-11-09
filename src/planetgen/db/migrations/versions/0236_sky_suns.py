"""The other suns in a planet's sky (GEN.1458).

Schema v236. `planets`, `moons` and `rogue_planets` gain `sky_suns` (VARCHAR(400), nullable): a planet's other suns as JSON `[[share, ux, uy, uz], ...]`, each other star of
its system's share of the mean light at the planet and its unit direction in the planet's equinox frame (`physics/sky_suns.py`). NULL for a lone star, a moon and a row made
before v236, which gets it when the galaxy is regenerated (a migration does not draw). A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0236"
down_revision = "0235"
branch_labels = None
depends_on = None

TABLES = ("planets", "moons", "rogue_planets")
"""tuple: The tables that gain the column."""


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        if not helpers.exists_column(connection, table, "sky_suns"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN sky_suns VARCHAR(400) NULL"))
