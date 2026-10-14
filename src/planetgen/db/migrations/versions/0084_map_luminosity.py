"""`phenomenon_scatter.map_luminosity_sol`: the luminosity a star of the row's mass would have on the map (DB.23).

Schema v84. The column is added empty; `planetgen plan --fill-map-luminosity` fills it for a galaxy scattered
before this version, and every new scatter writes it with the row. Skipped when the column is already there, so a
stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0084"
down_revision = "0083"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    present = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'phenomenon_scatter' AND COLUMN_NAME = 'map_luminosity_sol'")).scalar()
    if not present:
        op.execute("ALTER TABLE phenomenon_scatter ADD COLUMN map_luminosity_sol FLOAT NULL AFTER mass_solar")
