"""The scattered mass (MAP.165).

Schema v80. `phenomenon_scatter` gains `mass_solar`: the solar mass a scattered
black hole or neutron star was drawn with (NULL for the other kinds and for rows
scattered before v80, which the map sizes by class).
"""

import sqlalchemy as sa
from alembic import op

revision = "0080"
down_revision = "0079"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    present = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'phenomenon_scatter' AND COLUMN_NAME = 'mass_solar'")).scalar()
    if not present:
        op.execute("ALTER TABLE phenomenon_scatter ADD COLUMN mass_solar DOUBLE AFTER epoch_unix")
