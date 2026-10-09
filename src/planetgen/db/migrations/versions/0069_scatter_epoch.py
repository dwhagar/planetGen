"""The time a phenomenon scatter's positions hold at (GEN.137).

Schema v69. `phenomenon_scatter` gains `epoch_unix`: the database's orbit
epoch (`orbit_simulation_state.last_updated_at`, Unix seconds) when
`planetgen plan` drew the row, so a hypervelocity star's position at a later
time is `p0 + v (t - epoch_unix)` (`galaxy/straight_line.py`). Rows already
there keep NULL: their plan time is unknown, so they hold at whatever the
database's orbit epoch is.
"""

import sqlalchemy as sa
from alembic import op

revision = "0069"
down_revision = "0068"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    present = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.columns WHERE table_schema = DATABASE()"
        " AND table_name = 'phenomenon_scatter' AND column_name = 'epoch_unix'")).scalar()
    if not present:
        connection.execute(sa.text("ALTER TABLE phenomenon_scatter ADD COLUMN epoch_unix DOUBLE"))
