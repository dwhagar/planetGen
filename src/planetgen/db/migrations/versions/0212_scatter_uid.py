"""An Object ID on every scatter row (GEN.1254).

Schema v212. `bright_stars` and `phenomenon_scatter` gain `uid` (`BINARY(12)`, NULL until the scatter gives it) with a unique index each. The scatter writers fill it as they write
(`db/scatter_uids.py`); a charted system built round a scatter star, and a built scatter phenomenon, keep it. A row written before this step has none, and the regenerated scatter
replaces it. A column or index already there is skipped, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0212"
down_revision = "0211"
branch_labels = None
depends_on = None

_TABLES = (("bright_stars", "uq_bright_stars_uid"), ("phenomenon_scatter", "uq_phenomenon_scatter_uid"))
"""tuple: `(table, unique index)` of the two scatter tables."""


def upgrade():
    connection = op.get_bind()
    for table, index in _TABLES:
        helpers.add_column_instant(connection, table, "uid", "BINARY(12) NULL")
        helpers.add_index_online(connection, table, index, ["uid"], unique=True)


def downgrade():
    connection = op.get_bind()
    for table, index in _TABLES:
        helpers.drop_index_online(connection, table, index)
        helpers.drop_column_instant(connection, table, "uid")
