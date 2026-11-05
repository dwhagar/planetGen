"""A planet or moon keeps no description text (GEN.1212).

Schema v208. `planets` and `moons` lose `description`, the sentence of the body's class that generation copied from `tuning.PLANET_CLASSES` into every
row. The words are now read at view time from `physics/data/planet_class_text.json` (`physics/planet_class_text.py`), so one change to the
wording reaches every stored body. A column already gone is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0208"
down_revision = "0207"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        if helpers.exists_column(connection, table, "description"):
            connection.execute(sa.text(f"ALTER TABLE {table} DROP COLUMN description"))
