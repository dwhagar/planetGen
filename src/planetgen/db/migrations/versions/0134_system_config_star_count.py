"""A stored config keeps the number of stars it was forced to, the cap on it and the prevalence of five or more (GEN.267).

Schema v134. `system_configs` gains `stars` (TINYINT UNSIGNED, NULL: the forced star count, left to chance when NULL), `max_stars`
(TINYINT UNSIGNED, NOT NULL, default 12: the cap on the star count) and `prevalence_multi_star` (DOUBLE, NOT NULL, default 0: the
percent from the usual chance of five or more stars, as the other `prevalence_*` columns), each at the end of the table.

The defaults are what every config before this revision meant: no forced count, the cap of 12, the usual chance. So no row is touched
and nothing is backfilled, and a row made before the revision regenerates as it always did. A column already there is skipped, so a
stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0134"
down_revision = "0133"
branch_labels = None
depends_on = None

COLUMNS = (("stars", "TINYINT UNSIGNED NULL"), ("max_stars", "TINYINT UNSIGNED NOT NULL DEFAULT 12"),
           ("prevalence_multi_star", "DOUBLE NOT NULL DEFAULT 0"))
"""tuple: `(column, definition)` of each column added to `system_configs`, in table order."""


def upgrade():
    connection = op.get_bind()
    changes = [f"ADD COLUMN {column} {definition}" for column, definition in COLUMNS
               if not helpers.exists_column(connection, "system_configs", column)]
    if changes:
        connection.execute(sa.text(f"ALTER TABLE system_configs {', '.join(changes)}"))
