"""The density multiplier a sector was filled under and the override it came from (DB.60).

Schema v125. Two columns on `sectors`, added at the end of the table: `density_multiplier_used` (DECIMAL(10,4), NULL) is the
multiplier of the region density override (DB.38) that applied when the sector was filled, and `override_id` (BIGINT
UNSIGNED, NULL, the `id` of that row of `region_density_overrides`; no foreign key, since an override row is never deleted
and `sectors` is made before the table in a new galaxy) names that override, so a rebuilt galaxy can tell which setting a sector was filled under (docs/design/interface/map-ui-and-frontend-libraries.md
section 8). Both are NULL for a sector filled with no override, and for every sector saved before this revision. No row is
touched and nothing writes or reads the columns yet.

The columns are added in one statement with `ALGORITHM=COPY`, not `INSTANT`: `sectors` has a FULLTEXT index
(`ft_sectors_name`), and InnoDB refuses an instant or in-place change of a table that has one (error 1846), so the table
is rebuilt once. It holds one row per sector, not per object. A column already there is skipped, so a stopped run can be
started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0125"
down_revision = "0124"
branch_labels = None
depends_on = None

COLUMNS = (("density_multiplier_used", "DECIMAL(10,4) NULL"), ("override_id", "BIGINT UNSIGNED NULL"))
"""tuple: `(column, definition)`, added to `sectors` after `parity_dirty`."""


def upgrade():
    connection = op.get_bind()
    changes = [f"ADD COLUMN {column} {definition}" for column, definition in COLUMNS
               if not helpers.exists_column(connection, "sectors", column)]
    if changes:
        connection.execute(sa.text(f"ALTER TABLE sectors {', '.join(changes)}, ALGORITHM=COPY"))
