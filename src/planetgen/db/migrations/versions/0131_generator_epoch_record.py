"""The galaxy record and each run keep the generator epoch, and the record the fingerprint format and the battery (OPS.58).

Schema v131. `galaxy_shape` (one row, the code that planned the galaxy) gains `generator_epoch` (SMALLINT UNSIGNED), `fp_spec`
(SMALLINT UNSIGNED, the format of `util/fingerprint_encoding.py`), `battery_digest` (CHAR(64), the overall digest of the generator
battery, OPS.57) and `battery_cases` (JSON, `{case: full digest}`, so a check can name the first case that differs); `generation_runs`
gains `generator_epoch` (SMALLINT UNSIGNED), so a galaxy says which epoch each of its passes (a scatter, a fill, a redo) was made at.
All five are NULL on a galaxy and on runs from before this revision: the epoch they were made under was never recorded, and a check
says nothing of a NULL. No row is touched and nothing is backfilled.

`galaxy_shape` has one row and `generation_runs` one row a run, so each change is a single statement. A column already there is
skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0131"
down_revision = "0130"
branch_labels = None
depends_on = None

COLUMNS = {
    "galaxy_shape": (("generator_epoch", "SMALLINT UNSIGNED NULL"), ("fp_spec", "SMALLINT UNSIGNED NULL"),
                     ("battery_digest", "CHAR(64) NULL"), ("battery_cases", "JSON NULL")),
    "generation_runs": (("generator_epoch", "SMALLINT UNSIGNED NULL"),),
}
"""dict: `table -> ((column, definition), ...)`, added at the end of the table."""


def upgrade():
    connection = op.get_bind()
    for table, columns in COLUMNS.items():
        changes = [f"ADD COLUMN {column} {definition}" for column, definition in columns
                   if not helpers.exists_column(connection, table, column)]
        if changes:
            connection.execute(sa.text(f"ALTER TABLE {table} {', '.join(changes)}"))
