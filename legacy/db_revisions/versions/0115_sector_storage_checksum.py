"""The sector storage checksum and the parity dirty flag (DB.44).

Schema v115. Two columns on `sectors`, the two the rest of DB.9 reads: `storage_checksum` (BINARY(32), NULL) holds the
digest of a sector's stored rows once DB.9's code computes one, and `parity_dirty` (TINYINT NOT NULL, default 1) says the
sector's parity record is not current, so a sector that existed before this revision counts as not yet protected. No row
is touched and nothing reads the columns yet.

The columns are added in one statement with `ALGORITHM=COPY`, not `INSTANT`: `sectors` has a FULLTEXT index
(`ft_sectors_name`), and InnoDB refuses an instant or in-place change of a table that has one (error 1846), so the table
is rebuilt once. It holds one row per sector, not per object. A column already there is skipped, so a stopped run can be
started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0115"
down_revision = "0114"
branch_labels = None
depends_on = None

COLUMNS = (("storage_checksum", "BINARY(32) NULL"), ("parity_dirty", "TINYINT NOT NULL DEFAULT 1"))
"""tuple: `(column, definition)`, added to `sectors` after `max_mass_id`."""


def upgrade():
    connection = op.get_bind()
    missing = [f"ADD COLUMN {column} {definition}" for column, definition in COLUMNS
               if not helpers.exists_column(connection, "sectors", column)]
    if missing:
        connection.execute(sa.text(f"ALTER TABLE sectors {', '.join(missing)}, ALGORITHM=COPY"))
