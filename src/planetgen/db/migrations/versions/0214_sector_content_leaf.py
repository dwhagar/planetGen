"""The content part of a sector's leaf digest (TEST.602).

Schema v214. `sectors` gains `content_leaf` (BINARY(32), NULL): the leaf digest of the sector as it was saved (`fingerprint_tree.leaf_digest`, over
the sector in memory), kept apart from `leaf_digest`, which is now the leaf the fingerprint tree reads: the content part joined to a digest of the
rogue planets filed in the sector (`fingerprint_tree.resident_digest`) when it holds any. A rogue lands in a sector other than its donor's, often one
nobody has charted, and is written after that sector's leaf was taken, so the leaf has to be made again from a content part it can still find.
A row saved before v214 has `content_leaf` NULL and keeps the `leaf_digest` it has until it is saved again. The column is skipped when it is already
there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0214"
down_revision = "0213"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_column(connection, "sectors", "content_leaf"):
        connection.execute(sa.text("ALTER TABLE sectors ADD COLUMN content_leaf BINARY(32) NULL AFTER leaf_digest, ALGORITHM=COPY"))
