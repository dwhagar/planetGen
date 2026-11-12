"""The leaf digest of each sector and the shell and band tree of digests above it (DB.49).

Schema v120. `sectors.leaf_digest` (BINARY(32), NULL) is the fingerprint of one sector (GEN.136's canonical encoder, SHA-256),
stored at save time; NULL for a sector saved before it, which the fingerprint tree (GEN.305) fills as it meets them.
`sector_tree_nodes` holds the digests above the leaves (docs/design/generation/generation-determinism.md section 4.3): a band node
(`level` 1, one per shell and band) over the leaves of its wedges, and a shell node (`level` 2, one per shell, `band_index`
0) over its band nodes. The galaxy root is worked out on demand and not stored. `child_count` is the number of leaves under
a band node and of band nodes under a shell node, so that a region's digest can include the count; `dirty` is 1 while a
save or an edit has changed something under the node since `digest` was made, and a node never made has a NULL `digest`.

This revision adds the column and makes the table, and nothing else: no row is written, and no code reads either until
GEN.305. The column is added with `ALGORITHM=COPY`, not `INSTANT`: `sectors` has a FULLTEXT index (`ft_sectors_name`), and
InnoDB refuses an instant or in-place change of a table that has one (error 1846), so the table is rebuilt once. It holds one
row per sector, not per object. A step already done is left alone, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0120"
down_revision = "0119"
branch_labels = None
depends_on = None

CREATE_SECTOR_TREE_NODES = """
CREATE TABLE IF NOT EXISTS sector_tree_nodes (
    level         TINYINT UNSIGNED NOT NULL,
    shell_index   INT NOT NULL,
    band_index    INT NOT NULL DEFAULT 0,
    digest        BINARY(32) NULL,
    child_count   INT UNSIGNED NOT NULL DEFAULT 0,
    dirty         TINYINT(1) NOT NULL DEFAULT 1,
    updated_at    TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
    PRIMARY KEY (level, shell_index, band_index),
    KEY idx_sector_tree_nodes_dirty (dirty, level),
    CONSTRAINT chk_sector_tree_nodes_level CHECK (level IN (1, 2)),
    CONSTRAINT chk_sector_tree_nodes_shell_band CHECK (level = 1 OR band_index = 0),
    CONSTRAINT chk_sector_tree_nodes_dirty CHECK (dirty IN (0, 1))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
"""


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_column(connection, "sectors", "leaf_digest"):
        connection.execute(sa.text("ALTER TABLE sectors ADD COLUMN leaf_digest BINARY(32) NULL, ALGORITHM=COPY"))
    if not helpers.exists_table(connection, "sector_tree_nodes"):
        op.execute(CREATE_SECTOR_TREE_NODES)
