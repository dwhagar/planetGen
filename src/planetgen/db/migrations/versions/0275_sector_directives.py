"""The directive a sector was drawn under and how the draw came out (GEN.179).

Schema v275 (a placeholder: the train that merges this renumbers it to the next free revision, 0280 at the time of writing). `sector_directives` and `sector_directive_minimums` are new: a run given `--directive`
minimums (`generation/directives.py`, GEN.96) draws a sector again until it meets them, and keeps what the draw came out as. `sector_directives` has one row for each such sector: the digest of the directive, the status
(`met_naturally`, `met_after_k_attempts`, `met_forced` or `unmet`) and the attempts made; `sector_directive_minimums` one row for each minimum asked: its key, what was wanted and what the kept sector got. A sector
drawn with no directive has no rows, and nothing is backfilled (a migration does not draw). A table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0275"
down_revision = "0274"
branch_labels = None
depends_on = None

DIRECTIVES = """
CREATE TABLE sector_directives (
    sector_id         BIGINT UNSIGNED NOT NULL PRIMARY KEY,
    directive_digest  CHAR(16) NOT NULL,
    status            VARCHAR(24) NOT NULL CHECK (status IN ('met_naturally', 'met_after_k_attempts', 'met_forced', 'unmet')),
    attempts          SMALLINT UNSIGNED NOT NULL CHECK (attempts >= 1),
    CONSTRAINT fk_sector_directives_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
    KEY idx_sector_directives_status (status)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
"""

MINIMUMS = """
CREATE TABLE sector_directive_minimums (
    sector_id    BIGINT UNSIGNED NOT NULL,
    minimum_key  VARCHAR(24) NOT NULL,
    wanted       INT UNSIGNED NOT NULL,
    got          INT UNSIGNED NOT NULL,
    PRIMARY KEY (sector_id, minimum_key),
    CONSTRAINT fk_sector_directive_minimums_sector
        FOREIGN KEY (sector_id) REFERENCES sector_directives(sector_id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
"""


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_table(connection, "sector_directives"):
        connection.execute(sa.text(DIRECTIVES))
    if not helpers.exists_table(connection, "sector_directive_minimums"):
        connection.execute(sa.text(MINIMUMS))


def downgrade():
    connection = op.get_bind()
    for table in ("sector_directive_minimums", "sector_directives"):
        if helpers.exists_table(connection, table):
            connection.execute(sa.text(f"DROP TABLE {table}"))
