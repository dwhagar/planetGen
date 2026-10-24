"""The wiki links are gone: `sectors.wiki_url`, `star_systems.wikijs_url` and `star_systems.mediawiki_url` are dropped (OPS.69).

Schema v149. The wiki upload (and the hand-set wiki link on a sector) was taken out of the site, and nothing reads these three
columns any more. They are dropped; a column already gone is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0149"
down_revision = "0148"
branch_labels = None
depends_on = None

COLUMNS = (("sectors", "wiki_url"), ("star_systems", "wikijs_url"), ("star_systems", "mediawiki_url"))
"""tuple: The `(table, column)` pairs that go."""


def upgrade():
    connection = op.get_bind()
    for table, column in COLUMNS:
        if helpers.exists_column(connection, table, column):
            # Not INSTANT: `star_systems` has a FULLTEXT index, and InnoDB refuses an instant change beside one (as in 0143).
            connection.execute(sa.text(f"ALTER TABLE {table} DROP COLUMN {column}"))
