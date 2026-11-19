"""A planet that shares its giant's orbit at L4 or L5 (GEN.1151).

Schema v270. `planets` gains `lagrange_point` (VARCHAR(2), `L4` or `L5`) and `lagrange_host_id` (BIGINT UNSIGNED, the `planets.id` of the giant): a planet grown in the 1:1 zone of a giant
(`generation/trojan_planet.py`) orbits the star on the giant's own orbit, 60 degrees ahead or behind it. Both are NULL for every other planet and for a row made before v270; a migration
does not draw, and the regenerated galaxy holds the trojan planets. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0270"
down_revision = "0269"
branch_labels = None
depends_on = None

COLUMNS = (("lagrange_point", "VARCHAR(2)"), ("lagrange_host_id", "BIGINT UNSIGNED"))


def upgrade():
    connection = op.get_bind()
    for name, kind in COLUMNS:
        if not helpers.exists_column(connection, "planets", name):
            connection.execute(sa.text(f"ALTER TABLE planets ADD COLUMN {name} {kind} NULL"))

