"""The cell star ceiling defaults to 5,000 systems, from 50,000 (GEN.281).

Schema v100. A system takes about 20 ms to generate, so a cell of 50,000 took 14 to 17 minutes; Boss set the default to
5,000, about 100 seconds. `galaxy_shape.cell_star_ceiling` takes the new default, and a galaxy still at the old 50,000
is moved to it (one planned with another value keeps it). The column gains a CHECK (at least 1, as `plan` already
demands), which is also how the schema version is told. A step already done is left alone, so a stopped run can be
started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0100"
down_revision = "0099"
branch_labels = None
depends_on = None

CHECK = "chk_galaxy_shape_cell_star_ceiling"


def upgrade():
    connection = op.get_bind()
    op.execute("ALTER TABLE galaxy_shape ALTER COLUMN cell_star_ceiling SET DEFAULT 5000")
    op.execute("UPDATE galaxy_shape SET cell_star_ceiling = 5000 WHERE cell_star_ceiling = 50000")
    present = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.CHECK_CONSTRAINTS WHERE CONSTRAINT_SCHEMA = DATABASE()"
        " AND CONSTRAINT_NAME = :name"), {"name": CHECK}).scalar()
    if not present:
        op.execute(f"ALTER TABLE galaxy_shape ADD CONSTRAINT {CHECK} CHECK (cell_star_ceiling >= 1)")
