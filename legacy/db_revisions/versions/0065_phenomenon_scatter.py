"""The galaxy-wide phenomenon scatter (GEN.100).

Schema v64. `phenomenon_scatter` holds the black holes, neutron stars,
planetary nebulae, supernova remnants, hypervelocity stars and nucleus that
`planetgen plan` places before any sector is filled, and
`galaxy_shape.phenomenon_scatter_seed` records the seed it used. `schema.sql`
has already created the table by the time this runs; the column on an
existing `galaxy_shape` is the only thing to add. A galaxy planned before
this has no scatter (the seed stays NULL): its sectors keep rolling their own
phenomena until `planetgen plan` places them.
"""

import sqlalchemy as sa
from alembic import op

revision = "0065"
down_revision = "0064"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    present = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.columns WHERE table_schema = DATABASE()"
        " AND table_name = 'galaxy_shape' AND column_name = 'phenomenon_scatter_seed'")).scalar()
    if not present:
        connection.execute(sa.text("ALTER TABLE galaxy_shape ADD COLUMN phenomenon_scatter_seed BIGINT UNSIGNED"))
