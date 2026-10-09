"""The phenomenon scatter's mass cut (GEN.167).

Schema v71. `galaxy_shape` gains `phenomenon_min_mass_solar`: the lowest
mass, in solar masses, of the neutron stars and black holes the phenomenon
scatter placed (`planetgen plan --phenomenon-min-mass`). A sector fill draws
the ones below it itself (GEN.168). A galaxy scattered before v71 keeps NULL:
its scatter placed every mass, so its fills draw none.
"""

import sqlalchemy as sa
from alembic import op

revision = "0071"
down_revision = "0070"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    present = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.columns WHERE table_schema = DATABASE()"
        " AND table_name = 'galaxy_shape' AND column_name = 'phenomenon_min_mass_solar'")).scalar()
    if not present:
        connection.execute(sa.text(
            "ALTER TABLE galaxy_shape ADD COLUMN phenomenon_min_mass_solar DOUBLE AFTER phenomenon_scatter_seed"))
