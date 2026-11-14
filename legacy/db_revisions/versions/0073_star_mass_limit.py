"""The star scatter's mass limit (GEN.185).

Schema v73. `galaxy_shape` gains `bright_star_mass_limit_sol`: the mass, in
solar masses, from which the star scatter's mass pass placed every star
(the same limit as the phenomenon scatter's). A sector's fill and the
scatter's luminosity pass draw only lighter stars. A galaxy scattered before
v73 keeps NULL: its scatter had no mass pass, so nothing is excluded.
"""

import sqlalchemy as sa
from alembic import op

revision = "0073"
down_revision = "0072"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    present = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.columns WHERE table_schema = DATABASE()"
        " AND table_name = 'galaxy_shape' AND column_name = 'bright_star_mass_limit_sol'")).scalar()
    if not present:
        connection.execute(sa.text(
            "ALTER TABLE galaxy_shape ADD COLUMN bright_star_mass_limit_sol DOUBLE AFTER bright_star_seed"))
