"""A generated galaxy's young open clusters add 12 age-bin populations (GEN.244).

Schema v112. A young open cluster's stars all share its age, so each cluster is held in the one of twelve age bins
(`open_cluster_00` to `open_cluster_11`, a quarter dex wide from 1 Myr to 1 Gyr, `tuning.OPEN_CLUSTER_POPULATIONS`) its
age falls in. `bright_stars.population` takes the twelve new values; nothing else changes and no row is touched: the old
values stay valid. A step already done is left alone, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0112"
down_revision = "0111"
branch_labels = None
depends_on = None

POPULATIONS = ("young", "intermediate", "old", "bulge", "halo", "cluster", "satellite", "nuclear_cluster", "nuclear_disc",
               "satellite_young", "satellite_intermediate", *(f"open_cluster_{k:02d}" for k in range(12)))


def upgrade():
    connection = op.get_bind()
    clause = connection.execute(sa.text(
        "SELECT CHECK_CLAUSE FROM information_schema.CHECK_CONSTRAINTS WHERE CONSTRAINT_SCHEMA = DATABASE()"
        " AND CONSTRAINT_NAME = 'chk_bright_stars_population'")).scalar()
    if clause is not None and "open_cluster_00" in clause:
        return
    listed = ", ".join(f"'{name}'" for name in POPULATIONS)
    drop = "DROP CHECK chk_bright_stars_population, " if clause is not None else ""
    op.execute(f"ALTER TABLE bright_stars {drop}"
               f"ADD CONSTRAINT chk_bright_stars_population CHECK (population IN ({listed}))")
