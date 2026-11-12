"""A satellite galaxy's young and intermediate stars are populations of their own (GEN.526).

Schema v107. A satellite's stars are split by its star formation history into `satellite_young` (the last 0.1 Gy),
`satellite_intermediate` (the last 3 Gy before that) and `satellite` (the rest), as the disk's are
(`galaxy/satellite_history.py`). `bright_stars.population` takes the two new values; nothing else changes and no row is
touched: the old values stay valid. The column widens to VARCHAR(24) for the longer name (an in-place change that
copies no row). A step already done is left alone, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0107"
down_revision = "0106"
branch_labels = None
depends_on = None

POPULATIONS = ("young", "intermediate", "old", "bulge", "halo", "cluster", "satellite", "nuclear_cluster", "nuclear_disc",
               "satellite_young", "satellite_intermediate")


def upgrade():
    connection = op.get_bind()
    clause = connection.execute(sa.text(
        "SELECT CHECK_CLAUSE FROM information_schema.CHECK_CONSTRAINTS WHERE CONSTRAINT_SCHEMA = DATABASE()"
        " AND CONSTRAINT_NAME = 'chk_bright_stars_population'")).scalar()
    if clause is not None and "satellite_young" in clause:
        return
    listed = ", ".join(f"'{name}'" for name in POPULATIONS)
    drop = "DROP CHECK chk_bright_stars_population, " if clause is not None else ""
    op.execute(f"ALTER TABLE bright_stars MODIFY COLUMN population VARCHAR(24) NOT NULL, {drop}"
               f"ADD CONSTRAINT chk_bright_stars_population CHECK (population IN ({listed}))")
