"""The halo and cluster stellar populations (GEN.234).

Schema v94. The stellar halo replaces the flat density floor (`galaxy/density.py`), so a pre-placed bright star can
now belong to the `halo` population, and the globular clusters' stars to `cluster`. `bright_stars.population` was
limited to the four disk and bulge populations; it takes the two new ones. Nothing else changes, and nothing is
renumbered: the old values stay valid.
"""

import sqlalchemy as sa
from alembic import op

revision = "0094"
down_revision = "0093"
branch_labels = None
depends_on = None

POPULATIONS = ("young", "intermediate", "old", "bulge", "halo", "cluster")


def upgrade():
    connection = op.get_bind()
    clause = connection.execute(sa.text(
        "SELECT CHECK_CLAUSE FROM information_schema.CHECK_CONSTRAINTS WHERE CONSTRAINT_SCHEMA = DATABASE()"
        " AND CONSTRAINT_NAME = 'chk_bright_stars_population'")).scalar()
    if clause is not None and "halo" in clause:
        return
    listed = ", ".join(f"'{name}'" for name in POPULATIONS)
    drop = "DROP CHECK chk_bright_stars_population, " if clause is not None else ""
    op.execute(f"ALTER TABLE bright_stars {drop}ADD CONSTRAINT chk_bright_stars_population"
               f" CHECK (population IN ({listed}))")
