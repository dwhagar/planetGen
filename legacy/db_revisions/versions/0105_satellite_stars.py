"""A satellite galaxy's stars are a stellar population of their own (GEN.324).

Schema v105. The component pass (`generation/component_pass.py`) places the bright and massive stars of the globular
clusters and of the satellite galaxies (the Magellanic Clouds, the Sagittarius remnant, the dwarfs). A cluster's stars
were already the `cluster` population (v94); a satellite's are `satellite`, which `bright_stars.population` did not take.
It takes it now; nothing else changes and no row is touched: the old values stay valid. A step already done is left
alone, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0105"
down_revision = "0104"
branch_labels = None
depends_on = None

POPULATIONS = ("young", "intermediate", "old", "bulge", "halo", "cluster", "satellite", "nuclear_cluster", "nuclear_disc")


def upgrade():
    connection = op.get_bind()
    clause = connection.execute(sa.text(
        "SELECT CHECK_CLAUSE FROM information_schema.CHECK_CONSTRAINTS WHERE CONSTRAINT_SCHEMA = DATABASE()"
        " AND CONSTRAINT_NAME = 'chk_bright_stars_population'")).scalar()
    if clause is not None and "satellite" in clause:
        return
    listed = ", ".join(f"'{name}'" for name in POPULATIONS)
    drop = "DROP CHECK chk_bright_stars_population, " if clause is not None else ""
    op.execute(f"ALTER TABLE bright_stars {drop}ADD CONSTRAINT chk_bright_stars_population"
               f" CHECK (population IN ({listed}))")
