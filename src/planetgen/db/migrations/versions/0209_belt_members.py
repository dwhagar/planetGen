"""The stored members of a belt (GEN.1070).

Schema v209. A new table, `belt_members`: one row for each individual body a belt stores (the ledger's cut rule of `physics/belt_ledger.py`: every body above 1e19 kg and more of the next
heaviest until 90 percent of the belt's mass is held, at most 300), with its rank by mass (1 is the heaviest and keeps its number, so `<belt object ID>.<rank in hex>` names the body), mass,
radius, orbit (semi-major axis, eccentricity, inclination and the three angles), taxon, density, composition shares, status word, Soter's mu and the scatter tag (NULL until GEN.1092).

A belt saved before v209 has no members: they come from the belt's own draws when the galaxy is regenerated. The table is skipped when it is already there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0209"
down_revision = "0208"
branch_labels = None
depends_on = None

TABLE = """
    CREATE TABLE IF NOT EXISTS belt_members (
        id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
        belt_id           BIGINT UNSIGNED NOT NULL,
        member_rank       INT NOT NULL,
        mass_kg           DOUBLE NOT NULL,
        radius_km         DOUBLE NOT NULL,
        distance_km       DOUBLE NOT NULL,
        eccentricity      DOUBLE NOT NULL,
        inclination_deg   DOUBLE NOT NULL,
        ascending_node_deg DOUBLE NOT NULL,
        arg_periapsis_deg DOUBLE NOT NULL,
        mean_anomaly_deg  DOUBLE NOT NULL,
        taxon             VARCHAR(4) NOT NULL,
        density_g_cm3     DOUBLE NOT NULL,
        rock_fraction     DOUBLE NOT NULL,
        metal_fraction    DOUBLE NOT NULL,
        water_fraction    DOUBLE NOT NULL,
        ice_fraction      DOUBLE NOT NULL,
        status            VARCHAR(16) NOT NULL,
        soter_mu          DOUBLE NULL,
        scatter_tag       VARCHAR(24) NULL,
        CONSTRAINT fk_belt_members_belt
            FOREIGN KEY (belt_id) REFERENCES asteroid_belts(id) ON DELETE CASCADE,
        UNIQUE KEY uq_belt_members_belt_rank (belt_id, member_rank)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_table(connection, "belt_members"):
        connection.execute(sa.text(TABLE))
