"""The individual bodies of an interstellar asteroid field (GEN.1105).

Schema v226. A new table, `field_members`: one row for each body a field stores (the belt ledger's cut rule), with its rank (1 is the heaviest, `<field ID>.<rank hex>` names it), mass, radius, density, composition
shares and status. A field made before v226 has none until the galaxy is regenerated; a table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0226"
down_revision = "0225"
branch_labels = None
depends_on = None

TABLE = """
    CREATE TABLE IF NOT EXISTS field_members (
        id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
        field_id          BIGINT UNSIGNED NOT NULL,
        member_rank       INT NOT NULL,
        mass_kg           DOUBLE NOT NULL,
        radius_km         DOUBLE NOT NULL,
        density_g_cm3     DOUBLE NOT NULL,
        rock_fraction     DOUBLE NOT NULL,
        metal_fraction    DOUBLE NOT NULL,
        water_fraction    DOUBLE NOT NULL,
        ice_fraction      DOUBLE NOT NULL,
        status            VARCHAR(16) NOT NULL,
        CONSTRAINT fk_field_members_field
            FOREIGN KEY (field_id) REFERENCES asteroid_fields(id) ON DELETE CASCADE,
        UNIQUE KEY uq_field_members_field_rank (field_id, member_rank)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def upgrade():
    op.get_bind().execute(sa.text(TABLE))
