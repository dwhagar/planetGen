"""The per-species environment limits: one row per limit a species does not share with the human visitor (DB.30).

Schema v119. A species' environment limits are the human table (`planetgen.physics.protection_limits.LIMIT_ROWS`) with some
rows replaced (docs/design/habitability-index.md section 10). `species_limit_values` holds the replaced rows: the species,
the name of the limit and its value in the limit's own unit. A limit with no row is the human value, so a species with the
human limits has no rows, and a limit added to `LIMIT_ROWS` later needs no migration. The name is checked in code
(`species_limits.validate`), not by a constraint that would have to change with the table of limits.

This revision makes the table and nothing else; no code reads it yet. A table already there is left alone,
so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0119"
down_revision = "0118"
branch_labels = None
depends_on = None

TABLES = {
    "species_limit_values": """
CREATE TABLE IF NOT EXISTS species_limit_values (
    species_id   BIGINT UNSIGNED NOT NULL,
    limit_name   VARCHAR(48) NOT NULL,
    value        DOUBLE NOT NULL,
    PRIMARY KEY (species_id, limit_name),
    CONSTRAINT fk_species_limit_values_species
        FOREIGN KEY (species_id) REFERENCES species(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
""",
}


def upgrade():
    connection = op.get_bind()
    for name, create in TABLES.items():
        if not connection.execute(sa.text(
                "SELECT COUNT(*) FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()"
                " AND TABLE_NAME = :name"), {"name": name}).scalar():
            op.execute(create)
