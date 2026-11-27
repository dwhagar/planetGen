"""The coma of an active comet (GEN.455, part A).

Schema v275 (a placeholder: the train that merges this renumbers it to the next free revision). `comets` and `interstellar_comets` gain the four numbers of a coma at 1 AU lit distance
(`water_production_1au_per_s`, `coma_active_area_fraction`, `coma_outflow_1au_kms`, `nucleus_albedo`), and `comet_coma` and `interstellar_comet_coma` are new: one row for water and each volatile the
nucleus holds (`physics/comet_coma.py`): the molecules a second at 1 AU lit and the parent scale length at 1 AU lit (none for a gas the photolysis table lacks). An inactive comet, and a comet made
before v275, has none of it (NULL, no rows), and nothing is backfilled (a migration does not draw): the next regeneration makes them. A column or table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0275"
down_revision = "0274"
branch_labels = None
depends_on = None

COLUMNS = (
    ("water_production_1au_per_s", "DOUBLE NULL"),
    ("coma_active_area_fraction", "DOUBLE NULL"),
    ("coma_outflow_1au_kms", "DOUBLE NULL"),
    ("nucleus_albedo", "DOUBLE NULL"),
)
"""tuple: The columns `comets` and `interstellar_comets` gain, with their definitions."""

TABLE = """
CREATE TABLE {prefix}comet_coma (
    id                    BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    comet_id              BIGINT UNSIGNED NOT NULL,
    position              INT NOT NULL,
    species               VARCHAR(16) NOT NULL,
    production_1au_per_s  DOUBLE NOT NULL,
    scale_length_1au_km   DOUBLE NULL,
    CONSTRAINT fk_{prefix}comet_coma_comet
        FOREIGN KEY (comet_id) REFERENCES {owner}(id) ON DELETE CASCADE,
    UNIQUE KEY uq_{prefix}comet_coma (comet_id, position),
    KEY idx_{prefix}comet_coma_species (species)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
"""

OWNERS = (("", "comets"), ("interstellar_", "interstellar_comets"))
"""tuple: `(table prefix, owning table)` of the two coma tables."""


def upgrade():
    connection = op.get_bind()
    for _, owner in OWNERS:
        for name, definition in COLUMNS:
            if not helpers.exists_column(connection, owner, name):
                connection.execute(sa.text(f"ALTER TABLE {owner} ADD COLUMN {name} {definition}"))
    for prefix, owner in OWNERS:
        if not helpers.exists_table(connection, f"{prefix}comet_coma"):
            connection.execute(sa.text(TABLE.format(prefix=prefix, owner=owner)))


def downgrade():
    connection = op.get_bind()
    for prefix, _ in reversed(OWNERS):
        if helpers.exists_table(connection, f"{prefix}comet_coma"):
            connection.execute(sa.text(f"DROP TABLE {prefix}comet_coma"))
    for _, owner in reversed(OWNERS):
        for name, _ in reversed(COLUMNS):
            if helpers.exists_column(connection, owner, name):
                connection.execute(sa.text(f"ALTER TABLE {owner} DROP COLUMN {name}"))
