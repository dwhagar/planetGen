"""The tech-level columns of `species` (POP.9).

Schema v88. `species` gains the six tech indexes (`energy_index`, `materials_index`, `information_index`,
`medical_index`, `propulsion_index`, `defense_index`, each 0 to 7), `tech_level` (DECIMAL(3,2), indexed for
range search) and `tech_source` ('generated' or 'admin'). All NULL for species saved before v88. Each step is
skipped when it is already done, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0088"
down_revision = "0087"
branch_labels = None
depends_on = None

_INDEXES = ("energy_index", "materials_index", "information_index", "medical_index", "propulsion_index",
            "defense_index")


def _has_column(connection, column):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'species' AND COLUMN_NAME = :column_name"), {"column_name": column}).scalar()


def _has_index(connection, name):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.STATISTICS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'species' AND INDEX_NAME = :index_name"), {"index_name": name}).scalar()


def upgrade():
    connection = op.get_bind()
    after = "spacefaring"
    for column in _INDEXES:
        if not _has_column(connection, column):
            op.execute(f"ALTER TABLE species ADD COLUMN {column} TINYINT UNSIGNED NULL AFTER {after}"
                       f", ADD CONSTRAINT chk_species_{column} CHECK ({column} BETWEEN 0 AND 7)")
        after = column
    if not _has_column(connection, "tech_level"):
        op.execute(f"ALTER TABLE species ADD COLUMN tech_level DECIMAL(3,2) NULL AFTER {after}")
    if not _has_column(connection, "tech_source"):
        op.execute("ALTER TABLE species ADD COLUMN tech_source VARCHAR(16) NULL AFTER tech_level"
                   ", ADD CONSTRAINT chk_species_tech_source CHECK (tech_source IN ('generated', 'admin'))")
    if not _has_index(connection, "idx_species_tech_level"):
        op.execute("CREATE INDEX idx_species_tech_level ON species (tech_level)")
