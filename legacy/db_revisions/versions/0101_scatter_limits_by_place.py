"""The scatter limits by place as three settings of the galaxy (DB.54).

Schema v101. `galaxy_shape.core_min_luminosity_sol`
and `core_mass_limit_sol` (DOUBLE, NULL until set: the core follows the baseline's limits times 1.1,
`galaxy/scatter_limits.py`) and `halo_min_luminosity_sol` (DOUBLE NOT NULL, default 500 solar luminosities, at least
50), with one CHECK on all three. The halo has no mass limit of its own. A galaxy planned before this takes the
defaults, which change nothing until GEN.333 reads them. A step already done is left alone, so a stopped run can be
started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0101"
down_revision = "0100"
branch_labels = None
depends_on = None

CHECK = "chk_galaxy_shape_scatter_limits"
COLUMNS = (
    ("core_min_luminosity_sol", "DOUBLE NULL", "fill_luminosity_cap_sol"),
    ("core_mass_limit_sol", "DOUBLE NULL", "core_min_luminosity_sol"),
    ("halo_min_luminosity_sol", "DOUBLE NOT NULL DEFAULT 500", "core_mass_limit_sol"),
)


def _has_column(connection, name):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'galaxy_shape' AND COLUMN_NAME = :name"), {"name": name}).scalar())


def upgrade():
    connection = op.get_bind()
    for name, definition, after in COLUMNS:
        if not _has_column(connection, name):
            op.execute(f"ALTER TABLE galaxy_shape ADD COLUMN {name} {definition} AFTER {after}")
    present = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.CHECK_CONSTRAINTS WHERE CONSTRAINT_SCHEMA = DATABASE()"
        " AND CONSTRAINT_NAME = :name"), {"name": CHECK}).scalar()
    if not present:
        op.execute(
            f"ALTER TABLE galaxy_shape ADD CONSTRAINT {CHECK} CHECK (halo_min_luminosity_sol >= 50"
            " AND (core_min_luminosity_sol IS NULL OR core_min_luminosity_sol > 0)"
            " AND (core_mass_limit_sol IS NULL OR core_mass_limit_sol > 0))")
