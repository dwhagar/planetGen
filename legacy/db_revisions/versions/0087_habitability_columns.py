"""The Generation lane's columns on `planets` and `moons` (DB.26).

Schema v87. GEN.150 adds `age_myr`, `q_total` and `hosts` to `nebulae`. GEN.178 adds `induced_moment_a_m2` after `dipole_class`; GEN.148 adds `flux_variation` after
`energy_flux_w_m2`; GEN.207 (PHI-6) renames `phi4_chemistry` and `tier_chemistry` to `phi4_air` and `tier_air`
and adds `phi4_water`, `tier_water`, `phi4_suitability` and `tier_suitability`. New columns are NULL: rows
generated before v87 are unchanged until regenerated. Each step is skipped when it is already done, so a
stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0087"
down_revision = "0086"
branch_labels = None
depends_on = None

_RENAMES = (("phi4_chemistry", "phi4_air", "DOUBLE"), ("tier_chemistry", "tier_air", "TINYINT"))
_ADDS = (
    ("induced_moment_a_m2", "DOUBLE", "dipole_class"),
    ("flux_variation", "DOUBLE", "energy_flux_w_m2"),
    ("phi4_water", "DOUBLE", "phi4_air"),
    ("phi4_suitability", "DOUBLE", "phi4_water"),
    ("tier_water", "TINYINT", "tier_air"),
    ("tier_suitability", "TINYINT", "tier_water"),
)

_NEBULA_ADDS = (
    ("age_myr", "DOUBLE", "extinction_av"),
    ("q_total", "DOUBLE", "age_myr"),
    ("hosts", "VARCHAR(512)", "q_total"),
)


def _has_column(connection, table, column):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar()


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        for old, new, definition in _RENAMES:
            if _has_column(connection, table, old) and not _has_column(connection, table, new):
                op.execute(f"ALTER TABLE {table} CHANGE COLUMN {old} {new} {definition} NULL")
        for column, definition, after in _ADDS:
            if not _has_column(connection, table, column):
                op.execute(f"ALTER TABLE {table} ADD COLUMN {column} {definition} NULL AFTER {after}")
    for column, definition, after in _NEBULA_ADDS:
        if not _has_column(connection, "nebulae", column):
            op.execute(f"ALTER TABLE nebulae ADD COLUMN {column} {definition} NULL AFTER {after}")
