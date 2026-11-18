"""Surface radiation dose, UV and the galactic hazard (GEN.87).

Schema v75. `planets` and `moons` gain `surface_dose_msv_yr`,
`dose_gcr_msv_yr`, `dose_sep_msv_yr`, `dose_ground_msv_yr`,
`dose_helio_mult`, `uv_surface_index` and `ozone_loss_flag`; `stars` gains
`lethal_event_rate_per_gyr`. All NULL here: rows generated before v75 have
none (see `schema.sql`'s "v75" note).
"""

import sqlalchemy as sa
from alembic import op

revision = "0075"
down_revision = "0074"
branch_labels = None
depends_on = None

_BODY_COLUMNS = (
    ("surface_dose_msv_yr", "DOUBLE"), ("dose_gcr_msv_yr", "DOUBLE"), ("dose_sep_msv_yr", "DOUBLE"),
    ("dose_ground_msv_yr", "DOUBLE"), ("dose_helio_mult", "DOUBLE"), ("uv_surface_index", "DOUBLE"),
    ("ozone_loss_flag", "BOOLEAN"),
)


def _has_column(connection, table, column):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar()


def _add(connection, table, columns, after):
    if _has_column(connection, table, columns[0][0]):
        return
    clauses = []
    for column, definition in columns:
        clauses.append(f"ADD COLUMN {column} {definition} AFTER {after}")
        after = column
    op.execute(f"ALTER TABLE {table} " + ", ".join(clauses))


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        _add(connection, table, _BODY_COLUMNS, "phosphorus")
    _add(connection, "stars", (("lethal_event_rate_per_gyr", "DOUBLE"),), "xuv_fluence_j")
