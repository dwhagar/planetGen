"""Planets and moons store their mantle redox and their air (GEN.85).

Schema v70. `planets` and `moons` gain `mantle_redox`, `mantle_delta_iw`
and a partial pressure, kPa, for each of O2, CO2, CO, N2, Ar, H2, H2O, CH4,
H2S and SO2. All NULL here: rows generated before v70 have none (see
`schema.sql`'s "v70" note).
"""

import sqlalchemy as sa
from alembic import op

revision = "0070"
down_revision = "0069"
branch_labels = None
depends_on = None

_GASES = ("o2", "co2", "co", "n2", "ar", "h2", "h2o", "ch4", "h2s", "so2")

_COLUMNS = (
    ("mantle_redox", "VARCHAR(16) CHECK (mantle_redox IN ('reduced', 'intermediate', 'oxidized'))"),
    ("mantle_delta_iw", "DOUBLE"),
) + tuple((f"p_{gas}_kpa", "DOUBLE") for gas in _GASES)


def _has_column(connection, table, column):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar()


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        if _has_column(connection, table, "mantle_redox"):
            continue
        after = "axial_tilt_deg"
        clauses = []
        for column, definition in _COLUMNS:
            clauses.append(f"ADD COLUMN {column} {definition} AFTER {after}")
            after = column
        op.execute(f"ALTER TABLE {table} " + ", ".join(clauses))
