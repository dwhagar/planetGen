"""Stellar activity and planetary magnetic fields (GEN.86).

Schema v72. `stars` gains `log_lx_lbol`, `l_xuv_w`, `xuv_saturated`,
`flare_n33_per_yr`, `flare_alpha` and `xuv_fluence_j`; `planets` and `moons`
gain `magnetic_moment_a_m2`, `dipole_class`, `magnetopause_rp`,
`xuv_flux_earth`, `xuv_exposure_index` and `flare_irradiation_index`. All
NULL here: rows generated before v72 have none (see `schema.sql`'s "v72"
note).
"""

import sqlalchemy as sa
from alembic import op

revision = "0072"
down_revision = "0071"
branch_labels = None
depends_on = None

_STAR_COLUMNS = (
    ("log_lx_lbol", "DOUBLE"), ("l_xuv_w", "DOUBLE"), ("xuv_saturated", "BOOLEAN"),
    ("flare_n33_per_yr", "DOUBLE"), ("flare_alpha", "DOUBLE"), ("xuv_fluence_j", "DOUBLE"),
)
_BODY_COLUMNS = (
    ("magnetic_moment_a_m2", "DOUBLE"),
    ("dipole_class", "VARCHAR(16) CHECK (dipole_class IN ('none', 'weak', 'earth-like', 'strong', 'multipolar'))"),
    ("magnetopause_rp", "DOUBLE"), ("xuv_flux_earth", "DOUBLE"), ("xuv_exposure_index", "DOUBLE"),
    ("flare_irradiation_index", "DOUBLE"),
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
    _add(connection, "stars", _STAR_COLUMNS, "axial_tilt_deg")
    for table in ("planets", "moons"):
        _add(connection, table, _BODY_COLUMNS, "p_so2_kpa")
