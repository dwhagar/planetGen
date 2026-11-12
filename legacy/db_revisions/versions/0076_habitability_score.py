"""The habitability score of every planet and moon (GEN.89).

Schema v76. `planets` and `moons` gain `phi4` and the four domain scores
and tiers (`phi4_<domain>`, `tier_<domain>`), `phi_bio`, `phi_cpx`,
`phi_tech`, the four likelihoods `l_solv`, `l_chem`, `l_ener`, `l_rad`,
`equipment_tier`, `hab_note` and `energy_flux_w_m2`. All NULL here: rows
generated before v76 have none (see `schema.sql`'s "v76" note).
"""

import sqlalchemy as sa
from alembic import op

revision = "0076"
down_revision = "0075"
branch_labels = None
depends_on = None

_BODY_COLUMNS = (
    ("phi4", "DOUBLE"), ("phi4_pressure", "DOUBLE"), ("phi4_temperature", "DOUBLE"),
    ("phi4_chemistry", "DOUBLE"), ("phi4_radiation", "DOUBLE"),
    ("tier_pressure", "TINYINT"), ("tier_temperature", "TINYINT"),
    ("tier_chemistry", "TINYINT"), ("tier_radiation", "TINYINT"),
    ("phi_bio", "DOUBLE"), ("phi_cpx", "DOUBLE"), ("phi_tech", "DOUBLE"),
    ("l_solv", "DOUBLE"), ("l_chem", "DOUBLE"), ("l_ener", "DOUBLE"), ("l_rad", "DOUBLE"),
    ("equipment_tier", "TINYINT"), ("hab_note", "VARCHAR(255)"), ("energy_flux_w_m2", "DOUBLE"),
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
        _add(connection, table, _BODY_COLUMNS, "ozone_loss_flag")
