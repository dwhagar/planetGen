"""Hydrosphere and ocean chemistry (GEN.88).

Schema v74. `planets` and `moons` gain `water_mass_fraction`,
`hydrosphere`, `ocean_fraction`, `land_fraction`, `ocean_depth_km`,
`ice_shell_km`, `hp_ice_km`, `ocean_class`, `ocean_ph`, `water_activity`
and `phosphorus`; `rogue_planets` gains `hp_ice_km`. All NULL here: rows
generated before v74 have none (see `schema.sql`'s "v74" note).
"""

import sqlalchemy as sa
from alembic import op

revision = "0074"
down_revision = "0073"
branch_labels = None
depends_on = None

_BODY_COLUMNS = (
    ("water_mass_fraction", "DOUBLE"),
    ("hydrosphere", "VARCHAR(20) CHECK (hydrosphere IN ('dry', 'vapour', 'ice', 'ice-covered ocean', "
                    "'surface ocean', 'hycean'))"),
    ("ocean_fraction", "DOUBLE"), ("land_fraction", "DOUBLE"), ("ocean_depth_km", "DOUBLE"),
    ("ice_shell_km", "DOUBLE"), ("hp_ice_km", "DOUBLE"),
    ("ocean_class", "VARCHAR(16) CHECK (ocean_class IN ('ice-sealed', 'chloride brine', 'acid sulfate', "
                    "'soda', 'neutral'))"),
    ("ocean_ph", "DOUBLE"), ("water_activity", "DOUBLE"),
    ("phosphorus", "VARCHAR(8) CHECK (phosphorus IN ('high', 'limited', 'starved'))"),
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
        _add(connection, table, _BODY_COLUMNS, "flare_irradiation_index")
    _add(connection, "rogue_planets", (("hp_ice_km", "DOUBLE"),), "has_liquid_water")
