"""The bulk composition of every planet and moon (DB.69, GEN.428).

Schema v144. `planets` and `moons` gain the shares of a body's mass that the generator now draws first and derives the mass, radius and
density from: `metal_mass_fraction`, `rock_mass_fraction`, `ice_mass_fraction` (every ice but water) and `envelope_mass_fraction`, which add
up to 1 with the existing `water_mass_fraction` (GEN.88); the rock's `rock_feo_fraction`, `rock_al2o3_fraction` and `rock_mg_si_ratio`; the
envelope's `envelope_helium_fraction` and `envelope_z_fraction`; a small body's `porosity`; `carbon_rock` (silicon carbide and graphite in
place of silicate); and `composition_origin` ('drawn' for composition first, 'inverse' for one worked back from a mass and a radius).

All NULL here: a row made before v144 has no composition to read back (it needs its star's abundances and the draws, so the regenerated
galaxy of the generator epoch replaces it). Each column is skipped when it is already there, so a stopped run can be started again; the last
column marks the migration done.
"""

import sqlalchemy as sa
from alembic import op

revision = "0144"
down_revision = "0143"
branch_labels = None
depends_on = None

_COLUMNS = (
    ("metal_mass_fraction", "DOUBLE NULL"),
    ("rock_mass_fraction", "DOUBLE NULL"),
    ("ice_mass_fraction", "DOUBLE NULL"),
    ("envelope_mass_fraction", "DOUBLE NULL"),
    ("rock_feo_fraction", "DOUBLE NULL"),
    ("rock_al2o3_fraction", "DOUBLE NULL"),
    ("rock_mg_si_ratio", "DOUBLE NULL"),
    ("envelope_helium_fraction", "DOUBLE NULL"),
    ("envelope_z_fraction", "DOUBLE NULL"),
    ("porosity", "DOUBLE NULL"),
    ("carbon_rock", "BOOLEAN NULL"),
    ("composition_origin", "VARCHAR(8) NULL CHECK (composition_origin IN ('drawn', 'inverse'))"),
)


def _has_column(connection, table, column):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar()


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        for column, definition in _COLUMNS:
            if not _has_column(connection, table, column):
                op.execute(f"ALTER TABLE {table} ADD COLUMN {column} {definition}")
