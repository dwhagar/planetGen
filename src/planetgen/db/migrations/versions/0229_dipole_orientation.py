"""Where the magnetic dipole points and the surface field of a star (GEN.1350).

Schema v229. `planets`, `moons` and `rogue_planets` gain `dipole_tilt_deg` (DOUBLE, NULL: the angle between the dipole's axis and the spin axis, 0 to 90 degrees) and
`dipole_offset_class` (VARCHAR(8), NULL: `centred` or `offset`). `stars` gains `surface_dipole_field_gauss` (DOUBLE, NULL: the polar field of the large-scale dipole of a cool
dwarf or an Ap dwarf) and `neutron_stars` gains `magnetic_inclination_deg` (DOUBLE, NULL: the angle between the magnetic and the spin axis, 0 to 90 degrees). A row made before v229
has NULL in them and gets the values when the galaxy is regenerated (a migration does not draw; `physics/magnetic_geometry.py`). A column already there is skipped, so a stopped run
can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0229"
down_revision = "0228"
branch_labels = None
depends_on = None

COLUMNS = (
    ("planets", "dipole_tilt_deg", "DOUBLE NULL"), ("planets", "dipole_offset_class", "VARCHAR(8) NULL"),
    ("moons", "dipole_tilt_deg", "DOUBLE NULL"), ("moons", "dipole_offset_class", "VARCHAR(8) NULL"),
    ("rogue_planets", "dipole_tilt_deg", "DOUBLE NULL"), ("rogue_planets", "dipole_offset_class", "VARCHAR(8) NULL"),
    ("stars", "surface_dipole_field_gauss", "DOUBLE NULL"),
    ("neutron_stars", "magnetic_inclination_deg", "DOUBLE NULL"),
)
"""tuple: `(table, column, definition)` of every column the migration adds."""


def upgrade():
    connection = op.get_bind()
    for table, column, definition in COLUMNS:
        if not helpers.exists_column(connection, table, column):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} {definition}"))
