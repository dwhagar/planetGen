"""The ring system of a planet (GEN.732).

Schema v160. `planets` and `moons` gain `ring_kind` ('dusty', 'narrow' or 'broad'), `ring_origin` ('moon_erosion', 'disrupted_moon' or
'captured_body'), `ring_inner_km`, `ring_outer_km`, `ring_mass_kg`, `ring_optical_depth` (a mean over the area the ring covers) and
`ring_ice_fraction`: one summed ring system per body (`planetgen.physics.rings`).

All NULL here: the rings are drawn with the planet at the generator epoch, so a row made before v160 has none and the regenerated galaxy
replaces it. Each column is skipped when it is already there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0160"
down_revision = "0159"
branch_labels = None
depends_on = None

_COLUMNS = (
    ("ring_kind", "VARCHAR(8) NULL CHECK (ring_kind IN ('dusty', 'narrow', 'broad'))"),
    ("ring_origin", "VARCHAR(16) NULL CHECK (ring_origin IN ('moon_erosion', 'disrupted_moon', 'captured_body'))"),
    ("ring_inner_km", "DOUBLE NULL"),
    ("ring_outer_km", "DOUBLE NULL"),
    ("ring_mass_kg", "DOUBLE NULL"),
    ("ring_optical_depth", "DOUBLE NULL"),
    ("ring_ice_fraction", "DOUBLE NULL"),
)


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        changes = [f"ADD COLUMN {column} {definition}" for column, definition in _COLUMNS
                   if not helpers.exists_column(connection, table, column)]
        if changes:  # no ALGORITHM=INSTANT: both tables carry a FULLTEXT index, which refuses it
            connection.execute(sa.text(f"ALTER TABLE {table} {', '.join(changes)}"))
