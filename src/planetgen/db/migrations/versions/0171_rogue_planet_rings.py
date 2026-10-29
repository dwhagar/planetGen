"""The ring system of a rogue planet (GEN.836).

Schema v171. `rogue_planets` gains the seven ring columns `planets` and `moons` got at v160 (`ring_kind`, `ring_origin`, `ring_inner_km`,
`ring_outer_km`, `ring_mass_kg`, `ring_optical_depth`, `ring_ice_fraction`), drawn with the planet by `planetgen.physics.rings.apply_rogue`.

All NULL here: a row made before v169 has no ring and the regenerated galaxy replaces it. Each column is skipped when it is already there,
so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0171"
down_revision = "0170"
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
    changes = [f"ADD COLUMN {column} {definition}" for column, definition in _COLUMNS
               if not helpers.exists_column(connection, "rogue_planets", column)]
    if changes:
        connection.execute(sa.text(f"ALTER TABLE rogue_planets {', '.join(changes)}"))
