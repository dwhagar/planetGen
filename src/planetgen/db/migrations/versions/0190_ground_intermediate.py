"""Intermediate ground on every planet and moon (GEN.474) and the chondrite group of every belt and field (GEN.475).

Schema v190. `planets` and `moons` gain `ground_intermediate_fraction`, the areal share of the surface that is intermediate (andesitic) ground: the volcanic crust
of a body whose mantle melt is 52 to 63 percent silica (`physics/rock_types.py`, wired into `physics/surface_composition.py`).

Also `asteroid_belts` and `asteroid_fields` gain `chondrite_group` (VARCHAR(4)), the nearest meteorite group of the body's mean composition (`physics/meteorite_groups.py`, `none` far from every group).

NULL here: a row made before v190 has no such share, and its ground needs the draws of the generator epoch, so the regenerated galaxy replaces it. The column is
skipped when it is already there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0190"
down_revision = "0189"
branch_labels = None
depends_on = None

_COLUMNS = (
    ("planets", "ground_intermediate_fraction", "DOUBLE NULL"),
    ("moons", "ground_intermediate_fraction", "DOUBLE NULL"),
    ("asteroid_belts", "chondrite_group", "VARCHAR(4) NULL"),
    ("asteroid_fields", "chondrite_group", "VARCHAR(4) NULL"),
)
"""tuple: The `(table, column, definition)` added."""


def upgrade():
    connection = op.get_bind()
    for table, column, definition in _COLUMNS:
        if not helpers.exists_column(connection, table, column):    # a plain ALTER: `planets` and `moons` carry a FULLTEXT key, which refuses INSTANT
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} {definition}"))
