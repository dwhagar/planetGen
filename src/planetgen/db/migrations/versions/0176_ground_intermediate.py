"""Intermediate ground on every planet and moon (GEN.474).

Schema v176. `planets` and `moons` gain `ground_intermediate_fraction`, the areal share of the surface that is intermediate (andesitic) ground: the volcanic crust
of a body whose mantle melt is 52 to 63 percent silica (`physics/rock_types.py`, wired into `physics/surface_composition.py`).

NULL here: a row made before v176 has no such share, and its ground needs the draws of the generator epoch, so the regenerated galaxy replaces it. The column is
skipped when it is already there, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0176"
down_revision = "0175"
branch_labels = None
depends_on = None

_COLUMN = ("ground_intermediate_fraction", "DOUBLE NULL")
"""tuple: The `(column, definition)` added to both tables."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        helpers.add_column_instant(connection, table, *_COLUMN)
