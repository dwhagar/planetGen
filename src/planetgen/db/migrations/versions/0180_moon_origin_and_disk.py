"""A moon's origin and its giant parent's disk (DB.70).

Schema v180. `planets` and `moons` gain `moon_origin` (how a moon came to orbit its parent: `disk`, `impact`, `captured`, `ring` or `binary`,
`physics/moon_composition.py`) and the circumplanetary disk of a giant planet that its moons formed in: `moon_disk_lifetime_myr`,
`moon_disk_alpha_over_f`, `moon_disk_f_inf` and `moon_disk_ice_line_km`.

All NULL here: the origin and the disk are draws of the generator epoch, so a row made before v180 has nothing to read back (the regenerated
galaxy replaces it). Each column is skipped when it is already there, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0180"
down_revision = "0179"
branch_labels = None
depends_on = None

_COLUMNS = (
    ("moon_origin", "VARCHAR(8) NULL"),
    ("moon_disk_lifetime_myr", "DOUBLE NULL"),
    ("moon_disk_alpha_over_f", "DOUBLE NULL"),
    ("moon_disk_f_inf", "DOUBLE NULL"),
    ("moon_disk_ice_line_km", "DOUBLE NULL"),
)
"""tuple: The `(column, definition)` pairs added to both tables."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        for column, definition in _COLUMNS:
            if not helpers.exists_column(connection, table, column):
                op.execute(f"ALTER TABLE {table} ADD COLUMN {column} {definition}")
