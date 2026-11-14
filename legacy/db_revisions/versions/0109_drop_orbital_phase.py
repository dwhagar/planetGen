"""A planet or moon is its state vector: the orbital phase column goes (GEN.289).

Schema v109. The position, velocity, mass, epoch, parent and dynamical state of v102 hold a planet or moon whole, and the
nightly update moves them with `kepler.universal_step`, so `orbital_phase_deg` on `planets` and `moons` was a second
record of the same place on a circle. Position and phase were always written together, so nothing is lost.

The parent of a circumbinary planet (a close pair's, with no `star_id`) is its system: v102 left such rows pointing at a
star that does not exist. A column is dropped only where it is still there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0109"
down_revision = "0108"
branch_labels = None
depends_on = None


def _has_column(connection, table, column):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar())


def upgrade():
    connection = op.get_bind()
    op.execute("UPDATE planets SET parent_type = 'system', parent_id = star_system_id WHERE star_id IS NULL")
    for table in ("planets", "moons"):
        if _has_column(connection, table, "orbital_phase_deg"):
            op.execute(f"ALTER TABLE {table} DROP COLUMN orbital_phase_deg")
