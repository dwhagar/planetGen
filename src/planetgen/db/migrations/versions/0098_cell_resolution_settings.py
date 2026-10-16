"""The cell star ceiling and the fill luminosity cap, two settings of the galaxy (GEN.281).

Schema v98. `galaxy_shape.cell_star_ceiling` (DOUBLE, default 50,000 systems) and `galaxy_shape.fill_luminosity_cap_sol`
(DOUBLE, default 100 solar luminosities): a cell that expects more systems than the ceiling is filled only with the
stars at or above the cap (`galaxy/resolution.py`). The defaults are the code's own, so a galaxy planned before this
keeps the behaviour it had: no cell of the Milky Way's bar comes near 50,000. Columns already there are left alone, so
a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0098"
down_revision = "0097"
branch_labels = None
depends_on = None

COLUMNS = (
    ("cell_star_ceiling", "DOUBLE NOT NULL DEFAULT 50000", "reach_shell_index"),
    ("fill_luminosity_cap_sol", "DOUBLE NOT NULL DEFAULT 100", "cell_star_ceiling"),
)


def _has_column(connection, name):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'galaxy_shape' AND COLUMN_NAME = :name"), {"name": name}).scalar())


def upgrade():
    connection = op.get_bind()
    for name, definition, after in COLUMNS:
        if not _has_column(connection, name):
            op.execute(f"ALTER TABLE galaxy_shape ADD COLUMN {name} {definition} AFTER {after}")
