"""The nuclear star cluster and the nuclear stellar disc replace the Gaussian core (GEN.280).

Schema v99. `galaxy_shape.core_amplitude` (ADM.49's Gaussian on top of the bulge) is renamed `nuclear_amplitude`: how
much of the observed nuclear star cluster and nuclear stellar disc the centre holds, 1 being the observed pair scaled
to the bulge and 0 none (`galaxy/density.py`). A galaxy planned before this keeps its density, which had no such
structures, so every row is set to 0. `bright_stars.population` (VARCHAR(12), too short for `nuclear_cluster`) becomes
VARCHAR(16), an in-place change that copies no row, and takes the two new populations, `nuclear_cluster` and
`nuclear_disc`; the old values stay valid. A step already done is left alone, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0099"
down_revision = "0098"
branch_labels = None
depends_on = None

POPULATIONS = ("young", "intermediate", "old", "bulge", "halo", "cluster", "nuclear_cluster", "nuclear_disc")


def _has_column(connection, name):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'galaxy_shape' AND COLUMN_NAME = :name"), {"name": name}).scalar())


def upgrade():
    connection = op.get_bind()
    clause = connection.execute(sa.text(
        "SELECT CHECK_CLAUSE FROM information_schema.CHECK_CONSTRAINTS WHERE CONSTRAINT_SCHEMA = DATABASE()"
        " AND CONSTRAINT_NAME = 'chk_bright_stars_population'")).scalar()
    if clause is None or "nuclear_disc" not in clause:
        listed = ", ".join(f"'{name}'" for name in POPULATIONS)
        drop = "DROP CHECK chk_bright_stars_population, " if clause is not None else ""
        op.execute(f"ALTER TABLE bright_stars MODIFY COLUMN population VARCHAR(16) NOT NULL, {drop}"
                   f"ADD CONSTRAINT chk_bright_stars_population CHECK (population IN ({listed}))")
    if _has_column(connection, "core_amplitude"):
        op.execute("ALTER TABLE galaxy_shape CHANGE COLUMN core_amplitude nuclear_amplitude DOUBLE NOT NULL DEFAULT 0")
        op.execute("UPDATE galaxy_shape SET nuclear_amplitude = 0")
