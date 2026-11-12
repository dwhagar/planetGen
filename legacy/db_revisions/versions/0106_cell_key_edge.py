"""The cell key holds the last cube at the edge of its box (DB.65).

Schema v106. A coordinate just inside the box (`ABS(c) < 262144`) whose rounded cube number came to 2^21 (the last
double below 262,144 pc, and its neighbours) was put in cube 2^21 by the SQL, which the key's 21 bits cannot hold: the
weave dropped the bit and the point took the key of the opposite edge. Both the SQL (`cellkey.sql_cell_key`) and the
Python (`cellkey.cell_index`, `cell_keys`) now hold the cube number at 2^21 - 1.

Every generated `cell_key` column is recomputed with the new expression (a stored generated column's expression is
changed by a table rebuild, which also rebuilds its covering index). A column that already holds the new expression is
left alone, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner
from planetgen.galaxy import cellkey

revision = "0106"
down_revision = "0105"
branch_labels = None
depends_on = None

PC_COLUMNS = ("center_x_pc", "center_y_pc", "center_z_pc")
MPC_COLUMNS = ("position_x_mpc", "position_y_mpc", "position_z_mpc")
GALAXY_COLUMNS = ("galaxy_x_pc", "galaxy_y_pc", "galaxy_z_pc")

CELL_TABLES = (
    # table, the SQL for its x, y, z in parsecs
    ("black_holes", PC_COLUMNS),
    ("neutron_stars", PC_COLUMNS),
    ("nebulae", PC_COLUMNS),
    ("supernova_remnants", PC_COLUMNS),
    ("quasars", PC_COLUMNS),
    ("rogue_planets", PC_COLUMNS),
    ("interstellar_comets", PC_COLUMNS),
    ("asteroid_fields", PC_COLUMNS),
    ("facilities", PC_COLUMNS),
    ("nebula_field", PC_COLUMNS),
    ("galaxy_components", PC_COLUMNS),
    ("bright_stars", tuple(f"{c} / 1000.0" for c in MPC_COLUMNS)),
    ("phenomenon_scatter", tuple(f"{c} / 1000.0" for c in MPC_COLUMNS)),
    ("star_systems", GALAXY_COLUMNS),
)
"""tuple: `(table, SQL of its x, y, z in parsecs)` for every table that has a `cell_key`."""


def _expression(connection, table):
    return connection.execute(sa.text(
        "SELECT GENERATION_EXPRESSION FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = 'cell_key'"), {"table_name": table}).scalar()


def upgrade():
    connection = op.get_bind()
    for done, (table, expressions) in enumerate(CELL_TABLES):
        current = _expression(connection, table)
        if current is not None and "least(" in current.lower():
            continue
        alembic_runner.report_progress("Cell keys", done, len(CELL_TABLES))
        key = cellkey.sql_cell_key(*expressions)
        op.execute(f"ALTER TABLE {table} MODIFY COLUMN cell_key BIGINT UNSIGNED GENERATED ALWAYS AS ({key}) STORED")
    alembic_runner.report_progress("Cell keys", len(CELL_TABLES), len(CELL_TABLES))
