"""The scatter tables hold a band number as INT (DB.72).

Schema v108. A shell's bands run from `-(h + 1)` to `h + 1`, and from shell 20,861 outward `h + 1` is 32,768 or more
(`spherical.band_range`: +-102,942 at shell 65,535), which the SMALLINT `band_index` of `bright_stars`, `sector_stats`
and `phenomenon_scatter` cannot hold. The whole-sphere scatter stopped at 97 percent of its shells with "Out of range
value for column 'band_index'". The three columns become INT, as `sectors.band_index` and the block ladder already are,
so no column named for a band is narrower than INT.

A column changes in place by a table rebuild (its indexes with it) and no value is touched. A column that is already INT
is left alone, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner

revision = "0108"
down_revision = "0107"
branch_labels = None
depends_on = None

TABLES = ("bright_stars", "sector_stats", "phenomenon_scatter")


def _column_type(connection, table):
    return connection.execute(sa.text(
        "SELECT COLUMN_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = 'band_index'"), {"table_name": table}).scalar()


def upgrade():
    connection = op.get_bind()
    for done, table in enumerate(TABLES):
        if str(_column_type(connection, table)).lower() == "int":
            continue
        alembic_runner.report_progress("Band numbers", done, len(TABLES))
        op.execute(f"ALTER TABLE {table} MODIFY COLUMN band_index INT NOT NULL")
    alembic_runner.report_progress("Band numbers", len(TABLES), len(TABLES))
