"""The block ladder runs up to a root of three blocks, and the stored block aggregates are made again (MAP.258).

Schema v146. The ladder (`galaxy/drill.py`) was 243, 81, 27, 9 and 3 shells deep; the single-space map walks nodes
from a root down, so it now starts at 177,147 and every size is three times its child: 177,147, 59,049, 19,683, 6,561,
2,187, 729, then the old five. The root's one block shell is one row of three wedges and holds the whole grid (the
outermost shell is 65,535). A level's rows and wedges are the parent's count times a rounded ratio, so putting levels
above 243 moves the rows and wedges of every level below: a block `(243, shell, band, wedge)` is another block than it
was. So the stored aggregates of the old blocks mean nothing:

* `galaxy_block_stats.m` becomes INT UNSIGNED (it was SMALLINT UNSIGNED, which stops at 65,535), and the table is
  emptied and filled again from the generated sectors with a row for each of the eleven levels (the pass migration
  0092 made, a shell at a time);
* when the nebula block tables were filled (`nebula_block_state` holds its row), they are made again from the nebulae
  for the five levels they keep (`db.nebula_blocks.rebuild`).

A column that is already INT UNSIGNED is left alone and the tables are filled again either way, so a stopped run can
be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen import tuning
from planetgen.db import alembic_runner
from planetgen.db.migrations import helpers

revision = "0146"
down_revision = "0145"
branch_labels = None
depends_on = None

BATCH = 1000


def upgrade():
    connection = op.get_bind()
    if str(helpers.column_type(connection, "galaxy_block_stats", "m")).lower() != "int unsigned":
        op.execute("ALTER TABLE galaxy_block_stats MODIFY COLUMN m INT UNSIGNED NOT NULL")
    connection.execute(sa.text("DELETE FROM galaxy_block_stats"))
    _fill(connection)
    if connection.execute(sa.text("SELECT COUNT(*) FROM nebula_block_state")).scalar():
        from planetgen.db import nebula_blocks

        edge = connection.execute(sa.text("SELECT edge_pc FROM galaxy_shape WHERE id = 1")).scalar()
        nebula_blocks.rebuild(helpers.StoreConnection(connection),
                              float(edge) if edge else float(tuning.DEFAULT_SECTOR_EDGE_PC))


def _fill(connection):
    """Computes every block's sums from the generated sectors (`db.block_stats.accumulate`)."""
    from planetgen.db import block_stats

    shells = [row[0] for row in connection.execute(sa.text(
        "SELECT DISTINCT shell_index FROM sectors WHERE shell_index IS NOT NULL ORDER BY shell_index")).fetchall()]
    cells = {}
    for done, shell in enumerate(shells, start=1):
        rows = connection.execute(sa.text(
            "SELECT s.shell_index, s.band_index, s.wedge_index, st.actual_systems, st.expected_systems,"
            " st.actual_stars, st.mean_age_gy, st.total_luminosity_sol, st.habitable_worlds FROM (SELECT * FROM sectors WHERE"
            " shell_index = :shell) s LEFT JOIN sector_stats st ON st.shell_index = s.shell_index"
            " AND st.band_index = s.band_index AND st.wedge_index = s.wedge_index"), {"shell": shell})
        for r in rows.fetchall():
            block_stats.accumulate(cells, (r[0], r[1], r[2]), *r[3:])
        alembic_runner.report_progress("Block totals", done, len(shells))
    insert = sa.text(
        "INSERT INTO galaxy_block_stats (" + ", ".join(block_stats.COLUMNS) + ") VALUES ("
        + ", ".join(":" + column for column in block_stats.COLUMNS) + ")")
    rows = [dict(zip(block_stats.COLUMNS, row)) for row in block_stats.cell_rows(cells)]
    for start in range(0, len(rows), BATCH):
        connection.execute(insert, rows[start:start + BATCH])
