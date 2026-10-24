"""A generated sector and a block of sectors keep the number of habitable worlds they hold (MAP.143).

Schema v139. `sector_stats` gains `habitable_worlds` (INT, NULL until a sector is filled: the planets and moons of its systems with an
overall PHI-6 score that is Blue or Green, `habitability.POTENTIALLY_HABITABLE_MIN_SCORE`) and `galaxy_block_stats` gains `habitable_worlds` (BIGINT UNSIGNED, NOT NULL,
default 0: their sum over the block's generated sectors), so the Galaxy Map can colour a sector or block by how many it holds. Every
filled sector is counted from the scores its planets and moons already hold a shell at a time, and the block totals are rebuilt from them. A column already
there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner
from planetgen.db.migrations import helpers
from planetgen.physics import habitability

revision = "0139"
down_revision = "0138"
branch_labels = None
depends_on = None

BATCH = 1000


def upgrade():
    connection = op.get_bind()
    helpers.add_column_instant(connection, "sector_stats", "habitable_worlds", "INT NULL", after="total_luminosity_sol")
    helpers.add_column_instant(connection, "galaxy_block_stats", "habitable_worlds", "BIGINT UNSIGNED NOT NULL DEFAULT 0",
                               after="luminosity_sol")
    _count_sectors(connection)
    _rebuild_blocks(connection)


def _count_sectors(connection):
    """Counts the habitable planets and moons of every filled sector, a shell at a time."""
    floor = habitability.POTENTIALLY_HABITABLE_MIN_SCORE
    shells = [row[0] for row in connection.execute(sa.text(
        "SELECT DISTINCT shell_index FROM sector_stats WHERE bright_level_sol = 0 ORDER BY shell_index")).fetchall()]
    for done, shell in enumerate(shells, start=1):
        connection.execute(sa.text(
            "UPDATE sector_stats st JOIN sectors s ON s.shell_index = st.shell_index AND s.band_index = st.band_index"
            " AND s.wedge_index = st.wedge_index SET st.habitable_worlds = ("
            " (SELECT COUNT(*) FROM planets p JOIN star_systems ss ON ss.id = p.star_system_id"
            "  WHERE ss.sector_id = s.id AND p.phi4 >= :floor)"
            " + (SELECT COUNT(*) FROM moons m JOIN star_systems ss ON ss.id = m.star_system_id"
            "  WHERE ss.sector_id = s.id AND m.phi4 >= :floor))"
            " WHERE st.shell_index = :shell AND st.bright_level_sol = 0"), {"shell": shell, "floor": floor})
        alembic_runner.report_progress("Habitable worlds", done, len(shells))


def _rebuild_blocks(connection):
    """Computes every block's sums again from the generated sectors (`db.block_stats.accumulate`)."""
    from planetgen.db import block_stats

    shells = [row[0] for row in connection.execute(sa.text(
        "SELECT DISTINCT shell_index FROM sectors WHERE shell_index IS NOT NULL ORDER BY shell_index")).fetchall()]
    cells = {}
    for done, shell in enumerate(shells, start=1):
        rows = connection.execute(sa.text(
            "SELECT s.shell_index, s.band_index, s.wedge_index, st.actual_systems, st.expected_systems,"
            " st.actual_stars, st.mean_age_gy, st.total_luminosity_sol, st.habitable_worlds FROM (SELECT * FROM sectors"
            " WHERE shell_index = :shell) s LEFT JOIN sector_stats st ON st.shell_index = s.shell_index"
            " AND st.band_index = s.band_index AND st.wedge_index = s.wedge_index"), {"shell": shell})
        for r in rows.fetchall():
            block_stats.accumulate(cells, (r[0], r[1], r[2]), *r[3:])
        alembic_runner.report_progress("Block totals", done, len(shells))
    connection.execute(sa.text("DELETE FROM galaxy_block_stats"))
    insert = sa.text(
        "INSERT INTO galaxy_block_stats (" + ", ".join(block_stats.COLUMNS) + ") VALUES ("
        + ", ".join(":" + column for column in block_stats.COLUMNS) + ")")
    # `m` is SMALLINT UNSIGNED here; the levels above 243 of the ladder as it is now (MAP.258) do not fit, and
    # migration 0146 widens the column and fills the table again, with every level.
    rows = [dict(zip(block_stats.COLUMNS, row)) for row in block_stats.cell_rows(cells) if row[0] <= 65535]
    for start in range(0, len(rows), BATCH):
        connection.execute(insert, rows[start:start + BATCH])
