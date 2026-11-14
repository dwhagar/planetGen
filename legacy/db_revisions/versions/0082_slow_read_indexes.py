"""Indexes and a summary table for the pages that timed out on a big galaxy (PERF.68, PERF.70, PERF.74).

Schema v82. `phenomenon_scatter (kind, subtype, mass_solar)` serves the Galaxy
Map's coarse tiles; `star_systems (quadrant, name, id)` replaces the
quadrant-only index and `(is_binary, name, id)` is new, for the Systems list's
sorts; `sector_system_counts` is filled once here (one GROUP BY on the sector
index) and kept by `query.refresh_sector_system_counts`. Each step is skipped
when it is already there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner

revision = "0082"
down_revision = "0081"
branch_labels = None
depends_on = None

STEPS = 4


def _index_columns(connection, table, index):
    rows = connection.execute(sa.text(
        "SELECT COLUMN_NAME FROM information_schema.STATISTICS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :t AND INDEX_NAME = :i ORDER BY SEQ_IN_INDEX"), {"t": table, "i": index}).fetchall()
    return [row[0] for row in rows]


def _has_table(connection, table):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = :t"),
        {"t": table}).scalar()


def upgrade():
    connection = op.get_bind()
    alembic_runner.report_progress("steps", 0, STEPS)
    if _index_columns(connection, "phenomenon_scatter", "idx_phenomenon_scatter_class") != ["kind", "subtype", "mass_solar"]:
        op.execute("ALTER TABLE phenomenon_scatter ADD KEY idx_phenomenon_scatter_class (kind, subtype, mass_solar)")
    alembic_runner.report_progress("steps", 1, STEPS)
    if _index_columns(connection, "star_systems", "idx_star_systems_quadrant") != ["quadrant", "name", "id"]:
        op.execute("ALTER TABLE star_systems DROP KEY idx_star_systems_quadrant, "
                   "ADD KEY idx_star_systems_quadrant (quadrant, name, id)")
    alembic_runner.report_progress("steps", 2, STEPS)
    if not _index_columns(connection, "star_systems", "idx_star_systems_binary_name"):
        op.execute("ALTER TABLE star_systems ADD KEY idx_star_systems_binary_name (is_binary, name, id)")
    alembic_runner.report_progress("steps", 3, STEPS)
    if not _has_table(connection, "sector_system_counts"):
        op.execute("""
            CREATE TABLE sector_system_counts (
                sector_id     BIGINT UNSIGNED NOT NULL PRIMARY KEY,
                system_count  INT UNSIGNED NOT NULL,
                KEY idx_sector_system_counts_count (system_count, sector_id)
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""")
        op.execute("INSERT INTO sector_system_counts (sector_id, system_count) "
                   "SELECT sector_id, COUNT(*) FROM star_systems WHERE sector_id IS NOT NULL GROUP BY sector_id")
    alembic_runner.report_progress("steps", STEPS, STEPS)
