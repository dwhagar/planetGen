"""Take out the object timing samples and the timing histogram (PERF.197).

Schema v195. PERF.92 saved a repeatable 1-in-1,000 sample of single objects with their full timing (`object_timing_samples`) and a
24-bucket histogram for each kind (`object_timing_stats.histogram`). PERF.196 found that nothing read either: the admin Stats page
shows what the per-kind count, total, fastest and slowest give. The table and the column go, with the code that wrote them. A table
or column already gone is skipped, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0195"
down_revision = "0194"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    if helpers.exists_table(connection, "object_timing_samples"):
        op.execute("DROP TABLE object_timing_samples")
    helpers.drop_column_instant(connection, "object_timing_stats", "histogram")
