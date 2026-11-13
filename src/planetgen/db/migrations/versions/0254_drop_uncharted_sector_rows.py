"""Delete the placeholder sector rows of uncharted cells (OPS.154).

Schema v254, no schema change. MAP.402 made a `sectors` row for every cell that holds a waiting bright star or an unbuilt scattered phenomenon and has no cell in
the chart ledger: up to 1.5 billion rows, about 600 GB at full size. The Sectors list, the search and the sector pages read those cells from the block light parts
now (`db.query.virtual_cells`), so the rows are deleted (`db.uncharted_index.clear_unoccupied`: 5,000 at a time, each batch committed). A row a rogue planet or a
system refers to stays; it is that body's address. A second run finds nothing, so a stopped migration can be started again.
"""

from planetgen.db import alembic_runner
from planetgen.db.migrations import helpers

revision = "0254"
down_revision = "0253"
branch_labels = None
depends_on = None


class _Committing(helpers.StoreConnection):
    """The revision's connection that commits each batch, so a delete of this size is not one transaction."""

    def commit(self):
        self.connection.commit()


def upgrade():
    from alembic import op
    from planetgen.db import uncharted_index
    from planetgen.util import log

    conn = _Committing(op.get_bind())
    held = conn.execute("SELECT TABLE_ROWS AS n FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'sectors'").fetchone()
    total = int(held["n"] or 0) if held else 0
    deleted = uncharted_index.clear_unoccupied(conn, on_phase=log.normal, on_progress=lambda done: alembic_runner.report_progress(
        "Placeholder sector rows deleted", done, max(total, done)))
    log.normal(f"Migration 0254: {deleted:,} placeholder sector rows deleted.")
