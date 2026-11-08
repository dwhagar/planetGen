"""A sector row for every cell the scatters left something in (MAP.402).

Schema v227, no schema change. A scatter keeps only the address of the cell each bright star and scattered phenomenon lies in, so a cell nothing was
generated in had no `sectors` row unless a drifting rogue planet landed there. The migration makes the placeholder row (`db.uncharted_index.index`: a
placed sector with no cell in the chart ledger) for every cell that holds a waiting bright star or an unbuilt phenomenon and has none, so the Sectors
list, the search and the sector pages reach what the long-range scans detected; a scatter does the same when it ends. A second run adds nothing, so a
stopped migration can be started again.
"""

from planetgen.db.migrations import helpers

revision = "0227"
down_revision = "0226"
branch_labels = None
depends_on = None


def upgrade():
    from alembic import op
    from planetgen.db import uncharted_index

    uncharted_index.index(helpers.StoreConnection(op.get_bind()))
