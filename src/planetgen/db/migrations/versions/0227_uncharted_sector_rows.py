"""A sector row for every cell the scatters left something in (MAP.402), ended by OPS.154.

Schema v227, no schema change. The migration made the placeholder row (a placed sector with no cell in the chart ledger) for every cell that held a waiting bright
star or an unbuilt phenomenon. OPS.154 reads those cells from the block light parts instead and migration 0254 deletes the rows, so there is nothing left for this
one to make: a galaxy that reaches 0227 now is one that never had them, and 0254 would delete them straight away.
"""

revision = "0227"
down_revision = "0226"
branch_labels = None
depends_on = None


def upgrade():
    pass
