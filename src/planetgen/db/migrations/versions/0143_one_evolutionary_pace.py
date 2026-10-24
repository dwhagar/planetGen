"""A world's evolutionary pace is stored once: `evolutionary_speed` is gone.

Schema v143. `planets` and `moons` kept two draws of the same fact: `evolutionary_speed` (what the life sentence said) and
`evolutionary_pace` (what the timeline was built on), so a page could say "a normal pace" in one paragraph and "hyper-accelerated"
in the next. The generator now draws the pace once and stores it in `evolutionary_pace`, for every body with a life chemistry.

A body that has a speed and no pace (a moon, or a world outside the habitable zone) takes its speed as its pace: it is the one
pace its life sentence has always said. A body with both keeps its pace, which its stored ages and milestone were built on. Then
`evolutionary_speed` is dropped. A column already gone is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0143"
down_revision = "0142"
branch_labels = None
depends_on = None

TABLES = ("planets", "moons")
"""tuple: The tables that lose `evolutionary_speed`."""


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        if not helpers.exists_column(connection, table, "evolutionary_speed"):
            continue
        helpers.batched_update(
            connection, table, "evolutionary_pace = LOWER(evolutionary_speed)",
            "evolutionary_pace IS NULL AND LOWER(evolutionary_speed) IN ('fast', 'normal', 'slow')",
            label=f"{table} pace")
        # Not INSTANT: `planets` has a FULLTEXT index, and InnoDB refuses an instant change beside one (as in 0111).
        connection.execute(sa.text(f"ALTER TABLE {table} DROP COLUMN evolutionary_speed"))
