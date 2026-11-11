"""The stamp of the code that wrote each system's light and ultraviolet ledger rows (GEN.1498).

Schema v244. `system_summary` gains `ledger_stamp` (INT UNSIGNED, nullable): the stamp (`db/derived_stamp.py`, family `ledger`) of the code that wrote the system's light and ultraviolet
rows of the effect ledger. `refill-effect-ledger` takes the systems whose stamp is NULL or not the running code's and leaves the others alone. Every system already stored starts with no
stamp, so the first run after this writes the rows of every system once; a column already there is skipped, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0244"
down_revision = "0243"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    helpers.add_column_instant(connection, "system_summary", "ledger_stamp", "INT UNSIGNED NULL", after="habitability_stamp")
