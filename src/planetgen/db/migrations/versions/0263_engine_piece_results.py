"""The job engine keeps what each finished piece of a unit gave (OPS.219).

Schema v263. One table, `engine_piece_results`: one row for each piece of a job that returned something, in the job's own database, written in the same transaction as the cursor in `engine_items` and the
data the piece made. A unit that stops part way and is run again starts after its last stored cursor (`queue/engine.py`) and hands its `gather` the stored results beside the new ones, so the
unit's result is the same as if it had never stopped. `job_id` is the id of the job's row in the control database's `engine_jobs`; `piece_index` counts from 0 in the unit's own order; `result` is the
piece's JSON-able result. A table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0263"
down_revision = "0262"
branch_labels = None
depends_on = None

RESULTS = """
    CREATE TABLE IF NOT EXISTS engine_piece_results (
        job_id        BIGINT UNSIGNED NOT NULL,
        piece_index   INT UNSIGNED NOT NULL,
        result        MEDIUMTEXT NOT NULL,
        PRIMARY KEY (job_id, piece_index)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def upgrade():
    op.get_bind().execute(sa.text(RESULTS))
