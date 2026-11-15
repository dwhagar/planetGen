"""The job engine's per-galaxy ledger half (OPS.215).

Schema v261. One table, `engine_items`: one row per job that writes to this galaxy, holding the number of the attempt the supervisor handed it (`attempt`), the place it has reached (`cursor_key`) and the
time it finished (`done_at`). A job's own commits update its row `WHERE attempt = mine`, and the last one sets `done_at`, so a job's data and its done mark are one transaction (the fence of
`docs/design/operations/job-engine.md`). `job_id` is the id of the job's row in the control database's `engine_jobs`; the two databases are not joined by a key. Nothing writes the table until the supervisor
(OPS.207) and the fenced commit (OPS.223) are merged, and a galaxy made before v261 has no rows. A table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0261"
down_revision = "0260"
branch_labels = None
depends_on = None

ITEMS = """
    CREATE TABLE IF NOT EXISTS engine_items (
        job_id        BIGINT UNSIGNED NOT NULL PRIMARY KEY,
        attempt       INT UNSIGNED NOT NULL DEFAULT 0,
        cursor_key    VARCHAR(255) NULL,
        done_at       DATETIME(6) NULL,
        updated_at    DATETIME(6) NOT NULL
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def upgrade():
    op.get_bind().execute(sa.text(ITEMS))
