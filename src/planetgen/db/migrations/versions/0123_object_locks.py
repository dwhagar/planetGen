"""Lock marks: a standing keep mark on an object or one aspect of it (DB.37).

Schema v123. A lock is a standing keep mark (docs/design/interface/edit-everything.md): a regeneration or a plan keeps a locked object, or a locked aspect of it, without the admin marking it again, and a plan refuses to change a locked aspect unless it lifts the lock. `object_locks` holds one row per object and aspect: `object_uid` (the 96-bit object ID), `aspect` (a name from the aspect catalogue, or `*` for the whole object), who locked it (`locked_by_id`, the account's id in the control database, and `locked_by_name` as it was) and when, and an optional `reason`. An object has at most one lock per aspect, and a `*` lock beside aspect locks is allowed (the whole object holds anyway).

This revision makes the table and nothing else; no code reads it yet. A table already there is left alone,
so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0123"
down_revision = "0122"
branch_labels = None
depends_on = None

TABLES = {
    "object_locks": """
CREATE TABLE IF NOT EXISTS object_locks (
    object_uid       BINARY(12) NOT NULL,
    aspect           VARCHAR(64) NOT NULL DEFAULT '*',
    locked_by_id     BIGINT UNSIGNED NULL,
    locked_by_name   VARCHAR(255) NOT NULL,
    locked_at        TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    reason           VARCHAR(500) NULL,
    PRIMARY KEY (object_uid, aspect),
    KEY idx_object_locks_locked_at (locked_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
""",
}


def upgrade():
    connection = op.get_bind()
    for name, create in TABLES.items():
        if not connection.execute(sa.text(
                "SELECT COUNT(*) FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()"
                " AND TABLE_NAME = :name"), {"name": name}).scalar():
            op.execute(create)
