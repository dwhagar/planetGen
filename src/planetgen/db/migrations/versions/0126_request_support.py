"""Edit request support tables: request comments and notification records (DB.39).

Schema v126. Boss's 10.0 roles (docs/design/roles-and-requests.md): users submit edit requests, moderators approve user requests and admins approve moderator requests. The requester, the required approver role and the decision are columns of `change_plans` (DB.36), so a request is a plan; these two tables hold what the plan row does not. `request_comments` is the conversation on a request: the plan, who wrote it (`author_id`, the account's id in the control database, `author_name` and `author_role` as they were), the `body` and when. `request_notifications` is one row per notice owed to an account about a request: the plan, the `recipient_id`, the `kind` of notice (a comment, a decision) and when it was made and when it was `sent` (NULL until it goes). Both go with their plan.

This revision makes the tables and nothing else; no code reads them yet. A table already there is left alone,
so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0126"
down_revision = "0125"
branch_labels = None
depends_on = None

TABLES = {
    "request_comments": """
CREATE TABLE IF NOT EXISTS request_comments (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    plan_id       BIGINT UNSIGNED NOT NULL,
    author_id     BIGINT UNSIGNED NULL,
    author_name   VARCHAR(255) NOT NULL,
    author_role   VARCHAR(16) NOT NULL,
    body          TEXT NOT NULL,
    created_at    TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    KEY idx_request_comments_plan (plan_id, id),
    CONSTRAINT fk_request_comments_plan FOREIGN KEY (plan_id) REFERENCES change_plans(id) ON DELETE CASCADE,
    CONSTRAINT chk_request_comments_author_role CHECK (author_role IN ('user', 'moderator', 'admin'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
""",
    "request_notifications": """
CREATE TABLE IF NOT EXISTS request_notifications (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    plan_id       BIGINT UNSIGNED NOT NULL,
    recipient_id  BIGINT UNSIGNED NOT NULL,
    kind          VARCHAR(32) NOT NULL,
    created_at    TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    sent_at       TIMESTAMP(3) NULL,
    KEY idx_request_notifications_recipient (recipient_id, sent_at, id),
    KEY idx_request_notifications_plan (plan_id, id),
    CONSTRAINT fk_request_notifications_plan FOREIGN KEY (plan_id) REFERENCES change_plans(id) ON DELETE CASCADE
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
