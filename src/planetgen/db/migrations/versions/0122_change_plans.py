"""Change plans, their operations and the edit journal (DB.36).

Schema v122. The tables of the 8.6 headline (docs/design/interface/edit-everything.md): every admin change is a plan, an ordered list of
operations on stored objects, approved before anything is written, then applied as one job and journaled.

`change_plans` is one row per plan: who asked and in which role (`requester_id` is the account's id in the control database,
`requester_name` its name when asked, as the account can be removed), the request in the requester's words, the `state`
(draft, requested, proposed, approved, applied, undone, stale, rejected or withdrawn), the role that must approve it
(`approver_role`, NULL when the admin confirms their own plan), who approved or rejected it and when with a note, and the job
that applied it (`apply_job_id`, a control database `work_jobs` id). The requester, approver and state fields are here from the
start so that the 10.0 roles (a request is a plan) add only comments and notices (DB.39) and no rebuild.
`change_plan_ops` is one row per operation: the plan, its order, the operation `kind`, the target object ID, the values before
and after as JSON, and its status. `edit_journal` is the append-only list of every applied operation: plan, operation, object ID,
sector, aspect, before, after and the generator epoch the sector was under, indexed by object ID and by sector.

This revision makes the three tables and nothing else; no code reads them yet. A table already there is left alone, so a stopped
run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0122"
down_revision = "0121"
branch_labels = None
depends_on = None

TABLES = {
    "change_plans": """
CREATE TABLE IF NOT EXISTS change_plans (
    id               BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    requester_id     BIGINT UNSIGNED NULL,
    requester_name   VARCHAR(255) NOT NULL,
    requester_role   VARCHAR(16) NOT NULL,
    request_text     TEXT NULL,
    state            VARCHAR(16) NOT NULL DEFAULT 'draft',
    approver_role    VARCHAR(16) NULL,
    decided_by_id    BIGINT UNSIGNED NULL,
    decided_by_name  VARCHAR(255) NULL,
    decided_at       TIMESTAMP(3) NULL,
    decision_note    TEXT NULL,
    apply_job_id     BIGINT UNSIGNED NULL,
    created_at       TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    updated_at       TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
    KEY idx_change_plans_state (state, id),
    KEY idx_change_plans_requester (requester_id, id),
    CONSTRAINT chk_change_plans_state CHECK (state IN
        ('draft', 'requested', 'proposed', 'approved', 'applied', 'undone', 'stale', 'rejected', 'withdrawn')),
    CONSTRAINT chk_change_plans_requester_role CHECK (requester_role IN ('user', 'moderator', 'admin')),
    CONSTRAINT chk_change_plans_approver_role CHECK (approver_role IS NULL OR approver_role IN ('moderator', 'admin'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
""",
    "change_plan_ops": """
CREATE TABLE IF NOT EXISTS change_plan_ops (
    id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    plan_id     BIGINT UNSIGNED NOT NULL,
    op_order    INT NOT NULL,
    kind        VARCHAR(48) NOT NULL,
    target_uid  BINARY(12) NULL,
    before_json JSON NULL,
    after_json  JSON NULL,
    status      VARCHAR(16) NOT NULL DEFAULT 'pending',
    UNIQUE KEY uq_change_plan_ops_order (plan_id, op_order),
    KEY idx_change_plan_ops_target (target_uid),
    CONSTRAINT fk_change_plan_ops_plan FOREIGN KEY (plan_id) REFERENCES change_plans(id) ON DELETE CASCADE,
    CONSTRAINT chk_change_plan_ops_status CHECK (status IN ('pending', 'applied', 'undone', 'skipped', 'failed'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
""",
    "edit_journal": """
CREATE TABLE IF NOT EXISTS edit_journal (
    id               BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    plan_id          BIGINT UNSIGNED NOT NULL,
    op_id            BIGINT UNSIGNED NULL,
    object_uid       BINARY(12) NOT NULL,
    sector_id        BIGINT UNSIGNED NULL,
    aspect           VARCHAR(64) NOT NULL,
    before_json      JSON NULL,
    after_json       JSON NULL,
    generator_epoch  SMALLINT UNSIGNED NOT NULL,
    applied_at       TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    KEY idx_edit_journal_object (object_uid, id),
    KEY idx_edit_journal_sector (sector_id, id),
    KEY idx_edit_journal_plan (plan_id, id),
    CONSTRAINT fk_edit_journal_plan FOREIGN KEY (plan_id) REFERENCES change_plans(id)
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
