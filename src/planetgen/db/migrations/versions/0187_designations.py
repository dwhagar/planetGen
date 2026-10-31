"""The designations an object has stopped having, kept for ever (GEN.1138).

Schema v187. `designations` holds one row per retired designation: a retired object ID (`uid`, unique: an ID is never given twice) or a retired
name (`name`), the object that held it (`object_kind`, the table, and `object_row_id`, its `id`), the epoch it began in (NULL for a birth
designation) and the one it stopped in, why (`reason`), and `gone_epoch`, NULL while the object exists. The table is created empty: an object that
never changes its name or ID costs no row, and nothing is backfilled. A stopped run can be started again.
"""

from alembic import op

revision = "0187"
down_revision = "0186"
branch_labels = None
depends_on = None

_CREATE = """
CREATE TABLE IF NOT EXISTS designations (
    id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    designation_kind  VARCHAR(8) NOT NULL,
    uid               BINARY(12) NULL,
    name              VARCHAR(255) NULL,
    object_kind       VARCHAR(32) NOT NULL,
    object_row_id     BIGINT UNSIGNED NOT NULL,
    valid_from_epoch  SMALLINT UNSIGNED NULL,
    retired_epoch     SMALLINT UNSIGNED NOT NULL,
    reason            VARCHAR(16) NOT NULL,
    gone_epoch        SMALLINT UNSIGNED NULL,

    UNIQUE KEY uq_designations_uid (uid),
    KEY idx_designations_name (name),
    KEY idx_designations_object (object_kind, object_row_id, retired_epoch),
    CONSTRAINT chk_designations_kind CHECK (
        (designation_kind = 'uid' AND uid IS NOT NULL AND name IS NULL) OR (designation_kind = 'name' AND name IS NOT NULL AND uid IS NULL)),
    CONSTRAINT chk_designations_reason CHECK (reason IN ('collision bump', 'rename', 'regeneration', 'parent change', 'admin renumber', 'deleted'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
"""


def upgrade():
    op.execute(_CREATE)
