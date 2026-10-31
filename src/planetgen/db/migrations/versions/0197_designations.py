"""The designations an object stopped having (GEN.1138).

Schema v197. `designations` holds one row for each ID or name an object used to have: the retired text (`designation_uid` for an ID,
`designation_name` for a name), the object that held it (`object_kind`, the table, and `object_row_id`, its row number, which never
changes), the epoch it began in (NULL for the birth designation) and the epoch it stopped in, why (`reason`, five values), and, for an
object since deleted, the epoch it was deleted in (`gone_epoch`). A retired ID is reserved for ever: a unique index on `designation_uid`
and the objects' own unique `uid` columns together keep it from being given to any object again. A retired name is not reserved (the
collision bump hands a stored object's name to another by design), so `designation_name` has a plain index. No row is made for a
designation that changed before this table; nothing is backfilled. A table already there is kept, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0197"
down_revision = "0196"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    if helpers.exists_table(connection, "designations"):
        return
    op.execute("""
        CREATE TABLE designations (
            id                BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
            designation_kind  VARCHAR(4) NOT NULL,
            designation_uid   BINARY(12) NULL,
            designation_name  VARCHAR(255) NULL,
            object_kind       VARCHAR(24) NOT NULL,
            object_row_id     BIGINT UNSIGNED NOT NULL,
            valid_from_epoch  SMALLINT UNSIGNED NULL,
            retired_epoch     SMALLINT UNSIGNED NOT NULL,
            reason            VARCHAR(16) NOT NULL,
            gone_epoch        SMALLINT UNSIGNED NULL,
            retired_at        TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
            PRIMARY KEY (id),
            UNIQUE KEY uq_designations_uid (designation_uid),
            KEY idx_designations_name (designation_name, retired_epoch),
            KEY idx_designations_object (object_kind, object_row_id, retired_epoch),
            CONSTRAINT chk_designations_kind CHECK (designation_kind IN ('uid', 'name')),
            CONSTRAINT chk_designations_reason CHECK (reason IN ('collision bump', 'rename', 'regeneration', 'parent change', 'admin renumber')),
            CONSTRAINT chk_designations_shape CHECK (
                (designation_kind = 'uid' AND designation_uid IS NOT NULL AND designation_name IS NULL)
                OR (designation_kind = 'name' AND designation_name IS NOT NULL AND designation_uid IS NULL))
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""")
