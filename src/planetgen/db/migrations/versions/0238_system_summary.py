"""The stamp of the formula that scored each star system (GEN.1479).

Schema v238. `system_summary` has a narrow row for each star system: `habitability_stamp`, the stamp of the code that scored the system's planets and moons
(`db/derived_stamp.py`). `rescore-habitability` takes the systems whose stamp is NULL or not the running code's and leaves the others alone. Every system already stored
starts with no stamp, so the first run after this scores them all once (it is what the unmigrated run did every time); a table already there is kept and
the rows are not duplicated, so a stopped run can be started again. A side table because `star_systems` has a FULLTEXT key, and adding a column there copies every row.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0238"
down_revision = "0237"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_table(connection, "system_summary"):
        op.execute("""
            CREATE TABLE system_summary (
                star_system_id      BIGINT UNSIGNED NOT NULL PRIMARY KEY,
                habitability_stamp  INT UNSIGNED NULL,
                CONSTRAINT fk_system_summary_system
                    FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""")
