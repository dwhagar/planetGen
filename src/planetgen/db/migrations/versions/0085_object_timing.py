"""The `object_timing_stats` and `object_timing_samples` tables: how long each kind of object takes to generate (PERF.92).

Schema v85. The migration only creates the empty tables; the next run that saves sectors fills them. Each is
skipped when it is already there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0085"
down_revision = "0084"
branch_labels = None
depends_on = None

TABLES = {
    "object_timing_stats": """
        CREATE TABLE object_timing_stats (
            id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
            run_id        BIGINT UNSIGNED NOT NULL,
            sector_id     BIGINT UNSIGNED NULL,
            kind          VARCHAR(24) NOT NULL,
            object_count  INT UNSIGNED NOT NULL,
            total_us      BIGINT UNSIGNED NOT NULL,
            min_us        BIGINT UNSIGNED NOT NULL,
            max_us        BIGINT UNSIGNED NOT NULL,
            histogram     VARCHAR(255) NOT NULL,
            KEY idx_object_timing_stats_run (run_id, kind),
            KEY idx_object_timing_stats_sector (sector_id),
            CONSTRAINT fk_object_timing_stats_run FOREIGN KEY (run_id) REFERENCES generation_runs(id) ON DELETE CASCADE,
            CONSTRAINT fk_object_timing_stats_sector FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""",
    "object_timing_samples": """
        CREATE TABLE object_timing_samples (
            id         BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
            run_id     BIGINT UNSIGNED NOT NULL,
            sector_id  BIGINT UNSIGNED NULL,
            kind       VARCHAR(24) NOT NULL,
            ordinal    INT UNSIGNED NOT NULL,
            own_us     BIGINT UNSIGNED NOT NULL,
            total_us   BIGINT UNSIGNED NOT NULL,
            sampled    TINYINT(1) NOT NULL,
            KEY idx_object_timing_samples_run (run_id, kind),
            KEY idx_object_timing_samples_sector (sector_id),
            CONSTRAINT fk_object_timing_samples_run FOREIGN KEY (run_id) REFERENCES generation_runs(id) ON DELETE CASCADE,
            CONSTRAINT fk_object_timing_samples_sector FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""",
}


def upgrade():
    connection = op.get_bind()
    for name, statement in TABLES.items():
        exists = connection.execute(sa.text(
            "SELECT COUNT(*) FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = :name"),
            {"name": name}).scalar()
        if not exists:
            op.execute(statement)
