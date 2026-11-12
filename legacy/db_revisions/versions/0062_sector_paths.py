"""Sector paths: the spline knots of a body's path through its sector (GEN.123).

Schema v62. Adds `sector_paths` and `sector_path_knots`, as in `schema.sql`.
"""

from alembic import op

revision = "0062"
down_revision = "0061"
branch_labels = None
depends_on = None


def upgrade():
    op.execute(
        """
        CREATE TABLE IF NOT EXISTS sector_paths (
            id              BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
            sector_id       BIGINT UNSIGNED NOT NULL,
            object_table    VARCHAR(24) NOT NULL CHECK (object_table IN ('star_systems', 'rogue_planets', 'interstellar_comets')),
            object_id       BIGINT UNSIGNED NOT NULL,
            exited          TINYINT(1) NOT NULL,
            duration_years  DOUBLE NOT NULL,
            computed_at     TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

            UNIQUE (object_table, object_id),
            KEY idx_sector_paths_sector (sector_id),
            CONSTRAINT fk_sector_paths_sector
                FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
        """
    )
    op.execute(
        """
        CREATE TABLE IF NOT EXISTS sector_path_knots (
            id        BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
            path_id   BIGINT UNSIGNED NOT NULL,
            position  INT NOT NULL,
            t_years   DOUBLE NOT NULL,
            x_pc      DOUBLE NOT NULL,
            y_pc      DOUBLE NOT NULL,
            z_pc      DOUBLE NOT NULL,
            vx_kms    DOUBLE NOT NULL,
            vy_kms    DOUBLE NOT NULL,
            vz_kms    DOUBLE NOT NULL,

            UNIQUE (path_id, position),
            CONSTRAINT fk_sector_path_knots_path
                FOREIGN KEY (path_id) REFERENCES sector_paths(id) ON DELETE CASCADE
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
        """
    )
