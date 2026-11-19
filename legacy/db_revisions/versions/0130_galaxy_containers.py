"""Globular clusters and superbubbles stored as containers (DB.85).

Schema v130. The Location block of an item page (UX.144) has to say which globular cluster or superbubble a star is inside, and a cluster's or bubble's page has to list what it holds. Until now a globular cluster was a `galaxy_components` density with no way to ask who is inside it, and the superbubbles were not stored at all (they are drawn from the galaxy seed on demand, GEN.254). `galaxy_containers` stores each one as a container: its object ID (`uid`, a region serial of the sector holding its centre, see names in galaxy/object_uid.py), `kind` (`globular_cluster` or `superbubble`), `name`, the component it is (`component_id`, for a cluster) or the field address that drew it (`bubble_col`, `bubble_row`, `bubble_n`, for a bubble), its centre in the galaxy frame (`center_x_pc`, `center_y_pc`, `center_z_pc`), its `extent_pc` (a cluster's tidal radius, a bubble's radius) and the birth sector of its centre (`shell_index`, `band_index`, `wedge_index`). A plan writes the rows from the seed; nothing in them is fixed content.

This revision makes the table and nothing else; no code reads it yet. A table already there is left alone,
so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0130"
down_revision = "0129"
branch_labels = None
depends_on = None

TABLES = {
    "galaxy_containers": """
CREATE TABLE IF NOT EXISTS galaxy_containers (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    uid           BINARY(12) NOT NULL,
    kind          VARCHAR(16) NOT NULL,
    name          VARCHAR(255) NULL,
    component_id  BIGINT UNSIGNED NULL,
    bubble_col      INT NULL,
    bubble_row      INT NULL,
    bubble_n      SMALLINT UNSIGNED NULL,
    center_x_pc   DOUBLE NOT NULL,
    center_y_pc   DOUBLE NOT NULL,
    center_z_pc   DOUBLE NOT NULL,
    extent_pc     DOUBLE NOT NULL,
    shell_index   INT NOT NULL,
    band_index    INT NOT NULL,
    wedge_index   INT NOT NULL,
    created_at    TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    UNIQUE KEY uq_galaxy_containers_uid (uid),
    UNIQUE KEY uq_galaxy_containers_bubble (bubble_col, bubble_row, bubble_n),
    KEY idx_galaxy_containers_kind_x (kind, center_x_pc),
    KEY idx_galaxy_containers_birth (shell_index, band_index, wedge_index),
    KEY idx_galaxy_containers_component (component_id),
    CONSTRAINT fk_galaxy_containers_component FOREIGN KEY (component_id) REFERENCES galaxy_components(id) ON DELETE CASCADE,
    CONSTRAINT chk_galaxy_containers_kind CHECK (kind IN ('globular_cluster', 'superbubble')),
    CONSTRAINT chk_galaxy_containers_extent CHECK (extent_pc >= 1)
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
