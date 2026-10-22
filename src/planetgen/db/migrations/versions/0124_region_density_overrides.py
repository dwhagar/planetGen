"""Region density overrides: an admin's multiplier for the sectors of a region (DB.38).

Schema v124. The density edits of MAP.119 and ADM.57 (docs/design/edit-everything.md, docs/design/map-ui-and-frontend-libraries.md section 8) never edit the galaxy skeleton: the admin sets a multiplier for a region, applied to the density used when sectors there are filled. `region_density_overrides` is the append-only list of those settings: the region as a shape (`sphere`, `box`, `shell`, `band`, `wedge` or `sectors`) and its parameters as JSON (the region object of GEN.302), the `multiplier` (above zero; the 0.1 to 10 range is the interface's, so a typed value beyond it can carry a warning), who set it (`set_by_id`, the account's id in the control database, and `set_by_name` as it was), when, an optional `reason`, and `undone_at`, set once when the setting is undone (the only change a row ever gets).

This revision makes the table and nothing else; no code reads it yet. A table already there is left alone,
so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0124"
down_revision = "0123"
branch_labels = None
depends_on = None

TABLES = {
    "region_density_overrides": """
CREATE TABLE IF NOT EXISTS region_density_overrides (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    shape         VARCHAR(16) NOT NULL,
    params        JSON NOT NULL,
    multiplier    DECIMAL(10,4) NOT NULL,
    set_by_id     BIGINT UNSIGNED NULL,
    set_by_name   VARCHAR(255) NOT NULL,
    set_at        TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    reason        VARCHAR(500) NULL,
    undone_at     TIMESTAMP(3) NULL,
    KEY idx_region_density_overrides_active (undone_at, id),
    CONSTRAINT chk_region_density_overrides_shape CHECK (shape IN ('sphere', 'box', 'shell', 'band', 'wedge', 'sectors')),
    CONSTRAINT chk_region_density_overrides_multiplier CHECK (multiplier >= 0.0001)
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
