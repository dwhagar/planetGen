"""Ships and stations, user-owned or generated (DB.63).

Schema v129. The ships and stations of 9.3 (USR.10) need a place to live. `vessels` is one row per ship or station: its object ID (`uid`; the object-kind code for it is added to names/object_id.py with the code that writes the table), its `kind` (`ship` or `station`), `name`, `origin` (`generated` or `user`), `owner_id` and `owner_name` (the account's id in the control database and its name as it was; NULL for a generated one), the position it has now in the galaxy frame (`center_x_pc`, `center_y_pc`, `center_z_pc`), its point mass (`mass_kg`) and its movement vector (`velocity_x_kms`, `velocity_y_kms`, `velocity_z_kms`) at `epoch_unix`, the same fields as the SpatialPosition3D of GEN.294, the `course` (JSON, the route plan; NULL while it holds still) and what it was made under (`seed` and `version_key`, NULL for a user-made one).

This revision makes the table and nothing else; no code reads it yet. A table already there is left alone,
so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0129"
down_revision = "0128"
branch_labels = None
depends_on = None

TABLES = {
    "vessels": """
CREATE TABLE IF NOT EXISTS vessels (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    uid                  BINARY(12) NOT NULL,
    kind                 VARCHAR(8) NOT NULL,
    name                 VARCHAR(255) NOT NULL,
    origin               VARCHAR(12) NOT NULL,
    owner_id             BIGINT UNSIGNED NULL,
    owner_name           VARCHAR(255) NULL,
    center_x_pc          DOUBLE NULL,
    center_y_pc          DOUBLE NULL,
    center_z_pc          DOUBLE NULL,
    mass_kg              DOUBLE NULL,
    velocity_x_kms       DOUBLE NOT NULL DEFAULT 0,
    velocity_y_kms       DOUBLE NOT NULL DEFAULT 0,
    velocity_z_kms       DOUBLE NOT NULL DEFAULT 0,
    epoch_unix           DOUBLE NULL,
    course               JSON NULL,
    seed                 BIGINT UNSIGNED NULL,
    version_key          CHAR(22) NULL,
    created_at           TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    UNIQUE KEY uq_vessels_uid (uid),
    KEY idx_vessels_owner (owner_id, id),
    CONSTRAINT chk_vessels_kind CHECK (kind IN ('ship', 'station')),
    CONSTRAINT chk_vessels_origin CHECK (origin IN ('generated', 'user'))
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
