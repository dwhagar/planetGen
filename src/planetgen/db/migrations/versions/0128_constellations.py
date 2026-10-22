"""Stored constellations (DB.64).

Schema v128. The view from a planet (docs/design/sky-view.md, VIEW.3) groups the brighter stars of a system's sky into constellations, and a planet keeps the same ones each time. `constellations` stores them per system: the system's uid (`system_uid`), the constellation's number within that system (`constellation_index`), the object ID of its name (`name_id`, the `constellation` kind of VIEW.4, NULL until it is named), the member stars (`members`, JSON: the star object IDs) and the edges drawn between them (`edges`, JSON: pairs of members), and the seed and version key it was built under (`seed`, `version_key`) so it can be rebuilt and told stale.

This revision makes the table and nothing else; no code reads it yet. A table already there is left alone,
so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0128"
down_revision = "0127"
branch_labels = None
depends_on = None

TABLES = {
    "constellations": """
CREATE TABLE IF NOT EXISTS constellations (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    system_uid           BINARY(12) NOT NULL,
    constellation_index  SMALLINT UNSIGNED NOT NULL,
    name_id              BINARY(12) NULL,
    members              JSON NOT NULL,
    edges                JSON NOT NULL,
    seed                 BIGINT UNSIGNED NOT NULL,
    version_key          CHAR(22) NOT NULL,
    created_at           TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    UNIQUE KEY uq_constellations_system_index (system_uid, constellation_index),
    KEY idx_constellations_name (name_id)
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
