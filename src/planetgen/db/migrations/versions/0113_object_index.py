"""The object index table: one row per object of every kind the site lists (DB.79).

Schema v113. The lists are queries over the live tables (the phenomena list is a `UNION ALL` over about eight tables with a
`COUNT(*)` and an `ORDER BY` over the union; stars, planets and moons have no list at all and the phenomenon tables have no
name index). `object_index` is where PERF.140 keeps one row for each object of each kind the site lists, so a list is
one keyset read of one table and a name search is one FULLTEXT lookup. This revision makes the table and nothing else:
it is empty until PERF.140 fills it, and no row of any other table is touched.

A row is the object's kind (`system`, `star`, `planet`, `moon`, `belt`, `comet`, `facility` or a phenomenon kind) and its
96-bit object ID (`uid`, unique together: the primary key), the row id in its own table (`source_id`), its name, its sector
and system, the object ID of what it orbits or belongs to (`parent_uid`), its class or type (`object_class`), mass and
radius, whether it is placed on the galaxy, and when the row was indexed. The keyset keys put the kind first, then the
sort column, then the object ID, so a page after a row is one range read; `ft_object_index_name` is the FULLTEXT name
index. There is no foreign key: the builder keeps the table in step, and a key on a table this size costs every write.

A table already there is left alone, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0113"
down_revision = "0112"
branch_labels = None
depends_on = None

CREATE_OBJECT_INDEX = """
CREATE TABLE IF NOT EXISTS object_index (
    kind            VARCHAR(24) NOT NULL,
    uid             BINARY(12) NOT NULL,
    source_id       BIGINT UNSIGNED NOT NULL,
    name            VARCHAR(255) NOT NULL,
    sector_id       BIGINT UNSIGNED NULL,
    system_id       BIGINT UNSIGNED NULL,
    parent_uid      BINARY(12) NULL,
    object_class    VARCHAR(64) NULL,
    mass_kg         DOUBLE NULL,
    radius_km       DOUBLE NULL,
    placed          TINYINT(1) NOT NULL DEFAULT 1,
    indexed_at      TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    PRIMARY KEY (kind, uid),
    KEY idx_object_index_source (kind, source_id),
    KEY idx_object_index_name (kind, name, uid),
    KEY idx_object_index_class (kind, object_class, name, uid),
    KEY idx_object_index_mass (kind, mass_kg, uid),
    KEY idx_object_index_radius (kind, radius_km, uid),
    KEY idx_object_index_sector (sector_id, kind, name, uid),
    KEY idx_object_index_system (system_id, kind, uid),
    FULLTEXT KEY ft_object_index_name (name),
    CONSTRAINT chk_object_index_placed CHECK (placed IN (0, 1))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
"""


def upgrade():
    connection = op.get_bind()
    if connection.execute(sa.text(
            "SELECT COUNT(*) FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()"
            " AND TABLE_NAME = 'object_index'")).scalar():
        return
    op.execute(CREATE_OBJECT_INDEX)
