"""The nebula aggregates of the Galaxy Map's block ladder (MAP.172).

Schema v116. Three new tables: `nebula_block_stats` (what the nebulae add up to in a block of the ladder: count,
volume, centroid, bounding radius, the volume of each family, the dominant class and the cover), `nebula_block_classes`
(the volume of each class in a block, so one stored nebula can be added to a block) and `nebula_block_state` (the
highest nebula id folded in). They are filled from `nebula_field` and the stored nebulae the galaxy already holds
(`db.nebula_blocks.rebuild`), so the cover of a block no sector was generated in is known at once. A galaxy planned
before `nebula_field` (v83) has nothing to fill from; `planetgen fill-nebula-blocks` makes the tables once its
field is drawn. The fill is skipped once `nebula_block_state` holds its row (the mark that it was done), so a stopped run can
be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen import tuning
from planetgen.db.migrations import helpers

revision = "0116"
down_revision = "0115"
branch_labels = None
depends_on = None

TABLES = {
    "nebula_block_stats": """
        CREATE TABLE IF NOT EXISTS nebula_block_stats (
            m                         SMALLINT UNSIGNED NOT NULL,
            shell_index               INT NOT NULL,
            band_index                INT NOT NULL,
            wedge_index               INT NOT NULL,
            nebula_count              INT UNSIGNED NOT NULL,
            volume_pc3                DOUBLE NOT NULL,
            centroid_x_pc             DOUBLE NULL,
            centroid_y_pc             DOUBLE NULL,
            centroid_z_pc             DOUBLE NULL,
            bound_radius_pc           DOUBLE NULL,
            volume_diffuse_pc3        DOUBLE NOT NULL,
            volume_emission_pc3       DOUBLE NOT NULL,
            volume_reflection_pc3     DOUBLE NOT NULL,
            volume_planetary_pc3      DOUBLE NOT NULL,
            volume_dark_pc3           DOUBLE NOT NULL,
            volume_remnant_pc3        DOUBLE NOT NULL,
            dominant_class            CHAR(1) NULL,
            cover_volume_pc3          DOUBLE NOT NULL,
            cover                     DOUBLE NOT NULL,
            PRIMARY KEY (m, shell_index, band_index, wedge_index)
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""",
    "nebula_block_classes": """
        CREATE TABLE IF NOT EXISTS nebula_block_classes (
            m                         SMALLINT UNSIGNED NOT NULL,
            shell_index               INT NOT NULL,
            band_index                INT NOT NULL,
            wedge_index               INT NOT NULL,
            nebula_class              CHAR(1) NOT NULL,
            volume_pc3                DOUBLE NOT NULL,
            PRIMARY KEY (m, shell_index, band_index, wedge_index, nebula_class)
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""",
    "nebula_block_state": """
        CREATE TABLE IF NOT EXISTS nebula_block_state (
            id                        TINYINT UNSIGNED NOT NULL,
            nebulae_through_id        BIGINT UNSIGNED NOT NULL,
            PRIMARY KEY (id)
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""",
}


def upgrade():
    connection = op.get_bind()
    for statement in TABLES.values():    # `schema.sql` is applied first, so the tables are usually there already
        connection.execute(sa.text(statement))
    if connection.execute(sa.text("SELECT COUNT(*) FROM nebula_block_state")).scalar():
        return    # the state row is the mark of a fill that was done
    from planetgen.db import nebula_blocks

    edge = connection.execute(sa.text("SELECT edge_pc FROM galaxy_shape WHERE id = 1")).scalar()
    nebula_blocks.rebuild(helpers.StoreConnection(connection), float(edge) if edge else float(tuning.DEFAULT_SECTOR_EDGE_PC))
