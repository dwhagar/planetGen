"""The `galaxy_components` table: the discrete stellar components of a galaxy beyond the disk (GEN.159).

Schema v97. One table for every kind of component (the globular clusters now; the satellites and streams of GEN.237
later), each a density profile with a centre, velocity, mass, age and [Fe/H] and a `cell_key` for region queries. The
migration only creates the empty table; the plan step that draws the clusters from the galaxy seed comes with GEN.160,
so nothing writes rows yet. Skipped when the table is already there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0097"
down_revision = "0096"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    exists = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'galaxy_components'")).scalar()
    if exists:
        return
    op.execute("""
CREATE TABLE galaxy_components (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    kind                 VARCHAR(12) NOT NULL,
    name                 VARCHAR(255) NOT NULL,
    catalogue_id         VARCHAR(40),
    source               VARCHAR(12) NOT NULL,
    parent_component_id  BIGINT UNSIGNED,
    center_x_pc          DOUBLE NOT NULL,
    center_y_pc          DOUBLE NOT NULL,
    center_z_pc          DOUBLE NOT NULL,
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(center_x_pc) < 262144 AND ABS(center_y_pc) < 262144 AND ABS(center_z_pc) < 262144, ((((((((((((CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    KEY idx_galaxy_components_cell (cell_key, center_x_pc, center_y_pc, center_z_pc),
    velocity_x_kms       DOUBLE NOT NULL,
    velocity_y_kms       DOUBLE NOT NULL,
    velocity_z_kms       DOUBLE NOT NULL,
    epoch_unix           DOUBLE NOT NULL,
    mass_solar           DOUBLE NOT NULL,
    half_mass_radius_pc  DOUBLE NOT NULL,
    profile              VARCHAR(12) NOT NULL,
    king_w0              DOUBLE,
    tidal_radius_pc      DOUBLE,
    age_gy               DOUBLE NOT NULL,
    metallicity_fe_h     DOUBLE NOT NULL,
    core_collapsed       TINYINT(1) NOT NULL DEFAULT 0,
    seed                 BIGINT UNSIGNED NOT NULL,
    created_at           TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    UNIQUE KEY uq_galaxy_components_name (kind, name),
    KEY idx_galaxy_components_parent (parent_component_id),
    CONSTRAINT fk_galaxy_components_parent
        FOREIGN KEY (parent_component_id) REFERENCES galaxy_components(id) ON DELETE SET NULL,
    CONSTRAINT chk_galaxy_components_kind CHECK (kind IN ('cluster', 'satellite', 'stream')),
    CONSTRAINT chk_galaxy_components_source CHECK (source IN ('catalogue', 'generated')),
    CONSTRAINT chk_galaxy_components_profile CHECK (profile IN ('king', 'exponential', 'plummer')),
    CONSTRAINT chk_galaxy_components_core_collapsed CHECK (core_collapsed IN (0, 1))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""")
