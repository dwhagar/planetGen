"""The shape row carries the galaxy family and the numbers of its components (DB.66).

Schema v117. `galaxy_shape` gains `family` (`spiral`, `lenticular`, `elliptical`, `irregular` or `dwarf_spheroid`) and the
columns of the components a family may have that a spiral does not draw today: the bar's axes and angle, the thick disk's
amplitude, the Sersic spheroid (amplitude, index, effective radius, two axis ratios, three orientation angles) and the
irregular family's young clump share (docs/design/galaxy-types.md section 10). Every existing galaxy becomes a `spiral`
with no spheroid and no clump share, which changes no stored density. The bar and the thick disk of an existing galaxy
stay NULL, meaning "as the structure draw gives them from the seed", because they were never stored; a later family stores
its own, and an amplitude of 0 removes a component exactly. Nothing reads the columns yet (GEN.379 does). One statement
adds them all, so a stopped run leaves either none or all, and a step already done is left alone.
"""

import sqlalchemy as sa
from alembic import op

revision = "0117"
down_revision = "0116"
branch_labels = None
depends_on = None


def _has_column(connection, table, column):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar())


def upgrade():
    if _has_column(op.get_bind(), "galaxy_shape", "family"):
        return
    op.execute(
        "ALTER TABLE galaxy_shape"
        " ADD COLUMN family VARCHAR(16) NOT NULL DEFAULT 'spiral' AFTER pinned_stellar_mass_solar,"
        " ADD COLUMN bar_present TINYINT NULL AFTER family,"
        " ADD COLUMN bar_axis_ratio_b DOUBLE NULL AFTER bar_present,"
        " ADD COLUMN bar_axis_ratio_c DOUBLE NULL AFTER bar_axis_ratio_b,"
        " ADD COLUMN bar_angle_rad DOUBLE NULL AFTER bar_axis_ratio_c,"
        " ADD COLUMN thick_disk_amplitude DOUBLE NULL AFTER bar_angle_rad,"
        " ADD COLUMN spheroid_amplitude DOUBLE NOT NULL DEFAULT 0 AFTER thick_disk_amplitude,"
        " ADD COLUMN spheroid_sersic_n DOUBLE NOT NULL DEFAULT 0 AFTER spheroid_amplitude,"
        " ADD COLUMN spheroid_effective_radius_pc DOUBLE NOT NULL DEFAULT 0 AFTER spheroid_sersic_n,"
        " ADD COLUMN spheroid_axis_ratio_b DOUBLE NOT NULL DEFAULT 1 AFTER spheroid_effective_radius_pc,"
        " ADD COLUMN spheroid_axis_ratio_c DOUBLE NOT NULL DEFAULT 1 AFTER spheroid_axis_ratio_b,"
        " ADD COLUMN spheroid_orientation_1_rad DOUBLE NOT NULL DEFAULT 0 AFTER spheroid_axis_ratio_c,"
        " ADD COLUMN spheroid_orientation_2_rad DOUBLE NOT NULL DEFAULT 0 AFTER spheroid_orientation_1_rad,"
        " ADD COLUMN spheroid_orientation_3_rad DOUBLE NOT NULL DEFAULT 0 AFTER spheroid_orientation_2_rad,"
        " ADD COLUMN body_clump_share DOUBLE NOT NULL DEFAULT 0 AFTER spheroid_orientation_3_rad,"
        " ADD CONSTRAINT chk_galaxy_shape_family CHECK (family IN"
        " ('spiral', 'lenticular', 'elliptical', 'irregular', 'dwarf_spheroid')),"
        " ADD CONSTRAINT chk_galaxy_shape_bar_axis_ratio_b CHECK (bar_axis_ratio_b > 0 AND bar_axis_ratio_b <= 1),"
        " ADD CONSTRAINT chk_galaxy_shape_bar_axis_ratio_c CHECK (bar_axis_ratio_c > 0 AND bar_axis_ratio_c <= 1),"
        " ADD CONSTRAINT chk_galaxy_shape_thick_disk_amplitude CHECK (thick_disk_amplitude >= 0),"
        " ADD CONSTRAINT chk_galaxy_shape_spheroid_amplitude CHECK (spheroid_amplitude >= 0),"
        " ADD CONSTRAINT chk_galaxy_shape_spheroid_sersic_n CHECK (spheroid_sersic_n >= 0),"
        " ADD CONSTRAINT chk_galaxy_shape_spheroid_effective_radius CHECK (spheroid_effective_radius_pc >= 0),"
        " ADD CONSTRAINT chk_galaxy_shape_spheroid_axis_ratio_b CHECK (spheroid_axis_ratio_b > 0 AND spheroid_axis_ratio_b <= 1),"
        " ADD CONSTRAINT chk_galaxy_shape_spheroid_axis_ratio_c CHECK (spheroid_axis_ratio_c > 0 AND spheroid_axis_ratio_c <= 1),"
        " ADD CONSTRAINT chk_galaxy_shape_body_clump_share CHECK (body_clump_share BETWEEN 0 AND 1)")
