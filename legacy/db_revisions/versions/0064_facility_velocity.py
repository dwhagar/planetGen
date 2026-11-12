"""Stand-alone facilities store a velocity (GEN.125).

Schema v64. `facilities` gains `velocity_x_kms`/`_y_kms`/`_z_kms`, the
galactic-axes velocity a stand-alone facility has on the rotation curve at
its place, worked out here for the ones already stored (see `schema.sql`'s
"v64" note). Facilities on a body, in orbit or in a belt keep 0.
"""

import math

import sqlalchemy as sa
from alembic import op

from planetgen.galaxy.galactic_orbit import calculate_galactic_orbit
from planetgen.galaxy.system_position import galactic_velocity_ms
from planetgen.physics import constants

revision = "0064"
down_revision = "0063"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    has_column = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'facilities' AND COLUMN_NAME = 'velocity_x_kms'")).scalar()
    if not has_column:
        op.execute("ALTER TABLE facilities ADD COLUMN velocity_x_kms DOUBLE NOT NULL DEFAULT 0 AFTER galactic_radius_pc,"
                   " ADD COLUMN velocity_y_kms DOUBLE NOT NULL DEFAULT 0 AFTER velocity_x_kms,"
                   " ADD COLUMN velocity_z_kms DOUBLE NOT NULL DEFAULT 0 AFTER velocity_y_kms")
    rows = connection.execute(sa.text(
        "SELECT id, center_x_pc, center_y_pc, center_z_pc FROM facilities"
        " WHERE host_type = 'space' AND center_x_pc IS NOT NULL")).fetchall()
    for facility_id, x, y, z in rows:
        speed_kms, _period_gy = calculate_galactic_orbit(math.hypot(x, y) * constants.PARSEC_M / constants.LIGHTYEAR_M)
        velocity = galactic_velocity_ms((x, y, z), speed_kms)
        connection.execute(
            sa.text("UPDATE facilities SET velocity_x_kms = :vx, velocity_y_kms = :vy, velocity_z_kms = :vz,"
                    " modified_at = modified_at WHERE id = :id"),
            {"vx": velocity[0] / 1000.0, "vy": velocity[1] / 1000.0, "vz": velocity[2] / 1000.0, "id": facility_id})
