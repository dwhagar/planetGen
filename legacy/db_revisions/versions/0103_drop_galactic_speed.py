"""The scalar galactic speed is dropped; the movement vector is the one record (GEN.336).

Schema v103. Revision 0101 (DB.55) gave every free object a movement vector (`velocity_x/y/z_kms`, galactic axes) and
backfilled it from the scalar `galactic_orbital_speed_kms` as the tangent of the galaxy's rotation at the same speed.
The code now writes and reads the vector (a speed is its length: `store._speed_kms`), so the scalar goes from `stars`,
`black_holes`, `neutron_stars`, `nebulae`, `supernova_remnants`, `rogue_planets`, `interstellar_comets` and
`asteroid_fields`, and `binary_galactic_orbital_speed_kms` from `star_systems`. A star and a close pair take their
system's vector, so nothing is lost.

A row that has no vector yet (an object never put in the galaxy) keeps none. Each column is dropped only when it is
there, so a stopped run can be started again; dropping a column is an in-place change in MySQL 8 and copies no rows.
"""

import sqlalchemy as sa
from alembic import op

revision = "0103"
down_revision = "0102"
branch_labels = None
depends_on = None

DROPS = {
    "stars": "galactic_orbital_speed_kms",
    "black_holes": "galactic_orbital_speed_kms",
    "neutron_stars": "galactic_orbital_speed_kms",
    "nebulae": "galactic_orbital_speed_kms",
    "supernova_remnants": "galactic_orbital_speed_kms",
    "rogue_planets": "galactic_orbital_speed_kms",
    "interstellar_comets": "galactic_orbital_speed_kms",
    "asteroid_fields": "galactic_orbital_speed_kms",
    "star_systems": "binary_galactic_orbital_speed_kms",
}
"""dict: Table -> the scalar speed column dropped from it."""


def _has_column(connection, table, column):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar())


def upgrade():
    connection = op.get_bind()
    for table, column in DROPS.items():
        if _has_column(connection, table, column):
            op.execute(f"ALTER TABLE {table} DROP COLUMN {column}")
