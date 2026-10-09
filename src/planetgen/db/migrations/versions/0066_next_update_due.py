"""Every moving object gets its own update clock (GEN.106).

Schema v66. The objects the orbit update moves gain `epoch_unix` (when the
stored position holds) and an indexed `next_update_due` (when it will have
moved far enough to store again), both NULL here: the next update run
fills them in (see `schema.sql`'s "v65" note).
"""

import sqlalchemy as sa
from alembic import op

revision = "0066"
down_revision = "0065"
branch_labels = None
depends_on = None

_CLOCKS = (
    # (table, column the clock goes after, column prefix)
    ("star_systems", "velocity_z_kms", ""),
    ("star_systems", "binary_mutual_orbital_phase_deg", "binary_"),
    ("planets", "velocity_z_kms", ""),
    ("moons", "velocity_z_kms", ""),
    ("comets", "velocity_z_kms", ""),
    ("black_holes", "galactic_min_update_interval_years", ""),
    ("neutron_stars", "galactic_min_update_interval_years", ""),
    ("nebulae", "galactic_min_update_interval_years", ""),
    ("supernova_remnants", "galactic_min_update_interval_years", ""),
    ("rogue_planets", "galactic_min_update_interval_years", ""),
    ("interstellar_comets", "galactic_min_update_interval_years", ""),
    ("asteroid_fields", "galactic_min_update_interval_years", ""),
    ("facilities", "orbit_phase_deg", ""),
)


def _has_column(connection, table, column):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar()


def upgrade():
    connection = op.get_bind()
    for table, after, prefix in _CLOCKS:
        if _has_column(connection, table, f"{prefix}next_update_due"):
            continue
        op.execute(f"ALTER TABLE {table} ADD COLUMN {prefix}epoch_unix DOUBLE AFTER {after},"
                   f" ADD COLUMN {prefix}next_update_due DOUBLE AFTER {prefix}epoch_unix,"
                   f" ADD KEY idx_{table}_{prefix}next_update_due ({prefix}next_update_due)")
