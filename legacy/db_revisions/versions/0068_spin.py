"""Every rotating object stores its spin (GEN.104).

Schema v68. Stars, planets, moons, comets, rogue planets, interstellar
comets and standalone black holes and neutron stars gain `spin_axis_x/y/z`
and `axial_tilt_deg`; stars, comets, rogue planets and interstellar comets
also gain `rotation_period_hours`. All NULL here: rows generated before
v68 have no spin (see `schema.sql`'s "v68" note).
"""

import sqlalchemy as sa
from alembic import op

revision = "0068"
down_revision = "0067"
branch_labels = None
depends_on = None

_SPIN = (
    # (table, column the spin goes after, whether it gains rotation_period_hours)
    ("stars", "reflex_offset_z_km", True),
    ("planets", "rotation_period_hours", False),
    ("moons", "rotation_period_hours", False),
    ("comets", "next_update_due", True),
    ("black_holes", "next_update_due", False),
    ("neutron_stars", "next_update_due", False),
    ("rogue_planets", "next_update_due", True),
    ("interstellar_comets", "next_update_due", True),
)

_SPIN_COLUMNS = ("spin_axis_x", "spin_axis_y", "spin_axis_z", "axial_tilt_deg")


def _has_column(connection, table, column):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar()


def upgrade():
    connection = op.get_bind()
    for table, after, with_period in _SPIN:
        if _has_column(connection, table, "axial_tilt_deg"):
            continue
        columns = (("rotation_period_hours",) if with_period else ()) + _SPIN_COLUMNS
        clauses = []
        for column in columns:
            clauses.append(f"ADD COLUMN {column} DOUBLE AFTER {after}")
            after = column
        op.execute(f"ALTER TABLE {table} " + ", ".join(clauses))
