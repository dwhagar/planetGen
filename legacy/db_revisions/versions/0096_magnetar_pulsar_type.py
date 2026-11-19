"""The magnetar neutron star (GEN.153).

Schema v96. A neutron star younger than a magnetar's active life may be born a magnetar, and `neutron_stars.pulsar_type`
(`young`, `millisecond` or `non-pulsing`) takes a fourth value, `magnetar`. Its CHECK (`neutron_stars_chk_1`, the name
MySQL gave it) is replaced by one of the same name that lists all four; the quadrant CHECK beside it is not touched, and
no row is. The old values stay valid.
"""

import sqlalchemy as sa
from alembic import op

revision = "0096"
down_revision = "0095"
branch_labels = None
depends_on = None

PULSAR_TYPES = ("young", "millisecond", "magnetar", "non-pulsing")
NAME = "neutron_stars_chk_1"


def upgrade():
    connection = op.get_bind()
    clause = connection.execute(sa.text(
        "SELECT CHECK_CLAUSE FROM information_schema.CHECK_CONSTRAINTS WHERE CONSTRAINT_SCHEMA = DATABASE()"
        " AND CONSTRAINT_NAME = :name"), {"name": NAME}).scalar()
    if clause is not None and "magnetar" in clause:
        return
    listed = ", ".join(f"'{name}'" for name in PULSAR_TYPES)
    drop = f"DROP CHECK {NAME}, " if clause is not None else ""
    op.execute(f"ALTER TABLE neutron_stars {drop}ADD CONSTRAINT {NAME} CHECK (pulsar_type IN ({listed}))")
