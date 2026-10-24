"""A planet or moon keeps the tide heat flow its ocean's life can use (tides in the chemical-energy term of the habitability score).

Schema v152. `planets` and `moons` gain `tidal_flux_w_m2` (DOUBLE, NULL: the eccentricity-tide heat flow in W/m2 of an icy moon with a
frozen surface, NULL for a planet, which `habitability_world.energy_flux_w_m2` adds to the moon's own heat flow). The tide was drawn and thrown away
before, and the forced eccentricity behind it is not stored, so a body saved earlier keeps NULL and its stored energy flux
unchanged: only a moon generated now has the tide in it.

The column goes in with one ALTER for each table (no ALGORITHM=INSTANT: both carry a FULLTEXT index, which refuses it); a column already
there is left alone, so a stopped run starts again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0152"
down_revision = "0151"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        if not helpers.exists_column(connection, table, "tidal_flux_w_m2"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN tidal_flux_w_m2 DOUBLE NULL"))
