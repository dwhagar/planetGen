"""A planet or moon keeps the eccentricity of the orbit it was generated on (GEN.727).

Schema v161. `planets` and `moons` gain `eccentricity` (DOUBLE NOT NULL DEFAULT 0, the eccentricity of the orbit the body was made on;
`distance_km` stays its semi-major axis). A body saved earlier was made on a circle, so 0 is what it was; the regenerated galaxy of the
generator epoch of this change draws real eccentricities, inclinations and orbital phases.

The column goes in with one ALTER for each table (no ALGORITHM=INSTANT: both carry a FULLTEXT index, which refuses it); a column already
there is left alone, so a stopped run starts again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0161"
down_revision = "0160"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        if not helpers.exists_column(connection, table, "eccentricity"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN eccentricity DOUBLE NOT NULL DEFAULT 0"))
