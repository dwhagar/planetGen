"""Every kind of object may show a flavor line (GEN.1433).

Schema v244. `asteroid_belts`, `stars`, `comets`, `black_holes`, `neutron_stars`, `nebulae`, `supernova_remnants`, `quasars`, `interstellar_comets` and `asteroid_fields` gain `flavor_key` (VARCHAR(32), nullable):
the id of the flavor line the object shows, a line of its own scope of `physics/data/flavor_tree.json`. NULL for an object with no line and for a row made before v244, which gets it when the galaxy is regenerated (a
migration does not draw). A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0244"
down_revision = "0243"
branch_labels = None
depends_on = None

TABLES = ("asteroid_belts", "stars", "comets", "black_holes", "neutron_stars", "nebulae", "supernova_remnants", "quasars",
          "interstellar_comets", "asteroid_fields")


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        if not helpers.exists_column(connection, table, "flavor_key"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN flavor_key VARCHAR(32) NULL"))
