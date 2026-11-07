"""The Object ID is the stored name of an object with no star (GEN.1256).

Schema v221. No column changes. The objects named by their position ID (GEN.64) take their Object ID (GEN.171, `galaxy/object_uid.py`) as
the stored `name`, in its printed form (`FE81000A2B-0000005-000`): `rogue_planets`, `black_holes`, `neutron_stars`, `nebulae`,
`supernova_remnants`, `quasars`, `interstellar_comets` and `asteroid_fields`. The API shows these by their local name (`names/naming_key.py`)
and the galaxy's naming key, as before. A row with no `uid` keeps its name. Every step skips work already done, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0221"
down_revision = "0220"
branch_labels = None
depends_on = None

TABLES = ("rogue_planets", "black_holes", "neutron_stars", "nebulae", "supernova_remnants", "quasars", "interstellar_comets",
          "asteroid_fields")
"""tuple: The tables whose rows were named by a 19-digit position ID."""

PRINTED = "CONCAT(SUBSTR(HEX(uid), 1, 14), '-', SUBSTR(HEX(uid), 15, 7), '-', SUBSTR(HEX(uid), 22, 3))"
"""str: The SQL of an Object ID's printed form (`object_uid.format_id`): 14 hex digits of sector, 7 of serial and 3 of body."""


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        helpers.batched_update(connection, table, f"name = {PRINTED}, modified_at = modified_at",
                               "uid IS NOT NULL AND name REGEXP '^[0-9A-Fa-f]{19}$'", label=f"{table} name")
