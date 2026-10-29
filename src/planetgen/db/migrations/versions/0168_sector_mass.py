"""The mass of each generated sector, block and region (DB.108).

Schema v168. `sector_stats` gains `total_mass_kg` (the summed mass of the sector's stars, planets, moons, asteroid belts, comets,
black holes, neutron stars, quasars, rogue planets, interstellar comets and asteroid fields, NULL when a body of it has no stored
mass), and `galaxy_block_stats` and `region_stats` gain `mass_sum_kg`, the sum of their generated sectors' masses. Kilograms.

The sector masses are made from the stored bodies and the block and region rows again (`db.block_stats.rebuild`), so a database that
already holds sectors has every figure at once; the sectors that cannot be summed are named in the log and left NULL. Each column is
skipped when it is already there and the rebuild empties the tables before it fills them, so a stopped run can be started again.
"""

import logging

from alembic import op

from planetgen.db.migrations import helpers

revision = "0168"
down_revision = "0167"
branch_labels = None
depends_on = None

_COLUMNS = (
    ("sector_stats", "total_mass_kg", "DOUBLE NULL"),
    ("galaxy_block_stats", "mass_sum_kg", "DOUBLE NOT NULL DEFAULT 0"),
    ("region_stats", "mass_sum_kg", "DOUBLE NOT NULL DEFAULT 0"),
)
"""tuple: The `(table, column, definition)` added."""


def upgrade():
    connection = op.get_bind()
    for table, column, definition in _COLUMNS:
        helpers.add_column_instant(connection, table, column, definition)
    from planetgen.db import block_stats

    block_stats.rebuild(helpers.StoreConnection(connection), on_skipped=_warn)


def _warn(addresses):
    """Names the sectors that could not be summed."""
    logging.getLogger("alembic.runtime.migration").warning(
        "%d sectors hold a body with no stored mass and are left out of the mass sums: %s", len(addresses),
        ", ".join("-".join(str(part) for part in address) for address in addresses))
