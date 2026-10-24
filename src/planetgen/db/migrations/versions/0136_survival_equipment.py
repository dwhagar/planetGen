"""Survival equipment: required and recommended items replace the one equipment tier (GEN.655, DB.93).

Schema v136. `planets` and `moons` lose `equipment_tier` (0 ideal to 4 full life support) and gain the levels of the
survival equipment line of `physics/equipment.py` (`docs/design/survival-equipment-fields.md`), every one NULL for a gas
giant and for a body saved without a score: `eq_body` and `eq_body_kind` (Body Protection 0 to 6 and what set it),
`eq_body_rec` and `eq_body_rec_kind` (the recommended level, above the required one only), `eq_rad_suit` and
`eq_rad_suit_rec` (the radiation suit step: shielded garment, storm vest, hardened suit), `eq_mobility`, `eq_dust`,
`eq_gas_filter`, `eq_co2_scrub`, `eq_supply` and `eq_supply_rec`, `eq_shelter` and `eq_shelter_rec` (shielded shelter,
hardened habitat), `shelter_column_g_cm2` and `outside_hours` (DOUBLE).

The levels are a pure function of the stored conditions, so no generator epoch: `planetgen rescore-habitability`
(`store.refresh_habitability`) fills them for a galaxy made before this revision, as it fills every score the rules
change, and a body not re-scored shows no equipment line until then. The old tier is not kept anywhere (no
backward-compatibility layer). Each table is one `ALGORITHM=COPY` statement, not `INSTANT`: `planets` and `moons` have a FULLTEXT index, and
an instant change on such a table is refused (error 1846) when a column is dropped. A table already changed is skipped, so
a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0136"
down_revision = "0135"
branch_labels = None
depends_on = None

TINYINT_COLUMNS = (
    "eq_body", "eq_body_kind", "eq_body_rec", "eq_body_rec_kind", "eq_rad_suit", "eq_rad_suit_rec", "eq_mobility",
    "eq_dust", "eq_gas_filter", "eq_co2_scrub", "eq_supply", "eq_supply_rec", "eq_shelter", "eq_shelter_rec",
)
"""tuple: The TINYINT columns added to `planets` and `moons`."""

DOUBLE_COLUMNS = ("shelter_column_g_cm2", "outside_hours")
"""tuple: The DOUBLE columns added after them."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        changes = [f"ADD COLUMN {column} TINYINT NULL" for column in TINYINT_COLUMNS
                   if not helpers.exists_column(connection, table, column)]
        changes += [f"ADD COLUMN {column} DOUBLE NULL" for column in DOUBLE_COLUMNS
                    if not helpers.exists_column(connection, table, column)]
        if helpers.exists_column(connection, table, "equipment_tier"):
            changes.append("DROP COLUMN equipment_tier")
        if changes:
            connection.execute(sa.text(f"ALTER TABLE {table} {', '.join(changes)}, ALGORITHM=COPY"))
