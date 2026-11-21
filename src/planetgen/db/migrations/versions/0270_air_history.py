"""The record of a body's air history (GEN.1745 slice 4b).

Schema v270. `planets`, `moons` and `rogue_planets` gain `air_history` (MEDIUMTEXT, nullable): the record `physics/air_record.pack` makes of a rocky body's air history, as JSON, the steps of its history (age, surface
temperature, total pressure, the event that starts the step, the oxygen law, the partial pressure of every gas of 0.1 Pa or more) and its life ages. NULL for a giant, a body with an air under 500 Pa, a rogue planet
(no star, no age) and a row made before v270, which gets it when the galaxy is regenerated (a migration does not draw). A column already there is skipped, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0270"
down_revision = "0269"
branch_labels = None
depends_on = None

TABLES = ("planets", "moons", "rogue_planets")
"""tuple: The tables that gain the column."""


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        if not helpers.exists_column(connection, table, "air_history"):
            op.execute(f"ALTER TABLE {table} ADD COLUMN air_history MEDIUMTEXT NULL")


def downgrade():
    connection = op.get_bind()
    for table in reversed(TABLES):
        if helpers.exists_column(connection, table, "air_history"):
            op.execute(f"ALTER TABLE {table} DROP COLUMN air_history")
