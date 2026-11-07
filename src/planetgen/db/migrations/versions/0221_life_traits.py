"""The adaptation traits of a body's life (GEN.876).

Schema v221. `planets`, `moons` and `rogue_planets` gain `life_traits` (VARCHAR(48), NULL): the traits of the life that lives on the body, as text such as
`radiation:2,sulfide:3` (the stress and the grade, 2 tolerant or 3 extreme; `physics/adaptation_traits.py`), NULL for none. It is NULL here: a trait is drawn from
the stress a finished body has, in the life pass, so a row made before v221 gets its traits when the galaxy is regenerated. A column already there is skipped, so a
stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0221"
down_revision = "0220"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons", "rogue_planets"):
        if not helpers.exists_column(connection, table, "life_traits"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN life_traits VARCHAR(48) NULL AFTER life_complexity"))
