"""The belt an irregular moon was caught from (GEN.1085, part 2).

Schema v224. `irregular_moons` gains `belt_order` (INT) and `belt_rank` (INT), both nullable: for a body caught from an inner belt, the belt's number among the belts of the planet's star (nearest first), and, when the body
was one of the belt's stored members, its rank there (it keeps the designation `<belt ID>.<rank hex>`; the member's row stays in `belt_members` with the status `moon`). A row made before v224 has them NULL, which is a body of
the outer disk or of the planet's own zone; nothing is backfilled. Each column is skipped when it is already there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0224"
down_revision = "0223"
branch_labels = None
depends_on = None

COLUMNS = (("belt_order", "INT NULL"), ("belt_rank", "INT NULL"))


def upgrade():
    connection = op.get_bind()
    for name, definition in COLUMNS:
        if not helpers.exists_column(connection, "irregular_moons", name):
            connection.execute(sa.text(f"ALTER TABLE irregular_moons ADD COLUMN {name} {definition}, ALGORITHM=INSTANT"))
