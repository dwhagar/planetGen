"""A name kept when charting names a star (GEN.1257).

Schema v221. `designations.reason` takes a sixth value, `charted`: an uncharted star has a name its Object ID spells (a local name, GEN.1255) that is computed and never
stored, and when charting gives the star its own word-salad name the spelled name is kept as the first retired name of its system, so the "Formerly" line shows it and a
search for it finds the system. The CHECK (`chk_designations_reason`) is replaced by one that lists the six reasons. No row changes.
"""

import sqlalchemy as sa
from alembic import op

revision = "0221"
down_revision = "0220"
branch_labels = None
depends_on = None

REASONS = ("collision bump", "rename", "regeneration", "parent change", "admin renumber", "charted")
OLD_REASONS = REASONS[:-1]
NAME = "chk_designations_reason"


def _clause():
    return op.get_bind().execute(sa.text(
        "SELECT CHECK_CLAUSE FROM information_schema.CHECK_CONSTRAINTS WHERE CONSTRAINT_SCHEMA = DATABASE()"
        " AND CONSTRAINT_NAME = :name"), {"name": NAME}).scalar()


def _replace(reasons, clause):
    listed = ", ".join(f"'{name}'" for name in reasons)
    drop = f"DROP CHECK {NAME}, " if clause is not None else ""
    op.execute(f"ALTER TABLE designations {drop}ADD CONSTRAINT {NAME} CHECK (reason IN ({listed}))")


def upgrade():
    clause = _clause()
    if clause is None or "'charted'" not in clause:       # already as wanted otherwise, so a stopped run can be started again
        _replace(REASONS, clause)


def downgrade():
    op.execute("DELETE FROM designations WHERE reason = 'charted'")
    clause = _clause()
    if clause is None or "'charted'" in clause:
        _replace(OLD_REASONS, clause)
