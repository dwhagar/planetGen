"""A deleted object keeps its ID (GEN.1141).

Schema v265. `designations.reason` takes a seventh value, `deleted`: when an object is deleted its current ID is written as a retired one (so it is reserved for ever and never issued again, invariant 7 of
`docs/design/architecture/designation-history.md`) and the object's rows in `designations` are stamped with `gone_epoch`, so a link to it answers 410 and not 404 (`api/ids.py`). The CHECK
(`chk_designations_reason`) is replaced by one that lists the seven reasons. No row changes.
"""

import sqlalchemy as sa
from alembic import op

revision = "0265"
down_revision = "0264"
branch_labels = None
depends_on = None

REASONS = ("collision bump", "rename", "regeneration", "parent change", "admin renumber", "charted", "deleted")
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
    if clause is None or "'deleted'" not in clause:       # already as wanted otherwise, so a stopped run can be started again
        _replace(REASONS, clause)


def downgrade():
    op.execute("DELETE FROM designations WHERE reason = 'deleted'")
    _replace(OLD_REASONS, _clause())
