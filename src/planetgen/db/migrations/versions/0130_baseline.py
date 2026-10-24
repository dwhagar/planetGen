"""Baseline: the schema as of v130.

Schema v130. A database below v130 is refused (`SchemaTooOldError`); the
revisions before this one were deleted when Alembic's history was
re-baselined here (DB.98), and this one only marks where it starts, so it
changes nothing.
"""

revision = "0130"
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    pass
