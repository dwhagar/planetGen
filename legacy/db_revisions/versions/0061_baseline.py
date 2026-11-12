"""Baseline: the schema as of v61.

Schema v61. Databases up to v61 are created and upgraded by `schema.sql`
and the `_migrate_vN_to_vM` steps in `planetgen/db/store.py`; this
revision only marks where Alembic takes over, so it changes nothing.
"""

revision = "0061"
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    pass
