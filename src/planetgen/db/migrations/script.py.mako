"""${message}

Schema v${up_revision} (revision ${up_revision}, after ${down_revision | comma,n}).
"""

from alembic import op
import sqlalchemy as sa
${imports if imports else ""}

revision = ${repr(up_revision)}
down_revision = ${repr(down_revision)}
branch_labels = None
depends_on = None


def upgrade():
    ${upgrades if upgrades else "pass"}
