"""The cloud decks of a giant planet (GEN.1577).

Schema v269. `planets` gains `cloud_decks` (TEXT, nullable): the hot-condensate clouds of a giant (the base pressure in bar and temperature in K of each of chromium, manganese sulfide, sodium
sulfide, zinc sulfide and potassium chloride that forms along its profile, deepest first) and the ceiling of each trace gas of its star's gas in molecules per hydrogen atom, as JSON text
(`physics/hot_cloud_deck.py`). NULL for a body that is not a giant of a star and for a row made before v269; nothing is backfilled, a regenerated galaxy fills it. A column already there is skipped,
so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0269"
down_revision = "0268"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_column(connection, "planets", "cloud_decks"):
        connection.execute(sa.text("ALTER TABLE planets ADD COLUMN cloud_decks TEXT NULL AFTER sky_suns"))


def downgrade():
    connection = op.get_bind()
    if helpers.exists_column(connection, "planets", "cloud_decks"):
        connection.execute(sa.text("ALTER TABLE planets DROP COLUMN cloud_decks"))
