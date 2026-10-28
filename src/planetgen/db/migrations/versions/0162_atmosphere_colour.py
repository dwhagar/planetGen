"""The atmosphere colour of a planet or moon (DB.96).

Schema v162. `planets` and `moons` gain `atmosphere_colour` (INT UNSIGNED, nullable): the colour of the body's air as one sRGB integer,
`red * 65536 + green * 256 + blue` (`planetgen.physics.body_colour`). The cloud cover it was worked out with is `cloud_cover` (v138).

NULL here: the colour is worked out with the planet at the generator epoch, from its air, its star, its ozone and its aerosols, so a row made
before v162 has none and the regenerated galaxy replaces it. Each column is skipped when it is already there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0162"
down_revision = "0161"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        if not helpers.exists_column(connection, table, "atmosphere_colour"):
            # no ALGORITHM=INSTANT: both tables carry a FULLTEXT index, which refuses it
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN atmosphere_colour INT UNSIGNED NULL"))
