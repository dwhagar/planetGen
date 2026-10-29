"""The technology level of a civilization, on the Galaxy Map's civilization reason (MAP.344).

Schema v174. No column changes: the `civilization` rows of `system_interests` (v173) keep, as `best_score`, the technology level
(`species.tech_level`) of the most advanced species on a world the reason counts, set here from the species already stored; NULL for a
system whose civilizations have no species yet. The population pass sets it from now on (`planetgen.db.interests.refresh_technology`). The
map reads it as the band of the system (`tuning.TECH_BANDS`) and not as a habitability tier.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner, interests

revision = "0174"
down_revision = "0173"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    alembic_runner.report_progress("Civilization technology", 0, 1)
    connection.execute(sa.text(
        "UPDATE system_interests i SET i.best_score = (SELECT MAX(sp.tech_level) FROM species sp"
        " JOIN planets p ON p.id = sp.homeworld_planet_id"
        f" WHERE p.star_system_id = i.star_system_id AND {interests.REASONS[interests.CIVILIZATION]['planets']})"
        " WHERE i.reason = 'civilization'"))
    alembic_runner.report_progress("Civilization technology", 1, 1)
