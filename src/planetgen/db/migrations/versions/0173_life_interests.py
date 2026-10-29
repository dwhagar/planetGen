"""The life reasons a star system is worth a mark on the Galaxy Map (MAP.311).

Schema v173. No column changes: `system_interests` (v165) also holds the reasons `civilization` (a planet or moon whose timeline reached a
technological civilization, on a class that can carry life), `sentient` (a life score for sentient life of 1 or more) and `life` (any of the four
life scores of 1 or more), counted here from the planets and moons already stored. The rows keep the count and no score. A system that already has
a row for a reason is left as it is, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner, interests

revision = "0173"
down_revision = "0172"
branch_labels = None
depends_on = None

_REASONS = (interests.CIVILIZATION, interests.SENTIENT, interests.LIFE)
STEPS = 1 + 2 * len(_REASONS)


def upgrade():
    connection = op.get_bind()
    marks = ", ".join(f"'{reason}'" for reason in _REASONS)
    connection.execute(sa.text(f"DELETE FROM system_interests WHERE reason IN ({marks})"))
    done = 1
    alembic_runner.report_progress("Life interests", done, STEPS)
    for reason in _REASONS:
        for key, table in (("planets", "planets"), ("moons", "moons")):
            test = interests.REASONS[reason][key].replace("p.", "b.").replace("m.", "b.")
            connection.execute(sa.text(
                "INSERT INTO system_interests (star_system_id, reason, n, best_score, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc)"
                f" SELECT ss.id, '{reason}', LEAST(65535, c.n), NULL, ss.galaxy_x_pc, ss.galaxy_y_pc, ss.galaxy_z_pc FROM"
                f" (SELECT b.star_system_id AS star_system_id, COUNT(*) AS n FROM {table} b WHERE {test} GROUP BY b.star_system_id) c"
                " JOIN star_systems ss ON ss.id = c.star_system_id WHERE ss.galaxy_x_pc IS NOT NULL"
                " AND ss.galaxy_y_pc IS NOT NULL AND ss.galaxy_z_pc IS NOT NULL"
                " ON DUPLICATE KEY UPDATE n = LEAST(65535, system_interests.n + VALUES(n))"))
            done += 1
            alembic_runner.report_progress("Life interests", done, STEPS)
