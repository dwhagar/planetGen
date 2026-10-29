"""A third reason a star system is worth a mark on the Galaxy Map: an inhabited world (MAP.321).

Schema v167. No column changes: `system_interests` (v165) also holds the reason `inhabited`, counted here from the planets and moons already
stored the way the system page decides it: a world whose timeline reached a technological civilization (a class that can carry life only), or
that has a terrestrial colony facility on it. The row keeps the count and the best habitability score among them. A system that already has
a row for it is left as it is, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen import tuning
from planetgen.db import alembic_runner

revision = "0167"
down_revision = "0166"
branch_labels = None
depends_on = None

STEPS = 3

_CLASSES = ", ".join(f"'{letter}'" for letter in tuning.HABITABLE_PLANET_CLASSES)


def upgrade():
    connection = op.get_bind()
    connection.execute(sa.text("DELETE FROM system_interests WHERE reason = 'inhabited'"))
    alembic_runner.report_progress("Inhabited interests", 1, STEPS)
    for step, (table, host) in enumerate((("planets", "planet_id"), ("moons", "moon_id")), start=2):
        connection.execute(sa.text(
            "INSERT INTO system_interests (star_system_id, reason, n, best_score, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc)"
            " SELECT ss.id, 'inhabited', LEAST(65535, c.n), c.best, ss.galaxy_x_pc, ss.galaxy_y_pc, ss.galaxy_z_pc FROM"
            f" (SELECT b.star_system_id AS star_system_id, COUNT(*) AS n, MAX(b.phi4) AS best FROM {table} b"
            f" WHERE ((b.life_stage = 'technological_civilization' AND b.planet_class IN ({_CLASSES}))"
            f" OR EXISTS (SELECT 1 FROM facilities f WHERE f.{host} = b.id AND f.kind = 'colony' AND f.placement = 'terrestrial'))"
            " GROUP BY b.star_system_id) c"
            " JOIN star_systems ss ON ss.id = c.star_system_id WHERE ss.galaxy_x_pc IS NOT NULL"
            " AND ss.galaxy_y_pc IS NOT NULL AND ss.galaxy_z_pc IS NOT NULL"
            " ON DUPLICATE KEY UPDATE n = LEAST(65535, system_interests.n + VALUES(n)),"
            " best_score = GREATEST(COALESCE(system_interests.best_score, -1), COALESCE(VALUES(best_score), -1))"))
        alembic_runner.report_progress("Inhabited interests", step, STEPS)
