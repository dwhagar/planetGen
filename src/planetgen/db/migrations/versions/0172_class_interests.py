"""The planet classes as reasons a star system is worth a mark on the Galaxy Map (MAP.220).

Schema v172. No column changes: `system_interests` (v165) also holds one reason for each planet class, `class-A` to `class-V`
(`planetgen.db.interests.CLASS_REASONS`), counted here from the planets and moons already stored: a system has a row for each class it
has a planet or moon of, with how many. The rows have no score. A system that already has a row for a class is left as it is, so a stopped
run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen import tuning
from planetgen.db import alembic_runner

revision = "0172"
down_revision = "0171"
branch_labels = None
depends_on = None

STEPS = 3

_CLASSES = ", ".join(f"'{letter}'" for letter in sorted(tuning.PLANET_CLASSES))


def upgrade():
    connection = op.get_bind()
    connection.execute(sa.text("DELETE FROM system_interests WHERE reason LIKE 'class-%'"))
    alembic_runner.report_progress("Class interests", 1, STEPS)
    for step, table in enumerate(("planets", "moons"), start=2):
        connection.execute(sa.text(
            "INSERT INTO system_interests (star_system_id, reason, n, best_score, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc)"
            " SELECT ss.id, CONCAT('class-', c.planet_class), LEAST(65535, c.n), NULL, ss.galaxy_x_pc, ss.galaxy_y_pc, ss.galaxy_z_pc FROM"
            f" (SELECT b.star_system_id AS star_system_id, b.planet_class AS planet_class, COUNT(*) AS n FROM {table} b"
            f" WHERE b.planet_class IN ({_CLASSES}) GROUP BY b.star_system_id, b.planet_class) c"
            " JOIN star_systems ss ON ss.id = c.star_system_id WHERE ss.galaxy_x_pc IS NOT NULL"
            " AND ss.galaxy_y_pc IS NOT NULL AND ss.galaxy_z_pc IS NOT NULL"
            " ON DUPLICATE KEY UPDATE n = LEAST(65535, system_interests.n + VALUES(n))"))
        alembic_runner.report_progress("Class interests", step, STEPS)
