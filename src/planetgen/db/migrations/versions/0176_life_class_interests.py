"""The life classes as reasons a star system is worth a mark on the Galaxy Map (MAP.351).

Schema v176. No column changes: `system_interests` (v165) also holds one reason for each life class, `lifeclass-1a` to `lifeclass-9`
(`planetgen.db.interests.LIFE_CLASS_REASONS`), counted here from the planets and moons already stored: a system has a row for each life class
that lives on one of its planets or moons (`life_class` and `life_class_variant`), with how many. The rows have no score. A system that
already has a row for a class is left as it is, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner, interests
from planetgen.physics import life_classes

revision = "0176"
down_revision = "0175"
branch_labels = None
depends_on = None

STEPS = 3

_NUMBERS = ", ".join(str(number) for number in sorted({definition.number for definition in life_classes.CLASSES}))


def upgrade():
    connection = op.get_bind()
    connection.execute(sa.text(f"DELETE FROM system_interests WHERE reason LIKE '{interests.LIFE_CLASS_PREFIX}%'"))
    alembic_runner.report_progress("Life class interests", 1, STEPS)
    for step, table in enumerate(("planets", "moons"), start=2):
        connection.execute(sa.text(
            "INSERT INTO system_interests (star_system_id, reason, n, best_score, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc)"
            f" SELECT ss.id, CONCAT('{interests.LIFE_CLASS_PREFIX}', c.life_class, c.variant), LEAST(65535, c.n), NULL,"
            " ss.galaxy_x_pc, ss.galaxy_y_pc, ss.galaxy_z_pc FROM"
            " (SELECT b.star_system_id AS star_system_id, b.life_class AS life_class, COALESCE(b.life_class_variant, '') AS variant,"
            f" COUNT(*) AS n FROM {table} b WHERE b.life_class IN ({_NUMBERS})"
            " GROUP BY b.star_system_id, b.life_class, COALESCE(b.life_class_variant, '')) c"
            " JOIN star_systems ss ON ss.id = c.star_system_id WHERE ss.galaxy_x_pc IS NOT NULL"
            " AND ss.galaxy_y_pc IS NOT NULL AND ss.galaxy_z_pc IS NOT NULL"
            " ON DUPLICATE KEY UPDATE n = LEAST(65535, system_interests.n + VALUES(n))"))
        alembic_runner.report_progress("Life class interests", step, STEPS)
