"""The reasons a star system is worth a mark on the Galaxy Map (MAP.305, MAP.303).

Schema v165. `system_interests` replaces `habitable_systems` (v148): one row for each reason a star system is worth a mark, with how many
bodies give it, the largest score among them and the system's position, keyed on the position so the Galaxy Map's marker layer finds the
systems in a box without reading every planet. The reasons so far are a planet or moon rated Yellow or better overall on the PHI-6 score
(`habitable`; the old table held Blue and Green only) and a planet that is a gas giant (`gas-giant`). Every system already stored is counted
from the scores and classes its planets and moons hold, and the old table goes. A table already there is kept and filled again, so a
stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner
from planetgen.db.migrations import helpers
from planetgen.physics import habitability

revision = "0165"
down_revision = "0164"
branch_labels = None
depends_on = None

STEPS = 4


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_table(connection, "system_interests"):
        op.execute("""
            CREATE TABLE system_interests (
                star_system_id  BIGINT UNSIGNED NOT NULL,
                reason          VARCHAR(24) NOT NULL,
                n               SMALLINT UNSIGNED NOT NULL,
                best_score      DOUBLE NULL,
                galaxy_x_pc     DOUBLE NOT NULL,
                galaxy_y_pc     DOUBLE NOT NULL,
                galaxy_z_pc     DOUBLE NOT NULL,
                PRIMARY KEY (star_system_id, reason),
                KEY idx_system_interests_position (galaxy_x_pc, galaxy_y_pc, galaxy_z_pc),
                CONSTRAINT fk_system_interests_system
                    FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""")
    alembic_runner.report_progress("System interests", 1, STEPS)
    connection.execute(sa.text("DELETE FROM system_interests"))
    floor = habitability.YELLOW_FLOOR
    # habitable: a planet or moon rated Yellow or better overall
    for step, table in enumerate(("planets", "moons"), start=2):
        connection.execute(sa.text(
            "INSERT INTO system_interests (star_system_id, reason, n, best_score, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc)"
            " SELECT ss.id, 'habitable', LEAST(65535, c.n), c.best, ss.galaxy_x_pc, ss.galaxy_y_pc, ss.galaxy_z_pc FROM"
            f" (SELECT star_system_id, COUNT(*) AS n, MAX(phi4) AS best FROM {table} WHERE phi4 >= :floor GROUP BY star_system_id) c"
            " JOIN star_systems ss ON ss.id = c.star_system_id WHERE ss.galaxy_x_pc IS NOT NULL"
            " AND ss.galaxy_y_pc IS NOT NULL AND ss.galaxy_z_pc IS NOT NULL"
            " ON DUPLICATE KEY UPDATE n = LEAST(65535, system_interests.n + VALUES(n)),"
            " best_score = GREATEST(system_interests.best_score, VALUES(best_score))"),
            {"floor": floor})
        alembic_runner.report_progress("System interests", step, STEPS)
    # gas-giant: a planet whose body type is a gas giant's
    connection.execute(sa.text(
        "INSERT INTO system_interests (star_system_id, reason, n, best_score, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc)"
        " SELECT ss.id, 'gas-giant', LEAST(65535, c.n), NULL, ss.galaxy_x_pc, ss.galaxy_y_pc, ss.galaxy_z_pc FROM"
        " (SELECT star_system_id, COUNT(*) AS n FROM planets WHERE body_type = 'g' GROUP BY star_system_id) c"
        " JOIN star_systems ss ON ss.id = c.star_system_id WHERE ss.galaxy_x_pc IS NOT NULL"
        " AND ss.galaxy_y_pc IS NOT NULL AND ss.galaxy_z_pc IS NOT NULL"))
    connection.execute(sa.text("DROP TABLE IF EXISTS habitable_systems"))
    alembic_runner.report_progress("System interests", STEPS, STEPS)

