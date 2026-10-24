"""The table of the star systems that hold a habitable world (MAP.220).

Schema v148. `habitable_systems` has a row for each star system with a planet or moon whose overall PHI-6 score is Blue or Green
(`habitability.POTENTIALLY_HABITABLE_MIN_SCORE`): how many, and its position, with a key on the position so the Galaxy Map's marker layer
finds the ones inside a box without reading every planet. (`star_systems` has a FULLTEXT key, so adding a column there would copy the
whole table.) Every system already stored is counted from the scores its planets and moons hold. A table already there is kept and
filled again, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner
from planetgen.db.migrations import helpers
from planetgen.physics import habitability

revision = "0148"
down_revision = "0147"
branch_labels = None
depends_on = None

STEPS = 3


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_table(connection, "habitable_systems"):
        op.execute("""
            CREATE TABLE habitable_systems (
                star_system_id    BIGINT UNSIGNED NOT NULL PRIMARY KEY,
                habitable_worlds  SMALLINT UNSIGNED NOT NULL,
                galaxy_x_pc       DOUBLE NOT NULL,
                galaxy_y_pc       DOUBLE NOT NULL,
                galaxy_z_pc       DOUBLE NOT NULL,
                KEY idx_habitable_systems_position (galaxy_x_pc, galaxy_y_pc, galaxy_z_pc),
                CONSTRAINT fk_habitable_systems_system
                    FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""")
    alembic_runner.report_progress("Habitable systems", 1, STEPS)
    connection.execute(sa.text("DELETE FROM habitable_systems"))
    floor = habitability.POTENTIALLY_HABITABLE_MIN_SCORE
    for step, table in enumerate(("planets", "moons"), start=2):
        connection.execute(sa.text(
            "INSERT INTO habitable_systems (star_system_id, habitable_worlds, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc)"
            " SELECT ss.id, LEAST(65535, c.n), ss.galaxy_x_pc, ss.galaxy_y_pc, ss.galaxy_z_pc FROM"
            f" (SELECT star_system_id, COUNT(*) AS n FROM {table} WHERE phi4 >= :floor GROUP BY star_system_id) c"
            " JOIN star_systems ss ON ss.id = c.star_system_id WHERE ss.galaxy_x_pc IS NOT NULL"
            " AND ss.galaxy_y_pc IS NOT NULL AND ss.galaxy_z_pc IS NOT NULL"
            " ON DUPLICATE KEY UPDATE habitable_worlds = LEAST(65535, habitable_systems.habitable_worlds + VALUES(habitable_worlds))"),
            {"floor": floor})
        alembic_runner.report_progress("Habitable systems", step, STEPS)
