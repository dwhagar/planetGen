"""Two more reasons a star system is worth a mark on the Galaxy Map (MAP.321).

Schema v166. No column changes: `system_interests` (v165) also holds the reasons `uninhabitable` (a scored terrestrial planet or moon under
the Yellow floor, `habitability.YELLOW_FLOOR`) and `asteroid-field` (a system with an asteroid belt), counted here from the planets, moons and
belts already stored. A system that already has a row for either is left as it is, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner
from planetgen.physics import habitability

revision = "0166"
down_revision = "0165"
branch_labels = None
depends_on = None

STEPS = 4

_INSERT = ("INSERT INTO system_interests (star_system_id, reason, n, best_score, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc)"
           " SELECT ss.id, :reason, LEAST(65535, c.n), c.best, ss.galaxy_x_pc, ss.galaxy_y_pc, ss.galaxy_z_pc FROM ({source}) c"
           " JOIN star_systems ss ON ss.id = c.star_system_id WHERE ss.galaxy_x_pc IS NOT NULL"
           " AND ss.galaxy_y_pc IS NOT NULL AND ss.galaxy_z_pc IS NOT NULL"
           " ON DUPLICATE KEY UPDATE n = LEAST(65535, system_interests.n + VALUES(n)),"
           " best_score = GREATEST(COALESCE(system_interests.best_score, -1), COALESCE(VALUES(best_score), -1))")


def upgrade():
    connection = op.get_bind()
    floor = habitability.YELLOW_FLOOR
    connection.execute(sa.text("DELETE FROM system_interests WHERE reason IN ('uninhabitable', 'asteroid-field')"))
    alembic_runner.report_progress("More system interests", 1, STEPS)
    for step, (table, test) in enumerate((("planets", "body_type = 't' AND phi4 < :floor"), ("moons", "phi4 < :floor")), start=2):
        connection.execute(sa.text(_INSERT.format(
            source=f"SELECT star_system_id, COUNT(*) AS n, MAX(phi4) AS best FROM {table} WHERE {test} GROUP BY star_system_id")),
            {"reason": "uninhabitable", "floor": floor})
        alembic_runner.report_progress("More system interests", step, STEPS)
    connection.execute(sa.text(_INSERT.format(
        source="SELECT star_system_id, COUNT(*) AS n, NULL AS best FROM asteroid_belts GROUP BY star_system_id")),
        {"reason": "asteroid-field"})
    alembic_runner.report_progress("More system interests", STEPS, STEPS)
