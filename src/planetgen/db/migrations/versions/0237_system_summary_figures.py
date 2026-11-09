"""What each star system adds up to: its star, planet and moon masses, luminosities and body counts (GEN.1478).

Schema v237. `system_summary` (v236) gains `stars_mass_kg`, `system_mass_kg`, `brightest_star_luminosity_w`, `total_luminosity_w`, `star_count`, `planet_count`,
`moon_count`, `belt_count` and `comet_count` (`db/system_summary.py`), and every system already stored gets its figures, in runs of `BATCH` systems. A column or
row already there is kept and written again with the same values, so a stopped run can be started again. The table is narrow and has no FULLTEXT key, so the
columns are added in place.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner
from planetgen.db.migrations import helpers

revision = "0237"
down_revision = "0236"
branch_labels = None
depends_on = None

BATCH = 20000
"""int: Systems (by id) per statement."""

COLUMNS = (("stars_mass_kg", "DOUBLE NULL"), ("system_mass_kg", "DOUBLE NULL"), ("brightest_star_luminosity_w", "DOUBLE NULL"),
           ("total_luminosity_w", "DOUBLE NULL"), ("star_count", "SMALLINT UNSIGNED NULL"), ("planet_count", "SMALLINT UNSIGNED NULL"),
           ("moon_count", "SMALLINT UNSIGNED NULL"), ("belt_count", "SMALLINT UNSIGNED NULL"), ("comet_count", "SMALLINT UNSIGNED NULL"))

SELECT = (
    "SELECT ss.id,"
    " (SELECT SUM(b.mass_kg) FROM stars b WHERE b.star_system_id = ss.id),"
    " COALESCE((SELECT SUM(b.mass_kg) FROM stars b WHERE b.star_system_id = ss.id), 0)"
    " + COALESCE((SELECT SUM(b.mass_kg) FROM planets b WHERE b.star_system_id = ss.id), 0)"
    " + COALESCE((SELECT SUM(b.mass_kg) FROM moons b WHERE b.star_system_id = ss.id), 0),"
    " (SELECT MAX(b.luminosity_w) FROM stars b WHERE b.star_system_id = ss.id),"
    " (SELECT SUM(b.luminosity_w) FROM stars b WHERE b.star_system_id = ss.id),"
    " (SELECT COUNT(*) FROM stars b WHERE b.star_system_id = ss.id),"
    " (SELECT COUNT(*) FROM planets b WHERE b.star_system_id = ss.id),"
    " (SELECT COUNT(*) FROM moons b WHERE b.star_system_id = ss.id),"
    " (SELECT COUNT(*) FROM asteroid_belts b WHERE b.star_system_id = ss.id),"
    " (SELECT COUNT(*) FROM comets b WHERE b.star_system_id = ss.id)"
    " FROM star_systems ss WHERE ss.id BETWEEN :low AND :high")


def upgrade():
    connection = op.get_bind()
    for name, definition in COLUMNS:
        helpers.add_column_instant(connection, "system_summary", name, definition)
    top = connection.execute(sa.text("SELECT COALESCE(MAX(id), 0) FROM star_systems")).scalar()
    names = ", ".join(name for name, _ in COLUMNS)
    updates = ", ".join(f"{name} = VALUES({name})" for name, _ in COLUMNS)
    for low in range(1, top + 1, BATCH):
        connection.execute(sa.text(
            f"INSERT INTO system_summary (star_system_id, {names}) {SELECT} ON DUPLICATE KEY UPDATE {updates}"),
            {"low": low, "high": low + BATCH - 1})
        alembic_runner.report_progress("System figures", min(low + BATCH - 1, top), top)
