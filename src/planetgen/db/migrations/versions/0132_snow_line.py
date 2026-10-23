"""Each star's own snow line, and whether a planet or moon formed beyond it (GEN.458).

Schema v132. `stars` gains `snow_line_km` (DOUBLE, not null), the star's own snow line, drawn once around the 170 K blackbody position for its
luminosity (0.2 dex of scatter) and stored with it; `planets` and `moons` gain `born_beyond_snow_line` (BOOLEAN, nullable), whether the body formed
beyond its star's line (a moon takes its planet's). The generator no longer has a fixed 2.7 AU constant: the water share (GEN.428), the small-body
zones (GEN.442) and the moons (GEN.466) read the stored values.

A row made before this revision gets the blackbody position for its luminosity (no scatter: nothing was drawn for it) and the flag from its distance
against that line; a circumbinary planet (no `star_id`) is compared with the pair's line, the square root of the sum of the squares of its stars'. The
regenerated galaxy (the generator epoch of the 8.4 chain) replaces them with drawn values. Each step is skipped when it is already done, so a stopped
run can be started again; the `stars` column goes last and marks the migration done.
"""

import sqlalchemy as sa
from alembic import op

revision = "0132"
down_revision = "0131"
branch_labels = None
depends_on = None

# sqrt(L / (16 pi sigma)) / T^2 in metres for T = 170 K, then to km (the blackbody position `formation.snow_line_au` computes).
_BLACKBODY_KM = "SQRT(luminosity_w / (16 * PI() * 5.670374419e-8)) / (170 * 170) / 1000"


def _has_column(connection, table, column):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar()


def upgrade():
    connection = op.get_bind()
    if not _has_column(connection, "stars", "snow_line_km"):
        op.execute("ALTER TABLE stars ADD COLUMN snow_line_km DOUBLE NULL")
        op.execute(f"UPDATE stars SET snow_line_km = {_BLACKBODY_KM}")
    for table in ("planets", "moons"):
        if not _has_column(connection, table, "born_beyond_snow_line"):
            op.execute(f"ALTER TABLE {table} ADD COLUMN born_beyond_snow_line BOOLEAN NULL")
    op.execute(
        "UPDATE planets p"
        " LEFT JOIN stars st ON st.id = p.star_id"
        " LEFT JOIN (SELECT star_system_id, SQRT(SUM(snow_line_km * snow_line_km)) AS pair_km FROM stars GROUP BY star_system_id) pair"
        "   ON pair.star_system_id = p.star_system_id"
        " SET p.born_beyond_snow_line = p.distance_km >= COALESCE(st.snow_line_km, pair.pair_km)"
        " WHERE p.born_beyond_snow_line IS NULL AND p.body_type <> 'a'")
    op.execute(
        "UPDATE moons m JOIN planets p ON p.id = m.planet_id"
        " SET m.born_beyond_snow_line = p.born_beyond_snow_line WHERE m.born_beyond_snow_line IS NULL")
    op.execute("ALTER TABLE stars MODIFY COLUMN snow_line_km DOUBLE NOT NULL")
