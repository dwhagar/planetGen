"""Values in their own columns (DB.13): a run's command line, and the search indexes.

Schema v62. `generation_runs.arguments` (one JSON text) becomes one
`generation_run_arguments` row per argument, and the size and spectral
columns the search filters on get indexes. `schema.sql` has already created
the new table by the time this runs, so each step first checks what exists.
"""

import json

import sqlalchemy as sa
from alembic import op

revision = "0062"
down_revision = "0061"
branch_labels = None
depends_on = None

NEW_INDEXES = (
    ("star_systems", "idx_star_systems_quadrant", "quadrant"),
    ("stars", "idx_stars_radius_km", "radius_km"),
    ("stars", "idx_stars_star_type", "star_type"),
    ("planets", "idx_planets_radius_km", "radius_km"),
    ("moons", "idx_moons_radius_km", "radius_km"),
)
ARGUMENT_LENGTH = 1024


def _exists(connection, sql, **params):
    return connection.execute(sa.text(sql), params).scalar() > 0


def upgrade():
    connection = op.get_bind()
    connection.execute(sa.text(
        "CREATE TABLE IF NOT EXISTS generation_run_arguments ("
        " run_id BIGINT UNSIGNED NOT NULL, position INT NOT NULL, value VARCHAR(1024) NOT NULL,"
        " PRIMARY KEY (run_id, position),"
        " CONSTRAINT fk_generation_run_arguments_run FOREIGN KEY (run_id)"
        " REFERENCES generation_runs(id) ON DELETE CASCADE"
        ") ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"))
    if _exists(connection, "SELECT COUNT(*) FROM information_schema.columns WHERE table_schema = DATABASE()"
                           " AND table_name = 'generation_runs' AND column_name = 'arguments'"):
        for run_id, arguments in connection.execute(sa.text("SELECT id, arguments FROM generation_runs")).fetchall():
            try:
                values = [str(value) for value in json.loads(arguments)]
            except (TypeError, ValueError):
                values = []
            for position, value in enumerate(values):
                connection.execute(sa.text(
                    "INSERT IGNORE INTO generation_run_arguments (run_id, position, value)"
                    " VALUES (:run_id, :position, :value)"),
                    {"run_id": run_id, "position": position, "value": value[:ARGUMENT_LENGTH]})
        connection.execute(sa.text("ALTER TABLE generation_runs DROP COLUMN arguments"))
    for table, name, column in NEW_INDEXES:
        if not _exists(connection, "SELECT COUNT(*) FROM information_schema.statistics WHERE table_schema = DATABASE()"
                                   " AND table_name = :tbl AND index_name = :idx", tbl=table, idx=name):
            connection.execute(sa.text(f"CREATE INDEX {name} ON {table} ({column})"))
