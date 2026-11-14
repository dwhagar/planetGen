"""Structure columns for life and flavour on planets and moons (DB.50).

Schema v110. The evolutionary paragraphs and `flavor_text` say in prose what a body's life is and what its flavour line
is; the stage, the pace and the ages are only recoverable by reading the prose (`life_stage_from_paragraphs`, the
string searches in `decide_flavor_text`). This revision adds the structure that replaces them, on `planets` and
`moons`: `life_stage` (a `generation.evolution.MILESTONE_KEYS` key, NULL before any milestone), `evolutionary_pace`
(the timeline scale, `fast`, `normal` or `slow`), `life_milestone_age_gyr` and `life_system_age_gyr` (the age of the
milestone and the age of the system the timeline was written at, Gy) and `flavor_key` (which flavour line, replacing
`flavor_text`). Nothing is dropped and no row is touched: the columns are NULL until generation writes them (DB.51).
A column is added only where it is not there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0110"
down_revision = "0109"
branch_labels = None
depends_on = None

COLUMNS = (
    ("life_stage", "VARCHAR(32)"),
    ("evolutionary_pace", "VARCHAR(16)"),
    ("life_milestone_age_gyr", "DOUBLE"),
    ("life_system_age_gyr", "DOUBLE"),
    ("flavor_key", "VARCHAR(32)"),
)
"""tuple: `(column, type)`, added to `planets` and `moons` in this order, after `flavor_text_count`."""


def _has_column(connection, table, column):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar())


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        after = "flavor_text_count"
        for column, kind in COLUMNS:
            if not _has_column(connection, table, column):
                op.execute(f"ALTER TABLE {table} ADD COLUMN {column} {kind} NULL AFTER {after}")
            after = column
