"""The prose is gone: the evolutionary paragraph tables and the stored note and flavour text (DB.52).

Schema v111. The life structure of v110 (`life_stage`, `evolutionary_pace`, the ages and `flavor_key`) replaced the
evolutionary paragraphs and `flavor_text`, the tier columns replaced `hab_note`, and `system_configs.markdown` was no
longer read (DB.51). This revision drops `planet_evolutionary_paragraphs` and `moon_evolutionary_paragraphs`, `flavor_text`
and `hab_note` from `planets` and `moons`, and `system_configs.markdown`. `flavor_text_count` stays: the generator still
counts the line it draws. A table or column already gone is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0111"
down_revision = "0110"
branch_labels = None
depends_on = None

TABLES = ("planet_evolutionary_paragraphs", "moon_evolutionary_paragraphs")
COLUMNS = (("planets", "flavor_text"), ("planets", "hab_note"), ("moons", "flavor_text"), ("moons", "hab_note"),
           ("system_configs", "markdown"))
"""tuple: `(table, column)` dropped."""


def _has_column(connection, table, column):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar())


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        op.execute(f"DROP TABLE IF EXISTS {table}")
    for table, column in COLUMNS:
        if _has_column(connection, table, column):
            op.execute(f"ALTER TABLE {table} DROP COLUMN {column}")
