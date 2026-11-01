"""The chart ledger, and a sector row that cannot vanish under its residents (GEN.1127).

Schema v197. A new table, `chart_ledger`: one row for each sector cell that has been charted, with the run and the generator epoch that
charted it, when, and when the rogue top-up of its region run filled it (`topped_up_at`, NULL until it does). It replaces
`sector_stats.filled_at`, which this revision drops after copying every filled cell into the ledger (the run and the epoch come from the
cell's `sectors` row, the time from `filled_at`), so the charted state of a cell is a row of the ledger and there is no second record
to keep in step.

The `ON DELETE CASCADE` of `rogue_planets`, `interstellar_comets`, `supernova_remnants` and `facilities` on `sector_id` becomes
`RESTRICT`. `sector_id` is the sector that holds an object now (the nightly orbital update refiles it), not the one it was born in, so a
cascade deleted every resident that had drifted in. A sector row now cannot go while a resident points at it; deleting or regenerating a
sector removes only what was born there (`db/edits.delete_sector_birth_contents`). The foreign keys are changed with
`foreign_key_checks` off, in place and without a table copy: the rows already satisfy the rule, since it only stops a delete. A key that
already has the rule, or a table that does not exist yet, is skipped (the key is found by its column, whatever an old schema named it), and a key that was dropped and not yet added is added by the next run, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0197"
down_revision = "0196"
branch_labels = None
depends_on = None

TABLE = """
    CREATE TABLE IF NOT EXISTS chart_ledger (
        shell_index       INT NOT NULL,
        band_index        INT NOT NULL,
        wedge_index       INT NOT NULL,
        run_id            BIGINT UNSIGNED NULL,
        generator_epoch   SMALLINT UNSIGNED NULL,
        charted_at        TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
        topped_up_at      TIMESTAMP(3) NULL,
        PRIMARY KEY (shell_index, band_index, wedge_index),
        KEY idx_chart_ledger_run (run_id),
        CONSTRAINT fk_chart_ledger_run
            FOREIGN KEY (run_id) REFERENCES generation_runs(id) ON DELETE SET NULL
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""

BACKFILL = """
    INSERT IGNORE INTO chart_ledger (shell_index, band_index, wedge_index, run_id, generator_epoch, charted_at)
    SELECT st.shell_index, st.band_index, st.wedge_index, s.run_id, s.generator_epoch, st.filled_at
    FROM sector_stats st
    LEFT JOIN sectors s ON s.shell_index = st.shell_index AND s.band_index = st.band_index AND s.wedge_index = st.wedge_index
    WHERE st.filled_at IS NOT NULL"""

RESTRICTED = ("supernova_remnants", "rogue_planets", "interstellar_comets", "facilities")
"""tuple: The tables whose `sector_id` key stops cascading."""


def _sector_key(connection, table):
    """`(name, delete rule)` of the foreign key from `table.sector_id` to `sectors(id)`, or `None` when the table or the key is not there."""
    row = connection.execute(sa.text(
        "SELECT k.CONSTRAINT_NAME AS name, r.DELETE_RULE AS rule FROM information_schema.KEY_COLUMN_USAGE k"
        " JOIN information_schema.REFERENTIAL_CONSTRAINTS r ON r.CONSTRAINT_SCHEMA = k.CONSTRAINT_SCHEMA"
        " AND r.CONSTRAINT_NAME = k.CONSTRAINT_NAME AND r.TABLE_NAME = k.TABLE_NAME"
        " WHERE k.TABLE_SCHEMA = DATABASE() AND k.TABLE_NAME = :t AND k.COLUMN_NAME = 'sector_id'"
        " AND k.REFERENCED_TABLE_NAME = 'sectors'"), {"t": table}).first()
    return None if row is None else (row[0], row[1])


def upgrade():
    connection = op.get_bind()
    connection.execute(sa.text(TABLE))    # `schema.sql` is applied first, so the table is usually there already
    if helpers.exists_column(connection, "sector_stats", "filled_at"):
        connection.execute(sa.text(BACKFILL))
        helpers.drop_column_instant(connection, "sector_stats", "filled_at")
    connection.execute(sa.text("SET foreign_key_checks = 0"))
    try:
        for table in RESTRICTED:
            if not helpers.exists_column(connection, table, "sector_id"):
                continue
            key = _sector_key(connection, table)
            if key is not None and key[1] in ("RESTRICT", "NO ACTION"):
                continue
            if key is not None:
                connection.execute(sa.text(f"ALTER TABLE {table} DROP FOREIGN KEY {key[0]}, ALGORITHM=INPLACE, LOCK=NONE"))
            connection.execute(sa.text(
                f"ALTER TABLE {table} ADD CONSTRAINT fk_{table}_sector "
                "FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE RESTRICT, ALGORITHM=INPLACE, LOCK=NONE"))
    finally:
        connection.execute(sa.text("SET foreign_key_checks = 1"))


def downgrade():
    raise NotImplementedError("the migrations of this database only go forward")
