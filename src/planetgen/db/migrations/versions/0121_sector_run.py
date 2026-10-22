"""A sector records its generator epoch and run, not four version text columns (DB.16).

Schema v121. DB.7 (v67) put the code that made a sector on the sector as text: `version_key` (22 characters),
`planetgen_version`, `python_version` and `platform`, the same on every sector a run saves. This revision stores two small
columns instead and looks the text up from the run: `generator_epoch` (SMALLINT UNSIGNED, NULL; the epoch of
`util/generator_epoch.py` the code ran under, unknown for a sector saved before this revision) and `run_id` (BIGINT UNSIGNED,
NULL, a reference to `generation_runs.id`, which already holds the key, release, Python and platform of each run). A sector
with no recorded version (made before v67) keeps a NULL `run_id`.

A sector that has a version gets a run row to point at: each distinct set of the four values on the sectors becomes one
`generation_runs` row (command `migrated`, outcome `ok`, started and finished at the first such sector's `created_at`), and
the sectors of that set point at it, in batches. Then the four columns and their index (`idx_sectors_version_key`) are
dropped. Nothing else changes.

`sectors` has a FULLTEXT index (`ft_sectors_name`), and InnoDB refuses an instant or in-place change of a table that has
one (error 1846), so the two table changes (add the columns, drop the old ones) are each one `ALGORITHM=COPY` rebuild. It
holds one row per sector, not per object. Each step is skipped when it is already done, so a stopped run can be started
again: the sectors still without a `run_id` are the ones left to point.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0121"
down_revision = "0120"
branch_labels = None
depends_on = None

OLD_COLUMNS = ("version_key", "planetgen_version", "python_version", "platform")
"""tuple: The text columns that `generation_runs` already holds."""


def _point_sectors_at_runs(connection):
    groups = connection.execute(sa.text(
        "SELECT version_key, planetgen_version, python_version, platform, MIN(created_at) AS first_at FROM sectors"
        " WHERE version_key IS NOT NULL AND run_id IS NULL"
        " GROUP BY version_key, planetgen_version, python_version, platform")).mappings().all()
    for group in groups:
        values = {name: group[name] for name in OLD_COLUMNS}
        run_id = connection.execute(sa.text(
            "INSERT INTO generation_runs (command, version_key, planetgen_version, python_version, platform, started_at,"
            " finished_at, outcome) VALUES ('migrated', :version_key, COALESCE(:planetgen_version, ''),"
            " LEFT(COALESCE(:python_version, ''), 16), COALESCE(:platform, ''), :first_at, :first_at, 'ok')"),
            {**values, "first_at": group["first_at"]}).lastrowid
        connection.commit()
        helpers.batched_update(
            connection, "sectors", "run_id = :run_id",
            "run_id IS NULL AND version_key = :version_key AND planetgen_version <=> :planetgen_version"
            " AND python_version <=> :python_version AND platform <=> :platform",
            label="Pointing sectors at their run", params={**values, "run_id": run_id})


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_column(connection, "sectors", "run_id"):
        connection.execute(sa.text(
            "ALTER TABLE sectors ADD COLUMN generator_epoch SMALLINT UNSIGNED NULL, ADD COLUMN run_id BIGINT UNSIGNED NULL,"
            " ADD KEY idx_sectors_run (run_id),"
            " ADD CONSTRAINT fk_sectors_run FOREIGN KEY (run_id) REFERENCES generation_runs(id), ALGORITHM=COPY"))
    if helpers.exists_column(connection, "sectors", "version_key"):
        _point_sectors_at_runs(connection)
        drops = [f"DROP COLUMN {column}" for column in OLD_COLUMNS]
        if helpers.exists_index(connection, "sectors", "idx_sectors_version_key"):
            drops.insert(0, "DROP KEY idx_sectors_version_key")
        connection.execute(sa.text(f"ALTER TABLE sectors {', '.join(drops)}, ALGORITHM=COPY"))
