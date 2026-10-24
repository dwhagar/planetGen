"""A planet or moon keeps four life scores, Microbial, Flora, Fauna and Sentient, instead of the three habitability scores (GEN.649).

Schema v150. `planets` and `moons` gain `life_microbial`, `life_flora`, `life_fauna` and `life_sentient` (TINYINT UNSIGNED,
NULL for a body with no score: 0 is no life of that kind, 9 abundant) and lose `phi_bio`, `phi_cpx` and `phi_tech`. A row
saved at v149 gets its four scores from the two likelihoods it stored (`phi_bio`, `phi_cpx`) and the evolutionary path it
stored (`life_chemical`, `life_stage`, `evolutionary_pace` and the two ages), by `planetgen.physics.life_scores`; the human
operability score is dropped with nothing replacing it.

The columns go in with one ALTER for each table, rows are filled in pages of ids (a row with `life_microbial` set is done,
so a stopped run starts again where it stopped), and the three old columns are dropped last.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner
from planetgen.db.migrations import helpers
from planetgen.physics import habitability_world, life_scores

revision = "0150"
down_revision = "0149"
branch_labels = None
depends_on = None

OLD_COLUMNS = ("phi_bio", "phi_cpx", "phi_tech")
"""tuple: The columns dropped from `planets` and `moons`."""

PAGE = 5000
"""int: Rows read, worked out and written before a commit."""

_SELECT = """
    SELECT id, phi_bio, phi_cpx, life_chemical, life_stage, evolutionary_pace, life_milestone_age_gyr, life_system_age_gyr
    FROM {table}
    WHERE id > :last AND phi4 IS NOT NULL AND life_microbial IS NULL
    ORDER BY id LIMIT :n"""


def worked_out(row):
    """The four scores of one stored row."""
    native = habitability_world.native_life(row["life_chemical"], row["life_stage"], row["evolutionary_pace"] is not None)
    stage = life_scores.stage_of(native, row["life_stage"], row["evolutionary_pace"] is not None)
    return life_scores.scores(row["phi_bio"] or 0.0, row["phi_cpx"] or 0.0, stage, row["life_milestone_age_gyr"],
                              row["life_system_age_gyr"])


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        changes = [f"ADD COLUMN {column} TINYINT UNSIGNED NULL" for column in life_scores.FIELDS
                   if not helpers.exists_column(connection, table, column)]
        if changes:  # no ALGORITHM=INSTANT: both tables carry a FULLTEXT index, which refuses it
            connection.execute(sa.text(f"ALTER TABLE {table} {', '.join(changes)}"))
    for table in ("planets", "moons"):
        if not helpers.exists_column(connection, table, "phi_bio"):
            continue
        query = sa.text(_SELECT.format(table=table))
        total = connection.execute(sa.text(
            f"SELECT COUNT(*) FROM {table} WHERE phi4 IS NOT NULL AND life_microbial IS NULL")).scalar()
        done = last = 0
        assignments = ", ".join(f"{column} = :{column}" for column in life_scores.FIELDS)
        while True:
            rows = connection.execute(query, {"last": last, "n": PAGE}).mappings().all()
            if not rows:
                break
            connection.execute(sa.text(f"UPDATE {table} SET {assignments} WHERE id = :id"), [
                {"id": row["id"], **dict(zip(life_scores.FIELDS, worked_out(row)))} for row in rows])
            connection.commit()
            last = rows[-1]["id"]
            done += len(rows)
            alembic_runner.report_progress(f"Life scores of {table}", done, total)
    for table in ("planets", "moons"):
        old = [column for column in OLD_COLUMNS if helpers.exists_column(connection, table, column)]
        if old:
            connection.execute(sa.text(f"ALTER TABLE {table} {', '.join(f'DROP COLUMN {column}' for column in old)}"))
