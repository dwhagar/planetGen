"""A planet or moon keeps the whole milestone schedule of its evolutionary pace (DB.100).

Schema v156. `planets` and `moons` gain `life_age_abiogenesis_gyr`, `life_age_photosynthesis_gyr`, `life_age_complex_cells_gyr`,
`life_age_multicellularity_gyr` and `life_age_technological_civilization_gyr` (DOUBLE, Gyr, NULL for a body with no life
chemistry): the age at which its pace reaches each milestone, so the life paragraph and the timeline paragraph are written from
the row alone and a later tuning of the pace table moves no saved world. A row saved at v155 gets the ages of its stored
`evolutionary_pace` (the ages of `tuning.EVOLUTIONARY_TIMELINES`, one for each `tuning.LIFE_MILESTONES`); a row with no pace (or one the table does not know) stays NULL.

The columns go in with one ALTER for each table, rows are filled in pages of ids (a row with `life_age_abiogenesis_gyr` set is
done, so a stopped run starts again where it stopped).
"""

import sqlalchemy as sa
from alembic import op

from planetgen import tuning
from planetgen.db import alembic_runner
from planetgen.db.migrations import helpers

revision = "0156"
down_revision = "0155"
branch_labels = None
depends_on = None


def _ages(pace):
    """The ages (Gyr) of a pace of `tuning.EVOLUTIONARY_TIMELINES`, one for each milestone."""
    return tuple(tuning.EVOLUTIONARY_TIMELINES[pace][key] for key in tuning.LIFE_MILESTONES)


PAGE = 5000
"""int: Rows read, worked out and written before a commit."""

_FIRST = tuning.LIFE_SCHEDULE_FIELDS[0]

_SELECT = f"""
    SELECT id, evolutionary_pace FROM {{table}}
    WHERE id > :last AND evolutionary_pace IS NOT NULL AND {_FIRST} IS NULL
    ORDER BY id LIMIT :n"""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        changes = [f"ADD COLUMN {column} DOUBLE NULL" for column in tuning.LIFE_SCHEDULE_FIELDS
                   if not helpers.exists_column(connection, table, column)]
        if changes:  # no ALGORITHM=INSTANT: both tables carry a FULLTEXT index, which refuses it
            connection.execute(sa.text(f"ALTER TABLE {table} {', '.join(changes)}"))
    assignments = ", ".join(f"{column} = :{column}" for column in tuning.LIFE_SCHEDULE_FIELDS)
    for table in ("planets", "moons"):
        query = sa.text(_SELECT.format(table=table))
        total = connection.execute(sa.text(
            f"SELECT COUNT(*) FROM {table} WHERE evolutionary_pace IS NOT NULL AND {_FIRST} IS NULL")).scalar()
        done = last = 0
        while True:
            rows = connection.execute(query, {"last": last, "n": PAGE}).mappings().all()
            if not rows:
                break
            known = [row for row in rows if row["evolutionary_pace"] in tuning.EVOLUTIONARY_TIMELINES]  # an unknown pace stays NULL
            if known:
                connection.execute(sa.text(f"UPDATE {table} SET {assignments} WHERE id = :id"), [
                    {"id": row["id"], **dict(zip(tuning.LIFE_SCHEDULE_FIELDS, _ages(row["evolutionary_pace"])))}
                    for row in known])
            connection.commit()
            last = rows[-1]["id"]
            done += len(rows)
            alembic_runner.report_progress(f"Life schedule of {table}", done, total)
