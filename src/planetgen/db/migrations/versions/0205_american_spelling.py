"""American spelling of the stored values and one column (DOC.103).

Schema v205. The stored words that kept a British spelling after the visible text was fixed (DOC.102) take the American one, with no alias:
`planets.hydrosphere` and `moons.hydrosphere` `vapour` becomes `vapor`; `bright_stars.population` `nuclear_disc` becomes `nuclear_disk`;
`asteroid_belts.origin_cause` and the `co_causes` list `outer-disc` become `outer-disk`; and the column `ground_oxidised_fraction` of `planets`
and `moons` is renamed `ground_oxidized_fraction`. The CHECKs that list a renamed value (`planets_chk_5` and `moons_chk_5` for the hydrosphere,
`chk_bright_stars_population`) are replaced by ones of the same names that list the new words. A row keeps its other values; the generator draws the
same worlds, so no regeneration is needed. Every step skips work that is already done, so a stopped run can be started again.
"""

import re

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0205"
down_revision = "0204"
branch_labels = None
depends_on = None

VALUES = (
    ("planets", "hydrosphere", "vapour", "vapor"),
    ("moons", "hydrosphere", "vapour", "vapor"),
    ("bright_stars", "population", "nuclear_disc", "nuclear_disk"),
)
"""tuple: `(table, column, old value, new value)` of a stored word with a CHECK."""

COLUMNS = (
    ("planets", "ground_oxidised_fraction", "ground_oxidized_fraction"),
    ("moons", "ground_oxidised_fraction", "ground_oxidized_fraction"),
)
"""tuple: `(table, old column, new column)`."""


_WORDS = re.compile(r"\\?'([^'\\]+)\\?'")
"""re.Pattern: One quoted word of a CHECK clause as information_schema prints it (`_utf8mb4\\'dry\\'`)."""


def _check_names(connection, table, old):
    """The names and clauses of the CHECKs of `table` that still list `old`."""
    return connection.execute(sa.text(
        "SELECT c.CONSTRAINT_NAME, c.CHECK_CLAUSE FROM information_schema.CHECK_CONSTRAINTS c"
        " JOIN information_schema.TABLE_CONSTRAINTS t ON t.CONSTRAINT_SCHEMA = c.CONSTRAINT_SCHEMA AND t.CONSTRAINT_NAME = c.CONSTRAINT_NAME"
        " WHERE c.CONSTRAINT_SCHEMA = DATABASE() AND t.TABLE_NAME = :table AND t.CONSTRAINT_TYPE = 'CHECK' AND c.CHECK_CLAUSE LIKE :like"),
        {"table": table, "like": f"%{old}%"}).all()


def upgrade():
    connection = op.get_bind()
    for table, column, old, new in VALUES:
        for name, clause in _check_names(connection, table, old):
            op.execute(f"ALTER TABLE {table} DROP CHECK {name}")
            listed = ", ".join(f"'{new if word == old else word}'" for word in _WORDS.findall(clause))
            op.execute(f"ALTER TABLE {table} ADD CONSTRAINT {name} CHECK ({column} IN ({listed})) NOT ENFORCED")    # rows still old: enforced again below
        helpers.batched_update(connection, table, f"{column} = :new", f"{column} = :old", params={"old": old, "new": new},
                               label=f"{table} {column}")
        for name, _clause in connection.execute(sa.text(
                "SELECT c.CONSTRAINT_NAME, c.CHECK_CLAUSE FROM information_schema.CHECK_CONSTRAINTS c"
                " JOIN information_schema.TABLE_CONSTRAINTS t ON t.CONSTRAINT_SCHEMA = c.CONSTRAINT_SCHEMA AND t.CONSTRAINT_NAME = c.CONSTRAINT_NAME"
                " WHERE c.CONSTRAINT_SCHEMA = DATABASE() AND t.TABLE_NAME = :table AND t.CONSTRAINT_TYPE = 'CHECK' AND c.CHECK_CLAUSE LIKE :like"
                " AND t.ENFORCED = 'NO'"), {"table": table, "like": f"%{new}%"}).all():
            op.execute(f"ALTER TABLE {table} ALTER CHECK {name} ENFORCED")
    if helpers.exists_column(connection, "asteroid_belts", "origin_cause"):
        connection.execute(sa.text("UPDATE asteroid_belts SET origin_cause = 'outer-disk' WHERE origin_cause = 'outer-disc'"))
        connection.execute(sa.text("UPDATE asteroid_belts SET co_causes = REPLACE(co_causes, 'outer-disc', 'outer-disk') WHERE co_causes LIKE '%outer-disc%'"))
    if helpers.exists_column(connection, "asteroid_belts", "origin_sentence"):
        for column in ("origin_sentence", "origin_history"):    # the sentences of the belt origins (`the disc past 3 AU`), whole word only
            connection.execute(sa.text(
                f"UPDATE asteroid_belts SET {column} = REGEXP_REPLACE({column}, '\\\\b(d|D)isc(s?)\\\\b', '$1isk$2')"
                f" WHERE {column} REGEXP '\\\\b[dD]isc'"))
    for table, old, new in COLUMNS:
        if helpers.exists_column(connection, table, old) and not helpers.exists_column(connection, table, new):
            op.execute(f"ALTER TABLE {table} RENAME COLUMN {old} TO {new}, ALGORITHM=INPLACE")
