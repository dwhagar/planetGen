"""The life class switch (GEN.889, GEN.890): a body has a life class, not a life chemistry.

Schema v175. `planets` and `moons` lose `life_chemical` (and its index) and gain `life_class` (the number of the class, TINYINT UNSIGNED),
`life_class_variant` (its letter), `life_class_best` and `life_class_best_score`, `life_animal_grade` and `life_complexity`; the child tables
`planet_life_class_scores` and `moon_life_class_scores` (the classes a body supports) are new and empty; `species.life_chemical` becomes
`species.life_class` (the id of the class). A row made before v175 keeps the class its chemistry stood for (`CLASS_OF_CHEMISTRY`) and nothing
else: no scores, grade or complexity are worked out here, they come when the galaxy is regenerated. A column already there is skipped and
the chemistry is read only while it is still there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner
from planetgen.db.migrations import helpers

revision = "0175"
down_revision = "0174"
branch_labels = None
depends_on = None

CLASS_OF_CHEMISTRY = {
    "Chlorophyll a": ("1", "a"),
    "Blue-Optimized Porphyrins": ("1", "b"),
    "Bacteriochlorophylls": ("2", ""),
    "Retinal": ("2", ""),
    "Melanin": ("3", "b"),
    "Zinc-Bacteriochlorophyll": ("5", "b"),
}
"""dict: The class (number, letter) that each of the six chemistries of before v175 is folded into (docs/design/habitability/life-classes.md
section 3): the three pigments of Class 1 and Class 2 keep their light eaters, Retinal (a light supplement now) goes to the anoxygenic light
eaters it lived beside, Melanin (a radiation trait now) to the radiolytic hydrogen eaters and the acid pigment to the acid-world mats."""

COLUMNS = (
    ("life_class", "TINYINT UNSIGNED NULL"),
    ("life_class_variant", "CHAR(1) NULL"),
    ("life_class_best", "VARCHAR(2) NULL"),
    ("life_class_best_score", "TINYINT UNSIGNED NULL"),
    ("life_animal_grade", "TINYINT UNSIGNED NULL"),
    ("life_complexity", "DOUBLE NULL"),
)
STEPS = 7


def _case(number_part):
    whens = " ".join(f"WHEN '{name}' THEN {value!r}" for name, value in
                     ((n, c[0] if number_part else c[1]) for n, c in CLASS_OF_CHEMISTRY.items()))
    return f"CASE life_chemical {whens} ELSE NULL END"


def upgrade():
    connection = op.get_bind()
    done = 0
    for table in ("planets", "moons"):
        # no ALGORITHM=INSTANT: both tables carry a FULLTEXT index, which refuses it
        changes = [f"ADD COLUMN {column} {definition}" for column, definition in COLUMNS if not helpers.exists_column(connection, table, column)]
        if changes:
            connection.execute(sa.text(f"ALTER TABLE {table} {', '.join(changes)}"))
        done += 1
        alembic_runner.report_progress("Life classes", done, STEPS)
        if helpers.exists_column(connection, table, "life_chemical"):
            connection.execute(sa.text(
                f"UPDATE {table} SET life_class = {_case(True)}, life_class_variant = NULLIF({_case(False)}, '') WHERE life_chemical IS NOT NULL AND life_class IS NULL"))
            helpers.drop_index_online(connection, table, f"idx_{table}_life_chemical")
            connection.execute(sa.text(f"ALTER TABLE {table} DROP COLUMN life_chemical"))
        helpers.add_index_online(connection, table, f"idx_{table}_life_class", ["life_class"])
        done += 1
        alembic_runner.report_progress("Life classes", done, STEPS)
    if helpers.exists_column(connection, "species", "life_chemical"):
        if not helpers.exists_column(connection, "species", "life_class"):
            connection.execute(sa.text("ALTER TABLE species ADD COLUMN life_class VARCHAR(2) NULL"))
        whens = " ".join(f"WHEN '{name}' THEN '{number}{letter}'" for name, (number, letter) in CLASS_OF_CHEMISTRY.items())
        connection.execute(sa.text(f"UPDATE species SET life_class = CASE life_chemical {whens} ELSE NULL END WHERE life_class IS NULL"))
        connection.execute(sa.text("ALTER TABLE species DROP COLUMN life_chemical"))
    alembic_runner.report_progress("Life classes", STEPS, STEPS)
