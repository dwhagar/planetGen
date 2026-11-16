"""The birth mass of a star and the eccentricity a resonance holds (DB.132, DB.133).

Schema v268. `stars` gains `initial_mass_sol` (the mass the population model gave the star at birth, solar masses) and `phase_end_age_gy` (the age the star leaves the phase it was drawn in,
Gyr; NULL with an initial mass means it never does, a white dwarf's), the two numbers a star restored from the database lacked (`remnant_orbits` and `place_terms` read the birth mass, and the
age a star is drawn at stays below the phase end). `body_history` gains `held_eccentricity` (the eccentricity a resonance holds a planet or moon on, `_resonance_e`, the floor `settle_orbit` keeps
under the tide) and `forced_eccentricity` (the eccentricity the planets of a resonant chain force on a planet, the quantity oddity W5 reads), so a restored body that is settled or held again
starts from what it was held on. A galaxy made earlier has NULL in the new columns, which reads as 'not kept'; nothing is backfilled (the numbers of a body made earlier are stored nowhere). A
column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0268"
down_revision = "0267"
branch_labels = None
depends_on = None

COLUMNS = (
    ("stars", "initial_mass_sol", "DOUBLE NULL CHECK (initial_mass_sol IS NULL OR initial_mass_sol > 0)", None),
    ("stars", "phase_end_age_gy", "DOUBLE NULL CHECK (phase_end_age_gy IS NULL OR phase_end_age_gy > 0)", None),
    ("body_history", "held_eccentricity", "DOUBLE NULL CHECK (held_eccentricity IS NULL OR (held_eccentricity >= 0 AND held_eccentricity < 1))", "born_index"),
    ("body_history", "forced_eccentricity", "DOUBLE NULL CHECK (forced_eccentricity IS NULL OR (forced_eccentricity >= 0 AND forced_eccentricity < 1))", "held_eccentricity"),
)


def _has_column(table, name):
    return bool(op.get_bind().execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = :table AND COLUMN_NAME = :name"), {"table": table, "name": name}).scalar())


def upgrade():
    for table, name, definition, after in COLUMNS:
        if not _has_column(table, name):
            place = f" AFTER {after}" if after else ""
            op.execute(f"ALTER TABLE {table} ADD COLUMN {name} {definition}{place}")


def downgrade():
    for table, name, _definition, _after in reversed(COLUMNS):
        if _has_column(table, name):
            op.execute(f"ALTER TABLE {table} DROP COLUMN {name}")
