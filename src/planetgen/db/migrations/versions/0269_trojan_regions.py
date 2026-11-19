"""The trojan swarm of a giant, stored as belt rows with the Lagrange point and the libration of their bodies (GEN.1150).

Schema v269. `asteroid_belts` gains `lagrange_point` (VARCHAR(2), `L4` or `L5`), `libration_amplitude_deg`, `libration_phase_deg`, `mean_eccentricity` and `mean_inclination_deg` (DOUBLE): a trojan region
of a giant (`generation/trojan_region.py`) is a belt-like row whose parent is the giant (`parent_type` `planet`). All five are NULL for every other belt and for a row made before v269; a
migration does not draw, and the regenerated galaxy holds the swarms. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0269"
down_revision = "0268"
branch_labels = None
depends_on = None

COLUMNS = (("lagrange_point", "VARCHAR(2)"), ("libration_amplitude_deg", "DOUBLE"), ("libration_phase_deg", "DOUBLE"), ("mean_eccentricity", "DOUBLE"), ("mean_inclination_deg", "DOUBLE"))


def upgrade():
    connection = op.get_bind()
    for name, kind in COLUMNS:
        if not helpers.exists_column(connection, "asteroid_belts", name):
            connection.execute(sa.text(f"ALTER TABLE asteroid_belts ADD COLUMN {name} {kind} NULL"))
