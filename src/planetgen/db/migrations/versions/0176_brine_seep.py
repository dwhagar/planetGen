"""The brine seep, a liquid surface that is not a sea (GEN.908).

Schema v176. `planets.hydrosphere` and `moons.hydrosphere` take a seventh value, `brine seep`: a liquid surface whose ocean class is the
chloride brine (surface water under 2e-5 of the mass, GEN.762), which the water state used to call a `surface ocean`. The CHECKs
(`planets_chk_5` and `moons_chk_5`, the names MySQL gave them) are replaced by ones of the same names that list it, and a row stored as a
`surface ocean` with the `chloride brine` ocean class becomes a `brine seep`. No other column changes. The generator draws the state at the
generator epoch, so a row saved before it gets this label here and everything else from the regenerated galaxy.
"""

import sqlalchemy as sa
from alembic import op

revision = "0176"
down_revision = "0175"
branch_labels = None
depends_on = None

STATES = ("dry", "vapour", "ice", "ice-covered ocean", "surface ocean", "brine seep", "hycean")
CHECKS = (("planets", "planets_chk_5"), ("moons", "moons_chk_5"))


def upgrade():
    connection = op.get_bind()
    listed = ", ".join(f"'{state}'" for state in STATES)
    for table, name in CHECKS:
        clause = connection.execute(sa.text(
            "SELECT CHECK_CLAUSE FROM information_schema.CHECK_CONSTRAINTS WHERE CONSTRAINT_SCHEMA = DATABASE()"
            " AND CONSTRAINT_NAME = :name"), {"name": name}).scalar()
        if clause is None or "brine seep" not in clause:
            drop = f"DROP CHECK {name}, " if clause is not None else ""
            op.execute(f"ALTER TABLE {table} {drop}ADD CONSTRAINT {name} CHECK (hydrosphere IN ({listed}))")
        op.execute(f"UPDATE {table} SET hydrosphere = 'brine seep' WHERE hydrosphere = 'surface ocean' AND ocean_class = 'chloride brine'")
