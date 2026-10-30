"""The non-water ices of a body, one amount each (GEN.834).

Schema v182. `planets` and `moons` gain `ice_co2_mass_fraction`, `ice_ch4_mass_fraction`, `ice_nh3_mass_fraction`, `ice_n2_mass_fraction`,
`ice_so2_mass_fraction`, `ice_co_mass_fraction`, `ice_hcn_mass_fraction` and `ice_other_mass_fraction`: the shares of the whole body's mass
that each ice holds, adding up to the existing `ice_mass_fraction` (`physics/ice_species.py`).

All NULL here: the species are drawn with the body's bulk composition, so a row made before v182 gets them when the galaxy is regenerated.
Each column is skipped when it is already there, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0182"
down_revision = "0181"
branch_labels = None
depends_on = None

COLUMNS = (
    "ice_co2_mass_fraction", "ice_ch4_mass_fraction", "ice_nh3_mass_fraction", "ice_n2_mass_fraction", "ice_so2_mass_fraction",
    "ice_co_mass_fraction", "ice_hcn_mass_fraction", "ice_other_mass_fraction",
)
"""tuple: The DOUBLE columns added to `planets` and `moons`."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        for column in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                op.execute(f"ALTER TABLE {table} ADD COLUMN {column} DOUBLE NULL")
