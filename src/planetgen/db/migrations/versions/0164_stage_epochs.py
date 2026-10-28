"""The generator epoch each scatter and each species was made under (GEN.764).

Schema v164. A galaxy-wide scatter has no rows to carry an epoch, so `scatter_epochs` holds one row per kind of scatter drawn (`mass`, `luminosity`,
`band`, `phenomena`, `nebula-field`): the generator epoch and the run that drew it, and the row goes when the scatter is cleared. `species` gains
`generator_epoch` (SMALLINT UNSIGNED, NULL: unknown), the epoch the population pass made it under. A sector already has its own (v121) and the plan
its own (v131).

A galaxy made before v164 is given what it can be given: a kind of scatter that has rows gets the plan's epoch (the scatter is drawn right after
the plan; NULL when the plan has none), and a species gets the epoch of the sector its homeworld's system is in. Each step is skipped when it is
already done, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0164"
down_revision = "0163"
branch_labels = None
depends_on = None

# kind -> a statement that is true when the galaxy holds rows of that kind
_HELD = {
    "mass": "EXISTS (SELECT 1 FROM bright_stars WHERE initial_mass_sol >= g.bright_star_mass_limit_sol)",
    "luminosity": "EXISTS (SELECT 1 FROM bright_stars WHERE g.bright_star_mass_limit_sol IS NULL OR initial_mass_sol < g.bright_star_mass_limit_sol)",
    "phenomena": "EXISTS (SELECT 1 FROM phenomenon_scatter)",
    "nebula-field": "EXISTS (SELECT 1 FROM nebula_field)",
}

_CREATE = """
CREATE TABLE IF NOT EXISTS scatter_epochs (
    kind             VARCHAR(16) NOT NULL PRIMARY KEY,
    generator_epoch  SMALLINT UNSIGNED NULL,
    run_id           BIGINT UNSIGNED NULL,
    recorded_at      TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),

    CONSTRAINT chk_scatter_epochs_kind CHECK (kind IN ('mass', 'luminosity', 'band', 'phenomena', 'nebula-field')),
    CONSTRAINT fk_scatter_epochs_run FOREIGN KEY (run_id) REFERENCES generation_runs(id) ON DELETE SET NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
"""


def upgrade():
    connection = op.get_bind()
    op.execute(_CREATE)
    for kind, held in _HELD.items():
        connection.execute(sa.text(
            "INSERT IGNORE INTO scatter_epochs (kind, generator_epoch) "
            f"SELECT :kind, g.generator_epoch FROM galaxy_shape g WHERE g.id = 1 AND {held}"), {"kind": kind})
    connection.commit()
    helpers.add_column_instant(connection, "species", "generator_epoch", "SMALLINT UNSIGNED NULL", after="tech_source")
    helpers.batched_update(
        connection, "species", "generator_epoch = (SELECT s.generator_epoch FROM star_systems ss JOIN sectors s ON s.id = ss.sector_id"
        " WHERE ss.id = species.star_system_id)",
        "generator_epoch IS NULL AND EXISTS (SELECT 1 FROM star_systems ss JOIN sectors s ON s.id = ss.sector_id"
        " WHERE ss.id = species.star_system_id AND s.generator_epoch IS NOT NULL)", label="Species epochs")
