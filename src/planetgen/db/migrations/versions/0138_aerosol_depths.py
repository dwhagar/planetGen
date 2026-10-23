"""The haze, dust and cloud optical depths of a planet or moon (GEN.414, DB.68).

Schema v138. `planets` and `moons` gain `tau_haze`, `tau_dust`, `tau_cloud` (the optical depths at 550 nm), `cloud_cover` (the share of the sky under cloud, 0 to 1) and
`aerosol_draw` (the number in [0, 1) that picks the optics inside their brackets), all DOUBLE and nullable. The generator draws them from the body's own stream
(`physics/aerosols.py`), so they cannot be made from a stored row by SQL: a row saved before this revision keeps NULL, which `aerosols.of_body` reads as no aerosol,
and the regenerated galaxy (the generator epoch of this change) draws them. Each column is skipped when it is there, so a stopped run can be started again; the version
is detected from `planets.tau_cloud`.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0138"
down_revision = "0137"
branch_labels = None
depends_on = None

_COLUMNS = ("tau_haze", "tau_dust", "tau_cloud", "cloud_cover", "aerosol_draw")


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        for column in _COLUMNS:
            if not helpers.exists_column(connection, table, column):
                op.execute(f"ALTER TABLE {table} ADD COLUMN {column} DOUBLE NULL")
