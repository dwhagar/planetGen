"""The UV-C share of a star's blackbody (GEN.443).

Schema v145. `stars` gains `uv_c_share` (DOUBLE, nullable), the share of its blackbody's 200 to 280 nm light the star gives off: the Sun's
measured 0.256 is one value of the open bracket 0.1 to 1 of the radiation model (R17), and the generator draws it once per star from the
star's own stream (`physics/activity.py`), log-uniform. A remnant has none (NULL). The stored value cannot be made from a row by SQL, so a
row saved before this revision keeps NULL, which the radiation model reads as the blackbody; the regenerated galaxy (the generator epoch of this
change) draws it. A column already there is skipped, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0145"
down_revision = "0144"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_column(connection, "stars", "uv_c_share"):
        op.execute("ALTER TABLE stars ADD COLUMN uv_c_share DOUBLE NULL")
