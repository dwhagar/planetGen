"""Write the chart ledger row of every sector that was charted and has none (OPS.180).

Schema v247, no schema change. The Sectors list, the Galaxy page and the Generate page's "Galaxy at a glance" count a placed sector as charted when the chart ledger holds its cell (GEN.1127); a placed row whose
cell it does not hold reads as an uncharted placeholder. A sector that holds a star system or whose fill was recorded (`sector_stats.bright_level_sol` 0) is charted, so this writes the ledger row it lacks
(`db.chart_ledger.repair`) and says what it found: how many placed sectors there are, how many have a row, how many held content with none. A second run finds nothing.
"""

from planetgen.db.migrations import helpers

revision = "0253"
down_revision = "0252"
branch_labels = None
depends_on = None


def upgrade():
    from alembic import op
    from planetgen.db import chart_ledger
    from planetgen.util import log

    found = chart_ledger.repair(helpers.StoreConnection(op.get_bind()))
    log.normal(f"Chart ledger: {found['placed']:,} placed sectors, {found['charted']:,} with a ledger row ({found['ledger']:,} cells charted), "
               f"{found['missing']:,} held content with no row; {found['repaired']:,} rows written.")
