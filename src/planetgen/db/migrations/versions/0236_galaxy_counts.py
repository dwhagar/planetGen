"""The galaxy's counts, kept as they change (DB.121).

Schema v236. `sector_facets` holds, for each generated sector, its `total` counts (systems, stars, planets, moons, belts, habitable worlds) and the count of each value of each
facet the lists and the search page offer, keyed by the sector's address; `galaxy_totals` is the one row that is their sum plus the mass and luminosity of `sector_stats`
(`db.galaxy_counts`). Both are new tables: nothing is added to `star_systems`, whose FULLTEXT key makes any added column copy the whole table.

The rows are made from the stored bodies at once (`db.galaxy_counts.rebuild`), so a database that already holds sectors has every figure and its pages read them at the first
look. A table already there is skipped and the rebuild empties both before it fills them, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0236"
down_revision = "0235"
branch_labels = None
depends_on = None

_TABLES = (
    """CREATE TABLE IF NOT EXISTS sector_facets (
        shell_index   INT NOT NULL,
        band_index    INT NOT NULL,
        wedge_index   INT NOT NULL,
        facet         VARCHAR(16) NOT NULL,
        value         VARCHAR(40) NOT NULL,
        n             INT UNSIGNED NOT NULL,
        PRIMARY KEY (shell_index, band_index, wedge_index, facet, value),
        KEY idx_sector_facets_facet (facet, value, n)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""",
    """CREATE TABLE IF NOT EXISTS galaxy_totals (
        id                 TINYINT UNSIGNED NOT NULL PRIMARY KEY,
        systems            BIGINT NOT NULL DEFAULT 0,
        stars              BIGINT NOT NULL DEFAULT 0,
        planets            BIGINT NOT NULL DEFAULT 0,
        moons              BIGINT NOT NULL DEFAULT 0,
        belts              BIGINT NOT NULL DEFAULT 0,
        habitable_worlds   BIGINT NOT NULL DEFAULT 0,
        mass_kg            DOUBLE NOT NULL DEFAULT 0,
        luminosity_sol     DOUBLE NOT NULL DEFAULT 0,
        updated_at         TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""",
)


def upgrade():
    import sqlalchemy as sa

    connection = op.get_bind()
    for statement in _TABLES:
        connection.execute(sa.text(statement))
    from planetgen.db import galaxy_counts

    galaxy_counts.rebuild(helpers.StoreConnection(connection))
