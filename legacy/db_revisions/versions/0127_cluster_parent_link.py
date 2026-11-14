"""Unnamed clusters and the star-to-cluster parent link (DB.61).

Schema v127. `galaxy_components.name` becomes nullable, so a cluster that has no name of its own (a generated open cluster
nobody has named) can be a row; the unique key over `(kind, name)` stays, and MySQL lets any number of rows hold a NULL name,
so it still holds over the named rows. `stars` and `bright_stars` each gain `cluster_id` (BIGINT UNSIGNED, NULL), a
reference to `galaxy_components.id` with its own index (`idx_stars_cluster`, `idx_bright_stars_cluster`), so a star can point
to the cluster that is its parent; the link is cleared if the cluster row is removed. Every existing cluster keeps its name,
every star's link is NULL, and nothing writes or reads the link yet.

Each table is changed in one `ALGORITHM=COPY` statement: a foreign key cannot be added in place with checks on, and `stars`
has a FULLTEXT index (`ft_stars_name`), which InnoDB refuses to change instantly or in place (error 1846). So the two star
tables are rebuilt once, which takes time in proportion to their rows. A step already done is left alone, so a stopped run can
be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0127"
down_revision = "0126"
branch_labels = None
depends_on = None

TABLES = ("stars", "bright_stars")
"""tuple: the tables that gain `cluster_id`."""


def _is_nullable(connection, table, column):
    return connection.execute(sa.text(
        "SELECT IS_NULLABLE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :t AND COLUMN_NAME = :c"), {"t": table, "c": column}).scalar() == "YES"


def upgrade():
    connection = op.get_bind()
    if not _is_nullable(connection, "galaxy_components", "name"):
        connection.execute(sa.text("ALTER TABLE galaxy_components MODIFY name VARCHAR(255) NULL"))
    for table in TABLES:
        if not helpers.exists_column(connection, table, "cluster_id"):
            connection.execute(sa.text(
                f"ALTER TABLE {table} ADD COLUMN cluster_id BIGINT UNSIGNED NULL,"
                f" ADD KEY idx_{table}_cluster (cluster_id),"
                f" ADD CONSTRAINT fk_{table}_cluster FOREIGN KEY (cluster_id) REFERENCES galaxy_components(id)"
                " ON DELETE SET NULL, ALGORITHM=COPY"))
