"""Per-class totals of the phenomenon scatter, so the Phenomena table needn't scan it.

Schema v77. Adds `phenomenon_scatter_classes` (see `schema.sql`'s note) and
fills it from the scatter rows already there.
"""

from alembic import op

revision = "0077"
down_revision = "0076"
branch_labels = None
depends_on = None


def upgrade():
    op.execute(
        """
        CREATE TABLE IF NOT EXISTS phenomenon_scatter_classes (
            kind     VARCHAR(24) NOT NULL,
            subtype  VARCHAR(16) NOT NULL DEFAULT '',
            placed   BIGINT UNSIGNED NOT NULL DEFAULT 0,
            built    BIGINT UNSIGNED NOT NULL DEFAULT 0,
            PRIMARY KEY (kind, subtype)
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
        """
    )
    op.execute(
        "INSERT IGNORE INTO phenomenon_scatter_classes (kind, subtype, placed, built)"
        " SELECT kind, COALESCE(subtype, ''), COUNT(*), SUM(built_at IS NOT NULL)"
        " FROM phenomenon_scatter GROUP BY kind, COALESCE(subtype, '')"
    )


def downgrade():
    op.execute("DROP TABLE IF EXISTS phenomenon_scatter_classes")
