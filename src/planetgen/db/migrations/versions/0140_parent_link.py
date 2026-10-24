"""One parent link for every kind of object, and moons that can orbit a free planet (DB.73).

Schema v140. Every orbiter says what it orbits in the same two columns, `parent_type` and `parent_id`, which `planets` and
`moons` have had since v102: `star` (a `stars.id`), `planet` (a `planets.id`), `system` (the `star_systems.id`, for the
head of a system or a body orbiting the whole system), `node` (an `orbit_nodes.node_no` of the row's own system),
`rogue_planet` (a `rogue_planets.id`) or `free` (no parent, `parent_id` NULL). The column widens to VARCHAR(16) for the
longer names and its CHECK lists the six values. `stars`, `asteroid_belts`, `comets`, `black_holes`, `neutron_stars`,
`rogue_planets` and `interstellar_comets` gain both columns, backfilled from the columns that said it before: a star
hangs on its node (`parent_node_id`) or is the head of its system; a belt or comet orbits its star (`star_id`) or the
system; an anchored remnant orbits its star, any other remnant, rogue planet and interstellar comet is free. Every
table gets an index on `(parent_type, parent_id)`, so "what orbits this?" is one lookup for any kind.

`moons.planet_id` and `moons.star_system_id` stop being NOT NULL, so a moon can orbit a rogue planet: its parent is
`rogue_planet` and both columns are NULL. The new `moons.rogue_planet_id` (a foreign key, ON DELETE CASCADE) names the
rogue planet, so deleting it, or its sector, deletes its moons. A planet that leaves its star keeps its row and its ID and changes its
parent to `free`. Nothing else is touched.

Each step is skipped when it is already done, so a stopped run can be started again; the last column added
(`interstellar_comets.parent_id`) marks the revision done.
"""

import sqlalchemy as sa
from alembic import op

revision = "0140"
down_revision = "0139"
branch_labels = None
depends_on = None

PARENT_TYPES = ("star", "planet", "system", "node", "rogue_planet", "free")
"""tuple: What a parent link can name."""

OLD_PARENT_TYPES = ("star", "planet", "system")

# Table -> (default parent_type of a new row, backfill statements for the rows that exist).
NEW_LINKS = {
    "stars": ("system", (
        "UPDATE stars SET parent_type = 'node', parent_id = parent_node_id WHERE parent_node_id IS NOT NULL",
        "UPDATE stars SET parent_type = 'system', parent_id = star_system_id WHERE parent_node_id IS NULL")),
    "asteroid_belts": ("system", (
        "UPDATE asteroid_belts SET parent_type = 'star', parent_id = star_id WHERE star_id IS NOT NULL",
        "UPDATE asteroid_belts SET parent_type = 'system', parent_id = star_system_id WHERE star_id IS NULL")),
    "comets": ("system", (
        "UPDATE comets SET parent_type = 'star', parent_id = star_id WHERE star_id IS NOT NULL",
        "UPDATE comets SET parent_type = 'system', parent_id = star_system_id WHERE star_id IS NULL")),
    "black_holes": ("free", (
        "UPDATE black_holes SET parent_type = 'star', parent_id = star_id WHERE star_id IS NOT NULL",)),
    "neutron_stars": ("free", (
        "UPDATE neutron_stars SET parent_type = 'star', parent_id = star_id WHERE star_id IS NOT NULL",)),
    "rogue_planets": ("free", ()),
    "interstellar_comets": ("free", ()),
}
"""dict: The tables that gain the link, in the order they gain it (`interstellar_comets` last: it marks the revision)."""

WIDENED = ("planets", "moons")
"""tuple: The tables that already had the link."""


def _has_column(connection, table, column):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar())


def _has_index(connection, table, index):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.STATISTICS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND INDEX_NAME = :index_name"),
        {"table_name": table, "index_name": index}).scalar())


def _check_clause(connection, name):
    return connection.execute(sa.text(
        "SELECT CHECK_CLAUSE FROM information_schema.CHECK_CONSTRAINTS WHERE CONSTRAINT_SCHEMA = DATABASE()"
        " AND CONSTRAINT_NAME = :name"), {"name": name}).scalar()


def _set_check(connection, table):
    """Makes the table's `chk_<table>_parent_type` list `PARENT_TYPES`, dropping the old one first."""
    name = f"chk_{table}_parent_type"
    clause = _check_clause(connection, name)
    if clause is not None and "free" in clause:
        return
    drop = f"DROP CHECK {name}, " if clause is not None else ""
    listed = ", ".join(f"'{value}'" for value in PARENT_TYPES)
    op.execute(f"ALTER TABLE {table} {drop}ADD CONSTRAINT {name} CHECK (parent_type IN ({listed}))")


def _is_nullable(connection, table, column):
    return connection.execute(sa.text(
        "SELECT IS_NULLABLE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar() == "YES"


def _listed():
    return ", ".join(f"'{value}'" for value in PARENT_TYPES)


def upgrade():
    connection = op.get_bind()
    for table in WIDENED:
        default = "star" if table == "planets" else "planet"
        op.execute(f"ALTER TABLE {table} MODIFY COLUMN parent_type VARCHAR(16) NOT NULL DEFAULT '{default}'")
        _set_check(connection, table)
        if not _has_index(connection, table, f"idx_{table}_parent"):
            op.execute(f"ALTER TABLE {table} ADD KEY idx_{table}_parent (parent_type, parent_id)")
    for column in ("planet_id", "star_system_id"):
        if not _is_nullable(connection, "moons", column):
            op.execute(f"ALTER TABLE moons MODIFY COLUMN {column} BIGINT UNSIGNED NULL")
    if not _has_column(connection, "moons", "rogue_planet_id"):
        op.execute("ALTER TABLE moons ADD COLUMN rogue_planet_id BIGINT UNSIGNED NULL, "
                   "ADD KEY idx_moons_rogue_planet_id (rogue_planet_id), "
                   "ADD CONSTRAINT fk_moons_rogue_planet FOREIGN KEY (rogue_planet_id) "
                   "REFERENCES rogue_planets(id) ON DELETE CASCADE")
    for table, (default, backfill) in NEW_LINKS.items():
        if not _has_column(connection, table, "parent_id"):
            # One statement, so a stopped run never leaves a table with half the link.
            op.execute(f"ALTER TABLE {table} ADD COLUMN parent_type VARCHAR(16) NOT NULL DEFAULT '{default}', "
                       f"ADD COLUMN parent_id BIGINT UNSIGNED NULL, ADD KEY idx_{table}_parent (parent_type, parent_id), "
                       f"ADD CONSTRAINT chk_{table}_parent_type CHECK (parent_type IN ({_listed()}))")
        for statement in backfill:
            op.execute(statement)
