"""The members of a system's tree that are not stars (GEN.1013).

Schema v194. `dark_members` holds one row for each brown dwarf, neutron star, black hole or stripped helium star that hangs in a built
system's orbit tree without being one of its stars: where it hangs (`parent_node_no`, `slot`, the numbering of `orbit_nodes`), its
`kind`, its mass in solar masses and its radius in AU. A node's mass is the sum of what hangs below it, these included, so the
orbits and the gravity are the tree the system was drawn with. No row is made for a system already stored; a system gets its rows when
it is generated again. A table already there is kept, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0194"
down_revision = "0193"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    if helpers.exists_table(connection, "dark_members"):
        return
    op.execute("""
        CREATE TABLE dark_members (
            star_system_id  BIGINT UNSIGNED NOT NULL,
            member_no       SMALLINT NOT NULL,
            parent_node_no  SMALLINT NOT NULL,
            slot            SMALLINT NOT NULL,
            kind            VARCHAR(16) NOT NULL,
            mass_solar      DOUBLE NOT NULL,
            radius_au       DOUBLE NULL,
            PRIMARY KEY (star_system_id, member_no),
            CONSTRAINT fk_dark_members_star_system
                FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
            CONSTRAINT chk_dark_members_kind CHECK (kind IN ('brown-dwarf', 'neutron-star', 'black-hole', 'stripped-helium'))
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""")
