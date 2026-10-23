"""A star can be a member of an N-star system (GEN.630).

Schema v133. A system of more than two stars keeps a row in `stars` for every star: the pair a system has always had
(`primary`, `secondary`) and each further one (`member`), so `stars.role` takes a fourth value. Every star row of a
system with a tree of orbits also says which of its parent node's two children it is, in `slot` (0 or 1, the numbering of
`orbit_nodes.slot`), which `parent_node_id` alone could not: a node whose two children are both stars needs it. The CHECK
(`stars_chk_1`, the name MySQL gave it) is replaced by one of the same name that lists all four roles, and a pair already
stored gets its slots (the primary 0, the secondary 1) from the node it hangs on. No other row changes.
"""

import sqlalchemy as sa
from alembic import op

revision = "0133"
down_revision = "0132"
branch_labels = None
depends_on = None

ROLES = ("primary", "secondary", "single", "member")
NAME = "stars_chk_1"


def upgrade():
    connection = op.get_bind()
    clause = connection.execute(sa.text(
        "SELECT CHECK_CLAUSE FROM information_schema.CHECK_CONSTRAINTS WHERE CONSTRAINT_SCHEMA = DATABASE()"
        " AND CONSTRAINT_NAME = :name"), {"name": NAME}).scalar()
    if clause is None or "member" not in clause:
        listed = ", ".join(f"'{name}'" for name in ROLES)
        drop = f"DROP CHECK {NAME}, " if clause is not None else ""
        op.execute(f"ALTER TABLE stars {drop}ADD CONSTRAINT {NAME} CHECK (role IN ({listed}))")
    has_slot = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'stars'"
        " AND COLUMN_NAME = 'slot'")).scalar()
    if not has_slot:
        op.execute("ALTER TABLE stars ADD COLUMN slot TINYINT UNSIGNED NULL AFTER parent_node_id")
    op.execute("UPDATE stars SET slot = CASE role WHEN 'primary' THEN 0 ELSE 1 END"
               " WHERE slot IS NULL AND parent_node_id IS NOT NULL AND role IN ('primary', 'secondary')")
