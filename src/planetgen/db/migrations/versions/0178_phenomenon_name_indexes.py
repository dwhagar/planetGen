"""The Phenomena table pages by name in index order (PERF.205).

Schema v178. No column changes and no stored value changes: `asteroid_fields` and `interstellar_comets` get a `name` index like the other
six phenomenon tables already have (`idx_asteroid_fields_name`, `idx_interstellar_comets_name`), built online (`LOCK=NONE`); a stopped run can
be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0178"
down_revision = "0177"
branch_labels = None
depends_on = None

NEW = (("asteroid_fields", "idx_asteroid_fields_name"), ("interstellar_comets", "idx_interstellar_comets_name"))


def upgrade():
    connection = op.get_bind()
    for table, name in NEW:
        helpers.add_index_online(connection, table, name, ("name",))
