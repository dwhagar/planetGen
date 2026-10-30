"""The Species table's sorts read in index order (PERF.205).

Schema v177. No column changes and no stored value changes: `species` gets `idx_species_spacefaring_name (spacefaring, name)`,
`idx_species_tech_level_name (tech_level, name)` and `idx_species_era_name (era, name)`, and loses `idx_species_spacefaring` and
`idx_species_tech_level`, which are the first column of the new ones (the tech level search by range keeps using it). Online builds
(`LOCK=NONE`); a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0177"
down_revision = "0176"
branch_labels = None
depends_on = None

NEW = (
    ("idx_species_spacefaring_name", ("spacefaring", "name")),
    ("idx_species_tech_level_name", ("tech_level", "name")),
    ("idx_species_era_name", ("era", "name")),
)
OLD = ("idx_species_spacefaring", "idx_species_tech_level")


def upgrade():
    connection = op.get_bind()
    for name, columns in NEW:
        helpers.add_index_online(connection, "species", name, columns)
    for name in OLD:
        helpers.drop_index_online(connection, "species", name)
