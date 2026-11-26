"""The pair's orbit is its node's: the mutual-orbit columns of `star_systems` are gone (GEN.640 slice d).

Schema v280 (a placeholder: the train that merges this renumbers it to the next free revision). A system's own pair of stars has its orbit in its node of `orbit_nodes`, where
`orbit_trees.pair_places_km` and `generation.pair_orbit.PairOrbit` read it, and `store.advance_orbital_phases` turns that node (GEN.313, GEN.640 slice a); the close pair's host
(`generation.pair_host.PairHost`, slice c) keeps no copy of it either. The copy `star_systems` kept has no writer and no reader, so it is dropped: `binary_eccentricity`,
`binary_periapsis_km`, `binary_apoapsis_km`, `binary_mutual_orbital_period_years`, `binary_mutual_orbital_speed_kms`, `binary_mutual_orbital_inclination_deg`,
`binary_mutual_orbital_ascending_node_deg`, `binary_mutual_orbital_phase_deg`, `binary_epoch_unix`, `binary_next_update_due` (and its index), `binary_mutual_min_update_interval_years`,
`binary_mutual_position_x_km`/`y`/`z`, `binary_primary_position_x_km`/`y`/`z`, `binary_secondary_position_x_km`/`y`/`z` and `binary_secondary_mass_fraction`. `binary_separation_km` stays: the
sector lists, the maps and the facilities read the pair's axis from the row without the tree. A column already gone is skipped, so a stopped run can be started again. Nothing is
backfilled and nothing is lost: the node held the orbit all along.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0280"
down_revision = "0279"
branch_labels = None
depends_on = None

INDEX = "idx_star_systems_binary_next_update_due"
COLUMNS = (
    "binary_eccentricity", "binary_periapsis_km", "binary_apoapsis_km",
    "binary_mutual_orbital_period_years", "binary_mutual_orbital_speed_kms", "binary_mutual_orbital_inclination_deg",
    "binary_mutual_orbital_ascending_node_deg", "binary_mutual_orbital_phase_deg", "binary_epoch_unix", "binary_next_update_due",
    "binary_mutual_min_update_interval_years", "binary_mutual_position_x_km", "binary_mutual_position_y_km", "binary_mutual_position_z_km",
    "binary_primary_position_x_km", "binary_primary_position_y_km", "binary_primary_position_z_km",
    "binary_secondary_position_x_km", "binary_secondary_position_y_km", "binary_secondary_position_z_km",
    "binary_secondary_mass_fraction",
)
"""tuple: The columns of `star_systems` that held the pair's orbit, in the order of the table."""


def upgrade():
    connection = op.get_bind()
    if helpers.exists_index(connection, "star_systems", INDEX):
        helpers.drop_index_online(connection, "star_systems", INDEX)
    present = [column for column in COLUMNS if helpers.exists_column(connection, "star_systems", column)]
    if not present:
        return
    drops = ", ".join(f"DROP COLUMN {column}" for column in present)      # one statement: the columns together take one row version of the instant form
    try:
        helpers._with_retries(connection, f"ALTER TABLE star_systems {drops}, ALGORITHM=INSTANT")
    except Exception:  # a column check on a fresh schema can refuse the instant form; the plain drop takes the columns with it
        connection.execute(sa.text(f"ALTER TABLE star_systems {drops}"))


def downgrade():
    connection = op.get_bind()
    for column in COLUMNS:
        if not helpers.exists_column(connection, "star_systems", column):
            helpers.add_column_instant(connection, "star_systems", column, "DOUBLE NULL")
    if not helpers.exists_index(connection, "star_systems", INDEX):
        helpers.add_index_online(connection, "star_systems", INDEX, ("binary_next_update_due",))
