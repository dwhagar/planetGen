"""The stored composition of every nebula and supernova remnant (GEN.1276).

Schema v214. `nebulae` and `supernova_remnants` gain the columns of `physics/nebula_columns.py`: the gas's [O/H], the hydrogen, helium and metal mass fractions and those of C, N, O, Ne, S, Ar
and Fe, their gas-phase abundances (12 + log(X/H)), the ionized, atomic and molecular shares of the hydrogen, the dust-to-gas and dust-to-metals ratios, the share of the dust in each
grain kind and the progenitor's mass; a remnant adds the mass of its ejecta and of its swept-up gas and their hydrogen, helium and metal fractions. The rows already there are filled by the same
model the generator runs (`nebula_columns.values_for`, keyed on the object's ID (a field cloud, class M or N, on its place), so the numbers are the ones a regeneration of the seed gives), and their `composition` and `dominant_species` are
written again from the numbers. A column already there is skipped, and only a row with no composition is filled, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner
from planetgen.db.migrations import helpers

revision = "0214"
down_revision = "0213"
branch_labels = None
depends_on = None

_COMMON = (
    "gas_oh_dex", "mass_fraction_h", "mass_fraction_he", "mass_fraction_z", "mass_fraction_c", "mass_fraction_n", "mass_fraction_o", "mass_fraction_ne", "mass_fraction_s", "mass_fraction_ar",
    "mass_fraction_fe", "abundance_c", "abundance_n", "abundance_o", "abundance_ne", "abundance_s", "abundance_ar", "abundance_fe", "ionized_share", "atomic_share", "molecular_share",
    "dust_to_gas", "dust_to_metals", "grain_amorphous_olivine", "grain_amorphous_pyroxene", "grain_crystalline_silicate", "grain_iron", "grain_amorphous_carbon", "grain_graphite", "grain_pah", "grain_sic",
    "grain_ice", "progenitor_mass_sol",
)
_REMNANT = (
    "ejecta_mass_sol", "swept_mass_sol", "ejecta_mass_fraction_h", "ejecta_mass_fraction_he", "ejecta_mass_fraction_z", "swept_mass_fraction_h", "swept_mass_fraction_he", "swept_mass_fraction_z",
)
_COLUMNS = {"nebulae": _COMMON, "supernova_remnants": _COMMON + _REMNANT}
"""dict: Table -> the DOUBLE columns it gains."""

_READ = {
    "nebulae": "id, name, uid, nebula_class, density_cm3, radius_ly, galactic_radius_pc, center_x_pc, center_y_pc, center_z_pc, hosts",
    "supernova_remnants": "id, name, uid, remnant_class, progenitor_type, compact_remnant_kind, density_cm3, radius_ly, galactic_radius_pc",
}
BATCH = 500


def upgrade():
    connection = op.get_bind()
    for table, columns in _COLUMNS.items():
        for column in columns:
            helpers.add_column_instant(connection, table, column, "DOUBLE NULL")
    _fill(connection)


def _fill(connection):
    """Fills the composition of every row that has none, a batch of rows to a commit, and marks the sectors that hold them for the parity writer."""
    from planetgen.db import store
    from planetgen.physics import nebula_columns
    skeleton = store.get_galaxy_shape(helpers.StoreConnection(connection))
    seed, shape = (None, None) if skeleton is None else (skeleton.galaxy_seed, skeleton.shape)
    for table, columns in _COLUMNS.items():
        total = connection.execute(sa.text(f"SELECT COUNT(*) FROM {table} WHERE mass_fraction_h IS NULL")).scalar()
        done, last = 0, -1
        assignments = ", ".join(f"{column} = :{column}" for column in nebula_columns.COLUMNS[table]) + ", composition = :composition, dominant_species = :dominant_species"
        if table == "supernova_remnants":
            assignments = ", ".join(f"{column} = :{column}" for column in nebula_columns.COLUMNS[table]) + ", dominant_species = :dominant_species"
        while True:
            rows = connection.execute(sa.text(f"SELECT {_READ[table]} FROM {table} WHERE id > :last AND mass_fraction_h IS NULL ORDER BY id LIMIT :n"), {"last": last, "n": BATCH}).mappings().all()
            if not rows:
                break
            batch = []
            for row in rows:
                with nebula_columns.in_galaxy(seed, shape):
                    values = nebula_columns.values_for(table, dict(row))
                sentence, species = nebula_columns.sentences_for(values)
                batch.append({**values, "id": row["id"], "dominant_species": species[:255], **({"composition": sentence} if table == "nebulae" else {})})
            connection.execute(sa.text(f"UPDATE {table} SET {assignments} WHERE id = :id"), batch)
            connection.commit()
            last = rows[-1]["id"]
            done += len(rows)
            alembic_runner.report_progress(f"composition of {table}", done, total)
        connection.execute(sa.text(f"UPDATE sectors SET parity_dirty = 1 WHERE id IN (SELECT DISTINCT sector_id FROM {table} WHERE sector_id IS NOT NULL)"))
        connection.commit()
