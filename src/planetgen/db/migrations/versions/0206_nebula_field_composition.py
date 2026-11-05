"""The stored composition of the field clouds (GEN.1277).

Schema v206. `nebula_field` gains the columns of `nebula_columns.FIELD_COLUMNS` (the composition columns of `nebulae` without the progenitor mass), FLOAT. The rows already there are filled
from the galaxy's seed and shape: each cell is drawn again (`nebula_field.cell_draw`), which gives the cloud its own density, and `nebula_columns.field_values` works the composition out from it,
keyed on the cloud's cell and place (`field_key`), so a row reads the same as the `nebulae` row the cloud becomes when a sector stores it (GEN.1276 keys a class M or N nebula the same way). A
cloud the draw does not give back (a table planned under another seed) is filled from the density `nebula_composition.field_cloud_density_cm3` draws for its key. A column already there is
skipped, and only a row with no composition is filled, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner
from planetgen.db.migrations import helpers

revision = "0206"
down_revision = "0205"
branch_labels = None
depends_on = None

BATCH = 2000
CENTRE_MATCH_PC = 1e-6
"""float: How near a drawn cloud's centre must lie to the row's for the draw to be the row's cloud."""


def upgrade():
    from planetgen.physics import nebula_columns
    connection = op.get_bind()
    for column in nebula_columns.FIELD_COLUMNS:
        helpers.add_column_instant(connection, "nebula_field", column, "FLOAT NULL")
    _fill(connection, nebula_columns)


def _fill(connection, nebula_columns):
    """Fills the composition of every `nebula_field` row that has none, a batch of rows to a commit."""
    from planetgen.db import store
    from planetgen.galaxy import nebula_field
    from planetgen.physics import nebula_composition

    columns = nebula_columns.FIELD_COLUMNS
    total = connection.execute(sa.text("SELECT COUNT(*) FROM nebula_field WHERE mass_fraction_h IS NULL")).scalar()
    if not total:
        return
    skeleton = store.get_galaxy_shape(helpers.StoreConnection(connection))
    seed = None if skeleton is None else skeleton.galaxy_seed
    shape = None if skeleton is None else skeleton.shape
    statement = sa.text("UPDATE nebula_field SET " + ", ".join(f"{column} = :{column}" for column in columns)
                        + " WHERE cell_i = :cell_i AND cell_j = :cell_j AND cell_k = :cell_k AND n = :n")
    read = sa.text("SELECT cell_i, cell_j, cell_k, n, nebula_class, radius_ly, center_x_pc, center_y_pc, center_z_pc FROM nebula_field"
                   " WHERE mass_fraction_h IS NULL AND (cell_i, cell_j, cell_k, n) > (:i, :j, :k, :n) ORDER BY cell_i, cell_j, cell_k, n LIMIT :limit")
    last = (-2 ** 31, -2 ** 31, -2 ** 31, -1)
    done = 0
    while True:
        rows = connection.execute(read, {"i": last[0], "j": last[1], "k": last[2], "n": last[3], "limit": BATCH}).mappings().all()
        if not rows:
            break
        drawn = {}
        batch = []
        for row in rows:
            index = (row["cell_i"], row["cell_j"], row["cell_k"])
            center = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
            density = None
            if shape is not None:
                if index not in drawn:
                    drawn[index] = {n: (nebula, where) for n, nebula, where in nebula_field.cell_draw(seed, shape, index)}
                found = drawn[index].get(row["n"])
                if found is not None and found[0].nebula_class == row["nebula_class"] and all(abs(a - b) < CENTRE_MATCH_PC for a, b in zip(found[1], center)):
                    density = found[0].density_cm3
            if density is None:
                density = nebula_composition.field_cloud_density_cm3(row["nebula_class"], nebula_columns.field_key(center))
            with nebula_columns.in_galaxy(seed, shape):
                values = nebula_columns.field_values(row["nebula_class"], density, row["radius_ly"], center)
            batch.append({"cell_i": row["cell_i"], "cell_j": row["cell_j"], "cell_k": row["cell_k"], "n": row["n"], **dict(zip(columns, values))})
        connection.execute(statement, batch)
        connection.commit()
        last = (rows[-1]["cell_i"], rows[-1]["cell_j"], rows[-1]["cell_k"], rows[-1]["n"])
        done += len(rows)
        alembic_runner.report_progress("composition of the field clouds", done, total)
