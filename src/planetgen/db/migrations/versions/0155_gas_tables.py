"""The air of a planet or moon as child table rows (DB.71, GEN.424).

Schema v155. `planet_gases` and `moon_gases` hold one row per gas above zero: `gas`, its formula in the gas table, and `partial_pa`, its partial
pressure in Pa. The ten `p_<gas>_kpa` columns of `planets` and `moons` are copied into them (kPa times 1000, rows above zero only) and dropped. A
class's balance gas (helium, sodium, ammonia) was named in the text but never stored, so a row made before v155 has no row for it; the regenerated
galaxy of the generator epoch of this change stores it. Each step is skipped when it is already done, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0155"
down_revision = "0154"
branch_labels = None
depends_on = None

_GASES = (
    ("o2", "O2"), ("co2", "CO2"), ("co", "CO"), ("n2", "N2"), ("ar", "Ar"), ("h2", "H2"), ("h2o", "H2O"), ("ch4", "CH4"), ("h2s", "H2S"),
    ("so2", "SO2"),
)

_CREATE = """
CREATE TABLE IF NOT EXISTS {kind}_gases (
    id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    {owner}  BIGINT UNSIGNED NOT NULL,
    gas         VARCHAR(16) NOT NULL,
    partial_pa  DOUBLE NOT NULL CHECK (partial_pa > 0),

    CONSTRAINT fk_{kind}_gases_{kind}
        FOREIGN KEY ({owner}) REFERENCES {kind}s(id) ON DELETE CASCADE,
    UNIQUE KEY uq_{kind}_gases ({owner}, gas)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
"""


def upgrade():
    connection = op.get_bind()
    for kind in ("planet", "moon"):
        owner = f"{kind}_id"
        op.execute(_CREATE.format(kind=kind, owner=owner))
        for key, formula in _GASES:
            column = f"p_{key}_kpa"
            if helpers.exists_column(connection, f"{kind}s", column):
                connection.execute(sa.text(
                    f"INSERT IGNORE INTO {kind}_gases ({owner}, gas, partial_pa)"
                    f" SELECT id, :formula, {column} * 1000 FROM {kind}s WHERE {column} > 0"), {"formula": formula})
                connection.commit()
        gone = [f"p_{key}_kpa" for key, _ in _GASES if helpers.exists_column(connection, f"{kind}s", f"p_{key}_kpa")]
        if gone:
            # Not INSTANT: `planets` has a FULLTEXT index, and InnoDB refuses an instant change beside one (as in 0143).
            connection.execute(sa.text(f"ALTER TABLE {kind}s " + ", ".join(f"DROP COLUMN {column}" for column in gone)))
