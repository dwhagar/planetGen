"""The interior of every giant planet and moon (DB.89, GEN.642).

Schema v154. `planets` and `moons` gain the summary of the interior a giant is drawn with: `interior_kind` ('gas_giant', 'ice_giant' or
'brown_dwarf'), the compact core's and the dilute core's mass (`core_mass_earth`, `dilute_core_mass_earth`), the pressure at which the
hydrogen turns metallic (`metallization_pressure_gpa`), the molecular envelope's helium share against the natal share (`helium_depletion`)
and whether a helium-rich shell sits above the metallic envelope (`has_helium_shell`). The layers themselves, outside in, are rows of the
new `planet_interior_layers` and `moon_interior_layers`: the kind and the state, the pressure and the radius fraction at the top and the
bottom, and the mass of each material in Earth masses.

All NULL here, and no layer rows: the interior needs the draws of the generator epoch, so a row made before v150 has none and the
regenerated galaxy replaces it. Each column and table is skipped when it is already there, so a stopped run can be started again; the last
column marks the migration done.
"""

import sqlalchemy as sa
from alembic import op

revision = "0154"
down_revision = "0153"
branch_labels = None
depends_on = None

_INTERIOR_KINDS = ("gas_giant", "ice_giant", "brown_dwarf")
_KINDS = ", ".join(f"'{kind}'" for kind in _INTERIOR_KINDS)
_COLUMNS = (
    ("interior_kind", f"VARCHAR(16) NULL CHECK (interior_kind IN ({_KINDS}))"),
    ("core_mass_earth", "DOUBLE NULL"),
    ("dilute_core_mass_earth", "DOUBLE NULL"),
    ("metallization_pressure_gpa", "DOUBLE NULL"),
    ("helium_depletion", "DOUBLE NULL"),
    ("has_helium_shell", "BOOLEAN NULL"),
)
_LAYER_KINDS = ", ".join(f"'{kind}'" for kind in (
    "molecular_envelope", "metallic_envelope", "helium_gradient", "helium_shell", "dilute_core", "compact_core", "ice_mantle",
    "rock_core"))
_STATES = ", ".join(f"'{state}'" for state in ("molecular_fluid", "metallic_fluid", "ionic_fluid", "superionic", "molten"))
_MATERIALS = ("hydrogen", "helium", "iron", "rock", "water_liquid", "ice_x", "nh3_ice", "ch4_ice")

_LAYER_TABLE = """
    CREATE TABLE {table} (
        id                      BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
        {owner}_id              BIGINT UNSIGNED NOT NULL,
        layer_index             SMALLINT UNSIGNED NOT NULL,
        kind                    VARCHAR(24) NOT NULL CHECK (kind IN ({kinds})),
        state                   VARCHAR(16) NOT NULL CHECK (state IN ({states})),
        top_pressure_gpa        DOUBLE NOT NULL,
        bottom_pressure_gpa     DOUBLE NOT NULL,
        mass_earth              DOUBLE NOT NULL,
        top_radius_fraction     DOUBLE NOT NULL,
        bottom_radius_fraction  DOUBLE NOT NULL,
        {materials},
        UNIQUE KEY uq_{table}_index ({owner}_id, layer_index),
        CONSTRAINT fk_{table}_{owner}
            FOREIGN KEY ({owner}_id) REFERENCES {parent}(id) ON DELETE CASCADE
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def _has(connection, table, column=None):
    sql = ("SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = :table_name"
           + (" AND COLUMN_NAME = :column_name" if column else ""))
    return connection.execute(sa.text(sql), {"table_name": table, "column_name": column}).scalar()


def upgrade():
    connection = op.get_bind()
    for owner, parent in (("planet", "planets"), ("moon", "moons")):
        table = f"{owner}_interior_layers"
        if not _has(connection, table):
            op.execute(_LAYER_TABLE.format(
                table=table, owner=owner, parent=parent, kinds=_LAYER_KINDS, states=_STATES,
                materials=",\n        ".join(f"{material} DOUBLE NOT NULL" for material in _MATERIALS)))
    for table in ("planets", "moons"):
        for column, definition in _COLUMNS:
            if not _has(connection, table, column):
                op.execute(f"ALTER TABLE {table} ADD COLUMN {column} {definition}")
