"""A planet or moon keeps the Bond albedo and the greenhouse quantity its surface temperature is made of (GEN.412).

Schema v142. `planets` and `moons` gain `bond_albedo` (DOUBLE, NULL), `greenhouse_factor` (DOUBLE, NULL: G in
T_s = T_e (1 + G)^(1/4), 0 for a body with no air) and `albedo_basis` (TINYINT UNSIGNED, NULL: 0 recorded at generation,
1 the exact inverse of a stored row, 2 a seeded class draw). New rows carry the two numbers the generator drew; a row saved
before the revision gets them from `planetgen.physics.climate_draws` (docs/design/composition/atmosphere-colour.md section 8): its stored
scale height holds the temperature the body would have with no air, so both draws come back exactly. Where that has no
solution (no starlight, a temperature on the background floor, an albedo outside 0 to 1, a pressure above the barometric
law) the albedo is drawn from the class range, seeded by the table and id, and the greenhouse quantity is worked out from the
stored temperature, which never moves.

The three columns go in with one ALTER for each table, and rows are filled in pages of ids, so a stopped run starts again where
it stopped (a row with `albedo_basis` set is done). A row with no surface temperature or no starlight to work from keeps NULL.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner
from planetgen.db.migrations import helpers
from planetgen.physics import climate_draws

revision = "0142"
down_revision = "0141"
branch_labels = None
depends_on = None

COLUMNS = (("bond_albedo", "DOUBLE NULL"), ("greenhouse_factor", "DOUBLE NULL"), ("albedo_basis", "TINYINT UNSIGNED NULL"))
"""tuple: `(column, definition)` added to `planets` and `moons`, in table order."""

PAGE = 5000
"""int: Rows read, worked out and written before a commit."""

_SELECT = """
    SELECT b.id, b.planet_class, b.body_type, b.atmosphere, b.surface_temperature_k, b.scale_height_km,
           b.atm_molar_density, b.gravity_g, b.atm_density, b.atmospheric_pressure_pa,
           {orbit} AS orbit_km, COALESCE(st.luminosity_w, ss.binary_effective_luminosity_w) AS luminosity_w
    FROM {table} b {joins}
    WHERE b.id > :last AND b.albedo_basis IS NULL AND b.bond_albedo IS NULL
    ORDER BY b.id LIMIT :n"""

_JOINS = {
    "planets": ("b.distance_km",
                "LEFT JOIN stars st ON st.id = b.star_id LEFT JOIN star_systems ss ON ss.id = b.star_system_id"),
    "moons": ("pl.distance_km",
              "JOIN planets pl ON pl.id = b.planet_id LEFT JOIN stars st ON st.id = pl.star_id"
              " LEFT JOIN star_systems ss ON ss.id = pl.star_system_id"),
}
"""dict: Per table, the distance from the star and the joins that find its star's light."""


def worked_out(row, table):
    """`(albedo, greenhouse_factor, basis)` for one stored row (all `None` when there is no starlight or no temperature)."""
    luminosity, orbit_km = row["luminosity_w"], row["orbit_km"]
    flux = climate_draws.starlight_w_m2(luminosity, orbit_km) if luminosity and orbit_km else None
    if not row["surface_temperature_k"] or not flux:
        return None, None, None
    albedo, factor, _reason = climate_draws.from_stored(row, flux)
    if albedo is not None:
        return albedo, factor, climate_draws.BASIS_INVERSE
    albedo, factor = climate_draws.class_draw(row, flux, f"{table}:{row['id']}")
    return albedo, factor, climate_draws.BASIS_CLASS_DRAW


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        changes = [f"ADD COLUMN {column} {definition}" for column, definition in COLUMNS
                   if not helpers.exists_column(connection, table, column)]
        if changes:  # no ALGORITHM=INSTANT: both tables carry a FULLTEXT index, which refuses it
            connection.execute(sa.text(f"ALTER TABLE {table} {', '.join(changes)}"))
    for table in ("planets", "moons"):
        orbit, joins = _JOINS[table]
        query = sa.text(_SELECT.format(table=table, orbit=orbit, joins=joins))
        total = connection.execute(sa.text(
            f"SELECT COUNT(*) FROM {table} WHERE albedo_basis IS NULL AND bond_albedo IS NULL")).scalar()
        done = last = 0
        while True:
            rows = connection.execute(query, {"last": last, "n": PAGE}).mappings().all()
            if not rows:
                break
            updates = []
            for row in rows:
                albedo, factor, basis = worked_out(row, table)
                if basis is not None:
                    updates.append({"id": row["id"], "albedo": albedo, "factor": factor, "basis": basis})
            if updates:
                connection.execute(sa.text(
                    f"UPDATE {table} SET bond_albedo = :albedo, greenhouse_factor = :factor, albedo_basis = :basis"
                    " WHERE id = :id"), updates)
            connection.commit()
            last = rows[-1]["id"]
            done += len(rows)
            alembic_runner.report_progress(f"Albedo and greenhouse of {table}", done, total)
