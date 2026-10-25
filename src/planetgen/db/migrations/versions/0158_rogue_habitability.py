"""A rogue planet keeps the habitability score of its surface.

Schema v158. `rogue_planets` gain the score columns a star system's planet carries (`phi4` and its six domains, the six
tiers, the four life scores, the four likelihoods behind microbial life, the survival equipment levels with `shelter_column_g_cm2`
and `energy_flux_w_m2`; no regional temperatures, there is no star), all NULL for a gas giant or a
brown dwarf. A rogue planet saved before gets them from its stored surface conditions by
`planetgen.physics.rogue_habitability`: the score is a pure function of the surface regime, pressure, temperature, water and
heat the row already holds (nothing is drawn), so the values are the ones a new galaxy would give the same planet.

The columns go in with one ALTER, rows are filled in pages of ids (a row with `phi4` set is done, so a stopped run starts again
where it stopped).
"""

import types

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner
from planetgen.db.migrations import helpers
from planetgen.physics import rogue_habitability

revision = "0158"
down_revision = "0157"
branch_labels = None
depends_on = None

COLUMNS = (
    ("phi4", "DOUBLE"), ("phi4_pressure", "DOUBLE"), ("phi4_temperature", "DOUBLE"), ("phi4_air", "DOUBLE"),
    ("phi4_water", "DOUBLE"), ("phi4_radiation", "DOUBLE"), ("phi4_suitability", "DOUBLE"),
    ("tier_pressure", "TINYINT"), ("tier_temperature", "TINYINT"), ("tier_air", "TINYINT"), ("tier_water", "TINYINT"),
    ("tier_radiation", "TINYINT"), ("tier_suitability", "TINYINT"),
    ("life_microbial", "TINYINT UNSIGNED"), ("life_flora", "TINYINT UNSIGNED"), ("life_fauna", "TINYINT UNSIGNED"),
    ("life_sentient", "TINYINT UNSIGNED"),
    ("l_solv", "DOUBLE"), ("l_chem", "DOUBLE"), ("l_ener", "DOUBLE"), ("l_rad", "DOUBLE"),
    ("eq_body", "TINYINT"), ("eq_body_kind", "TINYINT"), ("eq_body_rec", "TINYINT"), ("eq_body_rec_kind", "TINYINT"),
    ("eq_rad_suit", "TINYINT"), ("eq_rad_suit_rec", "TINYINT"), ("eq_mobility", "TINYINT"), ("eq_dust", "TINYINT"),
    ("eq_gas_filter", "TINYINT"), ("eq_co2_scrub", "TINYINT"), ("eq_supply", "TINYINT"), ("eq_supply_rec", "TINYINT"),
    ("eq_shelter", "TINYINT"), ("eq_shelter_rec", "TINYINT"), ("shelter_column_g_cm2", "DOUBLE"),
    ("energy_flux_w_m2", "DOUBLE"),
)
"""tuple: `(column, definition)` added to `rogue_planets`, in table order."""

PAGE = 5000
"""int: Rows read, worked out and written before a commit."""

_SELECT = """
    SELECT id, planet_type, planet_class, mass_kg, radius_km, age_gy, internal_heat_flux_w_m2, surface_regime,
           surface_temperature_k, surface_pressure_pa, ice_shell_thickness_km, ocean_depth_km, has_liquid_water, hp_ice_km
    FROM rogue_planets
    WHERE id > :last AND phi4 IS NULL AND planet_type = 't' AND surface_temperature_k IS NOT NULL
    ORDER BY id LIMIT :n"""


def worked_out(row):
    """The `rogue_habitability.ROGUE_FIELDS` of one stored row, as a dict."""
    planet = types.SimpleNamespace(**dict(row))
    planet.has_liquid_water = bool(planet.has_liquid_water)
    rogue_habitability.generate(planet)
    return {name: getattr(planet, name) for name in rogue_habitability.ROGUE_FIELDS}


def upgrade():
    connection = op.get_bind()
    changes = [f"ADD COLUMN {column} {definition} NULL" for column, definition in COLUMNS
               if not helpers.exists_column(connection, "rogue_planets", column)]
    if changes:
        connection.execute(sa.text(f"ALTER TABLE rogue_planets {', '.join(changes)}"))
    query = sa.text(_SELECT)
    total = connection.execute(sa.text(
        "SELECT COUNT(*) FROM rogue_planets WHERE phi4 IS NULL AND planet_type = 't' AND surface_temperature_k IS NOT NULL")).scalar()
    assignments = ", ".join(f"{name} = :{name}" for name in rogue_habitability.ROGUE_FIELDS)
    done = last = 0
    while True:
        rows = connection.execute(query, {"last": last, "n": PAGE}).mappings().all()
        if not rows:
            break
        connection.execute(sa.text(f"UPDATE rogue_planets SET {assignments} WHERE id = :id"), [
            {"id": row["id"], **worked_out(row)} for row in rows])
        connection.commit()
        last = rows[-1]["id"]
        done += len(rows)
        alembic_runner.report_progress("Habitability score of rogue planets", done, total)
