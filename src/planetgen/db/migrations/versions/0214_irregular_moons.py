"""The irregular moons of a planet (GEN.1084).

Schema v214. Two new tables. `irregular_moons`: one row for each small body a planet caught (`generation/irregular_moons.py`; the stored ones are every body above 1e16 kg and more of the next heaviest until 90 percent of the
population's mass is held, at most 100) with its window, mass, radius, density, taxon, source distance, ice lost, composition shares and orbit. `irregular_sets`: one row for a planet that has any, with the expected mass of the bodies not stored, how many were above the floor, how many were lost and the temporary-capture flag.

These are a table, not columns on `planets`, because `planets` has a FULLTEXT key on `name` and MySQL refuses `ALGORITHM=INSTANT` for it and would rebuild the whole table.

A planet made before v214 has no irregular moons to read back (the generator epoch makes them), so the regenerated galaxy replaces it. A table already there is skipped, so a stopped run can be started
again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0214"
down_revision = "0213"
branch_labels = None
depends_on = None

SETS = """
    CREATE TABLE IF NOT EXISTS irregular_sets (
        id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
        planet_id         BIGINT UNSIGNED NOT NULL,
        remainder_kg      DOUBLE NOT NULL,
        above_floor       INT NOT NULL,
        lost              INT NOT NULL,
        temporary         TINYINT(1) NOT NULL,
        CONSTRAINT fk_irregular_sets_planet
            FOREIGN KEY (planet_id) REFERENCES planets(id) ON DELETE CASCADE,
        UNIQUE KEY uq_irregular_sets_planet (planet_id)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""

TABLE = """
    CREATE TABLE IF NOT EXISTS irregular_moons (
        id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
        planet_id         BIGINT UNSIGNED NOT NULL,
        moon_rank         INT NOT NULL,
        capture_window    VARCHAR(10) NOT NULL,
        mass_kg           DOUBLE NOT NULL,
        radius_km         DOUBLE NOT NULL,
        density_g_cm3     DOUBLE NOT NULL,
        taxon             VARCHAR(4) NOT NULL,
        source_distance_km DOUBLE NOT NULL,
        ice_loss          DOUBLE NOT NULL,
        rock_fraction     DOUBLE NOT NULL,
        metal_fraction    DOUBLE NOT NULL,
        water_fraction    DOUBLE NOT NULL,
        ice_fraction      DOUBLE NOT NULL,
        distance_km       DOUBLE NOT NULL,
        eccentricity      DOUBLE NOT NULL,
        inclination_deg   DOUBLE NOT NULL,
        ascending_node_deg DOUBLE NOT NULL,
        arg_periapsis_deg DOUBLE NOT NULL,
        mean_anomaly_deg  DOUBLE NOT NULL,
        CONSTRAINT fk_irregular_moons_planet
            FOREIGN KEY (planet_id) REFERENCES planets(id) ON DELETE CASCADE,
        UNIQUE KEY uq_irregular_moons_planet_rank (planet_id, moon_rank)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def upgrade():
    connection = op.get_bind()
    for name, statement in (("irregular_sets", SETS), ("irregular_moons", TABLE)):
        if not helpers.exists_table(connection, name):
            connection.execute(sa.text(statement))
