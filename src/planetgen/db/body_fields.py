# planetgen/db/body_fields.py

"""
The body columns (GEN.1658, the first slice of GEN.570): the columns of the three body tables, `planets`, `moons` and `rogue_planets`, declared once.

A planet, a moon and a rogue planet are one kind of object (`generation/planetary_object.PlanetaryObject`), and their tables share all but a few dozen of about 360 columns each.
Each column was typed out three times in `schema.sql`, once more in `column_units.py` for its unit, and a new field had to be named in ten places. Here a column is one row, `F(...)`:

* `name`, `type` (`DOUBLE`, `VARCHAR(16)`, `BIGINT UNSIGNED`), `nn` (`NOT NULL`) and `extra` (a `DEFAULT`, a `CHECK`, a key clause on the column) are the SQL;
* `unit` is the key of the unit registry (`physics.unit_registry`) or the word (`id`, `index`, `ratio`, `count` ...) the column is declared in, which `db/column_units.py` serves from
  (`units_by_group`) and `PlanetaryObject.FIELD_UNITS` will be served from (GEN.570); a column with no unit (text, a flag, a timestamp) has none;
* `c` is the comment written after it in `schema.sql`, `lead` the comment lines written above it and `post` the lines of a key written right after it (one, in `rogue_planets`);
* `t` is the tables that have it, as letters, `p` planets, `m` moons, `r` rogue planets (all three when left out).

Any of `type`, `nn`, `extra`, `unit`, `c`, `lead` and `post` can be a dict `{letters: value}` where the tables differ (a moon's `parent_type` defaults to `'planet'`, a planet's to `'star'`).
The columns are in the order a new database holds them; a table keeps its own columns in that order.

`scripts/build_body_tables.py` writes the three `CREATE TABLE` statements into `schema.sql` between the marker lines `-- BEGIN GENERATED body tables` and `-- END GENERATED body tables`
(`--check` fails when the file differs, `tests/test_body_fields.py` runs it); the keys, constraints and the engine are `CLAUSES` and `TRAILER` below. To add a column, add its row
here, add its Alembic revision, and run `python scripts/build_body_tables.py` and `python scripts/generate_db_models.py`.
"""

from dataclasses import dataclass

TABLES = {"p": "planets", "m": "moons", "r": "rogue_planets"}
"""dict: The letter of each body table."""
ALL = "pmr"


@dataclass(frozen=True)
class Field:
    """One column of the body tables (module note). A dict value is `{letters: value}` for the tables that differ."""

    name: str
    type: object
    nn: object = False
    extra: object = ""
    unit: object = None
    c: object = None
    lead: object = ()
    post: object = ()
    t: str = ALL

    def of(self, attribute, letter):
        """The value of `attribute` for the table `letter`: the value itself, or the entry of the dict whose letters hold it."""
        value = getattr(self, attribute)
        if isinstance(value, dict):
            for letters, entry in value.items():
                if letter in letters:
                    return entry
            raise KeyError((self.name, attribute, letter))
        return value


F = Field

FIELDS = (
    F('id', 'BIGINT UNSIGNED', extra='AUTO_INCREMENT PRIMARY KEY', unit='id'),
    F('sector_id', 'BIGINT UNSIGNED', unit='id', t='r'),
    F('planet_id', 'BIGINT UNSIGNED', unit='id', c='v140 (DB.73): NULL for a moon of a free planet', t='m'),
    F('star_system_id', 'BIGINT UNSIGNED', nn={'p': True, 'm': False}, unit='id', c={'p': None, 'm': 'v140 (DB.73): NULL for a moon of a free planet'}, t='pm'),
    F('star_id', 'BIGINT UNSIGNED', unit='id', lead={'p': ("-- The specific star this planet orbits, when that's a real stored", '-- `stars` row -- true for every single-star system and, as of v15,', "-- every 'wide' binary's planets too (each orbits one specific", "-- constituent star). Still NULL for a 'close' (P-type) binary's", '-- planets: systemData.py generates those against the merged', '-- `BinaryStarProxy`, never one individual constituent star, and the', '-- proxy is deliberately not stored as its own `stars` row (see', '-- star_systems.binary_* above). `star_system_id` above is always the', '-- reliable owning link regardless of star_id.'), 'm': ()}, t='pm'),
    F('orbital_index', 'INT', nn=True, unit='index', c={'p': None, 'm': "position in the parent planet's moons list"}, t='pm'),
    F('body_type', 'VARCHAR(4)', nn={'pm': True, 'r': False}, extra="CHECK (body_type IN ('t', 'g'))", c={'pm': None, 'r': 'v207 (GEN.1111)'}, lead={'pm': (), 'r': ('-- v207 (GEN.1111): a rogue planet is a planet with no star, so it carries every column a planet does except the orbit round a host (the columns below),', '-- and the system that ejected it, when, and how fast. The older columns named for the flat-rate rogue (planet_type, mass_bin, has_*) are NULL for these rows.')}),
    F('name', 'VARCHAR(255)', nn=True),
    F('planet_class', {'pm': 'VARCHAR(16)', 'r': 'VARCHAR(4)'}, lead={'pm': (), 'r': ('-- v47: its PLANET_CLASSES letter (GEN.8); NULL for a brown dwarf.',)}),
    F('mass_bin', 'VARCHAR(16)', nn=True, extra="DEFAULT 'terrestrial'", lead=('-- v37: the mass bin it was drawn from (program_constants.', "-- ROGUE_PLANET_MASS_BINS), or 'brown-dwarf' (13-80 Mjup)."), t='r'),
    F('distance_km', 'DOUBLE', nn=True, unit='km', c={'p': 'from the star', 'm': 'from the parent planet'}, t='pm'),
    F('radius_km', 'DOUBLE', nn=True, unit='km'),
    F('mass_kg', 'DOUBLE', nn=True, unit='kg'),
    F('volume_km3', 'DOUBLE', nn={'pm': True, 'r': False}, unit='km3', c={'p': 'derived, stored not recomputed', 'm': None, 'r': 'v207 (GEN.1111)'}),
    F('period_years', 'DOUBLE', nn=True, unit='yr', c={'p': 'derived, stored not recomputed', 'm': None}, t='pm'),
    F('zone', 'VARCHAR(4)', extra="CHECK (zone IN ('h', 'e', 'c'))", t='pm'),
    F('gravity_g', 'DOUBLE', unit='g0', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('surface_temperature_k', 'DOUBLE', unit='K'),
    F('surface_pressure_pa', 'DOUBLE', unit='Pa', t='r'),
    F('ice_shell_thickness_km', 'DOUBLE', unit='km', t='r'),
    F('density_g_cm3', 'DOUBLE', unit='g_cm3', c={'pm': None, 'r': 'v200 (GEN.1162): the mean density, g/cm3'}),
    F('atmosphere', 'TEXT', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('atm_density', 'DOUBLE', unit='kg_m3', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('atm_molar_density', 'DOUBLE', unit='kg/mol', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('atmospheric_pressure_pa', 'DOUBLE', unit='Pa', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('composition', 'TEXT', c={'p': 'descriptive string (contrast asteroid_belt_composition)', 'mr': None}),
    F('has_internal_heat', 'TINYINT(1)', extra='CHECK (has_internal_heat IN (0, 1))', unit='flag', c='v207: NULL for a row of the physics rogue', t='r'),
    F('age_gy', 'DOUBLE', unit='Gyr', lead=('-- v48: surface conditions (rogueSurface.rogue_surface_conditions).',), t='r'),
    F('scale_height_km', 'DOUBLE', unit='km', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('hill_radius_km', 'DOUBLE', unit='km', t='pm'),
    F('min_orbit_distance_km', 'DOUBLE', unit='km', t='pm'),
    F('habitable_zone_inner_km', 'DOUBLE', nn=True, unit='km', c={'p': 'copied from host star at generation time', 'm': None}, t='pm'),
    F('habitable_zone_outer_km', 'DOUBLE', nn=True, unit='km', t='pm'),
    F('flavor_text_count', 'INT', nn=True, extra='DEFAULT 0', unit='count', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('life_stage', 'VARCHAR(32)', c={'pm': 'v110 (DB.50): a MILESTONE_KEYS key, NULL before any milestone', 'r': 'v207 (GEN.1111)'}),
    F('evolutionary_pace', 'VARCHAR(16)', c={'pm': 'v110: the pace life evolves at (fast, normal, slow), drawn once, set with any life chemistry (v143)', 'r': 'v207 (GEN.1111)'}),
    F('life_milestone_age_gyr', 'DOUBLE', unit='Gyr', c={'pm': 'v110: the age of that milestone, Gy', 'r': 'v207 (GEN.1111)'}),
    F('life_system_age_gyr', 'DOUBLE', unit='Gyr', c={'pm': 'v110: the age of the system the timeline was written at, Gy', 'r': 'v207 (GEN.1111)'}),
    F('life_age_abiogenesis_gyr', 'DOUBLE', unit='Gyr', c={'pm': 'v156 (DB.100): the age the pace reaches abiogenesis, Gy', 'r': 'v207 (GEN.1111)'}),
    F('life_age_photosynthesis_gyr', 'DOUBLE', unit='Gyr', c={'pm': 'v156 (DB.100): the age the pace reaches photosynthesis, Gy', 'r': 'v207 (GEN.1111)'}),
    F('life_age_complex_cells_gyr', 'DOUBLE', unit='Gyr', c={'pm': 'v156 (DB.100): the age the pace reaches complex cells, Gy', 'r': 'v207 (GEN.1111)'}),
    F('life_age_multicellularity_gyr', 'DOUBLE', unit='Gyr', c={'pm': 'v156 (DB.100): the age the pace reaches multicellularity, Gy', 'r': 'v207 (GEN.1111)'}),
    F('life_age_technological_civilization_gyr', 'DOUBLE', unit='Gyr', c={'pm': 'v156 (DB.100): the age the pace reaches technological civilization, Gy', 'r': 'v207 (GEN.1111)'}),
    F('flavor_key', 'VARCHAR(32)', c={'pm': 'v110: which flavor line (replaces flavor_text)', 'r': 'v207 (GEN.1111)'}),
    F('orbital_inclination_deg', 'DOUBLE', nn=True, unit='deg', lead=('-- Orbital motion (v9, see header comment) -- inclination/ascending', '-- node are fixed at generation time. rotation_period_hours is a', '-- separate, static "day length" stat. There is no stored phase (v109,', '-- GEN.289): the body is its state vector below.'), t='pm'),
    F('orbital_ascending_node_deg', 'DOUBLE', nn=True, unit='deg', t='pm'),
    F('eccentricity', 'DOUBLE', nn=True, extra='DEFAULT 0', unit='orbit-e', c='v159 (GEN.727): of the orbit the body was generated on (distance_km is its semi-major axis)', t='pm'),
    F('position_x_km', 'DOUBLE', nn=True, unit='km', lead={'p': ('-- Position/speed (v11, see header comment) -- Cartesian position', "-- relative to this planet's orbital anchor (the star, or the", "-- BinaryStarProxy's combined center for a binary system), made from", '-- distance_km and the orbital-motion columns above and moved by', '-- planetgen.cli.orbits with kepler.universal_step, with the velocity', '-- below and epoch_unix as its state.', '-- orbital_speed_kms is constant around a circular orbit -- it only', '-- changes if distance_km/period_years do.'), 'm': ('-- Position/speed (v11, see header comment) -- Cartesian position', "-- relative to this moon's orbital anchor, its parent planet; made", '-- from distance_km and the orbital-motion columns above the same way', "-- planets.position_x/y/z_km are (see that table's comment).")}, t='pm'),
    F('position_y_km', 'DOUBLE', nn=True, unit='km', t='pm'),
    F('position_z_km', 'DOUBLE', nn=True, unit='km', t='pm'),
    F('orbital_speed_kms', 'DOUBLE', nn=True, unit='km_s', t='pm'),
    F('velocity_x_kms', 'DOUBLE', nn={'pm': True, 'r': False}, extra={'pm': 'DEFAULT 0', 'r': ''}, unit='km_s', c={'pm': None, 'r': 'v102 (DB.55): galactic-axes velocity, km/s'}, lead={'pm': ('-- v60: velocity relative to the orbital anchor, km/s, in step with the position (see the header comment\'s "v60" note).',), 'r': ()}),
    F('velocity_y_kms', 'DOUBLE', nn={'pm': True, 'r': False}, extra={'pm': 'DEFAULT 0', 'r': ''}, unit='km_s', c={'pm': None, 'r': 'v102 (DB.55): galactic-axes velocity, km/s'}),
    F('velocity_z_kms', 'DOUBLE', nn={'pm': True, 'r': False}, extra={'pm': 'DEFAULT 0', 'r': ''}, unit='km_s', c={'pm': None, 'r': 'v102 (DB.55): galactic-axes velocity, km/s'}, post={'pm': (), 'r': ('    KEY idx_rogue_planets_cell (cell_key, center_x_pc, center_y_pc, center_z_pc),',)}),
    F('galactic_radius_pc', 'DOUBLE', unit='pc', t='r'),
    F('quadrant', 'VARCHAR(4)', extra="CHECK (quadrant IN ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII'))", lead=("-- v41: the sector octant the center sits in (star_systems.quadrant's labels).",), t='r'),
    F('inside_nebula_id', 'BIGINT UNSIGNED', unit='id', lead=('', '-- v27: row timestamps -- see the header comment\'s "v27" note.', '-- v39: the innermost nebula or supernova remnant this sits inside', '-- (at most one set; _db.refresh_containment).'), t='r'),
    F('inside_remnant_id', 'BIGINT UNSIGNED', unit='id', t='r'),
    F('created_at', 'TIMESTAMP', nn=True, extra='DEFAULT CURRENT_TIMESTAMP', lead=('',), t='r'),
    F('modified_at', 'TIMESTAMP(3)', nn=True, extra='DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3)', t='r'),
    F('epoch_unix', 'DOUBLE', unit='unix_s', c='v66 (GEN.106), see header comment'),
    F('next_update_due', 'DOUBLE', unit='unix_s'),
    F('min_update_interval_years', 'DOUBLE', nn=True, unit='yr', c='v12, see header comment', t='pm'),
    F('rotation_period_hours', 'DOUBLE', nn={'pm': True, 'r': False}, unit='h', lead={'pm': (), 'r': ('-- v68 (GEN.104): the spin, see header comment.',)}),
    F('spin_axis_x', 'DOUBLE', unit='unit-vector', lead={'pm': ('-- v68 (GEN.104): the spin, see header comment.',), 'r': ()}),
    F('spin_axis_y', 'DOUBLE', unit='unit-vector'),
    F('spin_axis_z', 'DOUBLE', unit='unit-vector'),
    F('axial_tilt_deg', 'DOUBLE', unit='deg'),
    F('center_x_pc', 'DOUBLE', unit='pc', lead=('', '-- v28: this rogue planet\'s own galaxy-frame center -- see this file\'s "v28"', '-- header note. NULL together: never placed in the galaxy.'), t='r'),
    F('center_y_pc', 'DOUBLE', unit='pc', t='r'),
    F('center_z_pc', 'DOUBLE', unit='pc', t='r'),
    F('cell_key', 'BIGINT UNSIGNED', extra='GENERATED ALWAYS AS (IF(ABS(center_x_pc) < 262144 AND ABS(center_y_pc) < 262144 AND ABS(center_z_pc) < 262144, ((((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED', unit='key', lead=('-- v93 (region queries): where the object is now as one integer, the Morton key of its 0.25 pc cube; MySQL computes it,', '-- NULL outside the box of plus or minus 262,144 pc. `planetgen/galaxy/cellkey.py` has the same function.'), t='r'),
    F('mantle_redox', 'VARCHAR(16)', extra="CHECK (mantle_redox IN ('reduced', 'intermediate', 'oxidized'))", c={'pm': None, 'r': 'v207 (GEN.1111)'}, lead={'pm': ("-- v70 (GEN.85): mantle redox, see header comment (the air's partial pressures went to `planet_gases` and `moon_gases` in v155).",), 'r': ()}),
    F('mantle_delta_iw', 'DOUBLE', unit='dex', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('magnetic_moment_a_m2', 'DOUBLE', unit='A*m2', c={'pm': None, 'r': 'v207 (GEN.1111)'}, lead={'pm': ("-- v72 (GEN.86): magnetic field and the star's XUV and flares here, see header comment.",), 'r': ()}),
    F('dipole_class', 'VARCHAR(16)', extra="CHECK (dipole_class IN ('none', 'weak', 'earth-like', 'strong', 'multipolar'))", c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('dipole_tilt_deg', 'DOUBLE', unit='deg', c='v229 (GEN.1350): the angle between the dipole axis and the spin axis, 0 to 90 degrees, NULL without a field and for a row made before v229'),
    F('dipole_offset_class', 'VARCHAR(8)', c="v229 (GEN.1350): 'centered' or 'offset' (the dipole center over 0.15 radii from the center of the body), NULL without a field"),
    F('induced_moment_a_m2', 'DOUBLE', unit='A*m2', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('magnetopause_rp', 'DOUBLE', unit='R_planet', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('xuv_flux_earth', 'DOUBLE', unit='xuv-earth', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('xuv_exposure_index', 'DOUBLE', unit='xuv-ref', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('flare_irradiation_index', 'DOUBLE', unit='flares-per-au2', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('water_mass_fraction', 'DOUBLE', unit='frac', c={'pm': None, 'r': 'v207 (GEN.1111)'}, lead={'pm': ("-- v74 (GEN.88): water, land, ice and the ocean's chemistry, see header comment.",), 'r': ()}),
    F('hydrosphere', 'VARCHAR(20)', extra="CHECK (hydrosphere IN ('dry', 'vapor', 'ice', 'ice-covered ocean', 'surface ocean', 'brine seep', 'hycean'))", c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('ocean_fraction', 'DOUBLE', unit='frac', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('land_fraction', 'DOUBLE', unit='frac', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('ice_cover_fraction', 'DOUBLE', unit='frac', c={'pm': 'v249 (GEN.1567): the share of the surface the ice covers (an ice lid, or a body that is ice through), `ocean_fraction` is the liquid sea alone', 'r': 'v249 (GEN.1567)'}),
    F('ocean_depth_km', 'DOUBLE', unit='km'),
    F('has_liquid_water', 'TINYINT(1)', unit='flag', t='r'),
    F('ice_shell_km', 'DOUBLE', unit='km', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('hp_ice_km', 'DOUBLE', unit='km', c={'pm': None, 'r': 'v74 (GEN.88)'}),
    F('ocean_class', 'VARCHAR(16)', extra="CHECK (ocean_class IN ('ice-sealed', 'chloride brine', 'acid sulfate', 'soda', 'neutral'))", c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('ocean_ph', 'DOUBLE', unit='pH', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('water_activity', 'DOUBLE', unit='aw', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('phosphorus', 'VARCHAR(8)', extra="CHECK (phosphorus IN ('high', 'limited', 'starved'))", c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('surface_dose_msv_yr', 'DOUBLE', unit='mSv/yr', c={'pm': None, 'r': 'v207 (GEN.1111)'}, lead={'pm': ('-- v75 (GEN.87): surface dose, UV and ozone, see header comment.',), 'r': ()}),
    F('dose_gcr_msv_yr', 'DOUBLE', unit='mSv/yr', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('dose_sep_msv_yr', 'DOUBLE', unit='mSv/yr', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('dose_ground_msv_yr', 'DOUBLE', unit='mSv/yr', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('dose_belt_msv_yr', 'DOUBLE', unit='mSv/yr', c={'pm': "v187 (GEN.983): the dose a moon takes from its parent's trapped-particle belts, mSv/yr (0 for a planet)", 'r': 'v207 (GEN.1111)'}),
    F('dose_helio_mult', 'DOUBLE', unit='dose-mult', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('uv_surface_index', 'DOUBLE', unit='dna-mult', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('ozone_loss_flag', 'BOOLEAN', unit='flag', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('phi4', 'DOUBLE', unit='score1', c={'pm': None, 'r': 'v158 (GEN.701): the habitability score of the surface a visitor would find (rogue_habitability.py)'}, lead={'pm': ('-- v76 (GEN.89): the habitability score, see header comment.',), 'r': ()}),
    F('phi4_pressure', 'DOUBLE', unit='score1', c={'pm': None, 'r': 'v153'}),
    F('phi4_temperature', 'DOUBLE', unit='score1', c={'pm': None, 'r': 'v153'}),
    F('phi4_air', 'DOUBLE', unit='score1', c={'pm': None, 'r': 'v153'}),
    F('phi4_water', 'DOUBLE', unit='score1', c={'pm': None, 'r': 'v153'}),
    F('phi4_suitability', 'DOUBLE', unit='score1', c={'pm': None, 'r': 'v153'}),
    F('phi4_radiation', 'DOUBLE', unit='score1', c={'pm': None, 'r': 'v153'}),
    F('tier_pressure', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('tier_temperature', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('tier_air', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('tier_water', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('tier_suitability', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('tier_radiation', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('life_microbial', 'TINYINT UNSIGNED', unit='score9', c={'pm': 'v150 (GEN.649): 0 none to 9 abundant, from the evolutionary path and the surface', 'r': 'v153'}),
    F('life_flora', 'TINYINT UNSIGNED', unit='score9', c={'pm': 'v150', 'r': 'v153'}),
    F('life_fauna', 'TINYINT UNSIGNED', unit='score9', c={'pm': 'v150', 'r': 'v153'}),
    F('life_sentient', 'TINYINT UNSIGNED', unit='score9', c={'pm': 'v150', 'r': 'v153'}),
    F('l_solv', 'DOUBLE', unit='score1', c={'pm': None, 'r': 'v153'}),
    F('l_chem', 'DOUBLE', unit='score1', c={'pm': None, 'r': 'v153'}),
    F('l_ener', 'DOUBLE', unit='score1', c={'pm': None, 'r': 'v153'}),
    F('l_rad', 'DOUBLE', unit='score1', c={'pm': None, 'r': 'v153'}),
    F('eq_body', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('eq_body_kind', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('eq_body_rec', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('eq_body_rec_kind', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('eq_rad_suit', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('eq_rad_suit_rec', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('eq_mobility', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('eq_dust', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('eq_gas_filter', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('eq_co2_scrub', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('eq_supply', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('eq_supply_rec', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('eq_shelter', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('eq_shelter_rec', 'TINYINT', unit='index', c={'pm': None, 'r': 'v153'}),
    F('shelter_column_g_cm2', 'DOUBLE', unit='g_cm2', c={'pm': None, 'r': 'v153'}),
    F('outside_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band1_mean_k', 'DOUBLE', unit='K', c={'pm': 'v196 (GEN.1108): the climate of the body in seven bands of latitude (`physics/latitude_bands.py`): `band1` is the north pole, `band7` the south pole (a world locked to its star: the day center to the night center)', 'r': 'v207 (GEN.1111)'}),
    F('band1_swing_k', 'DOUBLE', unit='dK', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band1_day_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band1_night_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band1_longest_day_days', 'DOUBLE', unit='d', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band1_longest_night_days', 'DOUBLE', unit='d', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band1_state', 'VARCHAR(24)', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band2_mean_k', 'DOUBLE', unit='K', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band2_swing_k', 'DOUBLE', unit='dK', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band2_day_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band2_night_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band2_longest_day_days', 'DOUBLE', unit='d', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band2_longest_night_days', 'DOUBLE', unit='d', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band2_state', 'VARCHAR(24)', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band3_mean_k', 'DOUBLE', unit='K', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band3_swing_k', 'DOUBLE', unit='dK', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band3_day_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band3_night_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band3_longest_day_days', 'DOUBLE', unit='d', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band3_longest_night_days', 'DOUBLE', unit='d', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band3_state', 'VARCHAR(24)', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band4_mean_k', 'DOUBLE', unit='K', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band4_swing_k', 'DOUBLE', unit='dK', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band4_day_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band4_night_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band4_longest_day_days', 'DOUBLE', unit='d', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band4_longest_night_days', 'DOUBLE', unit='d', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band4_state', 'VARCHAR(24)', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band5_mean_k', 'DOUBLE', unit='K', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band5_swing_k', 'DOUBLE', unit='dK', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band5_day_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band5_night_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band5_longest_day_days', 'DOUBLE', unit='d', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band5_longest_night_days', 'DOUBLE', unit='d', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band5_state', 'VARCHAR(24)', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band6_mean_k', 'DOUBLE', unit='K', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band6_swing_k', 'DOUBLE', unit='dK', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band6_day_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band6_night_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band6_longest_day_days', 'DOUBLE', unit='d', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band6_longest_night_days', 'DOUBLE', unit='d', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band6_state', 'VARCHAR(24)', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band7_mean_k', 'DOUBLE', unit='K', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band7_swing_k', 'DOUBLE', unit='dK', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band7_day_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band7_night_hours', 'DOUBLE', unit='h', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band7_longest_day_days', 'DOUBLE', unit='d', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band7_longest_night_days', 'DOUBLE', unit='d', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('band7_state', 'VARCHAR(24)', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('temp_mean_swing_k', 'DOUBLE', unit='dK', c={'pm': 'v182 (GEN.848): the day and night swing of the whole surface, peak to peak, kelvin of difference', 'r': 'v207 (GEN.1111)'}),
    F('temp_livable_share', 'DOUBLE', unit='frac', c="v212 (GEN.1265): L, the share of the surface's area and time a person can be out in without a suit (Green band or better, wet bulb rule per sample), read off the temperature of every place through the year (`physics/temperature_coverage.py`), the Temperature score is made from it"),
    F('temp_livable_area', 'DOUBLE', unit='frac', c='v212: the share of the surface livable at least half of the time'),
    F('temp_livable_time', 'DOUBLE', unit='frac', c='v212: the mean share of the time that ground is livable'),
    F('temp_blue_share', 'DOUBLE', unit='frac', c='v212: the share of area and time in the Blue band (0 to 31 C)'),
    F('temp_green_share', 'DOUBLE', unit='frac', c='v212: the share in the Green band and not Blue (-20 to 45 C)'),
    F('temp_cold_share', 'DOUBLE', unit='frac', c='v212: the share below -20 C'),
    F('temp_hot_share', 'DOUBLE', unit='frac', c='v212: the share above 45 C'),
    F('temp_area_always', 'DOUBLE', unit='frac', c='v212: the share of the surface livable all the time'),
    F('temp_area_most', 'DOUBLE', unit='frac', c='v212: the share of the surface livable at least 90 percent of the time'),
    F('temp_livable_belt_only', 'TINYINT', unit='flag', c='v212: 1 when a world that turns once an orbit is livable only in its terminator belt'),
    F('temp_coldest_k', 'DOUBLE', unit='K', c='v212: the coldest sample anywhere, kelvin'),
    F('temp_coldest_lat_deg', 'DOUBLE', unit='deg', c='v212: where it is, degrees'),
    F('temp_coldest_lon_deg', 'DOUBLE', unit='deg', c='v212: where it is, degrees'),
    F('temp_coldest_day', 'DOUBLE', unit='d', c='v212: when it is, days from an equinox'),
    F('temp_hottest_k', 'DOUBLE', unit='K', c='v212: the hottest sample anywhere, kelvin'),
    F('temp_hottest_lat_deg', 'DOUBLE', unit='deg', c='v212: where it is, degrees'),
    F('temp_hottest_lon_deg', 'DOUBLE', unit='deg', c='v212: where it is, degrees'),
    F('temp_hottest_day', 'DOUBLE', unit='d', c='v212: when it is, days from an equinox'),
    F('temp_site_lat_deg', 'DOUBLE', unit='deg', c='v212: the best site, the cell livable the largest share of the time (ties go to the most moderate), latitude'),
    F('temp_site_lon_deg', 'DOUBLE', unit='deg', c='v212: its longitude, degrees'),
    F('temp_site_livable_time', 'DOUBLE', unit='frac', c='v212: the share of the time it is livable'),
    F('temp_site_coldest_k', 'DOUBLE', unit='K', c='v212: its coldest sample, kelvin'),
    F('temp_site_hottest_k', 'DOUBLE', unit='K', c='v212: its hottest sample, kelvin'),
    F('temp_site_mean_k', 'DOUBLE', unit='K', c='v212: its mean temperature, kelvin'),
    F('temp_site_longest_out_days', 'DOUBLE', unit='d', c='v212: its longest unbroken spell out of the livable band, days'),
    F('outside_hours_day', 'DOUBLE', unit='h', c={'pm': 'v187 (GEN.768): the hours of one local day a person can be outside in the listed equipment', 'r': 'v207 (GEN.1111)'}),
    F('outside_hours_limit', 'VARCHAR(16)', c={'pm': 'v187 (GEN.768): what sets them: daylight, dose, heat, cold, air, uv or vehicle', 'r': 'v207 (GEN.1111)'}),
    F('day_night_state', 'VARCHAR(24)', c={'pm': 'v188 (GEN.1097): the state of its day and night: day and night, tidally locked, spin-orbit resonance, non-rotating, polar day and night or no solid surface', 'r': 'v207 (GEN.1111)'}),
    F('night_hours', 'DOUBLE', unit='h', c={'pm': 'v188 (GEN.1097): hours of night in one solar day at the equator, NULL for a world that turns once an orbit', 'r': 'v207 (GEN.1111)'}),
    F('days_per_year', 'DOUBLE', unit='frac', c={'pm': 'v188 (GEN.1097): solar days in one orbit round the star, NULL for a world that turns once an orbit', 'r': 'v207 (GEN.1111)'}),
    F('roundness_ratio', 'DOUBLE', unit='round-ratio', c={'pm': 'v199 (GEN.1159): the mean radius over the round radius, 1 or more is round and less is not', 'r': 'v200 (GEN.1162): the mean radius over the round radius, 1 or more is round, NULL for a giant or a brown dwarf'}),
    F('yield_strength_pa', 'DOUBLE', unit='Pa', c={'pm': 'v199 (GEN.1159): the whole-body strength of the material mix of the body, Pa, drawn from the range of each material (GEN.1178)', 'r': 'v200 (GEN.1162): the strength of the rock and ice, Pa, NULL for a giant or a brown dwarf'}),
    F('binding_energy_j', 'DOUBLE', unit='J', c={'pm': 'v199 (GEN.1159): the energy that assembling the body releases, J', 'r': 'v200 (GEN.1162): the energy that assembling it releases, J'}),
    F('shape_temperature_k', 'DOUBLE', unit='K', c="v252 (GEN.1201): the temperature the body took its shape at, K (surface, accretion and radioactive heat, tides, capped at the rock's melting), the strength reads it"),
    F('planet_type', 'VARCHAR(12)', c={'pm': 'v206 (GEN.1219): the body type, from the envelope share of the mass', 'r': 'v207: the body type of GEN.1219: terrestrial, transitional or gas'}),
    F('class_variant', 'TINYINT', unit='class-digit', c={'pm': 'v206 (GEN.1220): the digit of the class code, see the header comment', 'r': 'v207 (GEN.1111)'}),
    F('relative_humidity', 'DOUBLE', unit='frac', c={'pm': "v212 (GEN.1221): the air's water vapor over saturation, 0.05 to 0.99", 'r': 'v212 (GEN.1221)'}),
    F('groundwater', 'VARCHAR(10)', c={'pm': 'v212 (GEN.1221): none, reachable or deep, the liquid water under a dry surface', 'r': 'v212 (GEN.1221)'}),
    F('class_match', 'DOUBLE', unit='frac', c='v218 (GEN.1207): the percent match of the class as a fraction, 0 when no class fits the body'),
    F('runner_up_class', 'CHAR(1)', c='v218 (GEN.1207): the letter of the next best class, NULL when only one passes the gates'),
    F('runner_up_match', 'DOUBLE', unit='frac', c='v218 (GEN.1207): its match'),
    F('class_fudge_values', 'VARCHAR(120)', c="v244 (GEN.1314): the values fudged to fit the class (the quantity names, comma separated, the temperatures count as 'temperature_k'), NULL for a full fit and for a body with no class"),
    F('class_fudge_amount', 'DOUBLE', unit='fudge-widths', c='v244 (GEN.1314): the total move of those values in widths of the fudge table, 0 for a full fit, NULL for a body with no class or a value it does not store'),
    F('oddities', 'VARCHAR(40)', c="v228 (GEN.1047): the ids of the oddities the body has, comma separated ('W10,G1'), NULL for none and for a row made before v221"),
    F('sky_suns', 'VARCHAR(400)', c="v236 (GEN.1458): the other suns in the planet's sky, JSON `[[share, ux, uy, uz], ...]` (each other star's share of the mean light and its unit direction in the planet's equinox frame, `physics/sky_suns.py`), NULL for a lone star and for a row made before v236"),
    F('light_factor', 'DOUBLE', unit='frac', c="v242 (GEN.1484): the share of its host's light a circumbinary planet keeps after its two stars eclipse each other (`physics/eclipses.py`), NULL for every other body and for a row made before v241"),
    F('obliquity_min_deg', 'DOUBLE', unit='deg', c="v247 (GEN.1535): the lowest tilt the planet's axis reaches over the age of its star (`physics/obliquity_variation.py`), degrees"),
    F('obliquity_max_deg', 'DOUBLE', unit='deg', c="v247 (GEN.1535): the highest tilt the planet's axis reaches over the age of its star (`physics/obliquity_variation.py`), degrees"),
    F('belt_order', 'INT', unit='index', c="v248 (GEN.1510): for a planet promoted from an asteroid belt, the belt's number among its star's belts (as `irregular_moons.belt_order`), NULL for every other planet and for a row made before v244", t='p'),
    F('belt_rank', 'INT', unit='index', c="v248 (GEN.1510): for a planet promoted from an asteroid belt, the rank it had among the belt's members (the number is not reused), NULL for every other planet", t='p'),
    F('energy_flux_w_m2', 'DOUBLE', unit='W_m2', c={'pm': None, 'r': 'v153'}),
    F('flux_variation', 'DOUBLE', unit='frac', c={'pm': None, 'r': 'v207 (GEN.1111)'}),
    F('reflex_offset_x_km', 'DOUBLE', unit='km', lead=('-- v20 (see header comment): this planet\'s own reflex-offset "wobble"', '-- from the combined pull of its own moons -- NULL/0 with no moons.'), t='p'),
    F('reflex_offset_y_km', 'DOUBLE', unit='km', t='p'),
    F('reflex_offset_z_km', 'DOUBLE', unit='km', t='p'),
    F('uid', 'BINARY(12)', lead=('-- v58: this object\'s unique ID, see the header comment\'s "v58" note. NULL for a row saved before v58.',)),
    F('parent_type', 'VARCHAR(16)', nn=True, extra={'p': "DEFAULT 'star'", 'm': "DEFAULT 'planet'", 'r': "DEFAULT 'free'"}, c={'pm': 'v102 (DB.55), v140 (DB.73): wider', 'r': 'v140 (DB.73): what it orbits, see the v140 note'}, lead={'pm': (), 'r': ('',)}),
    F('parent_id', 'BIGINT UNSIGNED', unit='id', c={'pm': 'v102 (DB.55)', 'r': 'v140 (DB.73)'}),
    F('dynamical_state', 'VARCHAR(8)', nn=True, extra="DEFAULT 'bound'", c='v102 (DB.55)', t='pm'),
    F('host_node_id', 'SMALLINT', unit='id', c='v104 (DB.53): the orbit node the planet orbits (NULL: its star)', t='p'),
    F('born_beyond_snow_line', 'BOOLEAN', unit='flag', c={'p': "v132 (GEN.458): formed beyond its star's snow line, NULL when the host has none to compare", 'm': "v132 (GEN.458): the parent planet's flag"}, t='pm'),
    F('ozone_column_du', 'DOUBLE', unit='DU', c={'pm': 'v135 (GEN.413): the ozone column, Dobson units', 'r': 'v207 (GEN.1111)'}),
    F('tau_haze', 'DOUBLE', unit='tau', c={'pm': 'v138 (GEN.414, DB.68): haze optical depth at 550 nm, NULL for a row saved before v138 or a body with no air', 'r': 'v207 (GEN.1111)'}),
    F('tau_dust', 'DOUBLE', unit='tau', c={'pm': 'v138 (GEN.414, DB.68): dust optical depth at 550 nm', 'r': 'v207 (GEN.1111)'}),
    F('tau_cloud', 'DOUBLE', unit='tau', c={'pm': 'v138 (GEN.414, DB.68): in-cloud optical depth at 550 nm (grey)', 'r': 'v207 (GEN.1111)'}),
    F('cloud_cover', 'DOUBLE', unit='frac', c={'pm': 'v138 (GEN.414, DB.68): share of the sky under cloud, 0 to 1', 'r': 'v207 (GEN.1111)'}),
    F('aerosol_draw', 'DOUBLE', unit='ratio', c={'pm': 'v138 (GEN.414, DB.68): the number in [0, 1) that picks the optics inside the brackets', 'r': 'v207 (GEN.1111)'}),
    F('rogue_planet_id', 'BIGINT UNSIGNED', unit='id', c='v140 (DB.73): the rogue planet a moon of a free planet orbits', t='m'),
    F('c_k0_ppm', 'DOUBLE', unit='ppm', c={'pm': 'v141 (GEN.459): potassium in the bulk silicate at formation, ppm by mass', 'r': 'v207 (GEN.1111)'}),
    F('c_th0_ppm', 'DOUBLE', unit='ppm', c={'pm': 'v141: thorium, the same', 'r': 'v207 (GEN.1111)'}),
    F('c_u0_ppm', 'DOUBLE', unit='ppm', c={'pm': 'v141: uranium, the same', 'r': 'v207 (GEN.1111)'}),
    F('al26_ratio', 'DOUBLE', unit='al26', c={'pm': 'v141: 26Al/27Al at formation', 'r': 'v207 (GEN.1111)'}),
    F('urey_ratio', 'DOUBLE', unit='frac', c={'pm': 'v141: the radiogenic share of the heat flowing out, today (an output of the thermal history since v146)', 'r': 'v207 (GEN.1111)'}),
    F('bond_albedo', 'DOUBLE', unit='frac', c={'pm': 'v142 (GEN.412): the Bond albedo drawn for the temperature', 'r': 'v207 (GEN.1111)'}),
    F('greenhouse_factor', 'DOUBLE', unit='greenhouse', c={'pm': 'v142: the greenhouse quantity G (T_s = T_e (1 + G)^(1/4)), 0 for no air', 'r': 'v207 (GEN.1111)'}),
    F('albedo_basis', 'TINYINT UNSIGNED', unit='index', c={'pm': 'v142: 0 recorded at generation, 1 exact inverse of a stored row, 2 seeded class draw', 'r': 'v207 (GEN.1111)'}),
    F('metal_mass_fraction', 'DOUBLE', unit='frac', c={'pm': "v144 (DB.69, GEN.428): the iron-rich core's share of the whole body's mass, light elements included", 'r': 'v207 (GEN.1111)'}),
    F('rock_mass_fraction', 'DOUBLE', unit='frac', c={'pm': 'v144: the silicate (or silicon carbide and graphite) share', 'r': 'v207 (GEN.1111)'}),
    F('ice_mass_fraction', 'DOUBLE', unit='frac', c={'pm': 'v144: every ice other than water, as one share', 'r': 'v207 (GEN.1111)'}),
    F('ice_co2_mass_fraction', 'DOUBLE', unit='frac', c={'pm': "v183 (GEN.834): the share of the whole body's mass that is carbon dioxide ice", 'r': 'v207 (GEN.1111)'}),
    F('ice_ch4_mass_fraction', 'DOUBLE', unit='frac', c={'pm': "v183 (GEN.834): the share of the whole body's mass that is methane ice", 'r': 'v207 (GEN.1111)'}),
    F('ice_nh3_mass_fraction', 'DOUBLE', unit='frac', c={'pm': "v183 (GEN.834): the share of the whole body's mass that is ammonia ice", 'r': 'v207 (GEN.1111)'}),
    F('ice_n2_mass_fraction', 'DOUBLE', unit='frac', c={'pm': "v183 (GEN.834): the share of the whole body's mass that is nitrogen ice", 'r': 'v207 (GEN.1111)'}),
    F('ice_so2_mass_fraction', 'DOUBLE', unit='frac', c={'pm': "v183 (GEN.834): the share of the whole body's mass that is sulfur dioxide ice", 'r': 'v207 (GEN.1111)'}),
    F('ice_co_mass_fraction', 'DOUBLE', unit='frac', c={'pm': "v183 (GEN.834): the share of the whole body's mass that is carbon monoxide ice", 'r': 'v207 (GEN.1111)'}),
    F('ice_hcn_mass_fraction', 'DOUBLE', unit='frac', c={'pm': "v183 (GEN.834): the share of the whole body's mass that is hydrogen cyanide ice", 'r': 'v207 (GEN.1111)'}),
    F('ice_other_mass_fraction', 'DOUBLE', unit='frac', c={'pm': "v183 (GEN.834): the share of the whole body's mass that is every other ice (ethane, methanol, formaldehyde, hydrogen sulfide) ice", 'r': 'v207 (GEN.1111)'}),
    F('envelope_mass_fraction', 'DOUBLE', unit='frac', c={'pm': 'v144: the hydrogen and helium share', 'r': 'v207 (GEN.1111)'}),
    F('rock_feo_fraction', 'DOUBLE', unit='frac', c={'pm': 'v144: FeO, a mass fraction of the rock', 'r': 'v207 (GEN.1111)'}),
    F('rock_al2o3_fraction', 'DOUBLE', unit='frac', c={'pm': 'v144: Al2O3, a mass fraction of the rock', 'r': 'v207 (GEN.1111)'}),
    F('rock_mg_si_ratio', 'DOUBLE', unit='mg-si', c={'pm': "v144: the rock's molar Mg/Si", 'r': 'v207 (GEN.1111)'}),
    F('rock_na2o_fraction', 'DOUBLE', unit='frac', c='v230 (GEN.1362): Na2O, a mass fraction of the rock'),
    F('rock_k2o_fraction', 'DOUBLE', unit='frac', c='v230 (GEN.1362): K2O, a mass fraction of the rock'),
    F('rock_tio2_fraction', 'DOUBLE', unit='frac', c='v230 (GEN.1362): TIO2, a mass fraction of the rock'),
    F('rock_cr2o3_fraction', 'DOUBLE', unit='frac', c='v230 (GEN.1362): CR2O3, a mass fraction of the rock'),
    F('rock_mno_fraction', 'DOUBLE', unit='frac', c='v230 (GEN.1362): MNO, a mass fraction of the rock'),
    F('rock_nio_fraction', 'DOUBLE', unit='frac', c='v230 (GEN.1362): NIO, a mass fraction of the rock'),
    F('rock_p2o5_fraction', 'DOUBLE', unit='frac', c='v230 (GEN.1362): P2O5, a mass fraction of the rock'),
    F('crust_p2o5_wt', 'DOUBLE', unit='pct', c="v233 (GEN.1367): The crust's P2O5, weight percent"),
    F('crust_apatite_share', 'DOUBLE', unit='crust-p-share', c="v233 (GEN.1367): Share of the crust's phosphorus in apatite"),
    F('crust_merrillite_share', 'DOUBLE', unit='crust-p-share', c='v233 (GEN.1367): Share in merrillite'),
    F('crust_schreibersite_share', 'DOUBLE', unit='crust-p-share', c='v233 (GEN.1367): Share in schreibersite'),
    F('ocean_phosphate_umol_l', 'DOUBLE', unit='umol/L', c="v233 (GEN.1367): The ocean's phosphate, umol/L"),
    F('crust_thickness_km', 'DOUBLE', unit='km', c='v234 (GEN.1380): the mean crust thickness, km'),
    F('craton_fraction', 'DOUBLE', unit='frac', c='v234 (GEN.1380): the share of the crust that is old thick craton, 0 to 1'),
    F('crust_sulfur_index', 'DOUBLE', unit='frac', c='v234 (GEN.1380): how much sulfur the crust offers a rising mafic magma, 0 to 1'),
    F('envelope_helium_fraction', 'DOUBLE', unit='frac', c={'pm': 'v144: Y, a mass fraction of the envelope', 'r': 'v207 (GEN.1111)'}),
    F('envelope_z_fraction', 'DOUBLE', unit='frac', c={'pm': 'v144: heavy elements, a mass fraction of the envelope', 'r': 'v207 (GEN.1111)'}),
    F('porosity', 'DOUBLE', unit='frac', c={'pm': "v144: the void share of a small body's volume", 'r': 'v207 (GEN.1111)'}),
    F('core_ni_fraction', 'DOUBLE', unit='frac', c='v231 (GEN.1364): Ni, a mass fraction of the core, NULL for a body with no metal and for a row made before v230'),
    F('core_light_fraction', 'DOUBLE', unit='frac', c='v231 (GEN.1364): the light elements together (S, Si, O, C, H), a mass fraction of the core, NULL for a body with no metal and for a row made before v230'),
    F('core_s_fraction', 'DOUBLE', unit='frac', c='v231 (GEN.1364): S, a mass fraction of the core, NULL for a body with no metal and for a row made before v230'),
    F('core_si_fraction', 'DOUBLE', unit='frac', c='v231 (GEN.1364): Si, a mass fraction of the core, NULL for a body with no metal and for a row made before v230'),
    F('core_o_fraction', 'DOUBLE', unit='frac', c='v231 (GEN.1364): O, a mass fraction of the core, NULL for a body with no metal and for a row made before v230'),
    F('core_c_fraction', 'DOUBLE', unit='frac', c='v231 (GEN.1364): C, a mass fraction of the core, NULL for a body with no metal and for a row made before v230'),
    F('core_h_fraction', 'DOUBLE', unit='frac', c='v231 (GEN.1364): H, a mass fraction of the core, NULL for a body with no metal and for a row made before v230'),
    F('carbon_rock', 'BOOLEAN', unit='flag', c={'pm': 'v144: the rock is silicon carbide and graphite, not silicate', 'r': 'v207 (GEN.1111)'}),
    F('composition_origin', 'VARCHAR(8)', extra="CHECK (composition_origin IN ('drawn', 'inverse'))", c={'pm': 'v144: drawn (composition first) or inverse (worked back from a mass and a radius)', 'r': 'v207 (GEN.1111)'}),
    F('ground_type', 'VARCHAR(16)', c={'pm': 'v159 (GEN.441): the ground type with the largest areal share (ocean, ice or a rock type)', 'r': 'v207 (GEN.1111)'}),
    F('ground_grain', 'VARCHAR(8)', c={'pm': 'v159: slab, crushed or powder', 'r': 'v207 (GEN.1111)'}),
    F('ground_ocean_fraction', 'DOUBLE', unit='frac', c={'pm': 'v159: the areal share of the surface that is ocean (the twelve shares add up to 1)', 'r': 'v207 (GEN.1111)'}),
    F('ground_ice_fraction', 'DOUBLE', unit='frac', c={'pm': 'v159: the areal share of the surface that is ice', 'r': 'v207 (GEN.1111)'}),
    F('ground_ultramafic_fraction', 'DOUBLE', unit='frac', c={'pm': 'v159: the areal share of the surface that is ultramafic', 'r': 'v207 (GEN.1111)'}),
    F('ground_anorthositic_fraction', 'DOUBLE', unit='frac', c={'pm': 'v159: the areal share of the surface that is anorthositic', 'r': 'v207 (GEN.1111)'}),
    F('ground_basaltic_fraction', 'DOUBLE', unit='frac', c={'pm': 'v159: the areal share of the surface that is basaltic', 'r': 'v207 (GEN.1111)'}),
    F('ground_felsic_fraction', 'DOUBLE', unit='frac', c={'pm': 'v159: the areal share of the surface that is felsic', 'r': 'v207 (GEN.1111)'}),
    F('ground_clay_fraction', 'DOUBLE', unit='frac', c={'pm': 'v159: the areal share of the surface that is clay', 'r': 'v207 (GEN.1111)'}),
    F('ground_oxidized_fraction', 'DOUBLE', unit='frac', c={'pm': 'v159: the areal share of the surface that is oxidized', 'r': 'v207 (GEN.1111)'}),
    F('ground_metal_fraction', 'DOUBLE', unit='frac', c={'pm': 'v159: the areal share of the surface that is metal', 'r': 'v207 (GEN.1111)'}),
    F('ground_carbon_fraction', 'DOUBLE', unit='frac', c={'pm': 'v159: the areal share of the surface that is carbon', 'r': 'v207 (GEN.1111)'}),
    F('ground_evaporite_fraction', 'DOUBLE', unit='frac', c={'pm': 'v159: the areal share of the surface that is evaporite', 'r': 'v207 (GEN.1111)'}),
    F('ground_sulfur_fraction', 'DOUBLE', unit='frac', c={'pm': 'v159: the areal share of the surface that is sulfur', 'r': 'v207 (GEN.1111)'}),
    F('ground_intermediate_fraction', 'DOUBLE', unit='frac', c={'pm': 'v190 (GEN.474): the areal share of the surface that is intermediate (andesitic) ground', 'r': 'v207 (GEN.1111)'}),
    F('surface_liquid', 'VARCHAR(32)', c={'pm': 'v188 (GEN.1035): the liquid the surface holds (an id of physics/surface_liquids.py, NULL: none)', 'r': 'v207 (GEN.1111)'}),
    F('surface_liquid_fraction', 'DOUBLE', unit='frac', c={'pm': 'v188: the share of the surface that liquid covers', 'r': 'v207 (GEN.1111)'}),
    F('surface_liquid_depth_km', 'DOUBLE', unit='km', c={'pm': 'v188: its mean depth over that share, km (NULL: not known)', 'r': 'v207 (GEN.1111)'}),
    F('surface_liquid_acid_wt', 'DOUBLE', unit='pct', c={'pm': 'v197 (GEN.1215): the strength of the liquid sulfuric acid the air makes, percent of its mass (NULL: none)', 'r': 'v207 (GEN.1111)'}),
    F('surface_liquid_acid_layer', 'VARCHAR(8)', c={'pm': 'v197: where that acid is, pool (at the ground) or cloud', 'r': 'v207 (GEN.1111)'}),
    F('surface_liquid_ammonia_wt', 'DOUBLE', unit='pct', c={'pm': "v197: the ammonia of a liquid water ocean, percent of the liquid's mass (NULL: no ocean or no ammonia)", 'r': 'v207 (GEN.1111)'}),
    F('surface_liquid_methane_mol', 'DOUBLE', unit='frac', c={'pm': 'v205 (GEN.1218b): the methane of a methane, ethane or nitrogen sea, mole fraction (NULL: another liquid)', 'r': 'v207 (GEN.1111)'}),
    F('surface_liquid_ethane_mol', 'DOUBLE', unit='frac', c={'pm': 'v212: its ethane, mole fraction', 'r': 'v207 (GEN.1111)'}),
    F('surface_liquid_nitrogen_mol', 'DOUBLE', unit='frac', c={'pm': 'v212: its dissolved nitrogen, mole fraction', 'r': 'v207 (GEN.1111)'}),
    F('surface_liquid_salt_wt', 'DOUBLE', unit='pct', c={'pm': "v212: the salt of a water sea as the sodium chloride that gives its water activity, percent of the solution's mass (NULL: not a water sea)", 'r': 'v207 (GEN.1111)'}),
    F('surface_liquid_2', 'VARCHAR(32)', c={'pm': 'v212: a second liquid beside the main one (an id of physics/surface_liquids.py, NULL: none)', 'r': 'v207 (GEN.1111)'}),
    F('surface_liquid_2_fraction', 'DOUBLE', unit='frac', c={'pm': 'v212: the share of the surface it covers', 'r': 'v207 (GEN.1111)'}),
    F('moon_origin', 'VARCHAR(8)', c={'pm': 'v180 (DB.70): how a moon came to orbit its parent (NULL for a planet and for a moon saved before)', 'r': 'v207 (GEN.1111)'}),
    F('moon_disk_lifetime_myr', 'DOUBLE', unit='Myr', c={'pm': "v180: a giant planet's circumplanetary disk: its lifetime, Myr (NULL for a body with no disk)", 'r': 'v207 (GEN.1111)'}),
    F('moon_disk_alpha_over_f', 'DOUBLE', unit='disk-alpha-f', c={'pm': "v180: the disk's alpha over f (the viscosity over the supply ratio of Canup and Ward)", 'r': 'v207 (GEN.1111)'}),
    F('moon_disk_f_inf', 'DOUBLE', unit='frac', c={'pm': "v180: the asymptotic ice fraction of the disk's ramp", 'r': 'v207 (GEN.1111)'}),
    F('moon_disk_ice_line_km', 'DOUBLE', unit='km', c={'pm': "v180: the distance of the disk's ice line from the giant, km", 'r': 'v207 (GEN.1111)'}),
    F('water_base_pressure_pa', 'DOUBLE', unit='Pa', c={'p': 'v151 (GEN.461): the pressure at the bottom of the water layer, Pa', 'm': 'v147 (GEN.461): the pressure at the bottom of the water layer, Pa', 'r': 'v207 (GEN.1111)'}),
    F('tidal_flux_w_m2', 'DOUBLE', unit='W_m2', c={'p': 'v152 (GEN.700): the eccentricity-tide heat flow of an icy moon under ice, W/m2 (NULL for a planet)', 'm': 'v152 (GEN.700): the eccentricity-tide heat flow of an icy moon under ice, W/m2', 'r': 'v207 (GEN.1111)'}),
    F('aquifer_depth_km', 'DOUBLE', unit='km', c={'pm': 'v157 (GEN.654, DB.92): the depth to the top of the water-bearing layer underground, km (NULL: none)', 'r': 'v207 (GEN.1111)'}),
    F('aquifer_thickness_km', 'DOUBLE', unit='km', c={'pm': 'v156: the thickness of that layer, km', 'r': 'v207 (GEN.1111)'}),
    F('mantle_eta0_pa_s', 'DOUBLE', unit='Pa*s', c={'pm': 'v153 (GEN.699): the viscosity prefactor of the mantle, eta = eta0 exp(Ea / R T), Pa s', 'r': 'v207 (GEN.1111)'}),
    F('mantle_activation_kj_mol', 'DOUBLE', unit='kJ/mol', c={'pm': 'v153: its viscosity activation energy Ea, kJ/mol', 'r': 'v207 (GEN.1111)'}),
    F('internal_heat_flux_w_m2', 'DOUBLE', unit='W_m2', c={'pm': "v153: the heat flowing out of the interior, radiogenic plus primordial, over the body's surface", 'r': None}),
    F('effective_temperature_k', 'DOUBLE', unit='K', t='r'),
    F('surface_regime', 'VARCHAR(24)', t='r'),
    F('interior_kind', 'VARCHAR(16)', extra="CHECK (interior_kind IN ('gas_giant', 'ice_giant', 'brown_dwarf'))", c={'pm': 'v154 (DB.89, GEN.642): the kind of interior a giant was drawn with', 'r': 'v207 (GEN.1111)'}),
    F('core_mass_earth', 'DOUBLE', unit='Mearth', c={'pm': "v154: the compact core's (an ice giant's rock core's) mass in Earth masses, 0 with none", 'r': 'v207 (GEN.1111)'}),
    F('dilute_core_mass_earth', 'DOUBLE', unit='Mearth', c={'pm': "v154: the dilute core's mass in Earth masses, 0 with none", 'r': 'v207 (GEN.1111)'}),
    F('metallization_pressure_gpa', 'DOUBLE', unit='GPa', c={'pm': 'v154: the drawn molecular-to-metallic pressure, GPa', 'r': 'v207 (GEN.1111)'}),
    F('helium_depletion', 'DOUBLE', unit='frac', c={'pm': "v154: the molecular envelope's helium share against the natal share, 1 with no helium rain", 'r': 'v207 (GEN.1111)'}),
    F('has_helium_shell', 'BOOLEAN', unit='flag', c={'pm': 'v154: a helium-rich shell above the metallic envelope', 'r': 'v207 (GEN.1111)'}),
    F('ring_kind', 'VARCHAR(8)', extra="CHECK (ring_kind IN ('dusty', 'narrow', 'broad'))", c={'pm': 'v160 (GEN.732): the kind of the ring system, NULL with none', 'r': None}, lead={'pm': (), 'r': ('-- v171 (GEN.836): the ring system, as on `planets` (v160); NULL for a planet with none.',)}),
    F('ring_origin', 'VARCHAR(16)', extra="CHECK (ring_origin IN ('moon_erosion', 'disrupted_moon', 'captured_body'))", c={'pm': 'v160: where its material came from', 'r': None}),
    F('ring_inner_km', 'DOUBLE', unit='km', c={'pm': "v160: the inner edge of the ring system, km from the body's center", 'r': None}),
    F('ring_outer_km', 'DOUBLE', unit='km', c={'pm': 'v160: its outer edge, km (inside the Roche limit of its own material and the innermost moon)', 'r': None}),
    F('ring_mass_kg', 'DOUBLE', unit='kg', c={'pm': 'v160: its mass, kg', 'r': None}),
    F('ring_optical_depth', 'DOUBLE', unit='tau', c={'pm': 'v160: its mean optical depth over the area it covers', 'r': None}),
    F('ring_ice_fraction', 'DOUBLE', unit='frac', c={'pm': 'v160: its water-ice share by mass', 'r': None}),
    F('ring_age_yr', 'DOUBLE', unit='yr', c={'pm': 'v192 (GEN.1057): years since the event that made it (a fallen moon or a passing comet), NULL for a resonance ring', 'r': 'v192 (GEN.1057): as on `planets`'}),
    F('galactic_orbital_period_gy', 'DOUBLE', nn=True, unit='Gyr', lead=('-- v17: always populated (a rogue planet is always standalone).',), t='r'),
    F('galactic_orbital_phase_deg', 'DOUBLE', nn=True, unit='deg', t='r'),
    F('galactic_min_update_interval_years', 'DOUBLE', nn=True, unit='yr', t='r'),
    F('atmosphere_colour', 'INT UNSIGNED', unit='rgb24', c={'p': "v162 (DB.96): the colour of the body's air, sRGB as red*65536 + green*256 + blue, NULL with no air", 'm': "v161 (DB.96): the colour of the body's air, sRGB as red*65536 + green*256 + blue, NULL with no air", 'r': 'v207 (GEN.1111)'}),
    F('lid_regime', 'VARCHAR(10)', c={'pm': 'v169 (GEN.721): mobile (plate tectonics) or stagnant (one rigid plate), drawn once, NULL for a gas giant', 'r': 'v207 (GEN.1111)'}),
    F('day_share_equator', 'DOUBLE', unit='frac', c={'pm': 'v170 (GEN.767): share of the time the star is up at the equator over a year, 0 to 1, NULL for a world that turns once an orbit', 'r': 'v207 (GEN.1111)'}),
    F('day_share_mean', 'DOUBLE', unit='frac', c={'pm': 'v170 (GEN.767): share of the time the star is up, averaged over the surface, 0 to 1', 'r': 'v207 (GEN.1111)'}),
    F('parent_shadow_share', 'DOUBLE', unit='frac', c={'pm': "v170 (GEN.767): a moon only, the share of its time in its planet's shadow, 0 to 1", 'r': 'v207 (GEN.1111)'}),
    F('solar_day_hours', 'DOUBLE', unit='h', c={'pm': 'v170 (GEN.767): hours from one noon to the next, NULL for a world that turns once an orbit', 'r': 'v207 (GEN.1111)'}),
    F('daylight_hours', 'DOUBLE', unit='h', c={'pm': 'v170 (GEN.767): hours of daylight in one solar day at the equator, NULL for a world that turns once an orbit', 'r': 'v207 (GEN.1111)'}),
    F('polar_circle_deg', 'DOUBLE', unit='deg', c={'pm': 'v170 (GEN.767): latitude beyond which the star stays up or down for at least a day, degrees', 'r': 'v207 (GEN.1111)'}),
    F('life_class', 'TINYINT UNSIGNED', unit='index', c={'pm': 'v175 (GEN.889): the number of the life class that lives on the body (physics/life_classes.py), NULL for none', 'r': 'v207 (GEN.1111)'}),
    F('life_class_variant', 'CHAR(1)', c={'pm': 'v175: its pigment or habitat letter (1a, 1b, 1c, 3a ...), NULL for a class with none', 'r': 'v207 (GEN.1111)'}),
    F('life_class_best', 'VARCHAR(2)', c={'pm': 'v175: the id of the best-supported class, NULL when no class is supported', 'r': 'v207 (GEN.1111)'}),
    F('life_class_best_score', 'TINYINT UNSIGNED', unit='score9', c={'pm': 'v175: its score, 0 to 9', 'r': 'v207 (GEN.1111)'}),
    F('life_animal_grade', 'TINYINT UNSIGNED', unit='grade3', c={'pm': 'v175 (GEN.890): 0 no animals, 1 sessile, 2 small and mobile, 3 large and active', 'r': 'v207 (GEN.1111)'}),
    F('life_complexity', 'DOUBLE', unit='complexity', c={'pm': 'v175: the complexity index of its life class, 1 for an Earth-like biosphere', 'r': 'v207 (GEN.1111)'}),
    F('life_traits', 'VARCHAR(48)', c={'pm': 'v223 (GEN.876): the adaptation traits of its life, text such as radiation:2,sulfide:3 (stress and grade), NULL for none', 'r': 'v223 (GEN.876)'}),
    F('cloud_chemistry', 'VARCHAR(48)', c='v245 (GEN.1425): the species of the Venus sulfur-cloud set its air gives it, ids such as V1,V2,V10 (`physics/sulfur_cloud_set.py`), NULL for none'),
    F('tide_range_spring_m', 'DOUBLE', unit='m', c={'pm': 'v191 (GEN.1034): open-ocean equilibrium spring tide range, metres, 0 under an ice cover, NULL without a sea', 'r': 'v207 (GEN.1111)'}),
    F('tide_range_neap_m', 'DOUBLE', unit='m', c={'pm': 'v191 (GEN.1034): the neap range, metres', 'r': 'v207 (GEN.1111)'}),
    F('ocean_liquid', 'VARCHAR(40)', c={'pm': 'v191 (GEN.1034): the liquid of the sea the tide is of (water until GEN.1035), NULL without a sea', 'r': 'v207 (GEN.1111)'}),
    F('donor_uid', 'BINARY(12)', c="v207 (GEN.1111): the ejecting system's uid, the event, the ejection, the note", t='r'),
    F('event_index', 'INT', unit='index', c="v207 (GEN.1111): the ejecting system's uid, the event, the ejection, the note", t='r'),
    F('ejection_delay_gy', 'DOUBLE', unit='Gyr', c="v207 (GEN.1111): the ejecting system's uid, the event, the ejection, the note", t='r'),
    F('since_ejection_gy', 'DOUBLE', unit='Gyr', c="v207 (GEN.1111): the ejecting system's uid, the event, the ejection, the note", t='r'),
    F('ejection_speed_kms', 'DOUBLE', unit='km_s', c="v207 (GEN.1111): the ejecting system's uid, the event, the ejection, the note", t='r'),
    F('rogue_note', 'TEXT', c="v207 (GEN.1111): the ejecting system's uid, the event, the ejection, the note", t='r'),
)
"""tuple: Every column of the body tables, in the order a new database holds them."""

ORDERS = {
    'r': (
        'id', 'sector_id', 'name', 'planet_type', 'planet_class', 'mass_bin',
        'mass_kg', 'radius_km', 'composition', 'has_internal_heat', 'age_gy', 'internal_heat_flux_w_m2',
        'effective_temperature_k', 'surface_regime', 'surface_temperature_k', 'surface_pressure_pa', 'ice_shell_thickness_km', 'ocean_depth_km',
        'has_liquid_water', 'hp_ice_km', 'ring_kind', 'ring_origin', 'ring_inner_km', 'ring_outer_km',
        'ring_mass_kg', 'ring_optical_depth', 'ring_ice_fraction', 'ring_age_yr', 'galactic_orbital_period_gy', 'galactic_orbital_phase_deg',
        'galactic_min_update_interval_years', 'epoch_unix', 'next_update_due', 'rotation_period_hours', 'spin_axis_x', 'spin_axis_y',
        'spin_axis_z', 'axial_tilt_deg', 'center_x_pc', 'center_y_pc', 'center_z_pc', 'cell_key',
        'velocity_x_kms', 'velocity_y_kms', 'velocity_z_kms', 'galactic_radius_pc', 'quadrant', 'inside_nebula_id',
        'inside_remnant_id', 'created_at', 'modified_at', 'uid', 'parent_type', 'parent_id',
        'phi4', 'phi4_pressure', 'phi4_temperature', 'phi4_air', 'phi4_water', 'phi4_radiation',
        'phi4_suitability', 'tier_pressure', 'tier_temperature', 'tier_air', 'tier_water', 'tier_radiation',
        'tier_suitability', 'life_microbial', 'life_flora', 'life_fauna', 'life_sentient', 'l_solv',
        'l_chem', 'l_ener', 'l_rad', 'eq_body', 'eq_body_kind', 'eq_body_rec',
        'eq_body_rec_kind', 'eq_rad_suit', 'eq_rad_suit_rec', 'eq_mobility', 'eq_dust', 'eq_gas_filter',
        'eq_co2_scrub', 'eq_supply', 'eq_supply_rec', 'eq_shelter', 'eq_shelter_rec', 'density_g_cm3',
        'roundness_ratio', 'yield_strength_pa', 'binding_energy_j', 'shelter_column_g_cm2', 'energy_flux_w_m2', 'body_type',
        'volume_km3', 'gravity_g', 'atmosphere', 'atm_density', 'atm_molar_density', 'atmospheric_pressure_pa',
        'scale_height_km', 'flavor_text_count', 'life_stage', 'evolutionary_pace', 'life_milestone_age_gyr', 'life_system_age_gyr',
        'life_age_abiogenesis_gyr', 'life_age_photosynthesis_gyr', 'life_age_complex_cells_gyr', 'life_age_multicellularity_gyr', 'life_age_technological_civilization_gyr', 'flavor_key',
        'mantle_redox', 'mantle_delta_iw', 'magnetic_moment_a_m2', 'dipole_class', 'dipole_tilt_deg', 'dipole_offset_class',
        'induced_moment_a_m2', 'magnetopause_rp', 'xuv_flux_earth', 'xuv_exposure_index', 'flare_irradiation_index', 'water_mass_fraction',
        'hydrosphere', 'ocean_fraction', 'land_fraction', 'ice_cover_fraction', 'shape_temperature_k', 'ice_shell_km',
        'ocean_class', 'ocean_ph', 'water_activity', 'phosphorus', 'surface_dose_msv_yr', 'dose_gcr_msv_yr',
        'dose_sep_msv_yr', 'dose_ground_msv_yr', 'dose_belt_msv_yr', 'dose_helio_mult', 'uv_surface_index', 'ozone_loss_flag',
        'outside_hours', 'band1_mean_k', 'band1_swing_k', 'band1_day_hours', 'band1_night_hours', 'band1_longest_day_days',
        'band1_longest_night_days', 'band1_state', 'band2_mean_k', 'band2_swing_k', 'band2_day_hours', 'band2_night_hours',
        'band2_longest_day_days', 'band2_longest_night_days', 'band2_state', 'band3_mean_k', 'band3_swing_k', 'band3_day_hours',
        'band3_night_hours', 'band3_longest_day_days', 'band3_longest_night_days', 'band3_state', 'band4_mean_k', 'band4_swing_k',
        'band4_day_hours', 'band4_night_hours', 'band4_longest_day_days', 'band4_longest_night_days', 'band4_state', 'band5_mean_k',
        'band5_swing_k', 'band5_day_hours', 'band5_night_hours', 'band5_longest_day_days', 'band5_longest_night_days', 'band5_state',
        'band6_mean_k', 'band6_swing_k', 'band6_day_hours', 'band6_night_hours', 'band6_longest_day_days', 'band6_longest_night_days',
        'band6_state', 'band7_mean_k', 'band7_swing_k', 'band7_day_hours', 'band7_night_hours', 'band7_longest_day_days',
        'band7_longest_night_days', 'band7_state', 'temp_mean_swing_k', 'temp_livable_share', 'temp_livable_area', 'temp_livable_time',
        'temp_blue_share', 'temp_green_share', 'temp_cold_share', 'temp_hot_share', 'temp_area_always', 'temp_area_most',
        'temp_livable_belt_only', 'temp_coldest_k', 'temp_coldest_lat_deg', 'temp_coldest_lon_deg', 'temp_coldest_day', 'temp_hottest_k',
        'temp_hottest_lat_deg', 'temp_hottest_lon_deg', 'temp_hottest_day', 'temp_site_lat_deg', 'temp_site_lon_deg', 'temp_site_livable_time',
        'temp_site_coldest_k', 'temp_site_hottest_k', 'temp_site_mean_k', 'temp_site_longest_out_days', 'outside_hours_day', 'outside_hours_limit',
        'day_night_state', 'night_hours', 'days_per_year', 'flux_variation', 'ozone_column_du', 'tau_haze',
        'tau_dust', 'tau_cloud', 'cloud_cover', 'aerosol_draw', 'c_k0_ppm', 'c_th0_ppm',
        'c_u0_ppm', 'al26_ratio', 'urey_ratio', 'bond_albedo', 'greenhouse_factor', 'albedo_basis',
        'metal_mass_fraction', 'rock_mass_fraction', 'ice_mass_fraction', 'ice_co2_mass_fraction', 'ice_ch4_mass_fraction', 'ice_nh3_mass_fraction',
        'ice_n2_mass_fraction', 'ice_so2_mass_fraction', 'ice_co_mass_fraction', 'ice_hcn_mass_fraction', 'ice_other_mass_fraction', 'envelope_mass_fraction',
        'rock_feo_fraction', 'rock_al2o3_fraction', 'rock_mg_si_ratio', 'rock_na2o_fraction', 'rock_k2o_fraction', 'rock_tio2_fraction',
        'rock_cr2o3_fraction', 'rock_mno_fraction', 'rock_nio_fraction', 'rock_p2o5_fraction', 'crust_p2o5_wt', 'crust_apatite_share',
        'crust_merrillite_share', 'crust_schreibersite_share', 'ocean_phosphate_umol_l', 'crust_thickness_km', 'craton_fraction', 'crust_sulfur_index',
        'envelope_helium_fraction', 'envelope_z_fraction', 'porosity', 'core_ni_fraction', 'core_light_fraction', 'core_s_fraction',
        'core_si_fraction', 'core_o_fraction', 'core_c_fraction', 'core_h_fraction', 'carbon_rock', 'composition_origin',
        'ground_type', 'ground_grain', 'ground_ocean_fraction', 'ground_ice_fraction', 'ground_ultramafic_fraction', 'ground_anorthositic_fraction',
        'ground_basaltic_fraction', 'ground_felsic_fraction', 'ground_clay_fraction', 'ground_oxidized_fraction', 'ground_metal_fraction', 'ground_carbon_fraction',
        'ground_evaporite_fraction', 'ground_sulfur_fraction', 'ground_intermediate_fraction', 'surface_liquid', 'surface_liquid_fraction', 'surface_liquid_depth_km',
        'surface_liquid_acid_wt', 'surface_liquid_acid_layer', 'surface_liquid_ammonia_wt', 'surface_liquid_methane_mol', 'surface_liquid_ethane_mol', 'surface_liquid_nitrogen_mol',
        'surface_liquid_salt_wt', 'surface_liquid_2', 'surface_liquid_2_fraction', 'class_variant', 'relative_humidity', 'groundwater',
        'class_match', 'runner_up_class', 'runner_up_match', 'class_fudge_values', 'class_fudge_amount', 'oddities',
        'sky_suns', 'light_factor', 'obliquity_min_deg', 'obliquity_max_deg', 'moon_origin', 'moon_disk_lifetime_myr',
        'moon_disk_alpha_over_f', 'moon_disk_f_inf', 'moon_disk_ice_line_km', 'water_base_pressure_pa', 'tidal_flux_w_m2', 'aquifer_depth_km',
        'aquifer_thickness_km', 'mantle_eta0_pa_s', 'mantle_activation_kj_mol', 'interior_kind', 'core_mass_earth', 'dilute_core_mass_earth',
        'metallization_pressure_gpa', 'helium_depletion', 'has_helium_shell', 'atmosphere_colour', 'lid_regime', 'day_share_equator',
        'day_share_mean', 'parent_shadow_share', 'solar_day_hours', 'daylight_hours', 'polar_circle_deg', 'life_class',
        'life_class_variant', 'life_class_best', 'life_class_best_score', 'life_animal_grade', 'life_complexity', 'life_traits',
        'cloud_chemistry', 'tide_range_spring_m', 'tide_range_neap_m', 'ocean_liquid', 'donor_uid', 'event_index',
        'ejection_delay_gy', 'since_ejection_gy', 'ejection_speed_kms', 'rogue_note',
    ),
}
"""dict: A table whose columns are not in the order of `FIELDS`, with its own order. `rogue_planets` keeps the order its database holds, because MySQL numbers a column's inline CHECK by its place (`rogue_planets_chk_3`) and an Alembic revision names the constraint by that number; reordering the columns of a new database would rename them."""

CLAUSES = {
    'p': (
        '',
        '    UNIQUE KEY uq_planets_uid (uid),',
        '    CONSTRAINT fk_planets_star_system',
        '        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,',
        '    CONSTRAINT fk_planets_star',
        '        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE SET NULL,',
        '    KEY idx_planets_star_system_id (star_system_id),',
        '    KEY idx_planets_next_update_due (next_update_due),',
        '    KEY idx_planets_star_id (star_id),',
        '    KEY idx_planets_name (name),',
        '    FULLTEXT KEY ft_planets_name (name),',
        '    KEY idx_planets_planet_class (planet_class),',
        '    KEY idx_planets_body_type (body_type),',
        '    KEY idx_planets_life_class (life_class),',
        '    KEY idx_planets_radius_km (radius_km),',
        "    CONSTRAINT chk_planets_parent_type CHECK (parent_type IN ('star', 'planet', 'system', 'node', 'rogue_planet', 'free')),",
        '    KEY idx_planets_parent (parent_type, parent_id),',
        "    CONSTRAINT chk_planets_dynamical_state CHECK (dynamical_state IN ('bound', 'unbound'))",
    ),
    'm': (
        '',
        '    UNIQUE KEY uq_moons_uid (uid),',
        '    CONSTRAINT fk_moons_planet',
        '        FOREIGN KEY (planet_id) REFERENCES planets(id) ON DELETE CASCADE,',
        '    CONSTRAINT fk_moons_star_system',
        '        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,',
        '    CONSTRAINT fk_moons_star',
        '        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE SET NULL,',
        '    CONSTRAINT fk_moons_rogue_planet',
        '        FOREIGN KEY (rogue_planet_id) REFERENCES rogue_planets(id) ON DELETE CASCADE,',
        '    KEY idx_moons_rogue_planet_id (rogue_planet_id),',
        '    KEY idx_moons_planet_id (planet_id),',
        '    KEY idx_moons_next_update_due (next_update_due),',
        '    KEY idx_moons_star_system_id (star_system_id),',
        '    KEY idx_moons_star_id (star_id),',
        '    KEY idx_moons_name (name),',
        '    FULLTEXT KEY ft_moons_name (name),',
        '    KEY idx_moons_planet_class (planet_class),',
        '    KEY idx_moons_body_type (body_type),',
        '    KEY idx_moons_life_class (life_class),',
        '    KEY idx_moons_radius_km (radius_km),',
        "    CONSTRAINT chk_moons_parent_type CHECK (parent_type IN ('star', 'planet', 'system', 'node', 'rogue_planet', 'free')),",
        '    KEY idx_moons_parent (parent_type, parent_id),',
        "    CONSTRAINT chk_moons_dynamical_state CHECK (dynamical_state IN ('bound', 'unbound'))",
    ),
    'r': (
        '',
        '    UNIQUE KEY uq_rogue_planets_uid (uid),',
        "    UNIQUE KEY uq_rogue_planets_donor_event (donor_uid, event_index),  -- v207 (GEN.1111): a donor's rogues are replaced, never duplicated",
        '    KEY idx_rogue_planets_parent (parent_type, parent_id),',
        "    CONSTRAINT chk_rogue_planets_parent_type CHECK (parent_type IN ('star', 'planet', 'system', 'node', 'rogue_planet', 'free')),",
        '    KEY idx_rogue_planets_name (name),  -- v40: name uniqueness lookups',
        "    KEY idx_rogue_planets_placed_name (((center_x_pc IS NOT NULL)), name),   -- v186 (PERF.211): the Phenomena table's sorts read in index order",
        "    KEY idx_rogue_planets_descriptor_name (((CASE WHEN mass_bin = 'brown-dwarf' THEN 'brown dwarf' WHEN planet_type = 'gas' THEN 'gas giant' ELSE 'terrestrial' END) COLLATE utf8mb4_unicode_ci), name),",
        '    KEY idx_rogue_planets_next_update_due (next_update_due),',
        '    -- v28: named explicitly -- see `nebulae`\'s identical "v18" CHECK comment.',
        '    CONSTRAINT chk_rogue_planets_placement CHECK (',
        '        (center_x_pc IS NULL) = (center_y_pc IS NULL) AND',
        '        (center_y_pc IS NULL) = (center_z_pc IS NULL) AND',
        '        (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)',
        '    ),',
        '',
        '    CONSTRAINT fk_rogue_planets_sector',
        '        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE RESTRICT,  -- v201 (GEN.1127): a sector row cannot go while a resident points at it',
        '    KEY idx_rogue_planets_sector_id (sector_id),',
        '    -- v28: see the header comment\'s "v28" note.',
        '    KEY idx_rogue_planets_galactic_radius_pc (galactic_radius_pc),',
        '    KEY idx_rogue_planets_center (center_x_pc, center_y_pc, center_z_pc),',
        '    -- v27: see the header comment\'s "v27" note.',
        '    KEY idx_rogue_planets_modified_at (modified_at),',
        '    -- v39: containment -- see the header comment\'s "v39" note.',
        '    CONSTRAINT fk_rogue_planets_inside_nebula',
        '        FOREIGN KEY (inside_nebula_id) REFERENCES nebulae(id) ON DELETE SET NULL,',
        '    CONSTRAINT fk_rogue_planets_inside_remnant',
        '        FOREIGN KEY (inside_remnant_id) REFERENCES supernova_remnants(id) ON DELETE SET NULL',
    ),
}
"""dict: The lines after the columns of each table: the keys, the constraints and what comes between them, as `schema.sql` holds them."""

TRAILER = ") ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;"
"""str: The end of each table's statement."""

BEGIN = "-- BEGIN GENERATED body tables (db/body_fields.py, scripts/build_body_tables.py): do not edit by hand"
END = "-- END GENERATED body tables"
"""str: The marker lines around the generated statements in `schema.sql`; each table is written between the two lines of its own (`BEGIN`/`END` plus the table's name)."""

NAME_WIDTH = 32
COMMENT_GAP = "  "


def fields_of(letter):
    """The columns of the table `letter` in table order."""
    fields = [field for field in FIELDS if letter in field.t]
    if letter in ORDERS:
        place = {name: index for index, name in enumerate(ORDERS[letter])}
        fields.sort(key=lambda field: place[field.name])
    return fields


def sql_of(field, letter):
    """The text after a column's name: its type, `NOT NULL` and the rest."""
    parts = [field.of("type", letter)]
    if field.of("nn", letter):
        parts.append("NOT NULL")
    extra = field.of("extra", letter)
    if extra:
        parts.append(extra)
    return " ".join(parts)


def table_sql(letter):
    """The `CREATE TABLE` statement of the table `letter`, as `schema.sql` holds it."""
    table = TABLES[letter]
    lines = [f"CREATE TABLE IF NOT EXISTS {table} ("]
    for field in fields_of(letter):
        lines.extend(("    " + line) if line else "" for line in field.of("lead", letter))
        text = sql_of(field, letter) + ","
        comment = field.of("c", letter)
        lines.append(f"    {field.name:<{NAME_WIDTH}} {text}" + (f"{COMMENT_GAP}-- {comment}" if comment else ""))
        lines.extend(field.of("post", letter))
    lines.extend(CLAUSES[letter])
    lines.append(TRAILER)
    return "\n".join(lines) + "\n"


def schema_block(letter):
    """The statement of the table `letter` with its marker lines, as `schema.sql` holds it."""
    table = TABLES[letter]
    return f"{BEGIN}: {table}\n{table_sql(letter)}{END}: {table}\n"


def units_by_group(table):
    """`{unit: columns}` of the table `table` in declaration order, the shape of `column_units.TABLES[table]`."""
    letter = next(key for key, name in TABLES.items() if name == table)
    groups = {}
    for field in fields_of(letter):
        unit = field.of("unit", letter)
        if unit is not None:
            groups.setdefault(unit, []).append(field.name)
    return {unit: tuple(names) for unit, names in groups.items()}


def field_units(table):
    """`{column: unit}` of the table `table`, the shape of a family's `FIELD_UNITS` (`generation/mass_object.py`)."""
    return {name: unit for unit, names in units_by_group(table).items() for name in names}
