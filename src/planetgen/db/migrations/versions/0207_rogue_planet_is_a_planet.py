"""A rogue planet is a planet with no star (GEN.1111).

Schema v207. `rogue_planets` gains every column a `planets` row has apart from the orbit round a host, so the one planet chain can write a rogue
(the same names and types; all NULL here, since the flat-rate rogue that wrote the rows before never had them and the regenerated galaxy replaces it),
and six columns of its own: `donor_uid` and `event_index` (the system that ejected it and which of that system's ejections it is; the pair is unique,
so saving a donor again replaces its rogues and never duplicates them), `ejection_delay_gy` (how long after the system formed it was ejected),
`since_ejection_gy` (how long ago), `ejection_speed_kms` (the kick) and `rogue_note` (the sentence that says where it came from).
Four child tables, `rogue_planet_reflection_spectrum`, `rogue_planet_interior_layers`, `rogue_planet_gases` and `rogue_planet_life_class_scores`,
hold what `planet_*` hold for a planet. The columns named for the flat-rate rogue (`planet_type`, `composition`, `has_internal_heat`, `has_moons`)
become NULLable: a row of the physics rogue leaves them NULL (the revision that removes the flat-rate generator drops them). A column or table
already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0207"
down_revision = "0206"
branch_labels = None
depends_on = None

COLUMNS = (
    ('body_type', "VARCHAR(4) CHECK (body_type IN ('t', 'g')) NULL"),
    ('volume_km3', 'DOUBLE NULL'),
    ('gravity_g', 'DOUBLE NULL'),
    ('atmosphere', 'TEXT NULL'),
    ('atm_density', 'DOUBLE NULL'),
    ('atm_molar_density', 'DOUBLE NULL'),
    ('atmospheric_pressure_pa', 'DOUBLE NULL'),
    ('scale_height_km', 'DOUBLE NULL'),
    ('flavor_text_count', 'INT NOT NULL DEFAULT 0'),
    ('life_stage', 'VARCHAR(32) NULL'),
    ('evolutionary_pace', 'VARCHAR(16) NULL'),
    ('life_milestone_age_gyr', 'DOUBLE NULL'),
    ('life_system_age_gyr', 'DOUBLE NULL'),
    ('life_age_abiogenesis_gyr', 'DOUBLE NULL'),
    ('life_age_photosynthesis_gyr', 'DOUBLE NULL'),
    ('life_age_complex_cells_gyr', 'DOUBLE NULL'),
    ('life_age_multicellularity_gyr', 'DOUBLE NULL'),
    ('life_age_technological_civilization_gyr', 'DOUBLE NULL'),
    ('flavor_key', 'VARCHAR(32) NULL'),
    ('mantle_redox', "VARCHAR(16) CHECK (mantle_redox IN ('reduced', 'intermediate', 'oxidized')) NULL"),
    ('mantle_delta_iw', 'DOUBLE NULL'),
    ('magnetic_moment_a_m2', 'DOUBLE NULL'),
    ('dipole_class', "VARCHAR(16) CHECK (dipole_class IN ('none', 'weak', 'earth-like', 'strong', 'multipolar')) NULL"),
    ('induced_moment_a_m2', 'DOUBLE NULL'),
    ('magnetopause_rp', 'DOUBLE NULL'),
    ('xuv_flux_earth', 'DOUBLE NULL'),
    ('xuv_exposure_index', 'DOUBLE NULL'),
    ('flare_irradiation_index', 'DOUBLE NULL'),
    ('water_mass_fraction', 'DOUBLE NULL'),
    ('hydrosphere', "VARCHAR(20) CHECK (hydrosphere IN ('dry', 'vapour', 'ice', 'ice-covered ocean', 'surface ocean', 'brine seep', 'hycean')) NULL"),
    ('ocean_fraction', 'DOUBLE NULL'),
    ('land_fraction', 'DOUBLE NULL'),
    ('ice_shell_km', 'DOUBLE NULL'),
    ('ocean_class', "VARCHAR(16) CHECK (ocean_class IN ('ice-sealed', 'chloride brine', 'acid sulfate', 'soda', 'neutral')) NULL"),
    ('ocean_ph', 'DOUBLE NULL'),
    ('water_activity', 'DOUBLE NULL'),
    ('phosphorus', "VARCHAR(8) CHECK (phosphorus IN ('high', 'limited', 'starved')) NULL"),
    ('surface_dose_msv_yr', 'DOUBLE NULL'),
    ('dose_gcr_msv_yr', 'DOUBLE NULL'),
    ('dose_sep_msv_yr', 'DOUBLE NULL'),
    ('dose_ground_msv_yr', 'DOUBLE NULL'),
    ('dose_belt_msv_yr', 'DOUBLE NULL'),
    ('dose_helio_mult', 'DOUBLE NULL'),
    ('uv_surface_index', 'DOUBLE NULL'),
    ('ozone_loss_flag', 'BOOLEAN NULL'),
    ('outside_hours', 'DOUBLE NULL'),
    ('band1_mean_k', 'DOUBLE NULL'),
    ('band1_swing_k', 'DOUBLE NULL'),
    ('band1_day_hours', 'DOUBLE NULL'),
    ('band1_night_hours', 'DOUBLE NULL'),
    ('band1_longest_day_days', 'DOUBLE NULL'),
    ('band1_longest_night_days', 'DOUBLE NULL'),
    ('band1_state', 'VARCHAR(24) NULL'),
    ('band2_mean_k', 'DOUBLE NULL'),
    ('band2_swing_k', 'DOUBLE NULL'),
    ('band2_day_hours', 'DOUBLE NULL'),
    ('band2_night_hours', 'DOUBLE NULL'),
    ('band2_longest_day_days', 'DOUBLE NULL'),
    ('band2_longest_night_days', 'DOUBLE NULL'),
    ('band2_state', 'VARCHAR(24) NULL'),
    ('band3_mean_k', 'DOUBLE NULL'),
    ('band3_swing_k', 'DOUBLE NULL'),
    ('band3_day_hours', 'DOUBLE NULL'),
    ('band3_night_hours', 'DOUBLE NULL'),
    ('band3_longest_day_days', 'DOUBLE NULL'),
    ('band3_longest_night_days', 'DOUBLE NULL'),
    ('band3_state', 'VARCHAR(24) NULL'),
    ('band4_mean_k', 'DOUBLE NULL'),
    ('band4_swing_k', 'DOUBLE NULL'),
    ('band4_day_hours', 'DOUBLE NULL'),
    ('band4_night_hours', 'DOUBLE NULL'),
    ('band4_longest_day_days', 'DOUBLE NULL'),
    ('band4_longest_night_days', 'DOUBLE NULL'),
    ('band4_state', 'VARCHAR(24) NULL'),
    ('band5_mean_k', 'DOUBLE NULL'),
    ('band5_swing_k', 'DOUBLE NULL'),
    ('band5_day_hours', 'DOUBLE NULL'),
    ('band5_night_hours', 'DOUBLE NULL'),
    ('band5_longest_day_days', 'DOUBLE NULL'),
    ('band5_longest_night_days', 'DOUBLE NULL'),
    ('band5_state', 'VARCHAR(24) NULL'),
    ('band6_mean_k', 'DOUBLE NULL'),
    ('band6_swing_k', 'DOUBLE NULL'),
    ('band6_day_hours', 'DOUBLE NULL'),
    ('band6_night_hours', 'DOUBLE NULL'),
    ('band6_longest_day_days', 'DOUBLE NULL'),
    ('band6_longest_night_days', 'DOUBLE NULL'),
    ('band6_state', 'VARCHAR(24) NULL'),
    ('band7_mean_k', 'DOUBLE NULL'),
    ('band7_swing_k', 'DOUBLE NULL'),
    ('band7_day_hours', 'DOUBLE NULL'),
    ('band7_night_hours', 'DOUBLE NULL'),
    ('band7_longest_day_days', 'DOUBLE NULL'),
    ('band7_longest_night_days', 'DOUBLE NULL'),
    ('band7_state', 'VARCHAR(24) NULL'),
    ('temp_mean_swing_k', 'DOUBLE NULL'),
    ('outside_hours_day', 'DOUBLE NULL'),
    ('outside_hours_limit', 'VARCHAR(16) NULL'),
    ('day_night_state', 'VARCHAR(24) NULL'),
    ('night_hours', 'DOUBLE NULL'),
    ('days_per_year', 'DOUBLE NULL'),
    ('flux_variation', 'DOUBLE NULL'),
    ('ozone_column_du', 'DOUBLE NULL'),
    ('tau_haze', 'DOUBLE NULL'),
    ('tau_dust', 'DOUBLE NULL'),
    ('tau_cloud', 'DOUBLE NULL'),
    ('cloud_cover', 'DOUBLE NULL'),
    ('aerosol_draw', 'DOUBLE NULL'),
    ('c_k0_ppm', 'DOUBLE NULL'),
    ('c_th0_ppm', 'DOUBLE NULL'),
    ('c_u0_ppm', 'DOUBLE NULL'),
    ('al26_ratio', 'DOUBLE NULL'),
    ('urey_ratio', 'DOUBLE NULL'),
    ('bond_albedo', 'DOUBLE NULL'),
    ('greenhouse_factor', 'DOUBLE NULL'),
    ('albedo_basis', 'TINYINT UNSIGNED NULL'),
    ('metal_mass_fraction', 'DOUBLE NULL'),
    ('rock_mass_fraction', 'DOUBLE NULL'),
    ('ice_mass_fraction', 'DOUBLE NULL'),
    ('ice_co2_mass_fraction', 'DOUBLE NULL'),
    ('ice_ch4_mass_fraction', 'DOUBLE NULL'),
    ('ice_nh3_mass_fraction', 'DOUBLE NULL'),
    ('ice_n2_mass_fraction', 'DOUBLE NULL'),
    ('ice_so2_mass_fraction', 'DOUBLE NULL'),
    ('ice_co_mass_fraction', 'DOUBLE NULL'),
    ('ice_hcn_mass_fraction', 'DOUBLE NULL'),
    ('ice_other_mass_fraction', 'DOUBLE NULL'),
    ('envelope_mass_fraction', 'DOUBLE NULL'),
    ('rock_feo_fraction', 'DOUBLE NULL'),
    ('rock_al2o3_fraction', 'DOUBLE NULL'),
    ('rock_mg_si_ratio', 'DOUBLE NULL'),
    ('envelope_helium_fraction', 'DOUBLE NULL'),
    ('envelope_z_fraction', 'DOUBLE NULL'),
    ('porosity', 'DOUBLE NULL'),
    ('carbon_rock', 'BOOLEAN NULL'),
    ('composition_origin', "VARCHAR(8) CHECK (composition_origin IN ('drawn', 'inverse')) NULL"),
    ('ground_type', 'VARCHAR(16) NULL'),
    ('ground_grain', 'VARCHAR(8) NULL'),
    ('ground_ocean_fraction', 'DOUBLE NULL'),
    ('ground_ice_fraction', 'DOUBLE NULL'),
    ('ground_ultramafic_fraction', 'DOUBLE NULL'),
    ('ground_anorthositic_fraction', 'DOUBLE NULL'),
    ('ground_basaltic_fraction', 'DOUBLE NULL'),
    ('ground_felsic_fraction', 'DOUBLE NULL'),
    ('ground_clay_fraction', 'DOUBLE NULL'),
    ('ground_oxidised_fraction', 'DOUBLE NULL'),
    ('ground_metal_fraction', 'DOUBLE NULL'),
    ('ground_carbon_fraction', 'DOUBLE NULL'),
    ('ground_evaporite_fraction', 'DOUBLE NULL'),
    ('ground_sulfur_fraction', 'DOUBLE NULL'),
    ('ground_intermediate_fraction', 'DOUBLE NULL'),
    ('surface_liquid', 'VARCHAR(32) NULL'),
    ('surface_liquid_fraction', 'DOUBLE NULL'),
    ('surface_liquid_depth_km', 'DOUBLE NULL'),
    ('surface_liquid_acid_wt', 'DOUBLE NULL'),
    ('surface_liquid_acid_layer', 'VARCHAR(8) NULL'),
    ('surface_liquid_ammonia_wt', 'DOUBLE NULL'),
    ('moon_origin', 'VARCHAR(8) NULL'),
    ('moon_disk_lifetime_myr', 'DOUBLE NULL'),
    ('moon_disk_alpha_over_f', 'DOUBLE NULL'),
    ('moon_disk_f_inf', 'DOUBLE NULL'),
    ('moon_disk_ice_line_km', 'DOUBLE NULL'),
    ('water_base_pressure_pa', 'DOUBLE NULL'),
    ('tidal_flux_w_m2', 'DOUBLE NULL'),
    ('aquifer_depth_km', 'DOUBLE NULL'),
    ('aquifer_thickness_km', 'DOUBLE NULL'),
    ('mantle_eta0_pa_s', 'DOUBLE NULL'),
    ('mantle_activation_kj_mol', 'DOUBLE NULL'),
    ('interior_kind', "VARCHAR(16) CHECK (interior_kind IN ('gas_giant', 'ice_giant', 'brown_dwarf')) NULL"),
    ('core_mass_earth', 'DOUBLE NULL'),
    ('dilute_core_mass_earth', 'DOUBLE NULL'),
    ('metallization_pressure_gpa', 'DOUBLE NULL'),
    ('helium_depletion', 'DOUBLE NULL'),
    ('has_helium_shell', 'BOOLEAN NULL'),
    ('atmosphere_colour', 'INT UNSIGNED NULL'),
    ('lid_regime', 'VARCHAR(10) NULL'),
    ('day_share_equator', 'DOUBLE NULL'),
    ('day_share_mean', 'DOUBLE NULL'),
    ('parent_shadow_share', 'DOUBLE NULL'),
    ('solar_day_hours', 'DOUBLE NULL'),
    ('daylight_hours', 'DOUBLE NULL'),
    ('polar_circle_deg', 'DOUBLE NULL'),
    ('life_class', 'TINYINT UNSIGNED NULL'),
    ('life_class_variant', 'CHAR(1) NULL'),
    ('life_class_best', 'VARCHAR(2) NULL'),
    ('life_class_best_score', 'TINYINT UNSIGNED NULL'),
    ('life_animal_grade', 'TINYINT UNSIGNED NULL'),
    ('life_complexity', 'DOUBLE NULL'),
    ('tide_range_spring_m', 'DOUBLE NULL'),
    ('tide_range_neap_m', 'DOUBLE NULL'),
    ('ocean_liquid', 'VARCHAR(40) NULL'),
    ('donor_uid', 'BINARY(12) NULL'),
    ('event_index', 'INT NULL'),
    ('ejection_delay_gy', 'DOUBLE NULL'),
    ('since_ejection_gy', 'DOUBLE NULL'),
    ('ejection_speed_kms', 'DOUBLE NULL'),
    ('rogue_note', 'TEXT NULL'),
)
"""tuple: The columns added to `rogue_planets`, each with its definition."""

RELAXED = (
    ("planet_type", "VARCHAR(4) NULL"),
    ("composition", "TEXT NULL"),
    ("has_internal_heat", "TINYINT(1) NULL"),
    ("has_moons", "TINYINT(1) NULL"),
)
"""tuple: The flat-rate rogue's NOT NULL columns that the physics rogue leaves empty."""

TABLES = (
    """
    CREATE TABLE IF NOT EXISTS rogue_planet_reflection_spectrum (
    id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    rogue_planet_id   BIGINT UNSIGNED NOT NULL,
    spectrum_type     VARCHAR(16) NOT NULL CHECK (spectrum_type IN ('visible', 'non_visible')),
    position          INT NOT NULL,
    value             VARCHAR(255) NOT NULL,

    CONSTRAINT fk_rogue_planet_reflection_spectrum_rogue_planet
        FOREIGN KEY (rogue_planet_id) REFERENCES rogue_planets(id) ON DELETE CASCADE,
    KEY idx_rogue_planet_reflection_spectrum_rogue_planet_id (rogue_planet_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""",
    """
    CREATE TABLE IF NOT EXISTS rogue_planet_interior_layers (
    id                      BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    rogue_planet_id         BIGINT UNSIGNED NOT NULL,
    layer_index             SMALLINT UNSIGNED NOT NULL,
    kind                    VARCHAR(24) NOT NULL CHECK (kind IN ('molecular_envelope', 'metallic_envelope', 'helium_gradient', 'helium_shell', 'dilute_core', 'compact_core', 'ice_mantle', 'rock_core')),
    state                   VARCHAR(16) NOT NULL CHECK (state IN ('molecular_fluid', 'metallic_fluid', 'ionic_fluid', 'superionic', 'molten')),
    top_pressure_gpa        DOUBLE NOT NULL,
    bottom_pressure_gpa     DOUBLE NOT NULL,
    mass_earth              DOUBLE NOT NULL,
    top_radius_fraction     DOUBLE NOT NULL,
    bottom_radius_fraction  DOUBLE NOT NULL,
    hydrogen                DOUBLE NOT NULL,
    helium                  DOUBLE NOT NULL,
    iron                    DOUBLE NOT NULL,
    rock                    DOUBLE NOT NULL,
    water_liquid            DOUBLE NOT NULL,
    ice_x                   DOUBLE NOT NULL,
    nh3_ice                 DOUBLE NOT NULL,
    ch4_ice                 DOUBLE NOT NULL,

    UNIQUE KEY uq_rogue_planet_interior_layers_index (rogue_planet_id, layer_index),
    CONSTRAINT fk_rogue_planet_interior_layers_rogue_planet
        FOREIGN KEY (rogue_planet_id) REFERENCES rogue_planets(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""",
    """
    CREATE TABLE IF NOT EXISTS rogue_planet_gases (
    id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    rogue_planet_id   BIGINT UNSIGNED NOT NULL,
    gas               VARCHAR(16) NOT NULL,
    partial_pa        DOUBLE NOT NULL CHECK (partial_pa > 0),

    CONSTRAINT fk_rogue_planet_gases_rogue_planet
        FOREIGN KEY (rogue_planet_id) REFERENCES rogue_planets(id) ON DELETE CASCADE,
    UNIQUE KEY uq_rogue_planet_gases (rogue_planet_id, gas)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""",
    """
    CREATE TABLE IF NOT EXISTS rogue_planet_life_class_scores (
    id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    rogue_planet_id   BIGINT UNSIGNED NOT NULL,
    life_class        VARCHAR(2) NOT NULL,
    score             TINYINT UNSIGNED NOT NULL CHECK (score BETWEEN 1 AND 9),

    CONSTRAINT fk_rogue_planet_life_class_scores_rogue_planet
        FOREIGN KEY (rogue_planet_id) REFERENCES rogue_planets(id) ON DELETE CASCADE,
    UNIQUE KEY uq_rogue_planet_life_class_scores (rogue_planet_id, life_class),
    KEY idx_rogue_planet_life_class_scores_class (life_class, score)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""",
)
"""tuple: The child tables (`schema.sql` applies them first on an existing database, so they are usually there already)."""


def upgrade():
    connection = op.get_bind()
    changes = [f"ADD COLUMN {column} {definition}" for column, definition in COLUMNS
               if not helpers.exists_column(connection, "rogue_planets", column)]
    if changes:
        connection.execute(sa.text(f"ALTER TABLE rogue_planets {', '.join(changes)}"))
    for column, definition in RELAXED:
        if helpers.exists_column(connection, "rogue_planets", column):
            connection.execute(sa.text(f"ALTER TABLE rogue_planets MODIFY {column} {definition}"))
    helpers.add_index_online(connection, "rogue_planets", "uq_rogue_planets_donor_event", ["donor_uid", "event_index"], unique=True)
    for statement in TABLES:
        connection.execute(sa.text(statement))


def downgrade():
    raise NotImplementedError("the migrations of this database only go forward")
