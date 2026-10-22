# planetgen/db/column_units.py

"""
The column registry (DB.76): every numeric column of the database declares the unit it is stored in, in one list.

A column's name carries its unit as a suffix (`radius_km`, `age_gy`, `position_x_mpc`), and until now that suffix was
the only declaration: 926 numeric columns, 422 of them with a unit in the name, and `point_columns` spelling the units
of a point again by hand. This module is the one place that declares them. `TABLES` lists, for each table, the columns
of each unit; `COLUMNS` is the same list by `(table, column)`. A unit is a key of the unit registry
(`physics.unit_registry`, GEN.568), which gives its kind and its factor to the base unit of its kind, or of
`LOCAL_UNITS` for the few the registry does not hold yet.

Rules (docs/design/unified-unit-quantity.md, Boss's cards Q3 and Q4):

* The stored unit of a column is chosen for precision, not forced to SI: a position keeps its frame's unit (galaxy pc,
  sector mpc, system km or au), because a double has a spacing that grows with the number.
* The suffix stays in the name and equals the declared unit: `SUFFIXES` maps each suffix to its unit and
  `tests/test_column_units.py` fails when a column's suffix and declaration disagree.
* No data moves. Values and units stay as stored; a column later unified is a normal migration.
* A column with no unit of its own is declared by what it is: `ratio` (a plain number), `count`, `index` (an ordinal
  or a score), `id`, `key`, `seed`, `flag`, `unit-vector` (one component of a unit vector), `row-defined` (the row says
  what it is in) or `unspecified` (a quantity whose code states no unit; each one waits for a source).
* The factors are the registry's (from `physics.constants`); a unit with no factor has no SI form (`dex`, `pH`, `mag`).

`point_columns.PointLayout` reads its position and mass units from here (`position_unit`, `mass_unit`).
"""

import dataclasses

from planetgen.physics import constants, unit_registry


@dataclasses.dataclass(frozen=True)
class LocalUnit:
    """
    A unit the unit registry (`physics.unit_registry`, GEN.568) does not hold yet, or a marker for a column that has no
    unit.

    Attributes:
        kind (str): What it measures.
        factor (float or None): The size of one of it in the SI unit of its kind, or `None` when it has none.
    """

    kind: str
    factor: float | None = None


_YEAR_S = constants.SECONDS_PER_YEAR

LOCAL_UNITS = {
    # kinds the registry has not got: volume, rates, dose, magnetic moment, logarithmic scales
    "km3": LocalUnit("volume", 1000.0 ** 3), "pc3": LocalUnit("volume", constants.PARSEC_M ** 3),
    "light-day": LocalUnit("length", constants.SPEED_OF_LIGHT_M_S * 86400.0),
    "1/s": LocalUnit("rate", 1.0), "1/yr": LocalUnit("rate", 1.0 / _YEAR_S),
    "1/Gy": LocalUnit("rate", 1.0 / (1e9 * _YEAR_S)),
    "Msun/yr": LocalUnit("mass rate", constants.SOLAR_MASS_TO_KG / _YEAR_S),
    "J/kg": LocalUnit("specific energy", 1.0),
    "mSv/yr": LocalUnit("dose rate", 1e-3 / _YEAR_S), "A*m2": LocalUnit("magnetic moment", 1.0),
    "xuv-earth": LocalUnit("flux"), "R_planet": LocalUnit("length"), "nebula-radius": LocalUnit("length"),
    "dex": LocalUnit("logarithmic"), "mag": LocalUnit("magnitude"), "pH": LocalUnit("logarithmic"),
    # no unit
    "ratio": LocalUnit("ratio"), "unit-vector": LocalUnit("ratio"), "count": LocalUnit("count"),
    "index": LocalUnit("index"), "flag": LocalUnit("flag"), "id": LocalUnit("id"), "key": LocalUnit("id"),
    "seed": LocalUnit("seed"), "row-defined": LocalUnit("row-defined"), "unspecified": LocalUnit("unspecified"),
}
"""dict: Unit key -> `LocalUnit`, for what `unit_registry.UNITS` does not hold. A unit that the registry gains moves
there, and its entry here goes."""


def unit_info(key):
    """The unit a key names, from the registry (`unit_registry.Unit`) or from `LOCAL_UNITS`; both have `kind` and
    `factor`. `KeyError` for a key in neither."""
    if key in unit_registry.UNITS:
        return unit_registry.UNITS[key]
    return LOCAL_UNITS[key]


SUFFIXES = (
    ("_j_per_kg", "J/kg"), ("_msv_yr", "mSv/yr"), ("_per_gyr", "1/Gy"), ("_per_yr", "1/yr"), ("_w_m2", "W_m2"), ("_a_m2", "A*m2"),
    ("_g_cm3", "g_cm3"), ("_light_days", "light-day"), ("_gauss", "G"), ("_myr", "Myr"), ("_ms", "ms"), ("_us", "us"), ("_mj", "Mjup"),
    ("_km3", "km3"), ("_pc3", "pc3"), ("_kms", "km_s"), ("_km", "km"), ("_pc", "pc"), ("_mpc", "mpc"), ("_au", "au"),
    ("_ly", "ly"), ("_deg", "deg"), ("_rad", "rad"), ("_years", "yr"), ("_yr", "yr"), ("_gyr", "Gyr"), ("_gy", "Gyr"),
    ("_unix", "unix_s"), ("_hours", "h"), ("_kpa", "kPa"), ("_pa", "Pa"), ("_kg", "kg"), ("_k", "K"), ("_w", "W"),
    ("_j", "J"),
)
"""tuple: A column-name suffix and the unit it stands for, longest-matching first where one is the tail of another. A
column ending in one of them must be declared in that unit (`unit_matching_suffix`). `_solar` and `_sol` are not here:
they stand for a mass or a luminosity by the rest of the name."""

SOLAR_SUFFIXES = ("_solar", "_sol")
"""tuple: The suffixes of a quantity in solar units: `Msun` when the name says mass, else `Lsun`."""

TABLES = {
    "asteroid_belt_composition": {
        "id": ("id", "belt_id"),
        "index": ("position",),
    },
    "asteroid_belts": {
        "id": ("id", "star_system_id", "star_id"),
        "index": ("orbital_index",),
        "km": ("distance_km", "lower_limit_km", "upper_limit_km", "position_x_km", "position_y_km", "position_z_km"),
        "km_s": ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
        "kg": ("mass_kg",),
        "unix_s": ("epoch_unix", "next_update_due"),
    },
    "asteroid_field_composition": {
        "id": ("id", "field_id"),
        "index": ("position",),
    },
    "asteroid_fields": {
        "id": ("id", "sector_id", "origin_event_id", "inside_nebula_id", "inside_remnant_id"),
        "ly": ("radius_ly",),
        "Gyr": ("galactic_orbital_period_gy",),
        "deg": ("galactic_orbital_phase_deg",),
        "yr": ("galactic_min_update_interval_years",),
        "unix_s": ("epoch_unix", "next_update_due"),
        "pc": ("center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc"),
        "key": ("cell_key",),
        "km_s": ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
        "kg": ("mass_kg",),
        "unspecified": ("formed_at",),
    },
    "black_holes": {
        "id": ("id", "star_id", "sector_id", "inside_nebula_id", "inside_remnant_id"),
        "Msun": ("mass_solar",),
        "km": ("event_horizon_radius_km",),
        "ratio": ("spin",),
        "flag": ("has_accretion_disk",),
        "K": ("temperature_k",),
        "W": ("luminosity_w",),
        "Gyr": ("age_gy", "galactic_orbital_period_gy"),
        "deg": ("galactic_orbital_phase_deg", "axial_tilt_deg"),
        "yr": ("galactic_min_update_interval_years",),
        "unix_s": ("epoch_unix", "next_update_due"),
        "unit-vector": ("spin_axis_x", "spin_axis_y", "spin_axis_z"),
        "pc": ("center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc"),
        "key": ("cell_key",),
        "km_s": ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
    },
    "bright_stars": {
        "id": ("id", "star_system_id"),
        "index": ("shell_index", "band_index", "wedge_index"),
        "mpc": ("position_x_mpc", "position_y_mpc", "position_z_mpc"),
        "key": ("cell_key",),
        "kg": ("mass_kg",),
        "km": ("radius_km",),
        "K": ("temperature_k",),
        "W": ("luminosity_w",),
        "Gyr": ("age_gy", "lifespan_gy", "phase_end_age_gy"),
        "Msun": ("initial_mass_sol", "system_mass_solar"),
        "seed": ("seed",),
        "flag": ("off_plane",),
        "km_s": ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
        "unix_s": ("epoch_unix", "next_update_due"),
        "count": ("member_count",),
        "dex": ("metallicity_fe_h",),
    },
    "change_plan_ops": {
        "id": ("id", "plan_id"),
        "index": ("op_order",),
    },
    "change_plans": {
        "id": ("id", "apply_job_id", "decided_by_id", "requester_id"),
    },
    "comet_composition": {
        "id": ("id", "comet_id"),
        "index": ("position",),
    },
    "comets": {
        "id": ("id", "star_system_id", "star_id"),
        "km": (
            "nucleus_diameter_km", "perihelion_distance_km", "distance_km", "position_x_km", "position_y_km",
            "position_z_km",
        ),
        "ratio": ("eccentricity",),
        "deg": ("inclination_deg", "arg_periapsis_deg", "ascending_node_deg", "mean_anomaly_deg", "axial_tilt_deg"),
        "yr": ("orbital_period_years", "min_update_interval_years"),
        "unspecified": ("parabolic_mean_anomaly",),
        "Msun": ("primary_mass_solar",),
        "flag": ("is_active",),
        "km_s": ("orbital_speed_kms", "velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
        "unix_s": ("epoch_unix", "next_update_due"),
        "h": ("rotation_period_hours",),
        "unit-vector": ("spin_axis_x", "spin_axis_y", "spin_axis_z"),
        "kg": ("mass_kg",),
    },
    "edit_journal": {
        "id": ("id", "op_id", "plan_id", "sector_id"),
        "index": ("generator_epoch",),
    },
    "facilities": {
        "id": (
            "id", "star_system_id", "star_id", "planet_id", "moon_id", "asteroid_belt_id", "asteroid_field_id",
            "sector_id", "type_id", "affiliation_polity_id",
        ),
        "pc": ("center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc"),
        "key": ("cell_key",),
        "kg": ("mass_kg",),
        "km_s": ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms", "orbital_speed_kms"),
        "km": ("orbit_distance_km",),
        "yr": ("orbit_period_years",),
        "deg": ("orbit_phase_deg",),
        "unix_s": ("epoch_unix", "next_update_due"),
        "index": ("rating_crime", "rating_housing", "rating_resources", "rating_maintenance", "rating_health"),
    },
    "facility_field_values": {
        "id": ("facility_id", "field_id"),
        "row-defined": ("value_int", "value_dec"),
        "flag": ("value_bool",),
    },
    "facility_type_fields": {
        "id": ("id", "type_id"),
        "flag": ("required",),
        "row-defined": ("min_value", "max_value"),
        "index": ("sort_order",),
    },
    "facility_types": {
        "id": ("id",),
        "index": ("sort_order",),
        "flag": ("active",),
    },
    "galaxy_block_stats": {
        "index": ("m", "shell_index", "band_index", "wedge_index"),
        "count": ("sector_count", "systems", "expected_systems", "stars"),
        "Gyr": ("age_sum",),
        "Lsun": ("luminosity_sol",),
    },
    "galaxy_components": {
        "id": ("id", "parent_component_id"),
        "pc": ("center_x_pc", "center_y_pc", "center_z_pc", "half_mass_radius_pc", "tidal_radius_pc"),
        "key": ("cell_key",),
        "km_s": ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
        "unix_s": ("epoch_unix",),
        "Msun": ("mass_solar",),
        "ratio": ("king_w0",),
        "Gyr": ("age_gy",),
        "dex": ("metallicity_fe_h",),
        "flag": ("core_collapsed",),
        "seed": ("seed",),
    },
    "galaxy_shape": {
        "id": ("id",),
        "pc": (
            "disk_scale_length_pc", "disk_scale_height_pc", "bulge_scale_radius_pc", "spiral_reference_radius_pc",
            "edge_pc", "spheroid_effective_radius_pc",
        ),
        "ratio": (
            "bulge_amplitude", "arm_amplitude", "arm_level", "nuclear_amplitude", "k_norm", "bar_axis_ratio_b",
            "bar_axis_ratio_c", "thick_disk_amplitude", "spheroid_amplitude", "spheroid_sersic_n",
            "spheroid_axis_ratio_b", "spheroid_axis_ratio_c", "body_clump_share", "density_ratio_avg",
        ),
        "count": ("arm_count", "expected_system_count_at_density_1", "cell_star_ceiling", "density_ratio_samples"),
        "rad": (
            "pitch_angle_rad", "spiral_reference_angle_rad", "bar_angle_rad", "spheroid_orientation_1_rad",
            "spheroid_orientation_2_rad", "spheroid_orientation_3_rad",
        ),
        "index": ("outer_shell_index", "reach_shell_index"),
        "Lsun": (
            "fill_luminosity_cap_sol", "core_min_luminosity_sol", "halo_min_luminosity_sol",
            "bright_star_min_luminosity_sol",
        ),
        "Msun": (
            "core_mass_limit_sol", "pinned_stellar_mass_solar", "bright_star_mass_limit_sol",
            "phenomenon_min_mass_solar",
        ),
        "flag": ("bar_present",),
        "seed": ("bright_star_seed", "phenomenon_scatter_seed"),
    },
    "galaxy_shape_pins": {
        "row-defined": ("value",),
    },
    "galaxy_shell": {
        "index": ("shell_index", "high_band_index"),
        "count": ("cell_count",),
    },
    "generation_run_arguments": {
        "id": ("run_id",),
        "index": ("position",),
    },
    "generation_runs": {
        "id": ("id",),
    },
    "id_blocks": {
        "id": ("next_id",),
    },
    "id_counters": {
        "id": ("next_value",),
    },
    "interstellar_comet_composition": {
        "id": ("id", "comet_id"),
        "index": ("position",),
    },
    "interstellar_comets": {
        "id": ("id", "sector_id", "inside_nebula_id", "inside_remnant_id"),
        "km": ("nucleus_diameter_km",),
        "km_s": ("velocity_kms", "velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
        "flag": ("is_active",),
        "Gyr": ("galactic_orbital_period_gy",),
        "deg": ("galactic_orbital_phase_deg", "axial_tilt_deg"),
        "yr": ("galactic_min_update_interval_years",),
        "unix_s": ("epoch_unix", "next_update_due"),
        "h": ("rotation_period_hours",),
        "unit-vector": ("spin_axis_x", "spin_axis_y", "spin_axis_z"),
        "pc": ("center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc"),
        "key": ("cell_key",),
        "kg": ("mass_kg",),
    },
    "moon_reflection_spectrum": {
        "id": ("id", "moon_id"),
        "index": ("position",),
    },
    "moons": {
        "id": ("id", "planet_id", "star_system_id", "star_id", "parent_id"),
        "index": (
            "orbital_index", "tier_pressure",
            "tier_temperature", "tier_air", "tier_water", "tier_suitability", "tier_radiation", "equipment_tier",
        ),
        "km": (
            "distance_km", "radius_km", "scale_height_km", "hill_radius_km", "min_orbit_distance_km",
            "habitable_zone_inner_km", "habitable_zone_outer_km", "position_x_km", "position_y_km", "position_z_km",
            "ocean_depth_km", "ice_shell_km", "hp_ice_km",
        ),
        "kg": ("mass_kg",),
        "km3": ("volume_km3",),
        "yr": ("period_years", "min_update_interval_years"),
        "g0": ("gravity_g",),
        "K": ("surface_temperature_k",),
        "g_cm3": ("density_g_cm3",),
        "unspecified": ("atm_density", "atm_molar_density"),
        "Pa": ("atmospheric_pressure_pa",),
        "count": ("flavor_text_count",),
        "Gyr": ("life_milestone_age_gyr", "life_system_age_gyr"),
        "deg": ("orbital_inclination_deg", "orbital_ascending_node_deg", "axial_tilt_deg"),
        "km_s": ("orbital_speed_kms", "velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
        "unix_s": ("epoch_unix", "next_update_due"),
        "h": ("rotation_period_hours",),
        "unit-vector": ("spin_axis_x", "spin_axis_y", "spin_axis_z"),
        "dex": ("mantle_delta_iw",),
        "kPa": (
            "p_o2_kpa", "p_co2_kpa", "p_co_kpa", "p_n2_kpa", "p_ar_kpa", "p_h2_kpa", "p_h2o_kpa", "p_ch4_kpa",
            "p_h2s_kpa", "p_so2_kpa",
        ),
        "A*m2": ("magnetic_moment_a_m2", "induced_moment_a_m2"),
        "R_planet": ("magnetopause_rp",),
        "xuv-earth": ("xuv_flux_earth",),
        "ratio": (
            "water_mass_fraction", "ocean_fraction", "land_fraction", "water_activity", "dose_helio_mult", "phi4",
            "phi4_pressure", "phi4_temperature", "phi4_air", "phi4_water", "phi4_suitability", "phi4_radiation",
            "phi_bio", "phi_cpx", "phi_tech", "l_solv", "l_chem", "l_ener", "flux_variation",
            "xuv_exposure_index", "flare_irradiation_index", "uv_surface_index",
        ),
        "pH": ("ocean_ph",),
        "mSv/yr": ("surface_dose_msv_yr", "dose_gcr_msv_yr", "dose_sep_msv_yr", "dose_ground_msv_yr"),
        "flag": ("ozone_loss_flag",),
        "rad": ("l_rad",),
        "W_m2": ("energy_flux_w_m2",),
    },
    "nearest_systems": {
        "id": ("id", "sector_id", "object_id", "star_system_id", "neighbor_system_id"),
        "index": ("neighbor_rank",),
        "pc": ("distance_pc",),
    },
    "nebula_block_classes": {
        "index": ("m", "shell_index", "band_index", "wedge_index"),
        "pc3": ("volume_pc3",),
    },
    "nebula_block_state": {
        "id": ("id", "nebulae_through_id"),
    },
    "nebula_block_stats": {
        "index": ("m", "shell_index", "band_index", "wedge_index"),
        "count": ("nebula_count",),
        "pc3": (
            "volume_pc3", "volume_diffuse_pc3", "volume_emission_pc3", "volume_reflection_pc3",
            "volume_planetary_pc3", "volume_dark_pc3", "volume_remnant_pc3", "cover_volume_pc3",
        ),
        "pc": ("centroid_x_pc", "centroid_y_pc", "centroid_z_pc", "bound_radius_pc"),
        "ratio": ("cover",),
    },
    "nebula_field": {
        "index": ("cell_i", "cell_k", "n"),
        "J": ("cell_j",),
        "ly": ("radius_ly",),
        "pc": ("center_x_pc", "center_y_pc", "center_z_pc"),
        "key": ("cell_key",),
        "flag": ("is_stored",),
    },
    "nebula_shape_balls": {
        "id": ("nebula_id",),
        "index": ("ball_index",),
        "nebula-radius": ("center_x", "center_y", "center_z", "radius"),
    },
    "nebulae": {
        "id": ("id", "sector_id", "inside_nebula_id", "inside_remnant_id"),
        "ly": ("radius_ly",),
        "g_cm3": ("density_cm3",),
        "K": ("temperature_k",),
        "mag": ("extinction_av",),
        "Myr": ("age_myr",),
        "1/s": ("q_total",),
        "Gyr": ("galactic_orbital_period_gy",),
        "deg": ("galactic_orbital_phase_deg",),
        "yr": ("galactic_min_update_interval_years",),
        "unix_s": ("epoch_unix", "next_update_due"),
        "pc": ("center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc"),
        "key": ("cell_key",),
        "km_s": ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
        "kg": ("mass_kg",),
        "ratio": (
            "shape_axis_x", "shape_axis_y", "shape_axis_z", "shape_warp_amplitude", "shape_warp_frequency",
            "shape_iso", "shape_scale",
        ),
        "rad": ("shape_angle_x", "shape_angle_y", "shape_angle_z"),
        "index": ("shape_warp_octaves",),
        "seed": ("shape_noise_seed",),
    },
    "neutron_stars": {
        "id": ("id", "star_id", "sector_id", "inside_nebula_id", "inside_remnant_id"),
        "Msun": ("mass_solar",),
        "km": ("radius_km",),
        "ms": ("spin_period_ms",),
        "G": ("magnetic_field_gauss",),
        "K": ("surface_temperature_k",),
        "W": ("luminosity_w",),
        "Gyr": ("age_gy", "galactic_orbital_period_gy"),
        "deg": ("galactic_orbital_phase_deg", "axial_tilt_deg"),
        "yr": ("galactic_min_update_interval_years",),
        "unix_s": ("epoch_unix", "next_update_due"),
        "unit-vector": ("spin_axis_x", "spin_axis_y", "spin_axis_z"),
        "pc": ("center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc"),
        "key": ("cell_key",),
        "km_s": ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
    },
    "object_index": {
        "id": ("source_id", "sector_id", "system_id"),
        "kg": ("mass_kg",),
        "km": ("radius_km",),
        "flag": ("placed",),
    },
    "object_locks": {
        "id": ("locked_by_id",),
    },
    "object_timing_samples": {
        "id": ("id", "run_id", "sector_id"),
        "index": ("ordinal",),
        "us": ("own_us", "total_us"),
        "flag": ("sampled",),
    },
    "object_timing_stats": {
        "id": ("id", "run_id", "sector_id"),
        "count": ("object_count",),
        "us": ("total_us", "min_us", "max_us"),
    },
    "orbit_nodes": {
        "id": ("star_system_id",),
        "index": ("node_no", "parent_node_no", "slot"),
        "Msun": ("mass_solar", "disk_mass_solar"),
        "km": ("separation_km",),
        "ratio": ("eccentricity",),
        "deg": ("inclination_deg", "ascending_node_deg", "arg_periapsis_deg", "mean_anomaly_deg"),
        "unix_s": ("epoch_unix",),
        "yr": ("period_years",),
        "flag": ("flags",),
        "au": ("disk_inner_au", "disk_outer_au"),
    },
    "orbit_simulation_state": {
        "id": ("id",),
    },
    "orbital_event_bodies": {
        "id": ("event_id", "sector_id"),
        "index": ("ordinal",),
        "kg": ("mass_kg",),
        "km": ("radius_km",),
        "pc": ("x_pc", "y_pc", "z_pc"),
        "km_s": ("vx_kms", "vy_kms", "vz_kms"),
        "flag": ("admin_edited",),
    },
    "orbital_events": {
        "id": ("id", "run_id"),
        "unix_s": ("event_time_unix", "detected_time_unix"),
        "km": ("separation_start_km", "separation_impact_km", "impact_parameter_km", "critical_impact_parameter_km",
               "new_radius_km"),
        "km_s": ("v_inf_kms", "v_esc_kms", "v_imp_kms", "catastrophic_velocity_kms"),
        "J/kg": ("q_r_j_per_kg", "q_star_rd_j_per_kg"),
        "ratio": ("erosion_x", "interacting_mass_fraction", "mass_error", "momentum_error", "conservation_tolerance",
                  "heat_share", "spin_breakup_ratio"),
        "s": ("window_s",),
        "J": ("impact_energy_j",),
        "kg": ("lost_mass_kg", "merged_mass_kg"),
        "Mjup": ("deuterium_limit_mj",),
        "Msun": ("hydrogen_limit_sol",),
        "count": ("influence_body_count", "sectors_dirtied", "suppressed_duplicates"),
        "flag": ("window_capped", "focusing_radius_used", "borderline", "already_overlapping", "numeric_warning"),
    },
    "phenomenon_scatter": {
        "id": ("id",),
        "index": ("shell_index", "band_index", "wedge_index"),
        "mpc": ("position_x_mpc", "position_y_mpc", "position_z_mpc"),
        "key": ("cell_key",),
        "km_s": ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
        "seed": ("seed",),
        "unix_s": ("epoch_unix", "next_update_due"),
        "Msun": ("mass_solar", "system_mass_solar"),
        "Lsun": ("map_luminosity_sol",),
        "count": ("member_count",),
    },
    "phenomenon_scatter_classes": {
        "count": ("placed", "built"),
    },
    "planet_reflection_spectrum": {
        "id": ("id", "planet_id"),
        "index": ("position",),
    },
    "planets": {
        "id": ("id", "star_system_id", "star_id", "parent_id", "host_node_id"),
        "index": (
            "orbital_index", "tier_pressure",
            "tier_temperature", "tier_air", "tier_water", "tier_suitability", "tier_radiation", "equipment_tier",
        ),
        "km": (
            "distance_km", "radius_km", "scale_height_km", "hill_radius_km", "min_orbit_distance_km",
            "habitable_zone_inner_km", "habitable_zone_outer_km", "position_x_km", "position_y_km", "position_z_km",
            "ocean_depth_km", "ice_shell_km", "hp_ice_km", "reflex_offset_x_km", "reflex_offset_y_km",
            "reflex_offset_z_km",
        ),
        "kg": ("mass_kg",),
        "km3": ("volume_km3",),
        "yr": ("period_years", "min_update_interval_years"),
        "g0": ("gravity_g",),
        "K": ("surface_temperature_k",),
        "g_cm3": ("density_g_cm3",),
        "unspecified": ("atm_density", "atm_molar_density"),
        "Pa": ("atmospheric_pressure_pa",),
        "count": ("flavor_text_count",),
        "Gyr": ("life_milestone_age_gyr", "life_system_age_gyr"),
        "deg": ("orbital_inclination_deg", "orbital_ascending_node_deg", "axial_tilt_deg"),
        "km_s": ("orbital_speed_kms", "velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
        "unix_s": ("epoch_unix", "next_update_due"),
        "h": ("rotation_period_hours",),
        "unit-vector": ("spin_axis_x", "spin_axis_y", "spin_axis_z"),
        "dex": ("mantle_delta_iw",),
        "kPa": (
            "p_o2_kpa", "p_co2_kpa", "p_co_kpa", "p_n2_kpa", "p_ar_kpa", "p_h2_kpa", "p_h2o_kpa", "p_ch4_kpa",
            "p_h2s_kpa", "p_so2_kpa",
        ),
        "A*m2": ("magnetic_moment_a_m2", "induced_moment_a_m2"),
        "R_planet": ("magnetopause_rp",),
        "xuv-earth": ("xuv_flux_earth",),
        "ratio": (
            "water_mass_fraction", "ocean_fraction", "land_fraction", "water_activity", "dose_helio_mult", "phi4",
            "phi4_pressure", "phi4_temperature", "phi4_air", "phi4_water", "phi4_suitability", "phi4_radiation",
            "phi_bio", "phi_cpx", "phi_tech", "l_solv", "l_chem", "l_ener", "flux_variation",
            "xuv_exposure_index", "flare_irradiation_index", "uv_surface_index",
        ),
        "pH": ("ocean_ph",),
        "mSv/yr": ("surface_dose_msv_yr", "dose_gcr_msv_yr", "dose_sep_msv_yr", "dose_ground_msv_yr"),
        "flag": ("ozone_loss_flag",),
        "rad": ("l_rad",),
        "W_m2": ("energy_flux_w_m2",),
    },
    "polities": {
        "id": ("id", "species_id", "capital_system_id"),
        "ly": ("reach_ly",),
    },
    "population_state": {
        "id": ("id", "scanned_planet_id"),
    },
    "quasars": {
        "id": ("id", "sector_id"),
        "Msun": ("black_hole_mass_solar",),
        "km": ("event_horizon_radius_km",),
        "ratio": ("eddington_ratio",),
        "W": ("luminosity_w",),
        "Msun/yr": ("accretion_rate_solar_per_year",),
        "light-day": ("broad_line_region_light_days",),
        "flag": ("is_radio_loud",),
        "ly": ("jet_length_ly",),
        "yr": ("active_age_years",),
        "pc": ("center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc"),
        "key": ("cell_key",),
        "km_s": ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
        "unix_s": ("epoch_unix", "next_update_due"),
    },
    "rogue_planets": {
        "id": ("id", "sector_id", "inside_nebula_id", "inside_remnant_id"),
        "kg": ("mass_kg",),
        "km": ("radius_km", "ice_shell_thickness_km", "ocean_depth_km", "hp_ice_km"),
        "flag": ("has_internal_heat", "has_moons", "has_liquid_water"),
        "Gyr": ("age_gy", "galactic_orbital_period_gy"),
        "W_m2": ("internal_heat_flux_w_m2",),
        "K": ("effective_temperature_k", "surface_temperature_k"),
        "Pa": ("surface_pressure_pa",),
        "deg": ("galactic_orbital_phase_deg", "axial_tilt_deg"),
        "yr": ("galactic_min_update_interval_years",),
        "unix_s": ("epoch_unix", "next_update_due"),
        "h": ("rotation_period_hours",),
        "unit-vector": ("spin_axis_x", "spin_axis_y", "spin_axis_z"),
        "pc": ("center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc"),
        "key": ("cell_key",),
        "km_s": ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
    },
    "scatter_members": {
        "index": ("source", "member_no", "parent_node_no", "slot"),
        "id": ("source_id",),
        "flag": ("flags",),
        "au": ("offset_x_au", "offset_y_au", "offset_z_au"),
        "km_s": ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
        "yr": ("decay_age_yr",),
        "Msun": ("mass_solar", "initial_mass_solar"),
        "seed": ("seed",),
        "unix_s": ("epoch_unix",),
    },
    "scatter_nodes": {
        "index": ("source", "node_no", "parent_node_no", "slot"),
        "id": ("source_id",),
        "Msun": ("mass_solar",),
        "km": ("separation_km",),
        "ratio": ("eccentricity",),
        "deg": ("inclination_deg", "ascending_node_deg", "arg_periapsis_deg", "mean_anomaly_deg"),
        "unix_s": ("epoch_unix",),
        "yr": ("period_years",),
        "flag": ("flags",),
    },
    "schema_migrations": {
        "index": ("version",),
    },
    "sector_name_registry": {
        "id": ("id", "first_sector_id"),
        "count": ("occurrence_count",),
    },
    "sector_path_knots": {
        "id": ("id", "path_id"),
        "index": ("position",),
        "yr": ("t_years",),
        "pc": ("x_pc", "y_pc", "z_pc"),
        "km_s": ("vx_kms", "vy_kms", "vz_kms"),
    },
    "sector_paths": {
        "id": ("id", "sector_id", "object_id"),
        "flag": ("exited",),
        "yr": ("duration_years",),
    },
    "sector_stats": {
        "index": ("shell_index", "band_index", "wedge_index"),
        "Lsun": ("bright_level_sol", "level_before_fill_sol", "mean_luminosity_sol", "total_luminosity_sol"),
        "Msun": ("bright_mass_sol",),
        "ratio": ("relative_density",),
        "count": ("expected_systems", "actual_systems", "actual_stars"),
        "K": ("mean_temperature_k",),
        "Gyr": ("mean_age_gy",),
    },
    "sector_system_counts": {
        "id": ("sector_id",),
        "count": ("system_count",),
    },
    "sector_tree_nodes": {
        "index": ("level", "shell_index", "band_index"),
        "count": ("child_count",),
        "flag": ("dirty",),
    },
    "sectors": {
        "id": ("id", "uid", "run_id"),
        "mpc": ("edge_mpc",),
        "pc": ("center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc"),
        "index": ("shell_index", "band_index", "wedge_index", "generator_epoch"),
        "kg": ("max_mass_kg",),
        "flag": ("parity_dirty",),
    },
    "species": {
        "id": ("id", "homeworld_planet_id", "star_system_id"),
        "yr": ("civilization_age_years",),
        "flag": ("spacefaring",),
        "index": (
            "energy_index", "materials_index", "information_index", "medical_index", "propulsion_index",
            "defense_index", "tech_level",
        ),
    },
    "species_limit_values": {
        "id": ("species_id",),
        "row-defined": ("value",),
    },
    "star_systems": {
        "id": ("id", "sector_id", "system_config_id", "inside_nebula_id", "inside_remnant_id"),
        "mpc": ("position_x_mpc", "position_y_mpc", "position_z_mpc"),
        "pc": ("galaxy_x_pc", "galaxy_y_pc", "galaxy_z_pc"),
        "key": ("cell_key",),
        "flag": ("is_binary",),
        "km": (
            "binary_separation_km", "binary_periapsis_km", "binary_apoapsis_km", "binary_radius_km",
            "binary_habitable_zone_inner_km", "binary_habitable_zone_outer_km", "binary_system_perimeter_km",
            "binary_heliosphere_radius_km", "binary_mutual_position_x_km", "binary_mutual_position_y_km",
            "binary_mutual_position_z_km", "binary_primary_position_x_km", "binary_primary_position_y_km",
            "binary_primary_position_z_km", "binary_secondary_position_x_km", "binary_secondary_position_y_km",
            "binary_secondary_position_z_km", "binary_planetary_wobble_x_km", "binary_planetary_wobble_y_km",
            "binary_planetary_wobble_z_km",
        ),
        "ratio": ("binary_eccentricity", "binary_secondary_mass_fraction"),
        "K": ("binary_temperature_k",),
        "kg": ("binary_effective_mass_kg",),
        "W": ("binary_effective_luminosity_w",),
        "Gyr": ("binary_age_gy", "binary_lifespan_gy", "binary_galactic_orbital_period_gy"),
        "deg": (
            "binary_galactic_orbital_phase_deg", "binary_mutual_orbital_inclination_deg",
            "binary_mutual_orbital_ascending_node_deg", "binary_mutual_orbital_phase_deg",
        ),
        "yr": (
            "binary_galactic_min_update_interval_years", "binary_mutual_orbital_period_years",
            "binary_mutual_min_update_interval_years",
        ),
        "km_s": (
            "binary_mutual_orbital_speed_kms", "runaway_speed_kms", "velocity_x_kms", "velocity_y_kms",
            "velocity_z_kms",
        ),
        "unix_s": ("binary_epoch_unix", "binary_next_update_due", "epoch_unix", "next_update_due"),
        "index": ("schema_version",),
    },
    "stars": {
        "id": ("id", "star_system_id", "parent_node_id"),
        "kg": ("mass_kg",),
        "km": (
            "radius_km", "habitable_zone_inner_km", "habitable_zone_outer_km", "system_perimeter_km",
            "heliosphere_radius_km", "wide_binary_a_crit_km", "reflex_offset_x_km", "reflex_offset_y_km",
            "reflex_offset_z_km",
        ),
        "K": ("temperature_k",),
        "W": ("luminosity_w", "l_xuv_w"),
        "Gyr": ("age_gy", "lifespan_gy", "galactic_orbital_period_gy"),
        "deg": ("galactic_orbital_phase_deg", "axial_tilt_deg"),
        "yr": ("galactic_min_update_interval_years", "decay_age_yr"),
        "h": ("rotation_period_hours",),
        "unit-vector": ("spin_axis_x", "spin_axis_y", "spin_axis_z"),
        "dex": ("log_lx_lbol", "metallicity_fe_h"),
        "flag": ("xuv_saturated", "flags"),
        "1/yr": ("flare_n33_per_yr",),
        "ratio": ("flare_alpha",),
        "J": ("xuv_fluence_j",),
        "1/Gy": ("lethal_event_rate_per_gyr",),
        "au": ("offset_x_au", "offset_y_au", "offset_z_au", "disk_inner_au", "disk_outer_au"),
        "km_s": ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
        "unix_s": ("epoch_unix",),
        "Msun": ("disk_mass_solar",),
    },
    "supernova_remnants": {
        "id": ("id", "sector_id", "compact_remnant_black_hole_id", "compact_remnant_neutron_star_id"),
        "yr": ("age_years", "galactic_min_update_interval_years"),
        "ly": ("radius_ly",),
        "g_cm3": ("density_cm3",),
        "K": ("temperature_k",),
        "mag": ("extinction_av",),
        "Gyr": ("galactic_orbital_period_gy",),
        "deg": ("galactic_orbital_phase_deg",),
        "unix_s": ("epoch_unix", "next_update_due"),
        "pc": ("center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc"),
        "key": ("cell_key",),
        "km_s": ("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"),
        "kg": ("mass_kg",),
    },
    "system_config_slots": {
        "id": ("id", "config_id"),
        "index": ("orbit_index",),
        "count": ("moons",),
    },
    "system_configs": {
        "id": ("id",),
        "flag": (
            "habitable_world", "asteroid_belt", "comets", "large_star", "moons", "planets", "intelligent_life",
            "binary_system", "wide_binary",
        ),
        "count": ("max_planets", "num_orbits"),
        "ratio": (
            "prevalence_habitable_world", "prevalence_asteroid_belt", "prevalence_comets", "prevalence_large_star",
            "prevalence_moons", "prevalence_max_planets", "prevalence_intelligent_life", "prevalence_binary_system",
            "prevalence_wide_binary", "prevalence_planets",
        ),
    },
    "system_name_registry": {
        "id": ("id", "first_star_system_id", "first_object_id"),
        "count": ("occurrence_count",),
        "index": ("diminutive_index",),
    },
    "system_owners": {
        "id": ("star_system_id", "polity_id"),
        "ly": ("distance_ly",),
    },
}
"""dict: Table -> unit -> the columns stored in it. Every numeric column of every table is on it once."""

COLUMNS = {(table, column): unit for table, by_unit in TABLES.items() for unit, columns in by_unit.items()
           for column in columns}
"""dict: `(table, column)` -> unit key."""


def unit_of(table, column):
    """The unit key a numeric column is stored in; `KeyError` for a column that is not declared."""
    return COLUMNS[(table, column)]


def kind_of(table, column):
    """What the column measures (the unit's `kind`)."""
    return unit_info(COLUMNS[(table, column)]).kind


def to_si(table, column):
    """The size of one stored unit of the column in the SI unit of its kind; `ValueError` when the unit has none."""
    unit = COLUMNS[(table, column)]
    factor = unit_info(unit).factor
    if factor is None:
        raise ValueError(f"{table}.{column} is stored in {unit!r}, which has no SI form")
    return factor


SUFFIX_EXCEPTIONS = {
    ("nebula_field", "cell_k"): "the third axis of a cell's index, not kelvin",
    ("orbital_events", "hydrogen_limit_sol"): "a mass (the hydrogen-burning limit in solar masses), though the name has no 'mass'",
}
"""dict: `(table, column)` -> why a name that ends in a unit's suffix is not in that unit."""


def unit_matching_suffix(column, table=None):
    """
    The unit the name of `column` says it is in, or `None` when the name carries none (or is an entry of
    `SUFFIX_EXCEPTIONS`): the first of `SUFFIXES` it ends with, and `Msun` or `Lsun` for a solar suffix (a name with
    `mass` in it is a mass).
    """
    if (table, column) in SUFFIX_EXCEPTIONS:
        return None
    for suffix, unit in SUFFIXES:
        if column.endswith(suffix):
            return unit
    for suffix in SOLAR_SUFFIXES:
        if column.endswith(suffix):
            return "Msun" if "mass" in column else "Lsun"
    return None
