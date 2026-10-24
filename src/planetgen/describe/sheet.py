# planetgen/describe/sheet.py

"""
The property sheet (UX.168): one registry of the stored properties each kind of body shows on its own page and
through `GET /api/objects/<ref>/sheet`, so the website, the API and the tests read one list.

A `PropertySpec` is a stored column, its label, the group it sits in and how it is shown. A star, planet, moon, belt,
comet or standalone phenomenon (UX.175) has an ordered tuple of them in `SPECS`. Every value goes through the site's one unit formatter: a
quantity with a ladder (a length, a mass, a time, a temperature, a pressure, a speed, a power, a gravity, a flux) is a
Document quantity (`describe.model.qty`), so the page shows it in the reader's unit system exactly as the description
does; any other number is written with its stored unit's symbol, never as a bare float. A column with no stored value
is a row that says "not stored", so an empty section means no data and not a bug.

The registry is held complete by `tests/test_property_sheet.py`: every column of a kind's table is either in `SPECS`
or in `NOT_SHOWN` with the reason, so a column a later change adds cannot be forgotten. A numeric column also has its
unit in `db/column_units.py`, and the unit decides how the sheet shows it.

Reads only: no generation code and no description text is touched.
"""

import dataclasses

from planetgen.db import column_units, models
from planetgen.describe import model
from planetgen.describe.quantity import inline_text
from planetgen.physics import constants, habitability_world, unit_registry
from planetgen.util.format import format_number

TABLES = {"star": "stars", "planet": "planets", "moon": "moons", "belt": "asteroid_belts", "comet": "comets",
          "nebula": "nebulae", "asteroid_field": "asteroid_fields", "black_hole": "black_holes",
          "neutron_star": "neutron_stars", "supernova_remnant": "supernova_remnants", "rogue_planet": "rogue_planets",
          "interstellar_comet": "interstellar_comets", "quasar": "quasars"}
"""dict: The table each kind of the sheet is stored in: a body of a system, or a phenomenon (UX.175)."""

GROUPS = {
    "identity": "Identity", "orbit": "Orbit", "physical": "Physical", "rotation": "Rotation",
    "composition": "Composition", "atmosphere": "Atmosphere", "climate": "Climate", "water": "Water",
    "magnetism": "Magnetic field", "radiation": "Radiation", "habitability": "Habitability", "life": "Life",
    "zones": "Zones and limits", "stellar": "Stellar activity", "galactic": "Galactic orbit", "disk": "Debris disk",
    "radiogenic": "Radiogenic heat", "gas": "Gas and light", "compact": "Compact object", "activity": "Activity",
}
"""dict: Every group key of a sheet and its heading, in the order the page shows them."""

NOT_STORED = "not stored"
"""str: What a row with no stored value says."""

REGISTRY_QUANTITIES = {
    "length": "length", "mass": "mass", "duration": "time", "temperature": "temperature", "pressure": "pressure",
    "speed": "speed", "power": "power", "acceleration": "gravity", "flux": "flux", "ratio": "ratio",
}
"""dict: A unit registry kind, to the Document quantity kind that shows it with its ladder. Any other kind is written
as a number and the unit's symbol."""

TEXT, FLAG, TIER, AUTO = "text", "flag", "tier", "auto"
"""str: How a column is shown besides a ladder quantity: as stored text, as Yes or No, as a PHI-6 colour, or from its
unit (`AUTO`)."""


@dataclasses.dataclass(frozen=True)
class PropertySpec:
    """
    One row of a sheet.

    Attributes:
        key (str): The stable name of the row in the API; the column's own name.
        label (str): What the page calls it.
        group (str): A key of `GROUPS`.
        quantity (str): `AUTO` (from the column's unit), `TEXT`, `FLAG`, `TIER`, or a Document quantity kind that
            overrides the unit's own (`radius` for a body's radius, `luminosity` with its ratio to the Sun).
        column (str): The stored column the value is read from.
        solar (bool): The row also carries its ratio to the Sun (a star's mass and luminosity).
    """

    key: str
    label: str
    group: str
    quantity: str = AUTO
    column: str = ""
    solar: bool = False

    def stored_column(self):
        """The column the value comes from."""
        return self.column or self.key


def _rows(group, *rows):
    """The specs of one group: each row is `(column, label)`, `(column, label, quantity)` or with `solar` as a fourth
    item."""
    return tuple(PropertySpec(row[0], row[1], group, row[2] if len(row) > 2 else AUTO, solar=len(row) > 3 and row[3])
                 for row in rows)


_STAR = (
    *_rows("identity", ("role", "Role", TEXT), ("star_type", "Spectral type", TEXT),
           ("yerkes_class", "Luminosity class", TEXT), ("state", "State", TEXT)),
    *_rows("physical", ("mass_kg", "Mass", AUTO, True), ("radius_km", "Radius", "radius"),
           ("temperature_k", "Temperature"), ("luminosity_w", "Luminosity", "luminosity", True),
           ("age_gy", "Age"), ("lifespan_gy", "Main-sequence lifetime"), ("metallicity_fe_h", "Metallicity [Fe/H]"),
           ("decay_age_yr", "Age since it stopped burning")),
    *_rows("zones", ("habitable_zone_inner_km", "Habitable zone, inner edge"),
           ("habitable_zone_outer_km", "Habitable zone, outer edge"), ("snow_line_km", "Snow line"),
           ("heliosphere_radius_km", "Heliosphere radius"), ("system_perimeter_km", "System perimeter"),
           ("wide_binary_a_crit_km", "Widest stable orbit around a pair")),
    *_rows("rotation", ("rotation_period_hours", "Rotation period"), ("axial_tilt_deg", "Axial tilt")),
    *_rows("stellar", ("log_lx_lbol", "X-ray to bolometric luminosity, log"), ("l_xuv_w", "X-ray and ultraviolet power"),
           ("xuv_saturated", "X-ray output saturated", FLAG), ("xuv_fluence_j", "X-ray and ultraviolet energy emitted"),
           ("uv_c_share", "Share of the ultraviolet that is UV-C"),
           ("flare_n33_per_yr", "Large flares a year"), ("flare_alpha", "Flare size distribution slope"),
           ("lethal_event_rate_per_gyr", "Lethal events per billion years")),
    *_rows("galactic", ("galactic_orbital_period_gy", "Galactic orbital period"),
           ("galactic_orbital_phase_deg", "Galactic orbital phase")),
    *_rows("disk", ("disk_kind", "Kind", TEXT), ("disk_mass_solar", "Mass"), ("disk_inner_au", "Inner edge"),
           ("disk_outer_au", "Outer edge")),
    *_rows("orbit", ("offset_x_au", "Offset from the system's centre, x"), ("offset_y_au", "Offset, y"),
           ("offset_z_au", "Offset, z")),
)

_BODY_GROUPS = (
    *_rows("identity", ("planet_class", "Class", TEXT), ("zone", "Zone", TEXT)),
    *_rows("orbit", ("distance_km", "Distance from its star"), ("period_years", "Orbital period"),
           ("orbital_speed_kms", "Orbital speed"), ("orbital_inclination_deg", "Inclination"),
           ("orbital_ascending_node_deg", "Longitude of the ascending node"), ("hill_radius_km", "Hill radius"),
           ("min_orbit_distance_km", "Closest stable orbit for a moon"), ("dynamical_state", "Dynamical state", TEXT),
           ("habitable_zone_inner_km", "Habitable zone, inner edge"),
           ("habitable_zone_outer_km", "Habitable zone, outer edge")),
    *_rows("physical", ("radius_km", "Radius", "radius"), ("mass_kg", "Mass"), ("volume_km3", "Volume"),
           ("density_g_cm3", "Density"), ("gravity_g", "Surface gravity")),
    *_rows("rotation", ("rotation_period_hours", "Rotation period"), ("axial_tilt_deg", "Axial tilt"),
           ("outside_hours", "Hours of daylight a person can be outside")),
    *_rows("composition", ("mantle_redox", "Mantle oxidation", TEXT), ("mantle_delta_iw", "Mantle oxidation, log below iron-wüstite"),
           ("phosphorus", "Phosphorus supply", TEXT), ("born_beyond_snow_line", "Formed beyond the snow line", FLAG),
           ("metal_mass_fraction", "Metal (iron) mass fraction"), ("rock_mass_fraction", "Rock mass fraction"),
           ("ice_mass_fraction", "Ice mass fraction"), ("envelope_mass_fraction", "Gas envelope mass fraction"),
           ("envelope_helium_fraction", "Helium share of the envelope"), ("envelope_z_fraction", "Heavy-element share of the envelope"),
           ("rock_feo_fraction", "Iron oxide share of the rock"), ("rock_al2o3_fraction", "Alumina share of the rock"),
           ("rock_mg_si_ratio", "Magnesium to silicon ratio of the rock"), ("porosity", "Porosity"),
           ("carbon_rock", "Carbon-rich rock", FLAG), ("composition_origin", "How the composition was found", TEXT)),
    *_rows("atmosphere", ("atmospheric_pressure_pa", "Surface pressure"), ("atm_density", "Air density"),
           ("atm_molar_density", "Mean molar mass of the air"), ("scale_height_km", "Scale height"),
           ("p_o2_kpa", "Oxygen partial pressure"), ("p_co2_kpa", "Carbon dioxide partial pressure"),
           ("p_co_kpa", "Carbon monoxide partial pressure"), ("p_n2_kpa", "Nitrogen partial pressure"),
           ("p_ar_kpa", "Argon partial pressure"), ("p_h2_kpa", "Hydrogen partial pressure"),
           ("p_h2o_kpa", "Water vapour partial pressure"), ("p_ch4_kpa", "Methane partial pressure"),
           ("p_h2s_kpa", "Hydrogen sulfide partial pressure"), ("p_so2_kpa", "Sulfur dioxide partial pressure"),
           ("tau_haze", "Haze optical depth"), ("tau_dust", "Dust optical depth"), ("tau_cloud", "Cloud optical depth"),
           ("cloud_cover", "Cloud cover"), ("ozone_column_du", "Ozone column"),
           ("ozone_loss_flag", "Ozone layer lost", FLAG)),
    *_rows("climate", ("bond_albedo", "Bond albedo"), ("greenhouse_factor", "Greenhouse factor"),
           ("surface_temperature_k", "Mean surface temperature"), ("temp_polar_k", "Polar temperature"),
           ("temp_north_k", "Northern temperature"), ("temp_equator_k", "Equatorial temperature"),
           ("temp_south_k", "Southern temperature"), ("temp_swing_k", "Day and night swing"),
           ("temp_day_side_k", "Day side temperature"), ("temp_terminator_k", "Terminator temperature"),
           ("temp_night_side_k", "Night side temperature"), ("energy_flux_w_m2", "Stellar energy received"),
           ("flux_variation", "Seasonal variation of that energy")),
    *_rows("water", ("hydrosphere", "Hydrosphere", TEXT), ("water_mass_fraction", "Water share of the mass"),
           ("ocean_fraction", "Surface under ocean"), ("land_fraction", "Surface above sea level"),
           ("ocean_depth_km", "Ocean depth"), ("ice_shell_km", "Ice shell thickness"),
           ("hp_ice_km", "High-pressure ice thickness"), ("ocean_class", "Ocean chemistry", TEXT),
           ("ocean_ph", "Ocean pH"), ("water_activity", "Water activity"),
           ("water_base_pressure_pa", "Pressure at the bottom of the water layer")),
    *_rows("magnetism", ("dipole_class", "Dipole field", TEXT), ("magnetic_moment_a_m2", "Magnetic moment"),
           ("induced_moment_a_m2", "Induced magnetic moment"), ("magnetopause_rp", "Magnetopause distance")),
    *_rows("radiation", ("surface_dose_msv_yr", "Surface dose"), ("dose_gcr_msv_yr", "Dose from cosmic rays"),
           ("dose_sep_msv_yr", "Dose from stellar particles"), ("dose_ground_msv_yr", "Dose from the ground"),
           ("dose_helio_mult", "Heliosphere shielding factor"), ("xuv_flux_earth", "X-ray and ultraviolet flux"),
           ("xuv_exposure_index", "X-ray and ultraviolet exposure index"),
           ("flare_irradiation_index", "Flare irradiation index"), ("uv_surface_index", "Surface ultraviolet index"),
           ("shelter_column_g_cm2", "Shelter column")),
    *_rows("habitability", ("phi4", "Habitability score (PHI-6)"), ("phi4_pressure", "Pressure score"),
           ("phi4_temperature", "Temperature score"), ("phi4_air", "Air score"), ("phi4_water", "Water score"),
           ("phi4_suitability", "Suitability score"), ("phi4_radiation", "Radiation score"),
           ("tier_pressure", "Pressure rating", TIER), ("tier_temperature", "Temperature rating", TIER),
           ("tier_air", "Air rating", TIER), ("tier_water", "Water rating", TIER),
           ("tier_suitability", "Suitability rating", TIER), ("tier_radiation", "Radiation rating", TIER)),
    *_rows("life", ("life_chemical", "Life chemistry", TEXT), ("evolutionary_pace", "Evolutionary pace", TEXT), ("life_stage", "Life stage", TEXT),
           ("life_milestone_age_gyr", "Age of that stage"), ("life_system_age_gyr", "Age of the system when it was written"),
           ("life_microbial", "Microbial life score"), ("life_flora", "Flora score"), ("life_fauna", "Fauna score"),
           ("life_sentient", "Sentient life score"),
           ("l_solv", "Solvent suitability"), ("l_chem", "Chemistry suitability"), ("l_ener", "Energy suitability"),
           ("l_rad", "Radiation tolerance")),
    *_rows("radiogenic", ("c_k0_ppm", "Potassium at formation"), ("c_th0_ppm", "Thorium at formation"),
           ("c_u0_ppm", "Uranium at formation"), ("al26_ratio", "Aluminium-26 ratio"), ("urey_ratio", "Urey ratio")),
)

_PLANET = _MOON = _BODY_GROUPS

_BELT = (
    *_rows("orbit", ("distance_km", "Distance from its star"), ("lower_limit_km", "Inner edge"),
           ("upper_limit_km", "Outer edge")),
    *_rows("physical", ("density", "Density", TEXT), ("mass_kg", "Mass")),
    *_rows("composition", ("composition_summary", "Composition", TEXT)),
)

_COMET = (
    *_rows("identity", ("orbit_type", "Orbit", TEXT), ("period_class", "Period class", TEXT),
           ("is_active", "Active", FLAG)),
    *_rows("orbit", ("perihelion_distance_km", "Closest approach to its star"), ("distance_km", "Distance from its star"),
           ("eccentricity", "Eccentricity"), ("inclination_deg", "Inclination"),
           ("arg_periapsis_deg", "Argument of periapsis"), ("ascending_node_deg", "Longitude of the ascending node"),
           ("mean_anomaly_deg", "Mean anomaly"), ("orbital_period_years", "Orbital period"),
           ("orbital_speed_kms", "Orbital speed"), ("primary_mass_solar", "Mass of the star it orbits")),
    *_rows("physical", ("nucleus_diameter_km", "Nucleus diameter"), ("mass_kg", "Mass")),
    *_rows("rotation", ("rotation_period_hours", "Rotation period"), ("axial_tilt_deg", "Axial tilt")),
    *_rows("composition", ("composition_summary", "Composition", TEXT)),
)

_GALACTIC_ORBIT = _rows("galactic", ("galactic_orbital_period_gy", "Galactic orbital period"),
                        ("galactic_orbital_phase_deg", "Galactic orbital phase"))
_DISTANCE = _rows("galactic", ("galactic_radius_pc", "Distance from the galactic centre"))

_NEBULA = (
    *_rows("identity", ("nebula_class", "Class", TEXT), ("nebula_type", "Type", TEXT),
           ("formation_cause", "Formation", TEXT)),
    *_rows("physical", ("radius_ly", "Radius"), ("mass_kg", "Mass")),
    *_rows("composition", ("composition", "Composition", TEXT), ("dominant_species", "Dominant gas", TEXT),
           ("hosts", "Stars it hosts", TEXT), ("age_myr", "Age of the stars"), ("q_total", "Ionizing photons")),
    *_rows("gas", ("density_cm3", "Density"), ("temperature_k", "Gas temperature"), ("extinction_av", "Extinction")),
    *_GALACTIC_ORBIT, *_DISTANCE,
)

_ASTEROID_FIELD = (
    *_rows("identity", ("field_class", "Class", TEXT), ("density", "Density", TEXT)),
    *_rows("physical", ("radius_ly", "Radius"), ("mass_kg", "Mass")),
    *_rows("composition", ("composition_family", "Composition family", TEXT),
           ("composition_summary", "Composition", TEXT)),
    *_GALACTIC_ORBIT, *_DISTANCE,
)

_BLACK_HOLE = (
    *_rows("identity", ("mass_class", "Class", TEXT)),
    *_rows("compact", ("mass_solar", "Mass"), ("event_horizon_radius_km", "Event horizon radius"), ("spin", "Spin"),
           ("has_accretion_disk", "Accretion disk", FLAG), ("temperature_k", "Temperature"),
           ("luminosity_w", "Luminosity"), ("age_gy", "Age")),
    *_rows("rotation", ("axial_tilt_deg", "Axial tilt")),
    *_GALACTIC_ORBIT, *_DISTANCE,
)

_NEUTRON_STAR = (
    *_rows("identity", ("pulsar_type", "Pulsar type", TEXT)),
    *_rows("compact", ("mass_solar", "Mass"), ("radius_km", "Radius"), ("magnetic_field_gauss", "Magnetic field"),
           ("surface_temperature_k", "Surface temperature"), ("luminosity_w", "Luminosity"), ("age_gy", "Age")),
    *_rows("rotation", ("spin_period_ms", "Spin period"), ("axial_tilt_deg", "Axial tilt")),
    *_GALACTIC_ORBIT, *_DISTANCE,
)

_SUPERNOVA_REMNANT = (
    *_rows("identity", ("remnant_class", "Class", TEXT), ("morphology", "Morphology", TEXT),
           ("progenitor_type", "Progenitor", TEXT), ("compact_remnant_kind", "Compact remnant left behind", TEXT)),
    *_rows("physical", ("age_years", "Age"), ("radius_ly", "Radius"), ("mass_kg", "Mass")),
    *_rows("gas", ("dominant_species", "Dominant gas", TEXT), ("density_cm3", "Density"),
           ("temperature_k", "Gas temperature"), ("extinction_av", "Extinction")),
    *_GALACTIC_ORBIT, *_DISTANCE,
)

_ROGUE_PLANET = (
    *_rows("identity", ("planet_type", "Type", TEXT), ("planet_class", "Class", TEXT),
           ("mass_bin", "Mass class", TEXT)),
    *_rows("physical", ("mass_kg", "Mass"), ("radius_km", "Radius", "radius"), ("age_gy", "Age"),
           ("has_moons", "Has moons", FLAG)),
    *_rows("composition", ("composition", "Composition", TEXT)),
    *_rows("climate", ("has_internal_heat", "Geologically active", FLAG),
           ("internal_heat_flux_w_m2", "Internal heat flow"), ("effective_temperature_k", "Effective temperature"),
           ("surface_regime", "Surface", TEXT), ("surface_temperature_k", "Surface temperature"),
           ("surface_pressure_pa", "Surface pressure")),
    *_rows("water", ("has_liquid_water", "Liquid water", FLAG), ("ice_shell_thickness_km", "Ice shell thickness"),
           ("ocean_depth_km", "Ocean depth"), ("hp_ice_km", "High-pressure ice thickness")),
    *_rows("rotation", ("rotation_period_hours", "Rotation period"), ("axial_tilt_deg", "Axial tilt")),
    *_GALACTIC_ORBIT, *_DISTANCE,
)

_INTERSTELLAR_COMET = (
    *_rows("identity", ("is_active", "Active", FLAG)),
    *_rows("physical", ("nucleus_diameter_km", "Nucleus diameter"), ("mass_kg", "Mass"),
           ("velocity_kms", "Speed")),
    *_rows("composition", ("composition_summary", "Composition", TEXT)),
    *_rows("rotation", ("rotation_period_hours", "Rotation period"), ("axial_tilt_deg", "Axial tilt")),
    *_GALACTIC_ORBIT, *_DISTANCE,
)

_QUASAR = (
    *_rows("compact", ("black_hole_mass_solar", "Black hole mass"),
           ("event_horizon_radius_km", "Event horizon radius")),
    *_rows("activity", ("luminosity_w", "Luminosity"), ("eddington_ratio", "Share of the Eddington limit"),
           ("accretion_rate_solar_per_year", "Gas swallowed"),
           ("broad_line_region_light_days", "Broad-line region radius"),
           ("is_radio_loud", "Radio-loud", FLAG), ("jet_length_ly", "Jet length"),
           ("active_age_years", "Active for")),
    *_DISTANCE,
)

SPECS = {"star": _STAR, "planet": _PLANET, "moon": _MOON, "belt": _BELT, "comet": _COMET, "nebula": _NEBULA,
         "asteroid_field": _ASTEROID_FIELD, "black_hole": _BLACK_HOLE, "neutron_star": _NEUTRON_STAR,
         "supernova_remnant": _SUPERNOVA_REMNANT, "rogue_planet": _ROGUE_PLANET,
         "interstellar_comet": _INTERSTELLAR_COMET, "quasar": _QUASAR}

"""dict: Each kind's ordered specs."""

_BOOKKEEPING = {
    "id": "the row key; the Object ID in the header is its public form",
    "uid": "the Object ID, shown in the header",
    "star_system_id": "the system, a link in the Location chain",
    "star_id": "the star it orbits, a link in the Location chain",
    "parent_type": "the parent kind; the Location chain is built from it",
    "parent_id": "the parent, a link in the Location chain",
    "epoch_unix": "the instant the orbit state was last updated; not a property of the body",
    "next_update_due": "when the orbit is next worked out; bookkeeping",
    "min_update_interval_years": "how often the orbit is worked out; bookkeeping",
    "galactic_min_update_interval_years": "how often the galactic orbit is worked out; bookkeeping",
    "name": "the page title",
    "orbital_index": "the order out from the star; the page's neighbours show it",
}
_VECTORS = "a component of the orbit state vector (position, velocity or spin axis); the Orbit and Rotation rows summarise it"
_EQUIPMENT = "a survival equipment level; the Required and Recommended lines of the description say it in words"


def _excuse(names, reason=None):
    """`{column: reason}` for each name, the reason given or the one `_BOOKKEEPING` holds."""
    return {name: reason or _BOOKKEEPING[name] for name in names}


_PHENOMENON_BOOKKEEPING = {
    "id": "the row key; the Object ID in the header is its public form",
    "uid": "the Object ID, shown in the header",
    "name": "the page title",
    "sector_id": "the sector, a link in the Location chain",
    "star_id": "the star of the system it belongs to, a link in the Location chain",
    "parent_type": "the parent kind; the Location chain is built from it",
    "parent_id": "the parent, a link in the Location chain",
    "inside_nebula_id": "the nebula it lies inside, a link in the Location chain",
    "inside_remnant_id": "the supernova remnant it lies inside, a link in the Location chain",
    "quadrant": "the octant, a badge in the page header",
    "cell_key": "the index cell it is filed under; bookkeeping",
    "created_at": "when the row was written; bookkeeping",
    "modified_at": "when the row last changed; bookkeeping",
    "epoch_unix": "the instant the galactic orbit was last updated; not a property of the object",
    "next_update_due": "when the galactic orbit is next worked out; bookkeeping",
    "galactic_min_update_interval_years": "how often the galactic orbit is worked out; bookkeeping",
    **{f"center_{axis}_pc": "its position in the galaxy frame; the Location block and the maps show it" for axis in "xyz"},
    **{f"velocity_{axis}_kms": "a component of its galactic motion; the galactic orbit rows summarise it" for axis in "xyz"},
    **{f"spin_axis_{axis}": _VECTORS for axis in "xyz"},
}
_SHAPE = "a parameter of the cloud's drawn 3D shape; bookkeeping for the 3D view"


def _phenomenon_not_shown(kind, **extra):
    """The bookkeeping columns of a phenomenon's table, each with why it is not a row, and the kind's own `extra`."""
    columns = models.metadata.tables[TABLES[kind]].columns
    return {**{name: reason for name, reason in _PHENOMENON_BOOKKEEPING.items() if name in columns}, **extra}


_BODY_NOT_SHOWN = {
    **_excuse(("id", "uid", "star_system_id", "star_id", "parent_type", "parent_id", "epoch_unix", "next_update_due",
               "min_update_interval_years", "name", "orbital_index")),
    "body_type": "terrestrial or giant, which the class says",
    "description": "the page's paragraph text",
    "atmosphere": "described in the paragraph text and the Atmosphere rows",
    "composition": "described in the paragraph text",
    "flavor_text_count": "how many flavour lines exist; bookkeeping",
    "flavor_key": "which flavour line the text uses; bookkeeping",
    "aerosol_draw": "the random draw behind the haze and dust depths; bookkeeping",
    "albedo_basis": "how the albedo and greenhouse draws were obtained (recorded, inverted or class draw); bookkeeping",
    **_excuse(("position_x_km", "position_y_km", "position_z_km", "velocity_x_kms", "velocity_y_kms", "velocity_z_kms",
               "spin_axis_x", "spin_axis_y", "spin_axis_z"), _VECTORS),
    **_excuse(("eq_body", "eq_body_kind", "eq_body_rec", "eq_body_rec_kind", "eq_rad_suit", "eq_rad_suit_rec", "eq_mobility",
               "eq_dust", "eq_gas_filter", "eq_co2_scrub", "eq_supply", "eq_supply_rec", "eq_shelter", "eq_shelter_rec"),
              _EQUIPMENT),
}

NOT_SHOWN = {
    "star": {
        **_excuse(("id", "uid", "star_system_id", "parent_type", "parent_id", "epoch_unix",
                   "galactic_min_update_interval_years", "name")),
        "parent_node_id": "the orbit tree node; bookkeeping",
        "slot": "the star's slot in its system; bookkeeping",
        "flags": "bit flags the generator keeps; each one is shown by the row it stands for",
        "cluster_id": "the cluster it belongs to, shown by the Location block",
        **_excuse(("reflex_offset_x_km", "reflex_offset_y_km", "reflex_offset_z_km", "spin_axis_x", "spin_axis_y",
                   "spin_axis_z", "velocity_x_kms", "velocity_y_kms", "velocity_z_kms"), _VECTORS),
    },
    "planet": {**_BODY_NOT_SHOWN, "host_node_id": "the orbit tree node; bookkeeping",
               **_excuse(("reflex_offset_x_km", "reflex_offset_y_km", "reflex_offset_z_km"), _VECTORS)},
    "moon": {**_BODY_NOT_SHOWN, "planet_id": "the planet it circles, a link in the Location chain",
             "rogue_planet_id": "the free planet it circles, a link in the Location chain"},
    "belt": {
        **_excuse(("id", "uid", "star_system_id", "star_id", "parent_type", "parent_id", "epoch_unix", "next_update_due",
                   "orbital_index")),
        **_excuse(("position_x_km", "position_y_km", "position_z_km", "velocity_x_kms", "velocity_y_kms",
                   "velocity_z_kms"), _VECTORS),
    },
    "comet": {
        **_excuse(("id", "uid", "star_system_id", "star_id", "parent_type", "parent_id", "epoch_unix", "next_update_due",
                   "min_update_interval_years", "name")),
        "parabolic_mean_anomaly": "the anomaly of an unbound orbit, which the Orbit rows summarise",
        **_excuse(("position_x_km", "position_y_km", "position_z_km", "velocity_x_kms", "velocity_y_kms", "velocity_z_kms",
                   "spin_axis_x", "spin_axis_y", "spin_axis_z"), _VECTORS),
    },
    "nebula": _phenomenon_not_shown(
        "nebula", **{name: _SHAPE for name in ("shape_axis_x", "shape_axis_y", "shape_axis_z", "shape_angle_x",
                                               "shape_angle_y", "shape_angle_z", "shape_warp_amplitude",
                                               "shape_warp_frequency", "shape_warp_octaves", "shape_noise_seed",
                                               "shape_iso", "shape_scale")}),
    "asteroid_field": _phenomenon_not_shown(
        "asteroid_field", formed_at="when the field formed, as a timestamp; bookkeeping",
        origin_event_id="the event that made it; bookkeeping"),
    "black_hole": _phenomenon_not_shown("black_hole"),
    "neutron_star": _phenomenon_not_shown("neutron_star"),
    "supernova_remnant": _phenomenon_not_shown(
        "supernova_remnant", compact_remnant_black_hole_id="the black hole it left behind, a link in its contents",
        compact_remnant_neutron_star_id="the neutron star it left behind, a link in its contents"),
    "rogue_planet": _phenomenon_not_shown("rogue_planet"),
    "interstellar_comet": _phenomenon_not_shown("interstellar_comet"),
    "quasar": _phenomenon_not_shown("quasar"),
}
"""dict: For each kind, the columns of its table the sheet does not show, each with why."""


def specs_for(kind):
    """The ordered specs of a kind; `ValueError` for a kind with no sheet."""
    if kind not in SPECS:
        raise ValueError(f"no property sheet for kind {kind!r}")
    return SPECS[kind]


_SYMBOLS = {
    "A*m2": " A·m²", "R_planet": " planet radii", "xuv-earth": "× Earth's", "1/yr": " a year", "1/Gy": " per billion years",
    "mSv/yr": " mSv/yr", "g/cm2": " g/cm²", "kg/mol": " kg/mol", "DU": " DU", "ppm": " ppm", "dex": " dex", "rad": " rad",
}
"""dict: How a unit the ladders do not cover is written after its number (with its own space, or none)."""


def _plain_number(value, unit):
    """A number the ladders do not cover, as text with its unit's symbol."""
    if unit in ("ratio", "unit-vector", "index", "count", "id", "key", "seed", "row-defined", "flag"):
        return format_number(value, ".4g")
    if unit in _SYMBOLS:
        symbol = _SYMBOLS[unit]
    else:
        symbol = unit_registry.UNITS[unit].symbol if unit in unit_registry.UNITS else unit
        symbol = symbol if symbol == "°" else f" {symbol}"
    return f"{format_number(value, '.4g')}{symbol}"


def _solar_ratio(kind, si):
    """A star's ratio to the Sun for its mass or luminosity."""
    return si / (constants.SOLAR_MASS_TO_KG if kind == "mass" else constants.SOLAR_LUMINOSITY)


def _value(spec, table, stored):
    """The row's value as an inline list (`describe.model`): a quantity, a number with its unit, or text."""
    if stored is None:
        return [NOT_STORED]
    quantity = spec.quantity
    if quantity == TEXT:
        text = str(stored)
        return [text[:1].upper() + text[1:]]
    if quantity == FLAG:
        return ["Yes" if stored else "No"]
    if quantity == TIER:
        names = habitability_world.TIER_NAMES
        return [names[int(stored)]] if 0 <= int(stored) < len(names) else [str(stored)]
    column = spec.stored_column()
    unit = column_units.COLUMNS[(table, column)]
    info = column_units.unit_info(unit)
    kind = quantity if quantity != AUTO else REGISTRY_QUANTITIES.get(info.kind)
    if kind is None or info.factor is None or (quantity == AUTO and unit not in unit_registry.UNITS):
        return [_plain_number(stored, unit)]
    value = float(stored)
    if spec.solar:
        si = value * info.factor + getattr(info, "offset", 0.0)
        return [model.qty(kind, value, unit, sol=_solar_ratio(kind, si))]
    return [model.qty(kind, value, unit)]


def sheet(conn, kind, object_id):
    """
    The property sheet of one body.

    Args:
        conn: An open, read-only connection.
        kind (str): A key of `SPECS`: `star`, `planet`, `moon`, `belt`, `comet` or a phenomenon kind.
        object_id (int): The row id.

    Returns:
        dict: `kind`, `groups`: the non-empty groups in `GROUPS` order, each `{"key", "label", "rows"}`; a row is
            `{"key", "label", "value" (the inline list), "text" (it as plain text), "stored" (whether the column held a
            value)}`.

    Raises:
        ValueError: For a kind with no sheet or a missing row.
    """
    specs = specs_for(kind)
    table = TABLES[kind]
    columns = ", ".join(sorted({spec.stored_column() for spec in specs}))
    row = conn.execute(f"SELECT {columns} FROM {table} WHERE id = ?", (object_id,)).fetchone()
    if row is None:
        raise ValueError(f"no {table} row with id {object_id}")
    by_group = {key: [] for key in GROUPS}
    for spec in specs:
        stored = row[spec.stored_column()]
        value = _value(spec, table, stored)
        by_group[spec.group].append({"key": spec.key, "label": spec.label, "value": value,
                                     "text": inline_text(value), "stored": stored is not None})
    groups = [{"key": key, "label": GROUPS[key], "rows": rows} for key, rows in by_group.items() if rows]
    return {"kind": kind, "groups": groups}
