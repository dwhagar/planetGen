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

from planetgen import tuning
from planetgen.datasets import fetch
from planetgen.db import column_units, models
from planetgen.describe import model
from planetgen.describe.quantity import inline_text
from planetgen.describe.sheet_sources import source_of
from planetgen.galaxy import scales
from planetgen.physics import (agb_carbon, chemically_peculiar, compaction, constants, daylight, gases, habitability_world, latitude_bands, outside_hours, regional_temperature,
                               stellar_abundances, stellar_composition, surface_liquid_world, unit_registry, units, variability)
from planetgen.population import facility_ratings
from planetgen.util.format import format_number, scientific_text

TABLES = {"star": "stars", "planet": "planets", "moon": "moons", "belt": "asteroid_belts", "comet": "comets",
          "nebula": "nebulae", "asteroid_field": "asteroid_fields", "black_hole": "black_holes",
          "neutron_star": "neutron_stars", "supernova_remnant": "supernova_remnants", "rogue_planet": "rogue_planets",
          "interstellar_comet": "interstellar_comets", "quasar": "quasars", "facility": "facilities",
          "globular_cluster": "galaxy_containers", "superbubble": "galaxy_containers"}
"""dict: The table each kind of the sheet is stored in: a body of a system, a phenomenon (UX.175), a facility (UX.183)
or a galaxy container (UX.184; a cluster and a superbubble share a table)."""

GROUPS = {
    "identity": "Identity", "orbit": "Orbit", "physical": "Physical", "rotation": "Rotation",
    "composition": "Composition", "atmosphere": "Atmosphere", "climate": "Climate", "water": "Water",
    "magnetism": "Magnetic field", "radiation": "Radiation", "habitability": "Habitability", "life": "Life",
    "zones": "Zones and limits", "stellar": "Stellar activity", "galactic": "Galactic orbit", "disk": "Disk",
    "radiogenic": "Radiogenic heat", "gas": "Gas and light", "compact": "Compact object", "activity": "Activity",
    "ratings": "Ratings",
}
"""dict: Every group key of a sheet and its heading, in the order the page shows them."""

SHOW_GROUPS = frozenset({"identity", "habitability", "life", "ratings"})
"""frozenset: The groups whose rows are listed on a page even when the value is none or 0 (UX.188): a world with no life
or a facility with no ratings set is itself the answer."""

_PLANET_ABSENT = frozenset({
    "atmospheric_pressure_pa", "outside_hours_day", "hydrosphere", "ocean_fraction", "water_mass_fraction",
    "dipole_class", "magnetic_moment_a_m2",
})
SHOW_WHEN_ABSENT = {
    "planet": _PLANET_ABSENT, "moon": _PLANET_ABSENT,
    "star": frozenset({"flare_n33_per_yr", "lethal_event_rate_per_gyr"}),
}
"""dict: For each kind, the row keys a page lists even when the value is none or 0, because the absence matters to a
captain (no air, no magnetic field, a quiet star); every other row is left out of a page when it holds nothing. The
API sheet always returns every row and says which are `absent` and which are `show_absent`."""

GLANCE = {
    "planet": ("planet_class", "zone", "surface_temperature_k", "flux_variation", "atmospheric_pressure_pa", "gravity_g",
               "surface_dose_msv_yr", "outside_hours_day", "rotation_period_hours", "daylight_hours", "day_share_equator", "distance_km", "period_years"),
    "moon": ("planet_class", "zone", "surface_temperature_k", "flux_variation", "atmospheric_pressure_pa", "gravity_g",
             "surface_dose_msv_yr", "outside_hours_day", "rotation_period_hours", "daylight_hours", "day_share_equator", "distance_km", "period_years"),
    "star": ("role", "star_type", "yerkes_class", "state", "mass_kg", "luminosity_w", "temperature_k",
             "habitable_zone_inner_km", "habitable_zone_outer_km", "snow_line_km", "flare_n33_per_yr",
             "lethal_event_rate_per_gyr", "heliosphere_radius_km", "system_perimeter_km"),
    "black_hole": ("mass_class", "mass_solar", "event_horizon_radius_km", "spin", "has_accretion_disk", "luminosity_w"),
    "neutron_star": ("pulsar_type", "mass_solar", "radius_km", "magnetic_field_gauss", "spin_period_ms"),
    "nebula": ("nebula_class", "nebula_type", "radius_ly", "mass_kg", "temperature_k", "density_cm3"),
    "supernova_remnant": ("remnant_class", "age_years", "radius_ly", "compact_remnant_kind", "temperature_k"),
    "rogue_planet": ("planet_class", "mass_kg", "radius_km", "surface_temperature_k", "surface_pressure_pa",
                     "has_liquid_water"),
    "asteroid_field": ("field_class", "radius_ly", "density", "mass_kg", "composition_summary"),
    "interstellar_comet": ("is_active", "nucleus_diameter_km", "velocity_kms", "composition_summary"),
    "quasar": ("luminosity_w", "black_hole_mass_solar", "jet_length_ly", "is_radio_loud", "event_horizon_radius_km"),
    "belt": ("distance_km", "density", "mass_kg", "composition_summary"),
    "comet": ("orbit_type", "is_active", "nucleus_diameter_km", "perihelion_distance_km", "orbital_period_years",
              "composition_summary"),
    "facility": ("kind", "placement", "host_type", "affiliation_name", "rating_crime", "rating_housing",
                 "rating_resources", "rating_maintenance", "rating_health"),
    "globular_cluster": ("extent_pc",),
    "superbubble": ("extent_pc",),
}
"""dict: For each kind, the row keys of its "At a glance" panel (UX.188), in a captain's order: what it is, can I
survive there, how long is a day and how far is it. Each is a row of the same sheet."""

SWING_FLAG_MIN = 0.75
"""float: The peak-to-peak swing of a body's starlight (`flux_variation`, about 4 times its orbit's eccentricity) from which the "At a glance"
panel flags it as very high (GEN.1135, which retired Class Q: a class that was only an orbit). 4 percent of bodies; 0.5 would flag 12 percent and the
old class line of 0.25 would flag 36 percent. A default until Boss gives a line."""

SWING_FLAG_LABEL = "Seasonal swing of starlight"
"""str: What the panel calls the flag row (the sheet's own row for the number is "Seasonal variation of that energy, peak to peak")."""

NOT_STORED = "not stored"
"""str: What a row with no stored value says."""

REGISTRY_QUANTITIES = {
    "length": "length", "mass": "mass", "duration": "time", "temperature": "temperature", "pressure": "pressure",
    "speed": "speed", "delta_temperature": "temperature difference", "fraction": "fraction", "power": "power", "acceleration": "gravity", "flux": "flux", "ratio": "ratio",
}
"""dict: A unit registry kind, to the Document quantity kind that shows it with its ladder. Any other kind is written
as a number and the unit's symbol."""

TEXT, FLAG, TIER, RATING, LIMIT, AUTO = "text", "flag", "tier", "rating", "limit", "auto"
"""str: How a column is shown besides a ladder quantity: as stored text, as Yes or No, as a PHI-6 colour, as a facility
rating (1 Green to 3 Red, `population/facility_ratings.py`), as what limits the hours a person can be outside (`LIMIT`, GEN.768), or from
its unit (`AUTO`)."""


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
           ("decay_age_yr", "Age since it stopped burning"),
           ("nuclear_arrival_age_gy", "Age it reaches hydrogen burning")),
    *_rows("zones", ("habitable_zone_inner_km", "Habitable zone, inner edge"),
           ("habitable_zone_outer_km", "Habitable zone, outer edge"), ("snow_line_km", "Snow line"),
           ("heliosphere_radius_km", "Heliosphere radius"), ("system_perimeter_km", "System perimeter"),
           ("wide_binary_a_crit_km", "Widest stable orbit around a pair")),
    *_rows("rotation", ("rotation_period_hours", "Rotation period"), ("axial_tilt_deg", "Axial tilt")),
    *_rows("stellar", ("log_lx_lbol", "X-ray to bolometric luminosity, log"), ("l_xuv_w", "X-ray and ultraviolet power"),
           ("xuv_saturated", "X-ray output saturated", FLAG), ("xuv_fluence_j", "X-ray and ultraviolet energy emitted"),
           ("uv_c_share", "UV-C power against its blackbody's"),
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
           ("orbital_speed_kms", "Orbital speed"), ("eccentricity", "Eccentricity"), ("orbital_inclination_deg", "Inclination"),
           ("orbital_ascending_node_deg", "Longitude of the ascending node"), ("hill_radius_km", "Hill radius"),
           ("min_orbit_distance_km", "Closest stable orbit for a moon"), ("dynamical_state", "Dynamical state", TEXT),
           ("habitable_zone_inner_km", "Habitable zone, inner edge"),
           ("habitable_zone_outer_km", "Habitable zone, outer edge"),
           ("moon_origin", "How a moon came to orbit its planet", TEXT)),
    *_rows("disk", ("moon_disk_lifetime_myr", "Moon-forming disk lifetime"),
           ("moon_disk_alpha_over_f", "Moon-forming disk turbulence over infall (alpha/f)"),
           ("moon_disk_f_inf", "Moon-forming disk ice fraction far past its ice line"),
           ("moon_disk_ice_line_km", "Moon-forming disk ice line")),
    *_rows("physical", ("radius_km", "Radius", "radius"), ("mass_kg", "Mass"), ("volume_km3", "Volume"),
           ("density_g_cm3", "Density"), ("gravity_g", "Surface gravity")),
    *_rows("rotation", ("rotation_period_hours", "Rotation period"), ("axial_tilt_deg", "Axial tilt"),
           ("outside_hours", "Hours a person can work outside per year (dose limit)"),
           ("outside_hours_day", "Hours a person can be outside in a day"), ("outside_hours_limit", "What limits those hours", LIMIT),
           ("solar_day_hours", "Length of the solar day"),
           ("day_night_state", "State of its day and night", TEXT), ("days_per_year", "Solar days in one year"),
           ("night_hours", "Hours of night in a solar day, equator"),
           ("daylight_hours", "Hours of daylight at the equator"), ("day_share_equator", "Share of the time the star is up, equator"),
           ("day_share_mean", "Share of the time the star is up, whole surface"),
           ("parent_shadow_share", "Share of the time in the planet's shadow"),
           ("polar_circle_deg", "Latitude of the polar circle")),
    *_rows("composition", ("mantle_redox", "Mantle oxidation", TEXT), ("mantle_delta_iw", "Mantle oxidation, log below iron-wüstite"),
           ("phosphorus", "Phosphorus supply", TEXT), ("born_beyond_snow_line", "Formed beyond the snow line", FLAG),
           ("metal_mass_fraction", "Metal (iron) mass fraction"), ("rock_mass_fraction", "Rock mass fraction"),
           ("ice_mass_fraction", "Ice mass fraction"), ("ice_co2_mass_fraction", "Carbon dioxide ice mass fraction"),
           ("ice_ch4_mass_fraction", "Methane ice mass fraction"), ("ice_nh3_mass_fraction", "Ammonia ice mass fraction"),
           ("ice_n2_mass_fraction", "Nitrogen ice mass fraction"), ("ice_so2_mass_fraction", "Sulfur dioxide ice mass fraction"),
           ("ice_co_mass_fraction", "Carbon monoxide ice mass fraction"), ("ice_hcn_mass_fraction", "Hydrogen cyanide ice mass fraction"),
           ("ice_other_mass_fraction", "Other ices mass fraction"), ("envelope_mass_fraction", "Gas envelope mass fraction"),
           ("envelope_helium_fraction", "Helium share of the envelope"), ("envelope_z_fraction", "Heavy-element share of the envelope"),
           ("rock_feo_fraction", "Iron oxide share of the rock"), ("rock_al2o3_fraction", "Alumina share of the rock"),
           ("rock_mg_si_ratio", "Magnesium to silicon ratio of the rock"), ("porosity", "Porosity"),
           ("carbon_rock", "Carbon-rich rock", FLAG)),
    *_rows("atmosphere", ("atmospheric_pressure_pa", "Surface pressure"), ("atm_density", "Air density"),
           ("atm_molar_density", "Mean molar mass of the air"), ("scale_height_km", "Scale height"),
           ("tau_haze", "Haze optical depth"), ("tau_dust", "Dust optical depth"), ("tau_cloud", "Cloud optical depth"),
           ("cloud_cover", "Cloud cover"), ("ozone_column_du", "Ozone column"),
           ("ozone_loss_flag", "Ozone layer lost", FLAG)),
    *_rows("climate", ("bond_albedo", "Bond albedo"), ("greenhouse_factor", "Greenhouse factor"),
           ("surface_temperature_k", "Mean surface temperature"),
           ("temp_mean_swing_k", "Day and night swing, whole surface (hottest to coldest)"), ("energy_flux_w_m2", "Stellar energy received"),
           ("flux_variation", "Seasonal variation of that energy, peak to peak")),
    *_rows("water", ("hydrosphere", "Hydrosphere", TEXT), ("water_mass_fraction", "Water share of the mass"),
           ("ocean_fraction", "Surface under ocean"), ("land_fraction", "Surface above sea level"),
           ("ocean_depth_km", "Ocean depth"), ("ice_shell_km", "Ice shell thickness"),
           ("surface_liquid", "Surface liquid", TEXT), ("surface_liquid_fraction", "Surface under that liquid"),
           ("surface_liquid_depth_km", "Depth of that liquid"),
           ("surface_liquid_acid_wt", "Sulfuric acid strength (liquid acid)"), ("surface_liquid_acid_layer", "Where the liquid acid is", TEXT),
           ("surface_liquid_ammonia_wt", "Ammonia in the ocean (share of the liquid)"),
           ("hp_ice_km", "High-pressure ice thickness"), ("ocean_class", "Ocean chemistry", TEXT),
           ("ocean_ph", "Ocean pH"), ("water_activity", "Water activity"),
           ("water_base_pressure_pa", "Pressure at the bottom of the water layer"),
           ("ocean_liquid", "Ocean liquid", TEXT),
           ("tide_range_spring_m", "Open-ocean spring tide range (equilibrium value; a coast can differ many times)"),
           ("tide_range_neap_m", "Open-ocean neap tide range (equilibrium value)"),
           ("tidal_flux_w_m2", "Tidal heat flow"), ("aquifer_depth_km", "Depth to underground water"),
           ("aquifer_thickness_km", "Thickness of the underground water")),
    *_rows("magnetism", ("dipole_class", "Dipole field", TEXT), ("magnetic_moment_a_m2", "Magnetic moment"),
           ("induced_moment_a_m2", "Induced magnetic moment"), ("magnetopause_rp", "Magnetopause distance")),
    *_rows("radiation", ("surface_dose_msv_yr", "Surface dose"), ("dose_gcr_msv_yr", "Dose from cosmic rays"),
           ("dose_sep_msv_yr", "Dose from stellar particles"), ("dose_ground_msv_yr", "Dose from the ground"),
           ("dose_belt_msv_yr", "Dose from the parent planet's radiation belts"),
           ("dose_helio_mult", "Cosmic-ray boost from a squeezed heliosphere"), ("xuv_flux_earth", "X-ray and ultraviolet flux"),
           ("xuv_exposure_index", "X-ray and ultraviolet exposure index"),
           ("flare_irradiation_index", "Flare irradiation index"), ("uv_surface_index", "Surface ultraviolet index"),
           ("shelter_column_g_cm2", "Shelter column")),
    *_rows("habitability", ("phi4", "Habitability score (PHI-6)"), ("phi4_pressure", "Pressure score"),
           ("phi4_temperature", "Temperature score"), ("phi4_air", "Air score"), ("phi4_water", "Water score"),
           ("phi4_suitability", "Suitability score"), ("phi4_radiation", "Radiation score"),
           ("tier_pressure", "Pressure rating", TIER), ("tier_temperature", "Temperature rating", TIER),
           ("tier_air", "Air rating", TIER), ("tier_water", "Water rating", TIER),
           ("tier_suitability", "Suitability rating", TIER), ("tier_radiation", "Radiation rating", TIER)),
    *_rows("life", ("life_class", "Life class number", TEXT), ("life_class_variant", "Life class letter", TEXT),
           ("life_class_best", "Best-supported life class", TEXT), ("life_class_best_score", "Best-supported class score"),
           ("life_animal_grade", "Animal grade"), ("life_complexity", "Complexity index"), ("evolutionary_pace", "Evolutionary pace", TEXT), ("life_stage", "Life stage", TEXT),
           ("life_milestone_age_gyr", "Age of that stage"), ("life_system_age_gyr", "Age of the system when it was written"),
           ("life_age_abiogenesis_gyr", "Age of abiogenesis"), ("life_age_photosynthesis_gyr", "Age of photosynthesis"),
           ("life_age_complex_cells_gyr", "Age of complex cells"), ("life_age_multicellularity_gyr", "Age of multicellularity"),
           ("life_age_technological_civilization_gyr", "Age of technological civilization"),
           ("life_microbial", "Microbial life score"), ("life_flora", "Flora score"), ("life_fauna", "Fauna score"),
           ("life_sentient", "Sentient life score"),
           ("l_solv", "Solvent suitability"), ("l_chem", "Chemistry suitability"), ("l_ener", "Energy suitability"),
           ("l_rad", "Radiation tolerance")),
    *_rows("radiogenic", ("c_k0_ppm", "Potassium at formation"), ("c_th0_ppm", "Thorium at formation"),
           ("c_u0_ppm", "Uranium at formation"), ("al26_ratio", "Aluminium-26 at formation"), ("urey_ratio", "Share of the interior heat that is radiogenic (Urey ratio)"),
           ("internal_heat_flux_w_m2", "Heat flowing out of the interior"), ("mantle_eta0_pa_s", "Mantle viscosity prefactor"),
           ("mantle_activation_kj_mol", "Mantle viscosity activation energy"), ("lid_regime", "Tectonic regime", TEXT)),
    *_rows("composition", ("interior_kind", "Kind of interior", TEXT), ("core_mass_earth", "Compact core mass"),
           ("dilute_core_mass_earth", "Dilute core mass"), ("metallization_pressure_gpa", "Pressure where hydrogen turns metallic"),
           ("helium_depletion", "Helium left in the molecular envelope, against its natal share"), ("has_helium_shell", "Helium-rich shell", FLAG)),
    *_rows("composition", ("ring_kind", "Kind of ring system", TEXT), ("ring_origin", "Where the ring came from", TEXT),
           ("ring_inner_km", "Ring inner edge"), ("ring_outer_km", "Ring outer edge"), ("ring_mass_kg", "Ring mass"),
           ("ring_optical_depth", "Ring optical depth"), ("ring_ice_fraction", "Ring ice share"), ("ring_age_yr", "Ring age")),
)

_PLANET = _MOON = _BODY_GROUPS

_BELT = (
    *_rows("orbit", ("distance_km", "Distance from its star"), ("lower_limit_km", "Inner edge"),
           ("upper_limit_km", "Outer edge")),
    *_rows("physical", ("density", "Density", TEXT), ("mass_kg", "Mass")),
    *_rows("composition", ("composition_summary", "Composition", TEXT), ("chondrite_group", "Nearest meteorite group", TEXT)),
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
           ("composition_summary", "Composition", TEXT), ("chondrite_group", "Nearest meteorite group", TEXT)),
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
    *_rows("habitability", ("phi4", "Habitability score (PHI-6)"), ("phi4_pressure", "Pressure score"),
           ("phi4_temperature", "Temperature score"), ("phi4_air", "Air score"), ("phi4_water", "Water score"),
           ("phi4_suitability", "Suitability score"), ("phi4_radiation", "Radiation score"),
           ("tier_pressure", "Pressure rating", TIER), ("tier_temperature", "Temperature rating", TIER),
           ("tier_air", "Air rating", TIER), ("tier_water", "Water rating", TIER),
           ("tier_suitability", "Suitability rating", TIER), ("tier_radiation", "Radiation rating", TIER),
           ("energy_flux_w_m2", "Chemical energy from its own heat")),
    *_rows("life", ("life_microbial", "Microbial life score"), ("life_flora", "Flora score"),
           ("life_fauna", "Fauna score"), ("life_sentient", "Sentient life score"),
           ("l_solv", "Solvent suitability"), ("l_chem", "Chemistry suitability"), ("l_ener", "Energy suitability"),
           ("l_rad", "Radiation tolerance")),
    *_rows("radiation", ("shelter_column_g_cm2", "Shelter column")),
    *_rows("composition", ("ring_kind", "Kind of ring system", TEXT), ("ring_origin", "Where the ring came from", TEXT),
           ("ring_inner_km", "Ring inner edge"), ("ring_outer_km", "Ring outer edge"), ("ring_mass_kg", "Ring mass"),
           ("ring_optical_depth", "Ring optical depth"), ("ring_ice_fraction", "Ring ice share"), ("ring_age_yr", "Ring age")),
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

_FACILITY = (
    *_rows("identity", ("kind", "Kind", TEXT), ("placement", "Placement", TEXT), ("host_type", "Built on", TEXT),
           ("affiliation_name", "Affiliation", TEXT)),
    *_rows("physical", ("mass_kg", "Mass")),
    *_rows("orbit", ("orbit_distance_km", "Orbit distance"), ("orbit_period_years", "Orbital period"),
           ("orbital_speed_kms", "Orbital speed"), ("orbit_phase_deg", "Orbital phase")),
    *_rows("ratings", ("rating_crime", "Crime", RATING), ("rating_housing", "Housing", RATING),
           ("rating_resources", "Resources", RATING), ("rating_maintenance", "Maintenance", RATING),
           ("rating_health", "Health", RATING), ("ratings_source", "Rated by", TEXT), ("ratings_note", "Rating note", TEXT)),
    *_DISTANCE,
)

_CONTAINER = (
    *_rows("physical", ("extent_pc", "Radius of its sphere")),
    *_rows("identity", ("shell_index", "Birth shell"), ("band_index", "Birth band"), ("wedge_index", "Birth wedge")),
)

SPECS = {"globular_cluster": _CONTAINER, "superbubble": _CONTAINER, "facility": _FACILITY, "star": _STAR, "planet": _PLANET, "moon": _MOON, "belt": _BELT, "comet": _COMET, "nebula": _NEBULA,
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
    "composition_origin": "how the generator found the composition; bookkeeping, not a fact about the world",
    "description": "the page's paragraph text",
    "atmosphere": "described in the paragraph text and the Atmosphere rows",
    "composition": "described in the paragraph text",
    "flavor_text_count": "how many flavour lines exist; bookkeeping",
    "flavor_key": "which flavour line the text uses; bookkeeping",
    "aerosol_draw": "the random draw behind the haze and dust depths; bookkeeping",
    "atmosphere_colour": "shown as the Atmosphere Colour swatch by the description builder, a colour rather than a number",
    "albedo_basis": "how the albedo and greenhouse draws were obtained (recorded, inverted or class draw); bookkeeping",
    **_excuse(("ground_type", "ground_grain", "ground_ocean_fraction", "ground_ice_fraction", "ground_ultramafic_fraction",
               "ground_anorthositic_fraction", "ground_basaltic_fraction", "ground_felsic_fraction", "ground_clay_fraction",
               "ground_oxidised_fraction", "ground_metal_fraction", "ground_carbon_fraction", "ground_evaporite_fraction",
               "ground_sulfur_fraction", "ground_intermediate_fraction"),
              "the stored mix of ground types; it feeds the ground dose and the surface rendering, and the page does "
              "not list it as a row yet"),
    **_excuse(regional_temperature.BAND_COLUMNS,
              "the seven bands of latitude (the rings of a world that turns once an orbit), side by side in the table of the bands: the annual mean, the day and night swing, "
              "the hours of day and night, the longest polar day and night and the state of each"),
    **_excuse(("position_x_km", "position_y_km", "position_z_km", "velocity_x_kms", "velocity_y_kms", "velocity_z_kms",
               "spin_axis_x", "spin_axis_y", "spin_axis_z"), _VECTORS),
    **_excuse(("eq_body", "eq_body_kind", "eq_body_rec", "eq_body_rec_kind", "eq_rad_suit", "eq_rad_suit_rec", "eq_mobility",
               "eq_dust", "eq_gas_filter", "eq_co2_scrub", "eq_supply", "eq_supply_rec", "eq_shelter", "eq_shelter_rec"),
              _EQUIPMENT),
}

def _container_not_shown():
    """The columns of `galaxy_containers` a cluster's or a superbubble's sheet leaves out."""
    return {
        **_excuse(("id", "uid", "name")),
        "kind": "the kind of the page; the header and the About text name it",
        "component_id": "the cluster's galaxy component; the About text gives its mass, age and metallicity",
        "bubble_col": "the superbubble's cell in the bubble grid; bookkeeping",
        "bubble_row": "the superbubble's cell in the bubble grid; bookkeeping",
        "bubble_n": "the superbubble's number in its cell; bookkeeping",
        "created_at": "when the row was written; bookkeeping",
        **{f"center_{axis}_pc": "its position in the galaxy frame; the Location block shows it" for axis in "xyz"},
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
    "rogue_planet": _phenomenon_not_shown(
        "rogue_planet", **_excuse(("eq_body", "eq_body_kind", "eq_body_rec", "eq_body_rec_kind", "eq_rad_suit", "eq_rad_suit_rec",
                                  "eq_mobility", "eq_dust", "eq_gas_filter", "eq_co2_scrub", "eq_supply", "eq_supply_rec",
                                  "eq_shelter", "eq_shelter_rec"), _EQUIPMENT)),
    "interstellar_comet": _phenomenon_not_shown("interstellar_comet"),
    "quasar": _phenomenon_not_shown("quasar"),
    "globular_cluster": _container_not_shown(),
    "superbubble": _container_not_shown(),
    "facility": {
        **_excuse(("id", "uid", "name", "star_system_id", "star_id", "epoch_unix", "next_update_due")),
        "planet_id": "the planet it stands on, a link in the Location chain",
        "moon_id": "the moon it stands on, a link in the Location chain",
        "asteroid_belt_id": "the belt it is built in, a link in the Location chain",
        "asteroid_field_id": "the asteroid field it is built in, a link in the Location chain",
        "sector_id": "the sector, a link in the Location chain",
        "type_id": "the facility type, named in the page's About text",
        "affiliation_polity_id": "the polity it belongs to; the Affiliation row names it",
        "description": "the page's paragraph text",
        "ratings_set_at": "when the ratings were set; bookkeeping",
        "cell_key": "the index cell it is filed under; bookkeeping",
        "created_at": "when the row was written; bookkeeping",
        "modified_at": "when the row last changed; bookkeeping",
        **{f"center_{axis}_pc": "its position in the galaxy frame; the Location block and the maps show it" for axis in "xyz"},
        **{f"velocity_{axis}_kms": "a component of its galactic motion; the orbit rows summarise it" for axis in "xyz"},
    },
}
"""dict: For each kind, the columns of its table the sheet does not show, each with why."""


def specs_for(kind):
    """The ordered specs of a kind; `ValueError` for a kind with no sheet."""
    if kind not in SPECS:
        raise ValueError(f"no property sheet for kind {kind!r}")
    return SPECS[kind]


_SYMBOLS = {
    "A*m2": " A·m²", "R_planet": " planet radii", "xuv-earth": "× Earth's", "1/yr": " a year", "1/Gy": " per billion years",
    "mSv/yr": " mSv/yr", "kg/mol": " kg/mol", "DU": " DU", "ppm": " ppm", "dex": " dex", "rad": " rad",
}
"""dict: How a unit the ladders do not cover is written after its number (with its own space, or none)."""


BARE_UNITS = ("ratio", "unit-vector", "index", "count", "id", "key", "seed", "row-defined", "flag")
"""tuple: Column units that say nothing about what a number is: a sheet row declared in one of them is a bug (GEN.771),
so `_plain_number` refuses it."""


def _plain_number(value, unit):
    """A number the ladders do not cover, as text with its unit's symbol, or the scale or noun a plain number's unit
    declares in `column_units.LOCAL_UNITS` (`prefix`, `suffix`; GEN.770)."""
    local = column_units.LOCAL_UNITS.get(unit)
    if local is not None and (local.prefix or local.suffix):
        number = scientific_text(value) if value and abs(value) < 1e-3 else format_number(value, ".4g")
        return f"{local.prefix}{number}{local.suffix}"
    if unit in BARE_UNITS:
        raise ValueError(f"a column of unit {unit!r} would print a bare number; declare it in a unit that says what it is "
                         f"(a share, a score of n, a scale or a noun) in db/column_units.py")
    if unit in _SYMBOLS:
        symbol = _SYMBOLS[unit]
    else:
        symbol = unit_registry.UNITS[unit].symbol if unit in unit_registry.UNITS else unit
        symbol = symbol if symbol == "°" else f" {symbol}"
    return f"{format_number(value, '.4g')}{symbol}"


NO_DAY = "the same through the day: a world that turns once an orbit has no day and night"
"""str: What the temperature row of a tidally locked world says (GEN.848): its places are fixed and their swing is 0."""


def _plus_minus(swing):
    """The swing as a clear range: the half of it either side, and the whole of it from the coldest hour to the hottest
    (GEN.795). A difference in kelvin, so it is written as one and never converted as a temperature."""
    return ["\u00b1", model.qty("temperature difference", swing / 2.0, "dK", companions=False), ", ",
            model.qty("temperature difference", swing, "dK", companions=False),
            " in all from the coldest hour to the hottest"]


def _range_value(spec, table, stored, row):
    """`_value`, with the day and night range of the whole surface (GEN.795, GEN.796): the mean swing as plus or minus, and on the mean
    temperature the coldest and the hottest hour with the stored mean exactly between them (a locked world says it has no day).
    `row` holds the stored columns of the body; the places are the bands, in the table of the bands."""
    value = _value(spec, table, stored)
    if stored is None or table not in ("planets", "moons"):
        return value
    column = spec.stored_column()
    if column == "temp_mean_swing_k":
        return _plus_minus(stored) if stored else value
    if column != "surface_temperature_k":
        return value
    swing = row.get("temp_mean_swing_k")
    span = regional_temperature.daily_range(stored, swing)
    if span is None:
        return value + [f" ({NO_DAY})"] if swing == 0 else value
    return value + [", from ", model.qty("temperature", span[0], "K", companions=False), " to ",
                    model.qty("temperature", span[1], "K", companions=False), " over a day"]


def _solar_ratio(kind, si):
    """A star's ratio to the Sun for its mass or luminosity."""
    return si / (constants.SOLAR_MASS_TO_KG if kind == "mass" else constants.SOLAR_LUMINOSITY)


def _value(spec, table, stored):
    """The row's value as an inline list (`describe.model`): a quantity, a number with its unit, or text."""
    if stored is None:
        return [NOT_STORED]
    quantity = spec.quantity
    if quantity == TEXT:
        text = str(stored).replace("_", " ")    # a stored key such as technological_civilization reads as words
        return [text[:1].upper() + text[1:]]
    if quantity == FLAG:
        return ["Yes" if stored else "No"]
    if quantity == RATING:
        return [facility_ratings.COLORS[int(stored) - 1].capitalize()] if 1 <= int(stored) <= 3 else [str(stored)]
    if quantity == LIMIT:
        return [outside_hours.LIMIT_WORDS.get(stored, str(stored))]
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


VACUUM_LIMIT_PA = 1e-3
"""float: A surface pressure below this is no atmosphere for a page: it says "None" and lists no gases (UX.203)."""

NO_WATER_MASS_FRACTION = 1e-6
"""float: A dry body holding less water than this share of its mass has no water for a page (UX.203)."""

NO_OCEAN_FRACTION = 1e-4
"""float: A dry body with less of its surface under ocean than this has no water for a page (UX.203)."""

LIFE_SCORES = ("life_microbial", "life_flora", "life_fauna", "life_sentient")
"""tuple: The four life scores of a body (0 to 9); a body with all four stored and 0 has no life."""

MISSING_SENTENCES = {
    "life": "No life has been found here.",
    "atmosphere": "It has no atmosphere.",
    "water": "It has effectively no water.",
}
"""dict: The sentence a page says for each section it leaves out (`missing`), in the order the page says them."""


def missing(kind, stored):
    """
    The sections of a body's page that hold only zeros (UX.203), as a list of keys of `MISSING_SENTENCES`: `life` (all
    four life scores stored and 0), `atmosphere` (a stored surface pressure under `VACUUM_LIMIT_PA`) and `water` (a
    dry hydrosphere with under `NO_WATER_MASS_FRACTION` of water by mass and under `NO_OCEAN_FRACTION` of the surface
    under ocean). A value that is not stored never counts as none. Only planets and moons have these sections.

    Args:
        kind (str): The kind of the sheet.
        stored (dict): The stored column values by sheet key.
    """
    if kind not in GAS_TABLES:
        return []
    found = []
    scores = [stored.get(key) for key in LIFE_SCORES]
    if all(score is not None for score in scores) and not any(scores):
        found.append("life")
    pressure = stored.get("atmospheric_pressure_pa")
    if pressure is not None and pressure < VACUUM_LIMIT_PA:
        found.append("atmosphere")
    if (stored.get("hydrosphere") == "dry" and (stored.get("water_mass_fraction") or 0) < NO_WATER_MASS_FRACTION
            and (stored.get("ocean_fraction") or 0) < NO_OCEAN_FRACTION):
        found.append("water")
    return found


GIANT_CLASSES = frozenset("IJT")
"""frozenset: The planet classes that are giants (ice giant, gas giant, gas dwarf): no surface, a deep envelope."""

NEGLIGIBLE_TIDE_ECCENTRICITY = 1e-12
"""float: The eccentricity below which a body's tide heat is taken as nil and is not stored (`interior_heat.update_heat`): the heat goes as the
square of it, so at 1e-12 it is 1e-24 of the heat at an eccentricity of 1, and the circularised orbits the generator draws are far smaller still (1e-170)."""
HP_ICE_MIN_BASE_PA = 6.2e8
"""float: The pressure at the base of a body's water, Pa, from which ice VI (the first high-pressure ice that `water_layer.column` counts in `hp_ice_km`)
forms: 620 MPa, where ice VI melts at 272 K (`water_layer.LIQUIDUS_K_MPA`)."""
POROSITY_MAX_MASS_KG = compaction.TAPER_END_KG
"""float: The mass, kg, below which a body keeps voids (1e21, the end of the compaction taper): `bulk_draw` stores a porosity for every body lighter than this and
none for a heavier one."""
AIRLESS_CLASSES = frozenset("CDS")
"""frozenset: The classes with no air (dead world, small icy body, airless super-Earth): no air column to describe."""

_GIANT_ONLY = ("interior_kind", "core_mass_earth", "dilute_core_mass_earth", "metallization_pressure_gpa",
               "helium_depletion", "has_helium_shell", "envelope_helium_fraction", "envelope_z_fraction")
_MEAN_SWING = ("temp_mean_swing_k",)
_SURFACE = (
    "hydrosphere", "ocean_fraction", "land_fraction", "ocean_depth_km", "ocean_class", "ocean_ph", "water_activity",
    "water_base_pressure_pa", "ice_shell_km", "hp_ice_km", "surface_liquid", "surface_liquid_fraction", "surface_liquid_depth_km", *surface_liquid_world.SOLUTE_FIELDS, "ocean_liquid", "tide_range_spring_m", "tide_range_neap_m",
    "mantle_redox", "mantle_delta_iw", "rock_feo_fraction", "rock_al2o3_fraction", "rock_mg_si_ratio",
    "ozone_column_du", "ozone_loss_flag", "surface_dose_msv_yr", "dose_gcr_msv_yr", "dose_sep_msv_yr",
    "dose_ground_msv_yr", "dose_belt_msv_yr", "dose_helio_mult", "uv_surface_index", "shelter_column_g_cm2", "outside_hours", "outside_hours_day", "outside_hours_limit",
    "phi4", "phi4_pressure", "phi4_temperature", "phi4_air", "phi4_water", "phi4_suitability", "phi4_radiation",
    "tier_pressure", "tier_temperature", "tier_air", "tier_water", "tier_suitability", "tier_radiation",
    "life_microbial", "life_flora", "life_fauna", "life_sentient", "l_solv", "l_chem", "l_ener", "l_rad",
    "c_k0_ppm", "c_th0_ppm", "c_u0_ppm", "al26_ratio", "urey_ratio", "internal_heat_flux_w_m2", "mantle_eta0_pa_s",
    "mantle_activation_kj_mol", *_MEAN_SWING,
    "life_class", "life_class_variant", "life_class_best", "life_class_best_score", "life_animal_grade", "life_complexity",
    "evolutionary_pace", "life_stage", "life_milestone_age_gyr", "life_system_age_gyr",
    "life_age_abiogenesis_gyr", "life_age_photosynthesis_gyr", "life_age_complex_cells_gyr",
    "life_age_multicellularity_gyr", "life_age_technological_civilization_gyr", "phosphorus", "magnetopause_rp",
)
_OCEAN = ("ocean_depth_km", "ocean_class", "ocean_ph", "water_activity", "phosphorus", "water_base_pressure_pa")
_TIDE_COLUMNS = ("ocean_liquid", "tide_range_spring_m", "tide_range_neap_m")
OCEAN_STATES = ("surface ocean", "brine seep", "ice-covered ocean", "hycean")
"""tuple: The hydrospheres that hold a liquid ocean, the ones the ocean columns describe."""
_LIFE_DETAIL = ("life_class", "life_class_variant", "life_class_best", "life_class_best_score", "life_animal_grade", "life_complexity",
                "evolutionary_pace", "life_stage", "life_milestone_age_gyr", "life_system_age_gyr",
                "life_age_abiogenesis_gyr", "life_age_photosynthesis_gyr", "life_age_complex_cells_gyr",
                "life_age_multicellularity_gyr", "life_age_technological_civilization_gyr")
_LIFE_AFTER_OXYGEN = ("life_age_complex_cells_gyr", "life_age_multicellularity_gyr", "life_age_technological_civilization_gyr")
_STAGES_AFTER_OXYGEN = tuple(tuning.LIFE_MILESTONES[tuning.LIFE_MILESTONES.index("complex_cells"):])
_DAYLIGHT = ("day_share_equator", "solar_day_hours", "daylight_hours", "night_hours", "days_per_year")
_RING_COLUMNS = ("ring_kind", "ring_origin", "ring_inner_km", "ring_outer_km", "ring_mass_kg", "ring_optical_depth", "ring_ice_fraction")
_AQUIFER_COLUMNS = ("aquifer_depth_km", "aquifer_thickness_km")
_SURFACE_LIQUID_COLUMNS = ("surface_liquid", "surface_liquid_fraction", "surface_liquid_depth_km")
ICE_SHELL_STATES = ("ice", "ice-covered ocean")
"""tuple: The hydrospheres with an ice shell (`ice_shell_km`)."""
_AIR = ("atm_density", "atm_molar_density", "scale_height_km", "tau_haze", "tau_dust", "tau_cloud", "cloud_cover")

def _locked(row):
    """Whether a body's stored rotation period is its orbital period (the regional temperatures are then three zones)."""
    if row is None or row.get("rotation_period_hours") is None or row.get("period_years") is None:
        return None
    return regional_temperature.is_locked(row["rotation_period_hours"], row["period_years"])


_DISK_COLUMNS = ("moon_disk_lifetime_myr", "moon_disk_alpha_over_f", "moon_disk_f_inf", "moon_disk_ice_line_km")
_BODY_RULES = [
    (_GIANT_ONLY, lambda klass, row: klass in GIANT_CLASSES and (row is None or row.get("body_type") != "t"),
     "only a giant has a deep envelope over a core, so these describe nothing on a rocky or icy world, nor on a gas dwarf that is a rocky body "
     "under a hydrogen air (class T)"),
    (_SURFACE, lambda klass, row: klass not in GIANT_CLASSES,
     "a giant has no surface: no ocean, rock chemistry, ground dose, habitability score, life or surface temperature"),
    (_AIR, lambda klass, row: klass not in AIRLESS_CLASSES and (row is None or (row.get("atmospheric_pressure_pa") or 0) > 0),
     "a dead world, a small icy body and an airless super-Earth, and any body whose air has gone to nothing, have no air column to give a density, scale height or haze"),
    (("energy_flux_w_m2", "flux_variation"), lambda klass, row: klass not in GIANT_CLASSES,
     "a giant is not scored for habitability, so it keeps no starlight or seasonal swing of its own (`insolation` and `habitability_world` leave them none); "
     "its moons, which are scored in its light, store the flux they receive and the swing"),
    (("lid_regime",), lambda klass, row: klass not in GIANT_CLASSES,
     "a giant has no lid: only a rocky or icy body has a crust over a mantle"),
    (("tidal_flux_w_m2",), lambda klass, row: klass not in GIANT_CLASSES,
     "the tide's heat is worked out for the interior of a rocky or icy body; a giant's own heat is the heat it radiates, which has its own row"),
    (("tidal_flux_w_m2",), lambda klass, row: row is None or (row.get("eccentricity") or 0.0) >= NEGLIGIBLE_TIDE_ECCENTRICITY,
     "a body on a circular orbit is not flexed, so the tide gives it no heat: the flow is stored only for an orbit with an eccentricity"),
    (("induced_moment_a_m2",), lambda klass, row: klass not in GIANT_CLASSES,
     "a giant has no induced moment: the model answers for a metal core and an ocean, not for a metallic envelope"),
    (_DAYLIGHT, lambda klass, row: _locked(row) is not True,
     "a tidally locked world keeps one face to its star, so it has no solar day, equatorial daylight or day share to give"),
    (("shelter_column_g_cm2",), lambda klass, row: row is None or (row.get("eq_shelter") or 0) > 0 or (row.get("eq_shelter_rec") or 0) > 0,
     "a shelter column is given only for a body whose radiation calls for a shelter, required or recommended"),
    (("outside_hours",), lambda klass, row: row is None or row.get("outside_hours") is not None,
     "the hours a person can work outside in a year are limited by dose only where the dose is a limit; a suit that allows a full working year outside has none"),
    (("ice_shell_km",), lambda klass, row: row is None or row.get("hydrosphere") in ICE_SHELL_STATES,
     "an ice shell covers a body whose water is ice or an ocean frozen over; a dry body, an open ocean and a vapour world have none"),
    (_RING_COLUMNS, lambda klass, row: row is None or any(row.get(key) is not None for key in _RING_COLUMNS),
     "only a body with a ring system has ring columns; a body without rings stores none"),
    (("ring_age_yr",), lambda klass, row: row is None or row.get("ring_age_yr") is not None,
     "a ring has an age only when a moon that fell or a comet that passed made it; the narrow ring of a small body's spin-orbit resonance has none"),
    (_SURFACE_LIQUID_COLUMNS, lambda klass, row: row is None or any(row.get(key) is not None for key in _SURFACE_LIQUID_COLUMNS),
     "only a body with a liquid at its surface (an ocean, a lava sea, a sea of a volatile) stores one; a dry, frozen or airless body has none"),
    (("surface_liquid_acid_wt", "surface_liquid_acid_layer"), lambda klass, row: row is None or row.get("surface_liquid_acid_wt") is not None,
     "only a body whose air holds sulfur dioxide, an oxidiser and some water makes liquid sulfuric acid, as a cloud layer or a pool; a body without that air has none"),
    (("surface_liquid_ammonia_wt",), lambda klass, row: row is None or row.get("surface_liquid_ammonia_wt") is not None,
     "only a body with a liquid water ocean and ice besides water holds ammonia in it; a dry or frozen body, or one with no other ice, has none"),
    (("surface_liquid_depth_km",), lambda klass, row: row is None or row.get("surface_liquid") not in surface_liquid_world.LAVA_IDS,
     "the depth of a lava sea is not known: the rules that give a molten surface say nothing of how deep it goes"),
    (_AQUIFER_COLUMNS, lambda klass, row: row is None or any(row.get(key) is not None for key in _AQUIFER_COLUMNS),
     "underground water needs a layer between the freezing and the boiling line of water under the crust, which a body either has or lacks; a body without one stores none"),
    (("hp_ice_km",), lambda klass, row: row is None or (row.get("hydrosphere") in ICE_SHELL_STATES
                                                         and (row.get("water_base_pressure_pa") or 0.0) >= HP_ICE_MIN_BASE_PA),
     "high-pressure ice (ice VI and VII) forms only where the water's weight reaches 620 MPa: a dry body, an open ocean, a vapour world and a "
     "body whose ice or ocean is too shallow for that pressure have none"),
    (("porosity",), lambda klass, row: row is None or (row.get("mass_kg") or 0.0) < POROSITY_MAX_MASS_KG,
     "only a body lighter than 1e21 kg is a heap with voids; above it gravity has closed them, so a planet or a large moon stores no porosity"),
    (_OCEAN, lambda klass, row: row is None or row.get("hydrosphere") in OCEAN_STATES,
     "an ocean's depth, chemistry, pH and water activity exist only where there is a liquid ocean"),
    (_TIDE_COLUMNS, lambda klass, row: row is None or row.get("hydrosphere") in ("surface ocean", "ice-covered ocean"),
     "an ocean tide is worked out for a sea with a free surface or one frozen over (zero); a brine seep, a deep hycean sea and a dry, ice or vapour world store none"),
    (_LIFE_DETAIL, lambda klass, row: row is None or any((row.get(key) or 0) > 0 for key in LIFE_SCORES),
     "the life chemistry, stage and milestone ages describe a living world; a lifeless one has none"),
    (_LIFE_AFTER_OXYGEN, lambda klass, row: row is None or row.get("life_stage") in _STAGES_AFTER_OXYGEN,
     "complex cells need free oxygen, which a world has only if its air holds it or its photosynthesis builds it; a world that never does has no age for the stages from complex cells on, and its stage stops before them"),
    (("magnetopause_rp",), lambda klass, row: row is None or (row.get("magnetic_moment_a_m2") or 0) > 0,
     "a world with no magnetic field has no magnetopause"),
    (_DISK_COLUMNS, lambda klass, row: row is None or any(row.get(key) is not None for key in _DISK_COLUMNS),
     "only a giant planet whose moons formed in a disk keeps that disk; a rocky body, a moon and a giant with no disk-born moons store none"),
]

_PLANET_RULES = [
    (("moon_origin",), lambda klass, row: False,
     "only a moon has an origin (disk, impact, captured, ring or binary); a planet orbits its star alone"),
    (("parent_shadow_share",), lambda klass, row: False,
     "only a moon sits in the shadow of its planet; a planet orbits its star alone"),
]

_STAR_RULES = [
    (("flare_alpha", "uv_c_share", "log_lx_lbol", "lifespan_gy"), lambda klass, row: row is None or row.get("yerkes_class") != "VII",
     "a white dwarf has no corona and no flares, so no X-ray share, flare law or UV-C draw, and it does not leave the stage it is in, "
     "so its lifespan has no end"),
    (("decay_age_yr",), lambda klass, row: row is None or row.get("yerkes_class") == "VII",
     "a star that is still burning has not stopped, so it has no age since it did; only a white dwarf, which cooled from the end of its "
     "progenitor's life, has one"),
    (("nuclear_arrival_age_gy",), lambda klass, row: row is None or row.get("nuclear_arrival_age_gy") is not None,
     "only a star still contracting toward the main sequence has an age at which hydrogen burning arrives; every other star is past it, or never had one"),
    (("disk_kind", "disk_inner_au", "disk_outer_au"), lambda klass, row: row is None or row.get("disk_kind") is not None,
     "most stars hold no disk: only a young star (class 0 to II), a main-sequence star with a cold debris belt (about one in seven "
     "to one in four) and the pair inside a planetary nebula do"),
    (("disk_mass_solar",), lambda klass, row: row is None or row.get("disk_kind") not in (None, "debris"),
     "the mass of a debris disk is not measured, only the light of its dust, so a star with a debris belt has no mass for it"),
    (("wide_binary_a_crit_km",), lambda klass, row: row is None or row.get("role") == "member" or row.get("binary_configuration") == "wide",
     "the widest stable orbit around a star is set by a companion that orbits beyond it: the two stars of a wide pair and every further member of a "
     "system have one, but a single star has no companion, and the stars of a close pair are orbited as one body with nothing outside them to set a limit"),
]

_COMET_RULES = [
    (("mean_anomaly_deg", "orbital_period_years", "period_class"), lambda klass, row: row is None or row.get("orbit_type") == "elliptical",
     "a comet on a parabolic path comes by once, so it has no period, period class or place on an orbit to give"),
]

APPLIES = {"planet": _BODY_RULES + _PLANET_RULES, "moon": _BODY_RULES, "star": _STAR_RULES, "comet": _COMET_RULES}
"""dict: Where a sheet column does not apply (UX.208), `{kind: [(keys, test, reason)]}`: each rule names the column keys, a
function `(class) -> bool` that is true for the classes the columns apply to, and the sentence that says why they do not
apply elsewhere. A column no rule names applies to every class of its kind, so one that is empty there is a gap in the
generation pathway, not a rule."""


OPEN_GAPS = {}
"""dict: The columns known to be empty where they apply (UX.208), `{(kind, key): item id}`: each is a gap in the generation
pathway, filed as the item named, and excused until that item is done. The fill audit (`tests/test_sheet_fill.py`) fails
for an empty applicable column that is not here, and the item's owner removes its entry when the column is filled."""


def applies(kind, klass, key, row=None):
    """Whether the column `key` applies to a body of `kind` and class `klass`: true unless a rule of `APPLIES` names it
    and rejects the class, or (for a rule that reads the body) the stored `row` (a mapping of its columns)."""
    for keys, test, _reason in APPLIES.get(kind, ()):
        if key in keys and not test(klass, row):
            return False
    return True


GAS_TABLES = {"planet": ("planet_gases", "planet_id"), "moon": ("moon_gases", "moon_id")}
"""dict: The child table and owner column that keep the air of a planet or a moon, one row per gas (DB.71)."""


def _gas_rows(conn, kind, object_id):
    """The rows of a body's air, one per gas it holds, the most abundant first: each is `(spec, partial pressure Pa)` with a spec
    whose value is read from the gas table."""
    table, owner = GAS_TABLES[kind]
    held = conn.execute(f"SELECT gas, partial_pa FROM {table} WHERE {owner} = ? ORDER BY partial_pa DESC, gas",
                        (object_id,)).fetchall()
    return [(PropertySpec(f"p_{row['gas']}", f"{gases.gas(row['gas']).name.capitalize()} partial pressure", "atmosphere", column="partial_pa"),
             row["partial_pa"]) for row in held]


ZONE_NAMES = {"h": "Hot Zone", "e": "Ecosphere / Habitable Zone", "c": "Cold Zone"}
"""dict: What a body's stored zone letter (`h`, `e`, `c`) is called on a page (UX.217)."""


def _zone_value(letter, by_key, specs, table):
    """The zone row of a planet or moon as an inline list (UX.217): the zone's name and the range of orbit distances it
    covers, from the star's habitable zone edges the body row stores: the Hot Zone lies inside the inner edge, the
    Ecosphere / Habitable Zone between the two, the Cold Zone beyond the outer edge."""
    name = ZONE_NAMES.get(str(letter).lower())
    if name is None:
        return None
    edge = {spec.key: spec for spec in specs}
    inner, outer = by_key.get("habitable_zone_inner_km"), by_key.get("habitable_zone_outer_km")
    if inner is None or outer is None:
        return [name]
    near, far = (_value(edge["habitable_zone_inner_km"], table, inner)[0], _value(edge["habitable_zone_outer_km"], table, outer)[0])
    letter = str(letter).lower()
    if letter == "h":
        return [f"{name}, closer than ", near]
    if letter == "c":
        return [f"{name}, farther than ", far]
    return [f"{name}, ", near, " to ", far, " from its star"]


_yr_to_h = unit_registry.converter("yr", "h")
_kg_to_msun = unit_registry.converter("kg", "Msun")

VARIABILITY_KEY = "variability"
"""str: The key of the star's derived row for the kind of variable it is (GEN.944, `physics.variability`); it is worked out from the stored
state, class, temperature and mass when the sheet is made and is not a column."""

CARBON_KEY = "carbon_star_kind"
"""str: The key of the star's derived row for a carbon star or an S star (GEN.1026, `physics.agb_carbon`); it is worked out from the stored state,
lifespan (which gives the initial mass), mass, [Fe/H] and age when the sheet is made and is not a column."""

PECULIAR_KEY = "chemically_peculiar_kind"
"""str: The key of the star's derived row for an Ap or Am star (GEN.1027, `physics.chemically_peculiar`); it is worked out from the stored class, mass,
radius, temperature and age when the sheet is made and is not a column."""

PECULIAR_FIELD_KEY = "chemically_peculiar_polar_field"
"""str: The key of an Ap star's derived row for its polar magnetic field (GEN.1027); an Am star has none."""

CHEMISTRY_ROWS = {
    "chemistry_branch": ("Disc population", "words"),
    "chemistry_alpha_fe": ("Alpha enhancement [alpha/Fe]", "dex"),
    "chemistry_metal_fraction": ("Metal mass fraction", "mass fraction"),
    "chemistry_dredge_up_helium": ("Helium added by the first dredge-up", "mass fraction"),
}
"""dict: The star's derived chemistry rows (UX.254, `physics.stellar_composition.star_makeup`) -> `(label, the unit the properties help page names)`: the
disc population, [alpha/Fe], the metal mass fraction and, for a red giant, the helium the first dredge-up added. They are worked out from the stored
class, mass, [Fe/H] and age when the sheet is made and are not columns."""

CHEMISTRY_KEYS = tuple(CHEMISTRY_ROWS)
"""tuple: The keys of `CHEMISTRY_ROWS`."""

_BRANCH_WORDS = {stellar_abundances.BRANCH_THIN: "Thin disc", stellar_abundances.BRANCH_HIGH_ALPHA: "Thick disc or halo (alpha-rich)"}


def _star_distance_ly(conn, system_id):
    """The distance from the galactic centre, light-years, the star's abundance draws are keyed on (its sector's, as the generator was given it),
    or the galaxy's default for a star in no placed sector."""
    found = None if system_id is None else conn.execute(
        "SELECT sectors.galactic_radius_pc FROM star_systems JOIN sectors ON sectors.id = star_systems.sector_id WHERE star_systems.id = ?",
        (system_id,)).fetchone()
    return scales.default_galactocentric_ly() if found is None or found["galactic_radius_pc"] is None else units.pc_to_ly(found["galactic_radius_pc"])


def _chemistry_rows(conn, row):
    """
    The chemistry rows of a star (UX.254): its disc population, [alpha/Fe], the metal mass fraction it formed with and, for a red giant, the helium the
    first dredge-up added. They come from `stellar_composition.star_makeup` on the stored class, mass, [Fe/H] and age, with the draws keyed as the
    generator keys its planets' host ratios (mass and the sector's distance from the centre, not the age), so nothing is stored. A star with no [Fe/H],
    one the composition does not cover or one whose downloaded data are missing has none.
    """
    if row["metallicity_fe_h"] is None or row["mass_kg"] is None or row["age_gy"] is None:
        return []
    mass_sol = _kg_to_msun(row["mass_kg"])
    key = stellar_composition.abundance_key(mass=row["mass_kg"], distance_ly=_star_distance_ly(conn, row["star_system_id"]))
    try:
        kind = stellar_composition.kind_of(row["yerkes_class"], mass_sol, row["state"], None)
        if kind not in stellar_composition.WITH_BIRTH_COMPOSITION:
            return []
        makeup = stellar_composition.star_makeup(kind, row["metallicity_fe_h"], row["age_gy"], mass_sol, key,
                                                 wolf_rayet_subtype=row["star_type"], teff_k=row["temperature_k"])
    except (stellar_composition.StarMakeupError, fetch.DataMissing):
        return []
    alpha = stellar_abundances.alpha_enhancement(makeup.birth)
    items = [("chemistry_branch", [_BRANCH_WORDS[makeup.branch]], None),
             ("chemistry_alpha_fe", [_plain_number(alpha, "dex")], alpha),
             ("chemistry_metal_fraction", [model.qty("fraction", makeup.birth.z)], makeup.birth.z)]
    if "dredge_up" in makeup.detail:
        rise = makeup.detail["dredge_up"].delta_y
        items.append(("chemistry_dredge_up_helium", [model.qty("fraction", rise)], rise))
    return [{"key": key, "label": CHEMISTRY_ROWS[key][0], "value": value, "text": inline_text(value), "raw": raw, "stored": True, "absent": False,
             "source": source_of("star", key), "show_absent": False} for key, value, raw in items]


ECLIPSE_KEY = "longest_eclipse_hours"
"""str: The key of the moon's derived row for the longest eclipse in its planet's shadow (UX.205); it is worked out from the stored
orbit and the planet's radius when the sheet is made and is not a column."""


def _longest_eclipse(conn, kind, object_id):
    """The longest time a moon is in its planet's shadow in one pass, hours (`daylight.longest_eclipse_hours`), or `None` for a moon that
    orbits no planet, is never in the shadow or lacks the orbit."""
    moon = conn.execute("SELECT parent_type, parent_id, distance_km, period_years, parent_shadow_share FROM moons WHERE id = ?",
                        (object_id,)).fetchone()
    if moon is None or moon["parent_type"] != "planet" or not moon["parent_shadow_share"]:
        return None
    parent = conn.execute("SELECT radius_km FROM planets WHERE id = ?", (moon["parent_id"],)).fetchone()
    if parent is None or not moon["distance_km"] or not moon["period_years"]:
        return None
    hours = daylight.longest_eclipse_hours(moon["distance_km"], parent["radius_km"], _yr_to_h(moon["period_years"]))
    return hours or None


BAND_ROW_LABELS = (("mean_k", "Annual mean temperature"), ("swing_k", "Day and night swing (hottest to coldest)"), ("day_hours", "Day"),
                   ("night_hours", "Night"), ("longest_day_days", "Longest day of the year"),
                   ("longest_night_days", "Longest night of the year"), ("state", "State"))
"""tuple: The rows of the table of the bands (GEN.1108), `regional_temperature.BAND_FIELDS` in order, each with its label."""

BAND_LABELS = {"north_polar": "North polar", "north_middle": "North middle", "north_tropic": "North tropics", "equator": "Equator",
               "south_tropic": "South tropics", "south_middle": "South middle", "south_polar": "South polar",
               "day_centre": "Day centre", "day_middle": "Day side", "day_edge": "Day edge", "terminator": "Terminator",
               "night_edge": "Night edge", "night_middle": "Night side", "night_centre": "Night centre"}
"""dict: What the table of the bands calls each band or ring (`latitude_bands.NAMES`, `RING_NAMES`)."""


def _band_cell(field, value):
    """One cell of the table of the bands as an inline list: the quantity, or the state word; a band with no such figure (a locked world has no hours) is empty."""
    if value is None:
        return []
    if field == "mean_k":
        return [model.qty("temperature", value, "K", companions=False)]
    if field == "swing_k":
        return ["no swing"] if not value else ["\u00b1", model.qty("temperature difference", value / 2.0, "dK", companions=False)]
    if field in ("day_hours", "night_hours"):
        return [model.qty("time", value, "h")]
    if field in ("longest_day_days", "longest_night_days"):
        return [model.qty("time", value, "d")]
    return [str(value).capitalize()]


def bands_table(row):
    """
    The seven bands of a body side by side (GEN.1108): `{"locked", "bands": [{"key", "label"}], "rows": [{"key", "label", "cells"}]}`, one row for each of
    `BAND_ROW_LABELS` that any band holds, one cell (an inline list, empty where the band has none) for each band from the north pole to the south pole (the day centre
    to the night centre for a world that turns once an orbit). `None` for a body that stores no bands.
    """
    if row is None or row.get(regional_temperature.column(1, "mean_k")) is None:
        return None
    locked = row.get(regional_temperature.column(1, "state")) in latitude_bands.RING_STATES
    names = latitude_bands.RING_NAMES if locked else latitude_bands.NAMES
    bands = [{"key": name, "label": BAND_LABELS[name]} for name in names]
    rows = []
    for field, label in BAND_ROW_LABELS:
        cells = [_band_cell(field, row.get(regional_temperature.column(number, field))) for number in range(1, regional_temperature.BAND_COUNT + 1)]
        if any(cells):
            rows.append({"key": field, "label": label, "cells": cells})
    return {"locked": locked, "bands": bands, "rows": rows}


def sheet(conn, kind, object_id):
    """
    The property sheet of one body.

    Args:
        conn: An open, read-only connection.
        kind (str): A key of `SPECS`: `star`, `planet`, `moon`, `belt`, `comet` or a phenomenon kind.
        object_id (int): The row id.

    Returns:
        dict: `kind`, `groups`: the non-empty groups in `GROUPS` order, each `{"key", "label", "rows"}`; a row is
            `{"key", "label", "value" (the inline list), "text" (it as plain text), "raw" (the stored number, `None` for text), "stored" (whether the column held a
            value), "absent" (none or 0), "source" (a one-line source for a physics row, UX.187, else `None`), "show_absent" (a page lists it even then)}`.

    Raises:
        ValueError: For a kind with no sheet or a missing row.
    """
    specs = specs_for(kind)
    table = TABLES[kind]
    wanted = {spec.stored_column() for spec in specs}
    if table in ("planets", "moons"):
        wanted.update(regional_temperature.BAND_COLUMNS)    # GEN.1108: the table of the bands
        wanted.add("temp_mean_swing_k")                      # the range on the mean temperature row
    if table == "stars":
        wanted.update(("lifespan_gy", "metallicity_fe_h", "age_gy", "star_system_id", "star_type"))    # GEN.1026: the carbon or S star row, UX.254: the chemistry rows
    columns = ", ".join(sorted(wanted))
    row = conn.execute(f"SELECT {columns} FROM {table} WHERE id = ?", (object_id,)).fetchone()
    if row is None:
        raise ValueError(f"no {table} row with id {object_id}")
    by_group = {key: [] for key in GROUPS}
    for spec in specs:
        stored = row[spec.stored_column()]
        value = _range_value(spec, table, stored, row)
        absent = stored is None or stored == 0 or stored == ""
        if spec.key == "zone" and kind in ("planet", "moon") and stored:
            value = _zone_value(stored, {s.key: row[s.stored_column()] for s in specs}, specs, table) or value
        by_group[spec.group].append({"key": spec.key, "label": spec.label, "value": value,
                                     "text": inline_text(value), "raw": stored if isinstance(stored, (int, float)) else None,
                                     "stored": stored is not None, "absent": absent, "source": source_of(kind, spec.key),
                                     "show_absent": spec.group in SHOW_GROUPS or spec.key in SHOW_WHEN_ABSENT.get(kind, ())})
    for spec, pascals in _gas_rows(conn, kind, object_id) if kind in GAS_TABLES else ():
        value = _value(spec, GAS_TABLES[kind][0], pascals)
        by_group[spec.group].append({"key": spec.key, "label": spec.label, "value": value, "text": inline_text(value), "raw": pascals, "stored": True,
                                     "absent": not pascals, "source": None, "show_absent": spec.key in SHOW_WHEN_ABSENT.get(kind, ())})
    eclipse = _longest_eclipse(conn, kind, object_id) if kind == "moon" else None
    if eclipse:
        value = [model.qty("time", eclipse, "h")]
        by_group["rotation"].append({"key": ECLIPSE_KEY, "label": "Longest eclipse in the planet's shadow", "value": value,
                                     "text": inline_text(value), "raw": eclipse, "stored": True, "absent": False, "source": None, "show_absent": False})
    variable = variability.kind_of(row["state"], row["yerkes_class"], row["temperature_k"],
                                   None if row["mass_kg"] is None else _kg_to_msun(row["mass_kg"])) if kind == "star" else None
    if variable:
        by_group["identity"].append({"key": VARIABILITY_KEY, "label": "Variability", "value": [variable], "text": variable, "raw": None,
                                     "stored": True, "absent": False, "source": None, "show_absent": False})
    carbon = None
    if kind == "star" and row["state"] == "asymptotic-giant" and row["metallicity_fe_h"] is not None:
        initial = agb_carbon.initial_mass_from_lifespan_sol(row["lifespan_gy"], row["metallicity_fe_h"])
        carbon = agb_carbon.kind_of(row["state"], initial, _kg_to_msun(row["mass_kg"]), row["metallicity_fe_h"], row["age_gy"])
    if carbon:
        by_group["identity"].append({"key": CARBON_KEY, "label": "Carbon or S star", "value": [carbon], "text": carbon, "raw": None,
                                     "stored": True, "absent": False, "source": None, "show_absent": False})
    peculiar = chemically_peculiar.of_star(row["yerkes_class"], row["mass_kg"], row["radius_km"], row["temperature_k"], row["age_gy"]) if kind == "star" else None
    if peculiar:
        by_group["identity"].append({"key": PECULIAR_KEY, "label": "Chemically peculiar star", "value": [peculiar.kind], "text": peculiar.kind, "raw": None,
                                     "stored": True, "absent": False, "source": None, "show_absent": False})
        if peculiar.field_gauss is not None:
            field = [model.qty("magnetic field", peculiar.field_gauss, "G")]
            by_group["physical"].append({"key": PECULIAR_FIELD_KEY, "label": "Polar magnetic field", "value": field, "text": inline_text(field),
                                         "raw": peculiar.field_gauss, "stored": True, "absent": False, "source": None, "show_absent": False})
    if kind == "star":
        by_group["composition"].extend(_chemistry_rows(conn, row))
    groups = [{"key": key, "label": GROUPS[key], "rows": rows} for key, rows in by_group.items() if rows]
    stored = {spec.key: row[spec.stored_column()] for spec in specs}
    result = {"kind": kind, "groups": groups, "missing": missing(kind, stored)}
    if table in ("planets", "moons"):
        result["bands"] = bands_table(row)
    return result
