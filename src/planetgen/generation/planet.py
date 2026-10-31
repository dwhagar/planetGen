# planetgen/generation/planet.py

"""
Planet Generation
==================

This module defines the `Planet` class, the representation of a single
planet or moon within a star system. `Planet` itself holds the object's
state and its presentation (`to_paragraph_list`); the logic that determines
what a planet physically *is* lives in `planetPhysics`, and the logic that
determines whether/how life arises on it lives in `life` — both
operate on `Planet` instances passed in as arguments, rather than as
methods on the class, to keep generation concerns separate from the object
itself.

A `Planet` is constructed *without* life data (see `planetPhysics`'s module
docstring); `StarSystem` applies `life.apply_life_data` to every
planet and moon once the whole system has been generated.

The `AsteroidBelt` class, a simpler representation of an asteroid belt, is
defined separately in `asteroidData`.
"""

import math

from planetgen import tuning
from planetgen.generation import habitable_request, prevalence
from planetgen.generation.config import SystemConfig
from planetgen.physics import aerosols, atmosphere, body_colour, bulk_composition, climate_draws, constants, daylight, giant_interior, habitability_world, hydrosphere, interior_heat, life_classes, magnetism, moon_composition, ocean_tides, planets as planetPhysics, radiation, radiogenic, rings, surface_composition
from planetgen.physics.position import axis_property, velocity_axis_property
from planetgen.util import draw
from planetgen.generation.planetary_object import PlanetaryObject
from planetgen.util.format import (
    format_number, format_pressure_pa, format_temperature_k,
)


DEFER_MOONS = "defer"
"""str: A `moon_count` that builds the planet bare; the caller adds its moons later with `Planet.add_moons` (a request that draws
the planet again until it is the one wanted does not generate moons for the ones it throws away)."""


class Planet(PlanetaryObject):
    """
    A class representing a single planet or moon and all of its properties.

    This class encapsulates the physical characteristics, orbital data, and
    atmospheric conditions of a celestial body. It can represent both
    planets and moons, and includes a list of any moons it has.

    The generation process is flexible, allowing for fully random creation or
    generation based on specified parameters like radius, mass, or class;
    this generation work is delegated to `planetPhysics.generate_planet_properties`
    and related functions in that module, called from `__init__`. Life
    chemistry and evolutionary data are NOT set during construction — see
    `life.apply_life_data`, applied later by `StarSystem`.

    A planet reads whichever of the host star's own properties it needs
    (luminosity, mass, spectral type) directly off the `star` object it
    holds a reference to, rather than duplicating them into separate
    `Planet` attributes -- `star` is the single source of truth for
    everything star-related.

    `body_type` (values `'t'`/`'g'`) is this class's half of a lightweight
    discriminator shared with `AsteroidBelt` (whose `body_type` is always
    `'a'`): `StarSystem` holds `Planet` and `AsteroidBelt` instances
    together in one `planets` list and switches on `.body_type` to tell
    them apart. `Star.type`/`BinaryStarProxy.type` is an unrelated,
    differently-shaped attribute (a full descriptive string like
    `"G2V Yellow Main Sequence Star"`, not a single-character code) that
    happens to share the same name on a different class.

    Attributes:
        system_config (SystemConfig): The shared SystemConfig object for the system.
        is_moon (bool): True if the object is a moon, False otherwise.
        moons (list): A list of `Planet` objects representing the moons of this planet.
        zone (str): The orbital zone of the planet ('h' for hot, 'e' for ecosphere, 'c' for cold).
        description (str): A text description of the planet's class.
        atm_molar_density (float): The molar density of the atmosphere.
        gravity (float): The surface gravity of the planet in g's.
        atm_density (float): The density of the atmosphere.
        surface_temperature (float): The surface temperature of the planet in Kelvin.
        density (float): The overall density of the planet in g/cm³.
        atmospheric_pressure (float): The atmospheric pressure at the surface in Pascals.
        mass (float): The mass of the planet in kilograms.
        atmosphere (str): The atmospheric composition.
        composition (str): The chemical composition of the planet, as a
                           descriptive string (contrast `AsteroidBelt.composition`,
                           a list of (component, concentration) tuples -- same
                           name, unrelated shape, on the two classes).
        radius (float): The radius of the planet in kilometers.
        planet_class (str): The classification of the planet (e.g., 'M', 'N').
        distance (float): The distance from the star in AU.
        body_type (str): 't' for terrestrial, 'g' for gas giant -- see the
                        class docstring for how this relates to `AsteroidBelt.body_type`
                        and the unrelated `Star.type`.
        scale_height (float): The atmospheric scale height in kilometers.
        hill_radius (float): The Hill radius of the planet in kilometers.
        min_orbit_distance (float): This body's own individual Hill
                                   radius x 5, in AU -- the outer bound
                                   for a *satellite*'s orbit around it
                                   (`planetPhysics.generate_moons`), and
                                   `StarSystem`'s own initial spacing
                                   estimate/belt-adjacent-planet spacing.
                                   Two real planets' minimum separation
                                   from each other uses their *mutual*
                                   Hill radius instead -- see
                                   `StarSystem._mutual_min_separation_au`
                                   -- not this single-body value.
        name (str): The generated name of the planet or moon.
        life_class, life_class_variant (int, str or None): The life class that lives on the body, its number and its pigment or
                             habitat letter (`physics.life_classes`; set by `life.apply_life_data`, None until then or with none).
        life_class_best, life_class_best_score (str, int or None): The best-supported class and its score, 0 to 9.
        life_class_scores (dict or None): `{class id: score}` for every class the body supports (a score of 1 or more).
        life_animal_grade (int or None): The grade of its animals, 0 none to 3 large and active (`physics.life_complexity`).
        life_complexity (float or None): The complexity index of its life class (`physics.life_complexity`).
        reflection_spectrum_visible (list): Visible-light reflection descriptors
                                            for any life present (set by
                                            `life.apply_life_data`, None until then).
        reflection_spectrum_non_visible (list): Non-visible-spectrum reflection
                                                descriptors for any life present
                                                (set by `life.apply_life_data`,
                                                None until then).
        habitable_zone (tuple): The inner and outer bounds of the host star's
                                habitable zone in AU (same value as `Star.habitable_zone`,
                                threaded down at construction time).
        star (Star): The Star object this planet orbits.
        volume (float): The planet's volume in km^3.
        period (float): The planet's orbital period in years.
        orbital_inclination_deg (float): The tilt of this body's orbital
                                         plane, in degrees, relative to its
                                         primary's reference plane.
        orbital_ascending_node_deg (float): The longitude of the ascending
                                            node, in degrees -- together
                                            with `orbital_inclination_deg`,
                                            fully orients the orbital
                                            plane in 3D.
        eccentricity (float): The eccentricity of the orbit it was generated on (GEN.727); `distance` is its semi-major axis.
                              The light it gets goes by the orbit's mean flux, `1 / sqrt(1 - e^2)` times that of a circle.
        velocity_x_kms, velocity_y_kms, velocity_z_kms (float): This body's
                                   velocity, km/s, relative to its orbital
                                   anchor (the same one `position_x/y/z`
                                   is measured from): the tangent of its
                                   orbit where it is,
                                   `orbital_speed_kms` long. Set with the
                                   position (`planetPhysics.
                                   update_orbital_position`). With the
                                   position, the mass and the epoch it is
                                   the body's state; the nightly update
                                   moves it (`physics.body_state`, GEN.289),
                                   and there is no stored phase.
        position_x, position_y, position_z (float): This body's current
                                   Cartesian position, in AU, relative to
                                   its orbital anchor -- the star (or, for
                                   a binary system, the `BinaryStarProxy`
                                   standing in for the system's combined
                                   center) for a planet, the parent planet
                                   for a moon. Made from `distance` and
                                   the orbital elements above (see
                                   `planetPhysics.update_orbital_position`);
                                   moved by the nightly update.
        orbital_speed_kms (float): This body's orbital speed where it was
                                   put, in km/s (the length of its
                                   velocity then; constant for a circle,
                                   and changed if `distance`/`period` are
                                   (see `StarSystem.validate_system`)).
        min_update_interval_years (float): The shortest `elapsed_years`
                                   worth calling `store.advance_orbital_phases`
                                   for (see
                                   `orbits.minimum_update_interval_years`).
                                   Not a narrative stat -- a guard value for
                                   the update path, derived purely from
                                   `period`.
        rotation_period_hours (float): This body's axial rotation period
                                       ("day length"), in hours.
        spin_axis_x/y/z (float): Unit spin axis in the orbit's frame (GEN.104,
                                 `physics.spin`); the spin vector is this
                                 times 2 pi / the rotation period.
        axial_tilt_deg (float): The angle, degrees, between that axis and the
                                orbit normal.
        mantle_redox (str or None): 'reduced', 'intermediate' or 'oxidized'
                                    (GEN.85, `physics.atmosphere`); `None`
                                    for a gas giant.
        mantle_delta_iw (float or None): The upper mantle's oxygen fugacity,
                                         log units above the iron-wustite buffer.
        gas_partials_pa (dict): The air: gas-table formula to partial pressure, Pa, only the gases above zero (stored in
                                the child tables `planet_gases` and `moon_gases`, DB.71; `physics.atmosphere`).
        life_stage (str or None): The latest `evolution.MILESTONE_KEYS` milestone the timeline reached (set by
                                  `life.apply_life_data`; `None` before any, or with no timeline).
        evolutionary_pace (str or None): How fast life evolves here, a `tuning.EVOLUTIONARY_TIMELINES` key (`fast`,
                                         `normal`, `slow`), drawn once by `life.apply_life_data`; `None` when the
                                         body has no life chemistry. Every sentence about the pace reads it.
        life_milestone_age_gyr, life_system_age_gyr (float or None): The milestone's age and the system's, Gy.
        life_age_abiogenesis_gyr ... life_age_technological_civilization_gyr (float or None): The age at which the
                                  world's pace reaches each milestone (`tuning.LIFE_SCHEDULE_FIELDS`, Gy; set by
                                  `life.apply_life_data` with the pace, `None` without one).
        flavor_key (str or None): The flavour line `life.decide_flavor_text` drew, as `category:index` into the
                                  `tuning` flavour lists (`flavor_text` reads it).
        flavor_text_count (int): The number of flavor texts added to this planet.
    """

    SERIALIZABLE_FIELDS = [
        "is_moon", "zone", "description", "atm_molar_density", "gravity",
        "atm_density", "surface_temperature", "density", "atmospheric_pressure",
        "mass", "atmosphere", "composition", "radius", "planet_class", "distance",
        "body_type", "scale_height", "hill_radius", "min_orbit_distance", "name",
        *life_classes.BODY_FIELDS, life_classes.SCORES_FIELD, "reflection_spectrum_visible",
        "reflection_spectrum_non_visible", "life_stage", "evolutionary_pace",
        "life_milestone_age_gyr", "life_system_age_gyr", *tuning.LIFE_SCHEDULE_FIELDS, "flavor_key",
        "flavor_text_count", "habitable_zone", "volume", "period",
        "orbital_inclination_deg", "orbital_ascending_node_deg", "eccentricity",
        "position_x", "position_y", "position_z", "orbital_speed_kms",
        "velocity_x_kms", "velocity_y_kms", "velocity_z_kms",
        "min_update_interval_years",
        "rotation_period_hours", "spin_axis_x", "spin_axis_y", "spin_axis_z", "axial_tilt_deg",
        *atmosphere.ATMOSPHERE_FIELDS,
        atmosphere.GAS_FIELD,
        *climate_draws.CLIMATE_FIELDS,
        *magnetism.BODY_FIELDS,
        *hydrosphere.HYDROSPHERE_FIELDS,
        *radiogenic.BODY_FIELDS,
        *interior_heat.BODY_FIELDS,
        *interior_heat.OUTPUT_FIELDS,
        *interior_heat.REGIME_FIELDS,
        *daylight.FIELDS,
        *ocean_tides.FIELDS,
        *bulk_composition.NEW_FIELDS,
        *giant_interior.BODY_FIELDS,
        giant_interior.LAYERS_ATTRIBUTE,
        *surface_composition.COLUMN_FIELDS,
        *moon_composition.BODY_FIELDS,
        *rings.BODY_FIELDS,
        *radiation.BODY_FIELDS,
        *habitability_world.BODY_FIELDS,
        *aerosols.BODY_FIELDS,
        *body_colour.BODY_FIELDS,
        "reflex_offset_x", "reflex_offset_y", "reflex_offset_z",
        "born_beyond_snow_line", "host_eccentricity",
    ]
    """
    Every attribute set by `__init__` (directly, or by the `planetPhysics`
    helpers it calls), excluding `star`/`system_config` (shared
    back-references, threaded into `from_dict` rather than serialized
    redundantly on every planet/moon) and `moons` (a nested list of further
    `Planet` dicts, handled separately in `to_dict`/`from_dict`). Covers
    moons too -- a moon is just a `Planet` with `is_moon=True`, and nesting
    is exactly 2 levels (a moon's own `moons` is always empty).
    """

    position_x = axis_property(0)
    position_y = axis_property(1)
    position_z = axis_property(2)
    velocity_x_kms = velocity_axis_property(0)
    velocity_y_kms = velocity_axis_property(1)
    velocity_z_kms = velocity_axis_property(2)
    """
    float: This body's offset from its orbital anchor, AU (the star, or a
    moon's parent planet): the "system" frame of `spatial`, its
    `SpatialPosition3D` (GEN.74), which `set_position_au` moves.
    """

    eccentricity = 0.0
    host_eccentricity = 0.0
    """float: The eccentricity of the orbit the body was generated on (GEN.727, `physics.orbit_elements`), beside `distance` and `period`; the
    orbit it is on now is whatever its state vector says (`orbit_from_vector`). A moon's `host_eccentricity` is its parent's, which the light
    it gets goes by. A circle by default: a body made without them is on one."""

    moon_origin = None
    moon_disk_system = None
    """str or None, DiskSystem or None: A moon's origin (`moon_composition.ORIGINS`, stored: DB.70) and its parent's circumplanetary disk while
    the moon is being generated (not stored itself; the giant keeps its numbers in `moon_composition.DISK_FIELDS`)."""

    reflex_offset_x = reflex_offset_y = reflex_offset_z = 0.0
    """
    float: Class-level defaults (schema v20) so a `Planet` reconstructed
    via `from_dict`'s `object.__new__`+`fields_from_dict` path still reads
    back `0.0` for a system saved before this field existed, rather than
    raising `AttributeError` -- same reasoning as `Star.a_crit_au`'s own
    class-level default. See this class's `__init__` for the full
    explanation of what this field means.
    """

    def __init__(self, system_config: SystemConfig, star, habitable_zone, distance,
                 radius=None, mass=None, zone_override=None, distance_override=None,
                 is_moon=False, moon_count=None, primary_mass_kg=None, moon_origin=None, moon_disk_system=None,
                 planned_objects=None, host_eccentricity=0.0, formation_distance_au=None):
        """
        Initializes a Planet object with its physical and orbital properties.

        Life chemistry and evolutionary data are intentionally not set here;
        call `life.apply_life_data(planet)` separately once the whole
        system's planets and moons have been generated (this is done
        automatically by `StarSystem.__init__`).

        Args:
            system_config (SystemConfig): The shared SystemConfig object for the system.
            star (Star): The Star object this planet orbits. Its luminosity,
                        mass, and spectral type are read directly off this
                        reference wherever needed, rather than being passed
                        in and duplicated as separate `Planet` attributes.
            habitable_zone (tuple): A tuple containing the inner and outer bounds
                              of the habitable zone in AU.
            distance (float): The planet's distance from its star in AU.
            radius (float, optional): The radius of the planet in kilometers.
            mass (float, optional): The mass of the planet in kilograms.
            zone_override (str, optional): A single character ('h', 'c', 'e') to
                                           override the planet's calculated zone.
            distance_override (float, optional): An override for the planet's
                                                 distance, used in specific calculations.
            is_moon (bool, optional): Flag indicating if the object is a moon.
            moon_count (int, optional): An exact number of moons to generate
                                        for this planet (0 for none). Ignored
                                        if `is_moon` is True. If None, moon
                                        generation falls back to
                                        `system_config.MOONS`/random chance.
                                        `DEFER_MOONS` builds none yet (`add_moons`).
            primary_mass_kg (float, optional): The mass (in kg) of the body
                                               this one actually orbits, for
                                               Kepler's third law. Defaults to
                                               `star.mass` -- correct for an
                                               ordinary planet, but a moon
                                               orbits its parent planet, not
                                               the star (`generate_moons`
                                               passes the parent planet's own
                                               `mass` here).
            formation_distance_au (float, optional): Where the planet formed, AU, when that is not where it orbits (a pulsar's disc planet
                                         formed in the volatile-poor debris of a disrupted star, so inside any snow line); it sets
                                         whether the planet is ice-rich or dry (`bulk_draw.origin_of`). Defaults to the distance.
            moon_origin (str, optional): For a moon, how it came to orbit its parent (`moon_composition.ORIGINS`), drawn by
                                         `generate_moons`: the orbit's inclination follows from it (GEN.714). A moon built
                                         without one keeps the old regular orbit.
            planned_objects (int, optional): How many planets and belts the system planned for this star: one means a single planet,
                                             whose eccentricity follows the broader law (GEN.727); `None` is taken as several.
            host_eccentricity (float, optional): For a moon, its parent planet's eccentricity (the light a moon gets goes by its
                                                 parent's orbit).
            moon_disk_system (DiskSystem, optional): For a disk moon, its giant parent's circumplanetary disk
                                                     (`moon_composition.draw_disk_system`), shared by that parent's moons.
        """
        self.system_config = system_config # Store SystemConfig
        self.is_moon = is_moon
        self.moon_origin = moon_origin
        self.moon_disk_system = moon_disk_system
        self.moons = []
        # This planet's own reflex-offset "wobble" (schema v20) from its
        # moons' combined pull -- 0.0 by default (no moons yet), recomputed
        # by `planetPhysics.generate_moons` once `self.moons` is populated
        # below. Always 0.0 for a moon itself (moons never generate their
        # own moons). See `Star.reflex_offset_x`'s docstring for the same
        # concept on the star side.
        self.reflex_offset_x = self.reflex_offset_y = self.reflex_offset_z = 0.0
        self.zone = None
        # GEN.458: whether it formed beyond its star's snow line; `StarSystem` sets it once the orbits are final
        # (a moon takes its planet's), and it stays None for a planet whose host has no snow line to compare.
        self.born_beyond_snow_line = None
        self.description = None
        self.atm_molar_density = None
        self.gravity = None
        self.atm_density = None
        self.surface_temperature = None
        self.density = None
        self.atmospheric_pressure = None
        self.mass = mass
        self.atmosphere = None
        self.composition = None
        self.radius = radius
        self.planet_class = None    # derived from the finished body (`planetPhysics.apply_class`)
        self.distance = distance
        # Where the body formed, in AU (a moon's is its planet's), which decides whether its water is ice (`bulk_draw.origin_of`).
        self.formation_distance_au = (formation_distance_au if formation_distance_au is not None
                                      else distance_override if distance_override is not None else distance)
        self.body_type = None
        self.scale_height = None
        self.hill_radius = None
        self.min_orbit_distance = None
        self.name = None # Set by the owning StarSystem (names.bodies.name_bodies) once orbits are final.
        for field in life_classes.BODY_FIELDS:
            setattr(self, field, None)
        setattr(self, life_classes.SCORES_FIELD, None)
        self.reflection_spectrum_visible = None
        self.reflection_spectrum_non_visible = None
        self.life_stage = None # Set later by life.apply_life_data, with the next three
        self.evolutionary_pace = None
        self.life_milestone_age_gyr = None
        self.life_system_age_gyr = None
        for field in tuning.LIFE_SCHEDULE_FIELDS:
            setattr(self, field, None) # The pace's milestone ages, set with the pace
        self.flavor_key = None # Set later by life.decide_flavor_text
        self.flavor_text_count = 0 # Initialize flavor text count for this planet
        self.orbital_inclination_deg = None
        self.orbital_ascending_node_deg = None
        self.eccentricity = 0.0
        self.host_eccentricity = host_eccentricity
        self.position_x = None
        self.position_y = None
        self.position_z = None
        self.velocity_x_kms = None
        self.velocity_y_kms = None
        self.velocity_z_kms = None
        self.orbital_speed_kms = None
        self.min_update_interval_years = None
        self.rotation_period_hours = None
        self.spin_axis_x = None
        self.spin_axis_y = None
        self.spin_axis_z = None
        self.axial_tilt_deg = None
        for field in atmosphere.ATMOSPHERE_FIELDS + magnetism.BODY_FIELDS + hydrosphere.HYDROSPHERE_FIELDS + radiogenic.BODY_FIELDS + interior_heat.BODY_FIELDS + interior_heat.OUTPUT_FIELDS + interior_heat.REGIME_FIELDS + daylight.FIELDS + ocean_tides.FIELDS + bulk_composition.NEW_FIELDS + rings.BODY_FIELDS + surface_composition.COLUMN_FIELDS + moon_composition.DISK_FIELDS + giant_interior.BODY_FIELDS + (giant_interior.LAYERS_ATTRIBUTE,) + radiation.BODY_FIELDS + habitability_world.BODY_FIELDS + aerosols.BODY_FIELDS + body_colour.BODY_FIELDS:
            setattr(self, field, None)
        setattr(self, atmosphere.GAS_FIELD, {})

        # From the star, should not be changed.
        self.habitable_zone = habitable_zone
        self.star = star # Store the Star object

        # The primary is the star for an ordinary planet, but a moon orbits
        # its parent planet -- primary_mass_kg (defaulted to star.mass)
        # carries that distinction in from generate_moons below. Its Hill
        # sphere (set inside generate_planet_properties) uses it (GEN.138).
        primary_mass_kg = primary_mass_kg if primary_mass_kg is not None else star.mass
        self.primary_mass_kg = primary_mass_kg
        # Calculate physical/orbital properties (no life data yet). A system that bars habitable worlds from its ecosphere
        # (`HABITABLE_WORLD` False) draws the body again where the finished body would have the conditions a habitable request asks for (GEN.916).
        given_mass, given_radius = self.mass, self.radius
        for _attempt in range(tuning.BARRED_CLASS_MAX_ATTEMPTS):
            self.mass, self.radius = given_mass, given_radius
            planetPhysics.generate_planet_properties(self, zone_override)  # sets self.volume (km^3) and self.mass
            self.period = planetPhysics.calculate_orbital_period_years(self.distance, primary_mass_kg)
            planetPhysics.calculate_surface_gravity(self)
            # GEN.727: the eccentricity first, since the light a planet gets goes by its orbit's mean flux
            planetPhysics.generate_orbit_shape(self, primary_mass_kg, planned_objects, distance_override)
            planetPhysics.calculate_atmospheric_conditions(self, distance_override)
            planetPhysics.apply_class(self)
            if not (system_config.HABITABLE_WORLD is False and self.zone == 'e' and habitable_request.meets(self)):
                break
            draw.redrawn("planet.barred_world")
        planetPhysics.generate_orbital_motion_properties(self, primary_mass_kg)
        self._age_drawn = star.age  # the age this body was made at; the star's rises once the planets are in (`planets.refinish_rocky`, GEN.992)

        if moon_count != DEFER_MOONS:
            self.add_moons(moon_count)

    def add_moons(self, moon_count=None):
        """
        Generates this planet's moons: an exact `moon_count` (0 for none) or, with `None`, whatever `system_config.MOONS` and the
        satellite draw allow. Called by `__init__`, and by a request that drew the planet bare (`DEFER_MOONS`) and keeps it.
        Never for a moon.
        """
        if self.is_moon:
            return
        if moon_count is not None:
            planetPhysics.generate_moons(self, moon_count=moon_count)
        elif self.system_config.MOONS is not False and (
                self.system_config.MOONS is True
                or self._draw_moons()):
            planetPhysics.generate_moons(self)
        # GEN.959: the moons' history with the planet's spin (outward drift, the fall of a moon into a ring, the planet's day), from its birth spin
        if self.moons and getattr(self, "_birth_spin", None) is not None:
            planetPhysics.generate_spin(self, self.star.mass, birth=self._birth_spin)
        # GEN.767: the day and night shares, the moons' with their parent's shadow (no draws)
        daylight.generate(self)
        for moon in self.moons:
            daylight.generate(moon, self)
        # GEN.1034: the open-ocean equilibrium tide of a world with a sea, from its star and its moons (no draws)
        ocean_tides.generate(self)
        for moon in self.moons:
            ocean_tides.generate(moon, self)

    @property
    def giant_kind(self):
        return self.body_type

    @property
    def parent_kind(self):
        return "planet" if self.is_moon else "star"

    def orbit_mu_au3_per_year2(self):
        """AU^3/yr^2 of what the planet or moon orbits, from its circular
        orbit's radius and period (Kepler's third law), so it agrees with
        the velocity `update_orbital_position` gives it."""
        return (2.0 * math.pi / self.period) ** 2 * self.distance ** 3

    def _draw_moons(self):
        """Whether a planet left to chance gets moons: an even chance, or for a small body the share of the hosting gate's
        `satellite-share` table for its population and size (GEN.761), moved by the run's moons prevalence (GEN.52)."""
        chance = planetPhysics.satellite_share(self)
        if prevalence.percent(self.system_config, "moons"):
            return draw.random() < prevalence.scaled_chance(self.system_config, "moons", 0.5 if chance is None else chance)
        return draw.randrange(2) == 1 if chance is None else draw.random() < chance

    def _encode(self, data):
        """`habitable_zone` is stored as a list (JSON has no tuple type); `moons` are serialized the same way,
        recursively."""
        data["habitable_zone"] = list(self.habitable_zone)
        data["moons"] = [moon.to_dict() for moon in self.moons]

    @classmethod
    def from_dict(cls, data, star, system_config):
        """
        Reconstructs a `Planet` (or moon) from a dict in the shape
        `to_dict()` produces, without re-running generation. Recurses into
        `moons`, threading the same `star`/`system_config` into each one.

        Args:
            data (dict): A dict in the shape `to_dict()` produces.
            star (Star or BinaryStarProxy): The system's star, threaded in
                                            rather than stored per-planet
                                            (see `StarSystem.from_dict`).
            system_config (SystemConfig): The system's shared config.

        Returns:
            Planet: The reconstructed planet or moon.
        """
        planet = cls._rebuild(data, system_config, star=star)
        planet.moons = [cls.from_dict(moon_data, star, system_config) for moon_data in data.get("moons", [])]
        return planet

    def _decode(self, data):
        self.habitable_zone = tuple(data["habitable_zone"])

    @property
    def life_class_id(self):
        """The id of the life class that lives on the body (`"1a"`, `"4c"`, `"6"` ...), `None` for none."""
        return life_classes.class_id(self.life_class, self.life_class_variant)

    @property
    def has_timeline(self):
        """Whether the body has an evolutionary timeline (the system's age at which `life_stage` was reached), as
        against a class that only has a pace."""
        return self.life_system_age_gyr is not None

    @property
    def flavor_text(self):
        """The flavour line `flavor_key` names (`life.flavor_line`), or `None`."""
        from planetgen.generation import life
        return life.flavor_line(self.flavor_key)

    def gas_giant_depth(self):
        """
        The numbers behind a gas giant's depth sentences, from its already-computed reference-level pressure
        (`self.atmospheric_pressure`): the depth relative to the standard 1-atmosphere level from the isothermal
        barometric law used to derive that pressure (see planetPhysics.calculate_atmospheric_conditions); the
        depth and temperature at which rising pressure would turn the atmosphere from gas to liquid metallic
        hydrogen, from extending `self.surface_temperature` down a dry adiabat (see
        `constants.ADIABATIC_INDEX_H2_HE` and `HYDROGEN_METALLIZATION_PRESSURE_PA`); and a core-pressure
        estimate from the analytic n=1 polytrope (Lane-Emden) solution to hydrostatic equilibrium. None of
        these are stitched into one continuous profile -- the adiabat only holds for the shallow gas layer,
        and the polytrope only describes the deep interior, with a poorly understood transition in between
        that neither approximation covers.

        Returns:
            dict: `depth_km` (below the 1 atmosphere level; negative above it), `doubling_km` (the depth over
                which pressure doubles), `transition_depth_km`, `transition_temperature_k` and
                `core_pressure_pa`; or `None` if the scale height, pressure or temperature aren't available.
        """
        if not self.scale_height or not self.atmospheric_pressure or not self.surface_temperature:
            return None

        reference_pressure_pa = 101300  # ~1 standard atmosphere, matching the "atmospheres" conversion above
        transition_pressure_pa = constants.HYDROGEN_METALLIZATION_PRESSURE_PA
        adiabatic_exponent = (constants.ADIABATIC_INDEX_H2_HE - 1) / constants.ADIABATIC_INDEX_H2_HE
        radius_m = self.radius * 1000
        return {
            "depth_km": self.scale_height * math.log(self.atmospheric_pressure / reference_pressure_pa),
            "doubling_km": self.scale_height * math.log(2),
            "transition_depth_km": self.scale_height * math.log(transition_pressure_pa / self.atmospheric_pressure),
            "transition_temperature_k": self.surface_temperature
            * (transition_pressure_pa / self.atmospheric_pressure) ** adiabatic_exponent,
            "core_pressure_pa": math.pi * constants.G * self.mass ** 2 / (8 * radius_m ** 4),
        }

    def _gas_giant_pressure_depth_text(self):
        """
        Short sentences putting the gas giant's reference-level pressure into physical context (see
        `gas_giant_depth` for what each figure is and where it holds).

        Returns:
            list[str]: Descriptive sentences, or an empty list if `gas_giant_depth` has nothing.
        """
        depth = self.gas_giant_depth()
        if depth is None:
            return []
        depth_km = depth["depth_km"]
        doubling_km = depth["doubling_km"]
        if depth_km > 0:
            depth_sentence = f"That's about {format_number(depth_km)} km below the 1 atmosphere level, where pressure roughly doubles every {format_number(doubling_km)} km."
        else:
            depth_sentence = f"That's about {format_number(abs(depth_km))} km above the 1 atmosphere level, where pressure roughly doubles every {format_number(doubling_km)} km with depth."

        return [
            depth_sentence,
            f"Around {format_number(depth['transition_depth_km'])} km down, rising pressure would compress the atmosphere into a liquid metallic state near {format_temperature_k(depth['transition_temperature_k'])}.",
            f"Deep in the interior, a simplified polytrope model estimates core pressure on the order of {format_pressure_pa(depth['core_pressure_pa'])}.",
        ]

    def habitability_world(self):
        """
        The `habitability.World` of a scored rocky body and the compact-object host it was rated against
        (`None` for an ordinary star), or `None` for a gas giant or a body saved without scores. A planet
        of a pulsar, neutron star or black hole is rated as a lethal dose, as its stored score was.
        """
        if self.body_type != "t" or getattr(self, "eq_body", None) is None:
            return None
        values = habitability_world.values_from_body(self)
        if values.get("surface_temperature_k") is None:
            return None
        host = habitability_world.compact_host_name(getattr(self, "star", None))
        return habitability_world.world_from(values, lethal_host=host), host

from planetgen.generation import object_timing  # noqa: E402 -- registers the generators below (PERF.92)
object_timing.track(Planet, "planet", kind_after=lambda planet: "moon" if planet.is_moon else "planet")
