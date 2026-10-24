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

from planetgen.generation import prevalence
from planetgen.generation.config import SystemConfig
from planetgen.physics import aerosols, atmosphere, climate_draws, constants, habitability_world, hydrosphere, magnetism, planets as planetPhysics, radiation, radiogenic
from planetgen.physics.position import axis_property, velocity_axis_property
from planetgen import tuning
from planetgen.util import draw
from planetgen.generation.planetary_object import PlanetaryObject
from planetgen.util.format import (
    format_number, to_paragraph, format_pressure_pa, format_temperature_k,
)


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
        life_chemical (str): The primary chemical basis for any potential life
                             (set by `life.apply_life_data`, None until then).
        evolutionary_speed (str): How fast life evolves here (e.g. 'fast',
                                  'normal', 'slow'); set by `life.apply_life_data`,
                                  None until then.
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
                                            fully orients the (circular)
                                            orbital plane in 3D.
        velocity_x_kms, velocity_y_kms, velocity_z_kms (float): This body's
                                   velocity, km/s, relative to its orbital
                                   anchor (the same one `position_x/y/z`
                                   is measured from): the tangent of its
                                   circular orbit where it is,
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
        orbital_speed_kms (float): This body's circular orbital speed, in
                                   km/s -- constant around a circular orbit,
                                   so (unlike position) unaffected by
                                   the body moving; only changes if
                                   `distance`/`period` do
                                   (see `StarSystem.validate_system`).
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
        p_o2_kpa ... p_so2_kpa (float): Each `physics.atmosphere.SPECIES`
                                        gas's partial pressure, kPa.
        life_stage (str or None): The latest `evolution.MILESTONE_KEYS` milestone the timeline reached (set by
                                  `life.apply_life_data`; `None` before any, or with no timeline).
        evolutionary_pace (str or None): The timeline scale (`fast`, `normal`, `slow`); `None` when the body has no
                                         timeline.
        life_milestone_age_gyr, life_system_age_gyr (float or None): The milestone's age and the system's, Gy.
        flavor_key (str or None): The flavour line `life.decide_flavor_text` drew, as `category:index` into the
                                  `tuning` flavour lists (`flavor_text` reads it).
        flavor_text_count (int): The number of flavor texts added to this planet.
    """

    SERIALIZABLE_FIELDS = [
        "is_moon", "zone", "description", "atm_molar_density", "gravity",
        "atm_density", "surface_temperature", "density", "atmospheric_pressure",
        "mass", "atmosphere", "composition", "radius", "planet_class", "distance",
        "body_type", "scale_height", "hill_radius", "min_orbit_distance", "name",
        "life_chemical", "evolutionary_speed", "reflection_spectrum_visible",
        "reflection_spectrum_non_visible", "life_stage", "evolutionary_pace",
        "life_milestone_age_gyr", "life_system_age_gyr", "flavor_key",
        "flavor_text_count", "habitable_zone", "volume", "period",
        "orbital_inclination_deg", "orbital_ascending_node_deg",
        "position_x", "position_y", "position_z", "orbital_speed_kms",
        "velocity_x_kms", "velocity_y_kms", "velocity_z_kms",
        "min_update_interval_years",
        "rotation_period_hours", "spin_axis_x", "spin_axis_y", "spin_axis_z", "axial_tilt_deg",
        *atmosphere.ATMOSPHERE_FIELDS,
        *climate_draws.CLIMATE_FIELDS,
        *magnetism.BODY_FIELDS,
        *hydrosphere.HYDROSPHERE_FIELDS,
        *radiogenic.BODY_FIELDS,
        *radiation.BODY_FIELDS,
        *habitability_world.BODY_FIELDS,
        *aerosols.BODY_FIELDS,
        "reflex_offset_x", "reflex_offset_y", "reflex_offset_z",
        "born_beyond_snow_line",
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
                 radius=None, planet_class=None, mass=None, zone_override=None, distance_override=None,
                 is_moon=False, moon_count=None, primary_mass_kg=None):
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
            planet_class (str, optional): The class of the planet (e.g., 'M', 'N').
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
            primary_mass_kg (float, optional): The mass (in kg) of the body
                                               this one actually orbits, for
                                               Kepler's third law. Defaults to
                                               `star.mass` -- correct for an
                                               ordinary planet, but a moon
                                               orbits its parent planet, not
                                               the star (`generate_moons`
                                               passes the parent planet's own
                                               `mass` here).
        """
        self.system_config = system_config # Store SystemConfig
        self.is_moon = is_moon
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
        self.planet_class = planet_class
        self.distance = distance
        self.body_type = None
        self.scale_height = None
        self.hill_radius = None
        self.min_orbit_distance = None
        self.name = None # Set by the owning StarSystem (names.bodies.name_bodies) once orbits are final.
        self.life_chemical = None
        self.evolutionary_speed = None
        self.reflection_spectrum_visible = None
        self.reflection_spectrum_non_visible = None
        self.life_stage = None # Set later by life.apply_life_data, with the next three
        self.evolutionary_pace = None
        self.life_milestone_age_gyr = None
        self.life_system_age_gyr = None
        self.flavor_key = None # Set later by life.decide_flavor_text
        self.flavor_text_count = 0 # Initialize flavor text count for this planet
        self.orbital_inclination_deg = None
        self.orbital_ascending_node_deg = None
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
        for field in atmosphere.ATMOSPHERE_FIELDS + magnetism.BODY_FIELDS + hydrosphere.HYDROSPHERE_FIELDS + radiogenic.BODY_FIELDS + radiation.BODY_FIELDS + habitability_world.BODY_FIELDS + aerosols.BODY_FIELDS:
            setattr(self, field, None)

        # From the star, should not be changed.
        self.habitable_zone = habitable_zone
        self.star = star # Store the Star object

        # The primary is the star for an ordinary planet, but a moon orbits
        # its parent planet -- primary_mass_kg (defaulted to star.mass)
        # carries that distinction in from generate_moons below. Its Hill
        # sphere (set inside generate_planet_properties) uses it (GEN.138).
        primary_mass_kg = primary_mass_kg if primary_mass_kg is not None else star.mass
        self.primary_mass_kg = primary_mass_kg
        # Calculate physical/orbital properties (no life data yet).
        planetPhysics.generate_planet_properties(self, zone_override)  # sets self.volume (km^3) and self.mass
        self.period = planetPhysics.calculate_orbital_period_years(self.distance, primary_mass_kg)
        planetPhysics.calculate_surface_gravity(self)
        planetPhysics.calculate_atmospheric_conditions(self, distance_override)
        planetPhysics.generate_orbital_motion_properties(self, primary_mass_kg)

        if not self.is_moon:
            if moon_count is not None:
                planetPhysics.generate_moons(self, moon_count=moon_count)
            elif self.system_config.MOONS is not False and (
                    self.system_config.MOONS is True
                    or self._draw_moons()):
                planetPhysics.generate_moons(self)

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
        """Whether a planet left to chance gets moons: an even chance, moved
        by the run's moons prevalence (GEN.52)."""
        if prevalence.percent(self.system_config, "moons"):
            return draw.random() < prevalence.scaled_chance(self.system_config, "moons", 0.5)
        return draw.randrange(2) == 1

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

    def _generate_life_and_flavor_paragraphs(self, object_type_desc, sentences):
        """
        Generates paragraphs related to life, evolution, and flavor text for the planet.

        Reads (but does not compute) `self.life_chemical`, `self.evolutionary_speed`,
        `self.evolutionary_data`, and `self.flavor_text`, all of which are
        set by `life.apply_life_data`/`life.decide_flavor_text`
        before this is ever called (`to_paragraph_list` is only invoked once
        a `StarSystem` is fully generated) -- this method is a pure read and
        may be called more than once without re-rolling or re-mutating
        `system_config.system_flavor_count`.

        Args:
            object_type_desc (str): Description of the object type ("planet" or "moon").
            sentences (list): A list of sentences to append life-related descriptions to.

        Returns:
            list: A list of generated paragraphs.
        """
        life_paragraphs = []

        if self.life_chemical:
            if self.evolutionary_speed:
                sentences.append(
                    f"The {object_type_desc}'s conditions are suitable for the development of life based on {self.life_chemical.lower()}, which is expected to evolve at a {self.evolutionary_speed} pace.")
            else:
                sentences.append(
                    f"The {object_type_desc}'s conditions are suitable for the development of life based on {self.life_chemical.lower()}.")
            if self.reflection_spectrum_visible:
                visible_spectrum = ", ".join(self.reflection_spectrum_visible)
                sentences.append(
                    f"The visible reflection spectrum of the life on this {object_type_desc} is characterized by {visible_spectrum.lower()}.")
            if self.reflection_spectrum_non_visible:
                non_visible_spectrum = ", ".join(self.reflection_spectrum_non_visible)
                sentences.append(f"In the non-visible spectrum, it exhibits {non_visible_spectrum.lower()}.")

        life_paragraphs.append(to_paragraph(sentences))

        # Add the evolutionary timeline data if available
        if self.evolutionary_data:
            life_paragraphs.extend(self.evolutionary_data)

        # Flavor text was already decided at generation time (life.decide_flavor_text);
        # this is a pure read so repeated renders don't re-roll or double-count it.
        if self.flavor_text:
            life_paragraphs.append(f"Sensors show {self.flavor_text}")

        return life_paragraphs

    @property
    def evolutionary_data(self):
        """The evolutionary timeline paragraph written from the stored structure (`evolution.timeline_paragraphs`),
        `[]` for a body with no timeline."""
        if self.evolutionary_pace is None:
            return []
        from planetgen.generation import evolution
        return evolution.timeline_paragraphs({
            "scale": self.evolutionary_pace, "stage": self.life_stage,
            "milestone_age_gyr": self.life_milestone_age_gyr, "system_age_gyr": self.life_system_age_gyr})

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
