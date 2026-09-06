# stellarObjects/starData.py

"""
Star Generation and Properties
==============================

This module defines the `Star` class, the central component for creating and
managing the properties of a star within a procedurally generated star system.
The class handles the generation of fundamental stellar characteristics such as
mass, radius, temperature, and luminosity based on astrophysical models and
probabilistic distributions.

The generation process can be random or guided by specific parameters like
spectral class. Once the core properties are established, the class calculates
derived attributes, including the habitable zone, the gravitational extent of
the system (Hill sphere), and the size of the star's stellar wind bubble
(heliosphere).
"""

import math
import random
import re

from .config import SystemConfig
from .names import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from . import physical_constants, program_constants
from .utils import (format_age_string, calculate_habitable_zone, calculate_hill_sphere,
                    format_length_km, format_relative_to_sol, generate_phoneme_salad_name,
                    get_star_evolutionary_profile, properties_to_string, reseed_rng)

class Star:
    """
    Represents a single star, encapsulating its physical and orbital properties.

    This class orchestrates the procedural generation of a star, from its
    spectral type and luminosity to its mass and radius. It serves as the
    gravitational and energetic center of a `StarSystem`, defining the context
    for all orbiting planets and other celestial bodies.

    The generation process can be random or guided by specific parameters, such
    as forcing a large star or specifying a precise star type (e.g., 'G2V').
    Once the core properties are established, it calculates derived
    attributes, including the habitable zone, the system's gravitational
    perimeter, and the heliosphere's radius.

    Note: `self.type` here is a full descriptive string (e.g.
    "G2V Yellow Main Sequence Star"), unrelated in shape and purpose to
    `Planet.body_type`/`AsteroidBelt.body_type` (single-character 't'/'g'/'a'
    discriminators used to tell those two classes apart within `StarSystem`'s
    mixed `planets` list) -- the two just happen to share a common English
    word as their name.
    """

    def _calculate_initial_star_age_and_lifespan(self):
        """
        Calculates the star's initial age and lifespan, using a method that
        depends on its Yerkes luminosity class.

        For a main-sequence (Yerkes 'V') star, spectral letter fully
        determines both current temperature and total lifespan, so lifespan
        is drawn from `STAR_EVOLUTION`'s fixed per-letter range.

        For any evolved/remnant star (giants, subgiants, bright giants,
        supergiants, hypergiants, subdwarfs), the current spectral letter
        reflects only present-day temperature -- an "M0III" red giant did
        not spend its life as an M-type star -- so `STAR_EVOLUTION` can't be
        used directly. Instead, the star's own already-generated mass is run
        through the same mass-luminosity-based main-sequence lifetime scaling
        used to derive that table, extended by the standard rule of thumb
        that ~90% of a star's total lifetime is spent on the main sequence;
        age is then drawn from the remaining post-main-sequence window (the
        star must already have completed its main-sequence phase to be
        observed in one of these classes).

        For white dwarfs, the lifespan is considered effectively infinite as
        they slowly cool, and age is instead drawn from a dedicated
        cooling-age range, since a white dwarf's spectral letter reflects its
        current temperature rather than a progenitor mass usable for either
        of the above.

        The initial age is randomly chosen to be a fraction of its total lifespan,
        ensuring the star is in a plausible evolutionary stage. This age can be
        further adjusted later by `adjust_age_for_planets` to accommodate the
        evolutionary requirements of its planets.

        Returns:
            tuple: A tuple containing:
                - age (float): The initial age of the star in billions of years (GY).
                - lifespan (float): The total expected lifespan of the star in
                                    billions of years (GY), or `float('inf')` for white dwarfs.
        """
        reseed_rng()

        if self.yerkes_class in ["VII", "D"]:
            # A white dwarf's spectral letter encodes its *current* surface
            # temperature, not its progenitor's mass/spectral class, so it
            # can't be looked up in STAR_EVOLUTION's main-sequence lifespan
            # table (that would e.g. treat an "F5VII" white dwarf as if it
            # were an F-type main-sequence star with a ~4-9 Gy lifespan). A
            # white dwarf's remaining "lifespan" is effectively infinite --
            # it just cools for trillions of years -- so age is drawn from a
            # dedicated cooling-age range instead of a lifespan fraction.
            lifespan = float('inf')
            min_age = program_constants.WHITE_DWARF_MIN_AGE_GY
            max_age = program_constants.WHITE_DWARF_MAX_AGE_GY
            if self.system_config.AGE == "old":
                min_age = min_age + (max_age - min_age) * program_constants.OLD_STAR_AGE_LIFESPAN_RATIO
            elif self.system_config.AGE == "young":
                max_age = max_age - (max_age - min_age) * (1 - program_constants.YOUNG_STAR_AGE_LIFESPAN_RATIO)
            age = random.uniform(min_age, max_age)
            return age, lifespan

        if self.yerkes_class != "V":
            # Evolved-star track (giants, subgiants, bright giants,
            # supergiants, hypergiants, subdwarfs): the current spectral
            # letter is just this star's present temperature, not a
            # main-sequence class it belongs to -- an "M0III" red giant did
            # not spend its whole life as an M-type star. Derive the star's
            # actual main-sequence lifespan from its own already-generated
            # mass instead, then extend it for the (much shorter)
            # post-main-sequence phase it's currently in.
            mass_sol = self.mass / physical_constants.SOLAR_MASS_TO_KG
            ms_lifespan = program_constants.SOLAR_MS_LIFESPAN_GY * mass_sol ** program_constants.MS_LIFESPAN_MASS_EXPONENT
            lifespan = ms_lifespan / program_constants.MS_LIFESPAN_FRACTION_OF_TOTAL

            # The star must already have completed its main-sequence phase to
            # be observed as an evolved class, so age is drawn from the
            # (comparatively short) post-main-sequence window only.
            min_age = ms_lifespan
            max_age = lifespan
            if self.system_config.AGE == "old":
                min_age = min_age + (max_age - min_age) * program_constants.OLD_STAR_AGE_LIFESPAN_RATIO
            elif self.system_config.AGE == "young":
                max_age = max_age - (max_age - min_age) * (1 - program_constants.YOUNG_STAR_AGE_LIFESPAN_RATIO)
            age = random.uniform(min_age, max_age)

            return age, lifespan

        spectral_class_char = self.type[0]
        star_info = program_constants.STAR_EVOLUTION.get(spectral_class_char, {})

        min_lifespan, max_lifespan = star_info["lifespan_gy"]
        lifespan = random.uniform(min_lifespan, max_lifespan)

        # MIN_INITIAL_STAR_AGE_GY is a sensible floor for long-lived stars (F/G/K/M),
        # but a fixed 100-million-year floor is nonsensical for a star whose entire
        # lifespan is shorter than that (O/B, and short-lived rolls of A) -- it would
        # make min_age exceed max_age below, silently discarding the young/old
        # request and any real age variance for those types. Scale it down for
        # short-lived stars so it never exceeds a small fraction of their own lifespan.
        min_age = min(program_constants.MIN_INITIAL_STAR_AGE_GY,
                      lifespan * program_constants.MIN_INITIAL_STAR_AGE_LIFESPAN_RATIO)
        max_age = lifespan * program_constants.MAX_INITIAL_STAR_AGE_LIFESPAN_RATIO

        # Values are calculated to be either in last ratio of the
        # lifespan (i.e. the last 1/3 of the star's life or the
        # first 1/3 of a star's life).
        if self.system_config.AGE == "old":
            min_age = min_age + lifespan * program_constants.OLD_STAR_AGE_LIFESPAN_RATIO
        elif self.system_config.AGE == "young":
            max_age = max_age - lifespan * (1 - program_constants.YOUNG_STAR_AGE_LIFESPAN_RATIO)

        age = random.uniform(min_age, max_age) # Ensure age is less than lifespan

        return age, lifespan

    def adjust_age_for_planets(self, planets):
        """
        Adjusts the star's age to ensure it is at least as old as the minimum
        age required for the most "evolved" planet in its system. This method
        should be called by the `StarSystem` class after planets have been generated.

        The star's age is a crucial factor in planetary evolution. This function
        iterates through all generated planets, identifies their minimum required
        ages based on their `planet_class` and the star's `supported_evolutionary_scales`,
        and then ensures the star's age meets or exceeds the highest of these
        minimum requirements. This prevents scenarios where a young star might
        host an "old" planet (e.g., a terraformed world requiring billions of years
        of development).

        If the star's lifespan is too short to accommodate the required
        planetary age even at its very end, the star's age is set to be near
        the end of its lifespan (the closest it can get), reflecting a system
        that has evolved as rapidly as physically possible for its type. If
        the requirement *is* reachable within the lifespan, `config.AGE`
        ("young"/"old") still biases where in the achievable window the age
        falls, the same way it biases the star's initial age in
        `_calculate_initial_star_age_and_lifespan` — a forced minimum age
        should narrow the age range, not silently discard a young/old
        request entirely.

        Args:
            planets (list): A list of `Planet` objects in the system.
        """
        reseed_rng()
        # Uses get_star_evolutionary_profile rather than a raw
        # STAR_EVOLUTION[spectral_class] lookup so giants/supergiants/white
        # dwarfs etc. (whose spectral letter reflects only current
        # temperature, not a main-sequence lifespan) are handled correctly.
        star_evolution_data = get_star_evolutionary_profile(self)
        supported_scales = star_evolution_data.get("supported_evolutionary_scales", [])

        min_required_age_for_system = 0.0

        for planet in planets:
            # Assuming planet has a 'planet_class' attribute
            if hasattr(planet, 'planet_class') and planet.planet_class in program_constants.PLANET_CLASSES:
                planet_class_data = program_constants.PLANET_CLASSES[planet.planet_class]
                age_ranges = planet_class_data.get("age_ranges", {})

                for scale in supported_scales:
                    if scale in age_ranges:
                        min_planet_age_for_scale, _ = age_ranges[scale]
                        min_required_age_for_system = max(min_required_age_for_system, min_planet_age_for_scale)

        # Ensure the star's age is at least the minimum required by its planets
        if self.age < min_required_age_for_system:
            if self.lifespan != float('inf'):
                max_reachable_age = self.lifespan * program_constants.MAX_PLANET_AGE_ADJUSTMENT_FACTOR
            else:
                max_reachable_age = min_required_age_for_system + program_constants.WHITE_DWARF_AGE_ADDITION_GY # Add 5 GY for WD if no upper bound

            if self.lifespan != float('inf') and max_reachable_age < min_required_age_for_system:
                # The requirement is physically unreachable within this star's lifespan
                # (e.g. a habitable world's minimum age exceeds a short-lived O/B/A star's
                # entire life). The best it can do is sit near the very end of its life --
                # there's no real "young" vs "old" choice left to honor here.
                self.age = random.uniform(self.lifespan * program_constants.UNREACHABLE_PLANET_AGE_MIN_LIFESPAN_RATIO, max_reachable_age)
            else:
                # The requirement is reachable within the lifespan: bias where in the
                # achievable [min_required_age_for_system, max_reachable_age] window the
                # age falls, mirroring the young/old bias used for the initial age roll.
                low, high = min_required_age_for_system, max_reachable_age
                if self.system_config.AGE == "old":
                    low = low + (high - low) * program_constants.OLD_STAR_AGE_LIFESPAN_RATIO
                elif self.system_config.AGE == "young":
                    high = high - (high - low) * (1 - program_constants.YOUNG_STAR_AGE_LIFESPAN_RATIO)
                self.age = random.uniform(low, high)

        # Ensure age doesn't exceed lifespan (unless lifespan is infinite)
        if self.lifespan != float('inf') and self.age >= self.lifespan:
            self.age = self.lifespan * program_constants.MAX_PLANET_AGE_ADJUSTMENT_FACTOR # Star is near end of life

    def calculate_system_perimeter(self):
        """
        Calculates the Hill sphere for the star system relative to the galaxy.

        The Hill sphere defines the gravitational sphere of influence of the star
        within the larger galactic environment. Any object orbiting the star
        beyond this radius would be more strongly perturbed by the galaxy's
        gravity than by the star's, making it unlikely to be gravitationally
        bound to the star system.

        This method uses the star's mass, the estimated mass of the Milky Way,
        and the star's distance from the galactic center to compute this boundary.
        The result is converted from meters to Astronomical Units (AU) for
        convenience and consistency with other system-scale measurements.

        Returns:
            float: The radius of the Hill sphere in Astronomical Units (AU).
        """
        galactic_center_dist_m = physical_constants.GALACTIC_CENTER_DISTANCE_LY * physical_constants.LY_TO_M
        hill_radius_m = calculate_hill_sphere(galactic_center_dist_m, self.mass, physical_constants.MILKY_WAY_MASS)
        return hill_radius_m / physical_constants.AU_TO_M

    @staticmethod
    def _calculate_heliosphere_radius_static(mass, luminosity, radius_km, star_type, yerkes_class):
        """
        Estimates the radius of the star's heliosphere (astrosphere).

        The heliosphere is a "bubble" in space created by the stellar wind, a
        stream of charged particles flowing from the star. The edge of this
        bubble, the heliopause, is where the stellar wind's pressure is balanced
        by the interstellar medium (ISM).

        This function employs a tiered approach to model the stellar wind's
        mass-loss rate and velocity, which are highly dependent on the star's
        evolutionary stage (Yerkes luminosity class) and spectral type. Three
        physical models are applied: the Nieuwenhuijzen & de Jager (1990)
        empirical mass-loss fit for every evolved star (giants through
        hypergiants), a radiation-driven-wind power law for hot O/B dwarfs,
        and an empirical coronal-wind scaling for cool main-sequence dwarfs
        (and white dwarfs) -- chosen to ensure physically plausible results
        across the stellar spectrum.

        The heliopause radius is then calculated based on the balance between
        the stellar wind's momentum flux and the pressure of the interstellar medium.

        Args:
            mass (float): The star's mass in kilograms.
            luminosity (float): The star's luminosity in Watts.
            radius_km (float): The star's radius in kilometers.
            star_type (str): The star's full spectral type string (e.g., 'G2V'),
                             whose first character (the spectral class) selects
                             between the hot-dwarf and cool-dwarf wind models
                             within the main-sequence tier.
            yerkes_class (str): The star's Yerkes luminosity class (e.g., 'V',
                                'III', '0'), which selects the wind/mass-loss model.

        Returns:
            float: The estimated radius of the heliosphere in Astronomical Units (AU).
        """
        # --- 1. Get Fundamental Stellar Properties (from arguments) ---
        radius_m = radius_km * physical_constants.KM_TO_M_FACTOR  # Convert radius from km to meters
        lum_sol = luminosity / physical_constants.SOLAR_LUMINOSITY
        mass_sol = mass / physical_constants.SOLAR_MASS_TO_KG
        radius_sol = radius_m / physical_constants.SOLAR_RADIUS_M
        spectral_class = star_type[0].upper() if star_type else None

        # --- 2. Calculate Mass-Loss Rate (M-dot) and Wind Velocity (v_inf) ---
        # A single formula for mass loss is insufficient. We use a tiered system based on
        # the star's Yerkes luminosity class (evolutionary stage) and spectral type.
        escape_velocity = math.sqrt(physical_constants.ESCAPE_VELOCITY_CONSTANT * physical_constants.G * mass / radius_m)


        # TIER 1: Evolved stars -- giants (III), subgiants (IV), bright giants
        # (II), supergiants (IB/IAB/IA/IA+), and hypergiants (0). These were
        # previously three separately hand-fit tiers (a hypergiant power-law,
        # plus Reimers' Law with two different eta values for "giants" vs
        # "supergiants") stitched together at Yerkes-class boundaries -- the
        # exact structure that let the old hypergiant constant drift ~9
        # orders of magnitude out of calibration unnoticed. Nieuwenhuijzen &
        # de Jager (1990) is a single empirical mass-loss fit spanning this
        # whole regime, so one formula and one continuous curve now covers
        # all of it (see physical_constants.py for the correction applied
        # above its known high-luminosity overestimation threshold).
        if yerkes_class in ["0", "IA+", "IA", "IAB", "IB", "II", "III", "IV"]:
            mass_loss_rate_smyr = (physical_constants.NDJ_MASS_LOSS_COEFFICIENT
                                   * (lum_sol**physical_constants.NDJ_LUMINOSITY_EXPONENT)
                                   * (mass_sol**physical_constants.NDJ_MASS_EXPONENT)
                                   * (radius_sol**physical_constants.NDJ_RADIUS_EXPONENT))
            if lum_sol > physical_constants.NDJ_HIGH_LUMINOSITY_THRESHOLD_LSUN:
                mass_loss_rate_smyr *= physical_constants.NDJ_HIGH_LUMINOSITY_CORRECTION_FACTOR

            if yerkes_class == "0":
                # Hypergiant winds are an extreme multiple of escape velocity.
                wind_velocity = physical_constants.HYPERGIANT_WIND_VELOCITY_FACTOR * escape_velocity
            else:
                # Wind velocity is a smaller fraction of escape velocity for these cooler giants/supergiants.
                wind_velocity = physical_constants.GIANT_WIND_VELOCITY_FACTOR * escape_velocity

        # TIER 2a: Hot main-sequence dwarfs (O and B, Class V/VI)
        # Unlike cool dwarfs, O/B stars have powerful radiation-driven winds (line-driven,
        # via the same physical mechanism as hypergiants/supergiants, just far less extreme):
        # mass loss scales steeply with luminosity rather than the coronal/magnetic-activity
        # scaling used for cool dwarfs below. Without this tier, a hot O-type star would be
        # (incorrectly) modeled with a weaker wind than a cool, dim M dwarf.
        elif spectral_class in ("O", "B") and yerkes_class in ["V", "VI"]:
            mass_loss_rate_smyr = physical_constants.OB_DWARF_MASS_LOSS_RATE_FACTOR * (lum_sol**physical_constants.OB_DWARF_MASS_LOSS_RATE_EXPONENT)
            wind_velocity = physical_constants.OB_DWARF_WIND_VELOCITY_FACTOR * escape_velocity

        # TIER 2b: Cool main-sequence dwarfs (A through M, Class V/VI) and White Dwarfs (VII)
        # For sun-like and cooler stars (plus stellar remnants), mass loss is very low. We
        # scale directly from the Sun's known properties for a stable, physically grounded result.
        else:
            # For cool main-sequence stars, mass loss is driven by magnetic activity, not
            # radiation pressure. The Reimers' law scaling (L*R/M) is inaccurate here.
            # We use a more modern, empirically-derived formula from Wood et al. (2005),
            # which relates mass-loss rate to surface X-ray flux.
            # M-dot ~ R^2 * F_x^1.34, and F_x ~ L^-0.5
            # This results in a scaling relationship of M-dot ~ R^2 * L^-0.5
            # We normalize this to the Sun's known mass-loss rate.
            # The mass term is added to account for gravitational binding.
            scaling_factor = (radius_sol**physical_constants.RADIUS_SOL_EXPONENT) * (lum_sol**physical_constants.LUMINOSITY_SOL_EXPONENT) * (mass_sol**physical_constants.MASS_SOL_EXPONENT)
            mass_loss_rate_smyr = physical_constants.SUN_MASS_LOSS_RATE_SOLAR_MASS_PER_YEAR * scaling_factor # 2e-14 is the Sun's M-dot in M_sol/yr

            # Scale wind velocity based on the star's escape velocity relative to the Sun's.
            wind_velocity = physical_constants.SOLAR_WIND_VELOCITY * (escape_velocity / physical_constants.SOLAR_ESCAPE_VELOCITY)

        # --- 3. Convert and Calculate Final Radius ---

        # Convert M-dot from (solar masses/year) to (kg/s)
        mass_loss_rate_kgs = mass_loss_rate_smyr * physical_constants.SOLAR_MASS_TO_KG / physical_constants.SECONDS_PER_YEAR

        # The heliopause radius is where the stellar wind's momentum flux balances the ISM pressure.
        # R = sqrt( (M-dot * v_inf) / (4 * pi * P_ism) )
        momentum_flux = max(mass_loss_rate_kgs * wind_velocity, physical_constants.MIN_MOMENTUM_FLUX) # Ensure momentum_flux is not zero or negative

        try:
            heliopause_radius_m = math.sqrt(momentum_flux / (physical_constants.FOUR_PI * physical_constants.ISM_PRESSURE))
        except ValueError:
            # This should not happen with the floor value, but as a final safety net.
            heliopause_radius_m = physical_constants.HELIOPAUSE_RADIUS_DEFAULT_M

        # Convert the final radius from meters to Astronomical Units (AU) for output.
        return heliopause_radius_m / physical_constants.AU_TO_M

    def calculate_heliosphere(self):
        """
        Estimates the radius of the star's heliosphere (astrosphere) by
        delegating to the `_calculate_heliosphere_radius_static` helper with
        this star's own mass, luminosity, radius, type, and Yerkes class.

        Returns:
            float: The estimated radius of the heliosphere in Astronomical Units (AU).
        """
        return Star._calculate_heliosphere_radius_static(self.mass, self.luminosity, self.radius, self.type, self.yerkes_class)

    def __init__(self, system_config: SystemConfig, name=None, _skip_property_init=False, **kwargs):
        """
        Initializes a Star object, generating its properties based on specified constraints.

        This constructor orchestrates the creation of a star. It can generate a
        star randomly or be guided by specific parameters provided through the
        `config` module. The `config.STAR_TYPE` parameter allows for the creation
        of a star with a specific spectral class, subclass, and Yerkes
        classification (e.g., 'G2V'). If `config.LARGE_STAR` is True, this
        influences the star's initial generation to favor larger, more massive
        stars.

        If a `name` is provided, it will be used for the star; otherwise, a
        random name will be generated using `generate_phoneme_salad_name`.

        Once the core properties (`type`, `radius`, `mass`, `temperature`, `luminosity`)
        are established by `generate_star()`, this constructor calculates derived
        properties such as the `age`, `lifespan`, `habitable_zone`,
        `system_perimeter` (Hill sphere), and `heliosphere_radius`.

        Args:
            system_config (SystemConfig): The shared SystemConfig object for the system.
            name (str, optional): The name to assign to the star. If None, a
                                  random name will be generated. Defaults to None.
            _skip_property_init (bool): If True, skips initialization of properties
                                        that might be managed by a subclass (e.g., BinaryStarProxy).
                                        Defaults to False.
        """
        self.system_config = system_config # Storing the SystemConfig instance
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)
        
        if not _skip_property_init:
            self.luminosity = None
            self.temperature = None
            self.mass = None
            self.radius = None
            self.type = None
            self.yerkes_class = None
            self.habitable_zone = None
            self.system_perimeter = None
            self.heliosphere_radius = None
            # Added mass_override for secondary star generation in binary systems
            self.generate_star(mass_override=kwargs.get('mass_override'))
            self.age, self.lifespan = self._calculate_initial_star_age_and_lifespan()
            self.habitable_zone = calculate_habitable_zone(self.luminosity)
            self.system_perimeter = self.calculate_system_perimeter()
            self.heliosphere_radius = self.calculate_heliosphere()

    def to_paragraph_list(self):
        """
        Generates a list of descriptive paragraphs about the star's properties.

        This method formats the star's physical characteristics (type, radius,
        mass, temperature, luminosity, habitable zone) into a structured block
        using `properties_to_string`, delegating the mass/luminosity/radius
        value formatting to `format_relative_to_sol`/`format_length_km`. It
        also constructs a sentence describing the star's age and expected
        lifespan (via `format_age_string`), incorporating evolutionary notes
        from `STAR_EVOLUTION` if available.

        The output is designed to be human-readable and can be formatted either
        as Wikitext or Markdown based on the `config.MARKDOWN` flag.

        Returns:
            list: A list of strings, where each string represents a paragraph
                  or a formatted data block describing the star.
        """
        paragraphs = []

        if round(self.habitable_zone[0], program_constants.ROUND_HABITABLE_ZONE_AU) == round(self.habitable_zone[1], program_constants.ROUND_HABITABLE_ZONE_AU):
            hab_lower = str(round(self.habitable_zone[0], program_constants.ROUND_HABITABLE_ZONE_AU_SMALL))
            hab_upper = str(round(self.habitable_zone[1], program_constants.ROUND_HABITABLE_ZONE_AU_SMALL))
        else:
            if self.habitable_zone[0] < program_constants.PERCENT_SOL_THRESHOLD_LOW:
                hab_lower = str(round(self.habitable_zone[0], program_constants.ROUND_HABITABLE_ZONE_AU_SMALL))
            else:
                hab_lower = str(round(self.habitable_zone[0], program_constants.ROUND_HABITABLE_ZONE_AU))

            if self.habitable_zone[1] < program_constants.PERCENT_SOL_THRESHOLD_LOW:
                hab_upper = str(round(self.habitable_zone[1], program_constants.ROUND_HABITABLE_ZONE_AU_SMALL))
            else:
                hab_upper = str(round(self.habitable_zone[1], program_constants.ROUND_HABITABLE_ZONE_AU))

        mass_string = format_relative_to_sol(self.system_config, self.mass, physical_constants.SOLAR_MASS_TO_KG, "kg", low_percent_precision=2)
        lum_string = format_relative_to_sol(self.system_config, self.luminosity, physical_constants.SOLAR_LUMINOSITY, "W", low_percent_precision=4)

        radius_string = format_length_km(self.system_config, self.radius, program_constants.RADIUS_KM_SCIENTIFIC_NOTATION_THRESHOLD,
                                         program_constants.ROUND_RADIUS_KM, program_constants.SCIENTIFIC_NOTATION_DECIMAL_PLACES)

        star_properties = {
            "type": self.type,
            "radius": radius_string,
            "mass": mass_string,
            "temp": f"{self.temperature} K",
            "lum": lum_string,
            "hab": f"Between {hab_lower} and {hab_upper} AU",
            "loc": self.name # Adding the star's name as location
        }

        markdown_key_map = {
            "type": "Type",
            "radius": "Radius",
            "mass": "Mass",
            "temp": "Temperature",
            "lum": "Luminosity",
            "hab": "Habitable Zone",
            "loc": "Location"
        }
        
        star_block = properties_to_string(self.system_config, star_properties, "Star Data", markdown_key_map=markdown_key_map)
        paragraphs.append(star_block)

        # Construct the age and evolutionary notes sentence
        age_str = format_age_string(self.age)

        if self.lifespan == float('inf'):
            age_sentence_base = f"The star is approximately {age_str} old and is now a white dwarf, which will cool for trillions of years"
        else:
            lifespan_str = format_age_string(self.lifespan)
            age_sentence_base = f"The star is approximately {age_str} old, with an expected lifespan of {lifespan_str}"

        full_age_and_notes_sentence = age_sentence_base
        # STAR_EVOLUTION's notes describe a *main-sequence* star's habitable-zone
        # lifespan window, keyed by spectral letter. For a giant, supergiant,
        # subdwarf, or white dwarf, that letter reflects only the star's current
        # temperature/color -- not a main-sequence evolutionary stage -- so the
        # notes (and their quoted lifespan figures) would be meaningless or
        # self-contradictory (e.g. a white dwarf whose age already exceeds the
        # "main-sequence lifespan" quoted for its current temperature class).
        if self.yerkes_class == "V":
            spectral_class_char = self.type[0]
            star_info = program_constants.STAR_EVOLUTION.get(spectral_class_char, {})
            if "evolutionary_constraint_notes" in star_info:
                # Append notes directly as they are a continuation, ensuring a space and period at the end
                full_age_and_notes_sentence += " and " + star_info["evolutionary_constraint_notes"]
        full_age_and_notes_sentence += "." # Ensure the entire combined sentence ends with a period

        paragraphs.append(full_age_and_notes_sentence)

        return paragraphs

    def __str__(self):
        """
        Returns a complete string representation of the star, suitable for display.

        This method compiles all the descriptive paragraphs generated by
        `to_paragraph_list()` and joins them with double newlines to create
        a well-formatted, multi-paragraph string.

        Returns:
            str: A formatted string describing the star and its properties.
        """
        return '\n\n'.join(self.to_paragraph_list())

    def generate_star(self, mass_override=None):
        """
        Generates the physical properties of the star based on a hierarchical,
        physically-constrained procedural model.

        The generation process is fundamentally different for user-specified stars
        (via `config.STAR_TYPE`) versus randomly generated ones.

        For specified stars (when `config.STAR_TYPE` is set):
        1.  **Parsing Input**: The method first parses the `config.STAR_TYPE` string
            (e.g., 'G2V') to extract the spectral class (O, B, A, F, G, K, M),
            subclass (0-9), and Yerkes luminosity class (e.g., V for Main Sequence).
            It raises a `ValueError` if the format is invalid.
        2.  **Temperature Calculation**: The star's surface temperature is calculated
            based on its spectral class and subclass, interpolating within the
            `TEMP_RANGES` data.
        3.  **Luminosity Determination**: A valid luminosity range is established
            by finding the *intersection* of the allowed luminosity ranges for
            both the spectral class (`SPECTRAL_LUMINOSITY_RANGES`) and the
            Yerkes class (`YERKES_LUMINOSITY_RANGES`). A random luminosity is
            then chosen from this overlapping range. Special handling is included
            for white dwarfs, which can have temporarily high luminosities.
        4.  **Mass Calculation**: The star's mass is determined based on the
            physical constraints defined for its Yerkes class in `YERKES_MASS_CONSTRAINTS`.
            For main-sequence stars, a mass-luminosity relation is used. For white
            dwarfs, an inverse mass-radius relation is considered, with mass
            constrained by the Chandrasekhar limit. For giants and supergiants,
            mass is randomly selected within their class-specific bounds.
        5.  **Radius Calculation**: The star's radius is calculated last. For white
            dwarfs, an inverse mass-radius relationship is applied. For all other
            star types, the Stefan-Boltzmann law is used, deriving radius from
            luminosity and temperature.

        For random stars (when `config.STAR_TYPE` is None):
        1.  **Spectral Class Selection**: A spectral class is chosen probabilistically
            based on their prevalence in the galactic population.
            `self.system_config.LARGE_STAR` being True biases this selection
            towards hotter, more massive stars.
        2.  **Luminosity Generation**: A luminosity is randomly chosen from the
            typical range for the selected spectral class.
        3.  **Yerkes Class Determination**: The Yerkes luminosity class (e.g.,
            Supergiant, Main Sequence, White Dwarf) is *derived* from the
            generated luminosity, using thresholds defined in `YERKES_LUMINOSITY_RANGES`.
        4.  **Temperature and Subclass Calculation**: The star's surface temperature
            is randomly selected within the `TEMP_RANGES` for its spectral class,
            and the subclass is then derived from this temperature.
        5.  **Mass and Radius Calculation**: Similar to the specified star path,
            mass is determined based on Yerkes class constraints and luminosity,
            and radius is calculated using the appropriate physical laws.

        This comprehensive approach ensures that all generated stars, whether
        user-specified or randomly created, adhere to the fundamental physical
        relationships between their core properties, resulting in astrophysically
        plausible stellar objects.
        """
        reseed_rng()
        yerkes_lookup = {
            "0": "Hypergiant", "IA+": "Luminous Supergiant", "IA": "Supergiant",
            "IAB": "Intermediate-size Luminous Supergiant", "IB": "Less Luminous Supergiant",
            "II": "Bright Giant", "III": "Giant", "IV": "Subgiant", "V": "Main Sequence",
            "VI": "Subdwarf", "VII": "White Dwarf", "D": "White Dwarf" # D is an alias for VII
        }

        if self.system_config.STAR_TYPE:
            # --- GENERATE STAR FROM SPECIFIED TYPE ---

            # 1. Parse and validate the specified star type string.
            match = re.match(r"([OBAFGKM])([0-9])(IA\+|IAB|VII|III|IA|IB|II|IV|VI|0|V|D)", self.system_config.STAR_TYPE.upper())
            if not match:
                raise ValueError("Invalid star type format. Expected format is e.g., G2V.")
            spectral_class, subclass_str, yerkes_class_str = match.groups()
            subclass = int(subclass_str)
            self.yerkes_class = yerkes_class_str
            yerkes_type = yerkes_lookup[yerkes_class_str]

            # 2. Calculate Temperature from spectral class and subclass.
            min_temp, max_temp = physical_constants.TEMP_RANGES[spectral_class]
            temp_range_size = max_temp - min_temp
            temperature = min_temp + (physical_constants.SUBCLASS_MAX_VALUE - subclass) * (temp_range_size / physical_constants.SUBCLASS_MAX_VALUE)
            temperature = int(round(temperature, program_constants.ROUND_TEMPERATURE_NEAREST_HUNDRED))

            # 3. Determine a physically valid Luminosity.
            # Find the overlapping luminosity range between the spectral and Yerkes classes.
            spec_min_lum, spec_max_lum = physical_constants.SPECTRAL_LUMINOSITY_RANGES[spectral_class]
            yerkes_min_lum, yerkes_max_lum = physical_constants.YERKES_LUMINOSITY_RANGES[yerkes_class_str]
            
            # Special case for hot, young white dwarfs, which can be temporarily very luminous.
            if yerkes_class_str in ["VII", "D"] and spectral_class in ["O", "B"]:
                yerkes_min_lum, yerkes_max_lum = 0.1, 100

            valid_min_lum = max(spec_min_lum, yerkes_min_lum)
            valid_max_lum = min(spec_max_lum, yerkes_max_lum)
            
            if valid_min_lum > valid_max_lum:
                # If there's no overlap, the combination is physically unlikely.
                # We will prioritize the Yerkes class range, as it's the primary
                # driver of luminosity for evolved stars.
                valid_min_lum, valid_max_lum = yerkes_min_lum, yerkes_max_lum

            luminosity = random.uniform(valid_min_lum, valid_max_lum)

        else:
            # --- GENERATE STAR RANDOMLY ---

            # 1. Generate Spectral Class based on galactic population.
            if self.system_config.LARGE_STAR:
                spectral_probabilities = program_constants.SPECTRAL_PROBABILITIES_LARGE_STAR
            else:
                spectral_probabilities = program_constants.SPECTRAL_PROBABILITIES_NORMAL
            spectral_class = random.choices(list(spectral_probabilities.keys()), weights=spectral_probabilities.values(), k=1)[0]

            # 2. Generate Luminosity from the spectral class's typical range.
            min_luminosity, max_luminosity = physical_constants.SPECTRAL_LUMINOSITY_RANGES[spectral_class]
            luminosity = random.uniform(min_luminosity, max_luminosity)

            # 3. Determine Yerkes Class from the resulting luminosity.
            if luminosity > physical_constants.YERKES_LUMINOSITY_RANGES["0"][0]:
                self.yerkes_class, yerkes_type = "0", "Hypergiant"
            elif luminosity > physical_constants.YERKES_LUMINOSITY_RANGES["IA"][0]:
                self.yerkes_class, yerkes_type = "IA", "Supergiant"
            elif luminosity > physical_constants.YERKES_LUMINOSITY_RANGES["IAB"][0]:
                self.yerkes_class, yerkes_type = "IAB", "Intermediate-size Luminous Supergiant"
            elif luminosity > physical_constants.YERKES_LUMINOSITY_RANGES["IB"][0]:
                self.yerkes_class, yerkes_type = "IB", "Less Luminous Supergiant"
            elif luminosity > physical_constants.YERKES_LUMINOSITY_RANGES["II"][0]:
                self.yerkes_class, yerkes_type = "II", "Bright Giant"
            elif luminosity > physical_constants.YERKES_LUMINOSITY_RANGES["III"][0]:
                self.yerkes_class, yerkes_type = "III", "Giant"
            elif luminosity > physical_constants.YERKES_LUMINOSITY_RANGES["IV"][0]:
                self.yerkes_class, yerkes_type = "IV", "Subgiant"
            elif luminosity > physical_constants.SPECTRAL_LUMINOSITY_RANGES["M"][0]: # Check against dimmest main sequence
                self.yerkes_class, yerkes_type = "V", "Main Sequence"
            else:
                self.yerkes_class, yerkes_type = "VII", "White Dwarf"

            # 4. Calculate Temperature and Subclass.
            min_temp, max_temp = physical_constants.TEMP_RANGES[spectral_class]
            temperature = int(round(random.uniform(min_temp, max_temp), program_constants.ROUND_TEMPERATURE_NEAREST_HUNDRED))
            temp_range_size = max_temp - min_temp
            subclass = physical_constants.SUBCLASS_MAX_VALUE - round((temperature - min_temp) / temp_range_size * physical_constants.SUBCLASS_MAX_VALUE)

        # --- CALCULATE FINAL PROPERTIES (COMMON TO BOTH PATHS) ---

        # 5. Calculate Mass based on Yerkes class constraints.
        min_mass, max_mass = physical_constants.YERKES_MASS_CONSTRAINTS[self.yerkes_class]
        
        if self.yerkes_class == "V":
            # For main-sequence stars, the mass-luminosity relation is strong.
            if mass_override:
                mass_sol = mass_override / physical_constants.SOLAR_MASS_TO_KG
            else:
                mass_sol = luminosity ** (1 / physical_constants.MAIN_SEQUENCE_MASS_LUMINOSITY_EXPONENT)
        elif self.yerkes_class in ["VII", "D"]:
            # For white dwarfs, mass is tightly constrained. Hotter (younger) ones
            # are typically more massive, closer to the Chandrasekhar limit.
            if spectral_class in ["O", "B"]:
                mass_sol = random.uniform(physical_constants.HOT_WHITE_DWARF_MIN_MASS_SOL, physical_constants.CHANDRASEKHAR_LIMIT_SOL)
            else:
                mass_sol = random.uniform(min_mass, physical_constants.COOL_WHITE_DWARF_MAX_MASS_SOL)
            if mass_override: # If mass is overridden for a WD, ensure it's within limits
                mass_sol = max(min(mass_override / physical_constants.SOLAR_MASS_TO_KG, physical_constants.CHANDRASEKHAR_LIMIT_SOL), min_mass)
        else:
            # For giants and supergiants, mass is less predictable from luminosity alone.
            # We choose a random mass within the physically allowed range for the class.
            if mass_override:
                mass_sol = mass_override / physical_constants.SOLAR_MASS_TO_KG
            else:
                mass_sol = random.uniform(min_mass, max_mass)


        # Ensure the calculated mass is within the absolute physical bounds for its class.
        mass_sol = max(min(mass_sol, max_mass), min_mass)
        mass = mass_sol * physical_constants.SOLAR_MASS_TO_KG

        # 6. Calculate Radius based on the star's type.
        if self.yerkes_class in ["VII", "D"]:
            # White dwarf radius follows an inverse mass-radius relationship.
            # R ∝ M^(-1/3). A 1 solar mass WD is the base.
            radius = physical_constants.WHITE_DWARF_BASE_RADIUS_KM * (mass_sol ** physical_constants.WHITE_DWARF_MASS_RADIUS_EXPONENT)
        else:
            # For all other stars, radius is calculated from luminosity and temperature
            # using the Stefan-Boltzmann law.
            luminosity_watts = luminosity * physical_constants.SOLAR_LUMINOSITY
            radius = math.sqrt(luminosity_watts / (physical_constants.FOUR_PI * physical_constants.STEFAN_BOLTZMANN_CONSTANT * temperature ** 4)) / physical_constants.KM_TO_M_FACTOR

        # 7. Set final star properties.
        color_descriptions = physical_constants.SPECTRAL_CLASS_COLORS
        star_type_str = f"{spectral_class}{subclass}{self.yerkes_class} {color_descriptions[spectral_class]} {yerkes_type} Star"

        self.type = star_type_str
        self.radius = radius
        self.mass = mass
        self.temperature = temperature
        self.luminosity = luminosity * physical_constants.SOLAR_LUMINOSITY