# planetgen/generation/star.py

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
import re

from planetgen.generation import forced_star, helium_core_draw
from planetgen.generation.config import SystemConfig
from planetgen.physics import constants
from planetgen.physics import (activity, age_metallicity, debris_disk, formation, imf, metallicity, pre_main_sequence,
                              spin, stellar_state)
from planetgen.util import draw
from planetgen.util import log
from planetgen.physics.stellar_evolution import (YERKES_CLASS_NAMES, birth_fe_h, evolve_star, sample_living_star,
                                                  sample_star_age_gy, star_params)
from planetgen.galaxy import scales
from planetgen.generation.stellar_object import StellarObject
from planetgen.physics.orbits import calculate_habitable_zone
from planetgen.util.checks import finite_domain
from planetgen.physics import unit_registry as units

_percm3_to_perm3 = units.converter("cm-3", "m-3")
_km_to_m = units.converter("km", "m")
_m_to_au = units.converter("m", "au")
_m_to_km = units.converter("m", "km")

STAR_TYPE_PATTERN = re.compile(r"([OBAFGKM])([0-9])(IA\+|IAB|VII|III|IA|IB|II|IV|VI|0|V|D)")
"""A forced `SystemConfig.STAR_TYPE` (uppercased): spectral class, subclass
digit, Yerkes class. Always matched with `fullmatch` -- a prefix match
accepted 'G2Vjunk' as G2V and parsed 'G10V' as a G1 hypergiant ('0')."""


def cloud_pressure_pa(density_cm3, temperature_k):
    """
    The pressure a nebula or supernova remnant's gas puts on a star's
    heliopause: ram pressure from the star moving through it at
    `constants.STAR_CLOUD_RELATIVE_SPEED_MS` plus thermal pressure.

    Args:
        density_cm3 (float): The cloud's hydrogen density (nH), per cm^3.
        temperature_k (float): The cloud's gas temperature in kelvin.

    Returns:
        float: The pressure in pascals.
    """
    n_m3 = _percm3_to_perm3(max(density_cm3 or 0.0, 0.0))
    rho = n_m3 * constants.ISM_MASS_PER_HYDROGEN * constants.HYDROGEN_ATOM_MASS_KG
    ram = rho * constants.STAR_CLOUD_RELATIVE_SPEED_MS ** 2
    thermal = n_m3 * constants.BOLTZMANN * max(temperature_k or 0.0, 0.0)
    return ram + thermal


def compressed_heliosphere_radius(radius_au, density_cm3, temperature_k):
    """
    A heliopause radius (computed against the open interstellar medium,
    `constants.ISM_PRESSURE`) pressed in by the cloud the star sits
    inside. The heliopause sits where the wind's momentum flux balances the
    outside pressure, so R scales as P^-1/2: a dense cold cloud (nH ~3,000
    cm^-3) shrinks the Sun's from ~85 AU to well under 1 AU. A cloud thinner
    than the open medium leaves it as it is.

    Args:
        radius_au (float): The open-space heliopause radius, AU.
        density_cm3 (float): The cloud's hydrogen density (nH), per cm^3.
        temperature_k (float): The cloud's gas temperature in kelvin.

    Returns:
        float: The compressed radius, AU.
    """
    pressure = max(cloud_pressure_pa(density_cm3, temperature_k), constants.ISM_PRESSURE)
    return radius_au * math.sqrt(constants.ISM_PRESSURE / pressure)


class Star(StellarObject):
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

    @staticmethod
    @finite_domain()
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
        radius_m = _km_to_m(radius_km)
        lum_sol = luminosity / constants.SOLAR_LUMINOSITY
        mass_sol = _kg_to_msun(mass)
        radius_sol = radius_m / constants.SOLAR_RADIUS_M
        spectral_class = star_type[0].upper() if star_type else None

        # --- 2. Calculate Mass-Loss Rate (M-dot) and Wind Velocity (v_inf) ---
        # A single formula for mass loss is insufficient. We use a tiered system based on
        # the star's Yerkes luminosity class (evolutionary stage) and spectral type.
        escape_velocity = math.sqrt(constants.ESCAPE_VELOCITY_CONSTANT * constants.G * mass / radius_m)


        # TIER 1: Evolved stars -- giants (III), subgiants (IV), bright giants
        # (II), supergiants (IB/IAB/IA/IA+), and hypergiants (0). These were
        # previously three separately hand-fit tiers (a hypergiant power-law,
        # plus Reimers' Law with two different eta values for "giants" vs
        # "supergiants") stitched together at Yerkes-class boundaries -- the
        # exact structure that let the old hypergiant constant drift ~9
        # orders of magnitude out of calibration unnoticed. Nieuwenhuijzen &
        # de Jager (1990) is a single empirical mass-loss fit spanning this
        # whole regime, so one formula and one continuous curve now covers
        # all of it (see constants.py for the correction applied
        # above its known high-luminosity overestimation threshold).
        if yerkes_class in ["0", "IA+", "IA", "IAB", "IB", "II", "III", "IV"]:
            mass_loss_rate_smyr = (constants.NDJ_MASS_LOSS_COEFFICIENT
                                   * (lum_sol**constants.NDJ_LUMINOSITY_EXPONENT)
                                   * (mass_sol**constants.NDJ_MASS_EXPONENT)
                                   * (radius_sol**constants.NDJ_RADIUS_EXPONENT))
            if lum_sol > constants.NDJ_HIGH_LUMINOSITY_THRESHOLD_LSUN:
                mass_loss_rate_smyr *= constants.NDJ_HIGH_LUMINOSITY_CORRECTION_FACTOR

            if yerkes_class == "0":
                # Hypergiant winds are an extreme multiple of escape velocity.
                wind_velocity = constants.HYPERGIANT_WIND_VELOCITY_FACTOR * escape_velocity
            else:
                # Wind velocity is a smaller fraction of escape velocity for these cooler giants/supergiants.
                wind_velocity = constants.GIANT_WIND_VELOCITY_FACTOR * escape_velocity

        # TIER 2a: Hot main-sequence dwarfs (O and B, Class V/VI)
        # Unlike cool dwarfs, O/B stars have powerful radiation-driven winds (line-driven,
        # via the same physical mechanism as hypergiants/supergiants, just far less extreme):
        # mass loss scales steeply with luminosity rather than the coronal/magnetic-activity
        # scaling used for cool dwarfs below. Without this tier, a hot O-type star would be
        # (incorrectly) modeled with a weaker wind than a cool, dim M dwarf.
        elif spectral_class in ("O", "B") and yerkes_class in ["V", "VI"]:
            mass_loss_rate_smyr = constants.OB_DWARF_MASS_LOSS_RATE_FACTOR * (lum_sol**constants.OB_DWARF_MASS_LOSS_RATE_EXPONENT)
            wind_velocity = constants.OB_DWARF_WIND_VELOCITY_FACTOR * escape_velocity

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
            scaling_factor = (radius_sol**constants.RADIUS_SOL_EXPONENT) * (lum_sol**constants.LUMINOSITY_SOL_EXPONENT) * (mass_sol**constants.MASS_SOL_EXPONENT)
            mass_loss_rate_smyr = constants.SUN_MASS_LOSS_RATE_SOLAR_MASS_PER_YEAR * scaling_factor # 2e-14 is the Sun's M-dot in M_sol/yr

            # Scale wind velocity based on the star's escape velocity relative to the Sun's.
            wind_velocity = constants.SOLAR_WIND_VELOCITY * (escape_velocity / constants.SOLAR_ESCAPE_VELOCITY)

        # --- 3. Convert and Calculate Final Radius ---

        # Convert M-dot from (solar masses/year) to (kg/s)
        mass_loss_rate_kgs = _msun_to_kg(mass_loss_rate_smyr) / constants.SECONDS_PER_YEAR

        # The heliopause radius is where the stellar wind's momentum flux balances the ISM pressure.
        # R = sqrt( (M-dot * v_inf) / (4 * pi * P_ism) )
        momentum_flux = max(mass_loss_rate_kgs * wind_velocity, constants.MIN_MOMENTUM_FLUX) # Ensure momentum_flux is not zero or negative

        try:
            heliopause_radius_m = math.sqrt(momentum_flux / (constants.FOUR_PI * constants.ISM_PRESSURE))
        except ValueError:
            # This should not happen with the floor value, but as a final safety net.
            heliopause_radius_m = constants.HELIOPAUSE_RADIUS_DEFAULT_M

        # Convert the final radius from meters to Astronomical Units (AU) for output.
        return _m_to_au(heliopause_radius_m)

    def calculate_heliosphere(self):
        """
        Estimates the radius of the star's heliosphere (astrosphere) by
        delegating to the `_calculate_heliosphere_radius_static` helper with
        this star's own mass, luminosity, radius, type, and Yerkes class.

        Returns:
            float: The estimated radius of the heliosphere in Astronomical Units (AU).
        """
        return Star._calculate_heliosphere_radius_static(self.mass, self.luminosity, self.radius, self.type, self.yerkes_class)

    def __init__(self, system_config: SystemConfig, name=None, **kwargs):
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

        Keyword Args:
            galactic_center_dist_ly (float, optional): This star's actual
                distance from the galactic center, in light-years, threaded
                into `calculate_system_perimeter` -- see that method's
                docstring. `None` (the default) uses the radius that holds half of the galaxy's disc stars (`galaxy.scales`).
            galactic_orbital_phase_deg (float, optional): This star's
                current angular position around its galactic orbit, in
                degrees -- the value `planetgen.db.store.advance_orbital_phases`
                advances over time, the same role `orbital_phase_deg` plays
                for a planet/moon's own orbit (see
                `physics.planets.generate_orbital_motion_properties`).
                `StarSystem.__init__` rolls this once per system and passes
                the *same* value to the primary star, the secondary star
                (if any), and the `BinaryStarProxy` -- a binary pair's
                negligible internal separation (AU) next to its galactic
                orbit radius (light-years) means both stars move around the
                galaxy together, not independently, so they must share one
                phase rather than each rolling their own. `None` (the
                default) rolls a fresh random value in `[0, 360)` -- the
                correct behavior for a standalone `Star` with no system
                threading this through.
        """
        super().__init__(system_config, name=name, galactic_center_dist_ly=kwargs.get('galactic_center_dist_ly'))

        self.luminosity = None
        self.temperature = None
        self.mass = None
        self.radius = None
        self.type = None
        self.yerkes_class = None
        self.habitable_zone = None
        self.system_perimeter = None
        self.heliosphere_radius = None
        self.galactic_orbital_speed_kms = None
        self.galactic_orbital_period_gy = None
        self.galactic_orbital_phase_deg = None
        self.galactic_min_update_interval_years = None
        self.generate_star(initial_mass_sol=kwargs.get('initial_mass_sol'),
                           age_gy=kwargs.get('age_gy'), fe_h=kwargs.get('metallicity_fe_h'))
        self._finish_init(kwargs)

    def _finish_init(self, kwargs):
        """The derived properties every generated star gets once its core
        ones exist: age and lifespan (fixed by the model with the rest),
        habitable zone, Hill sphere, heliosphere and galactic
        orbit."""
        self.age, self.lifespan = self._model_age_lifespan
        self.habitable_zone = calculate_habitable_zone(self.luminosity, self.temperature)
        self.system_perimeter = self.calculate_system_perimeter(self.galactic_center_dist_ly)
        self.heliosphere_radius = self.calculate_heliosphere()
        self.set_galactic_motion(self.galactic_center_dist_ly, kwargs.get('galactic_orbital_phase_deg'))
        # GEN.104: the spin, last so nothing above draws differently.
        # GEN.846: a stored star keeps the period its row has, so the row and the star agree.
        self.rotation_period_hours = self._model_rotation_hours
        if self.rotation_period_hours is None:
            self.rotation_period_hours = spin.star_rotation_period_hours(
                self.mass, self.radius, self.temperature, self.age, self.yerkes_class)
        spin.set_spin(self, spin.GALACTIC_POLE, spin.rayleigh_tilt_deg())
        # GEN.86: coronal X-rays, XUV and flares, from mass and age.
        activity.generate_activity(self)
        # GEN.158: [Fe/H], from its own stream so no older draw moves.
        self.metallicity_fe_h = kwargs.get('metallicity_fe_h')
        if self.metallicity_fe_h is None:
            # A horizontal-branch star's [Fe/H] was drawn with its colour (GEN.286).
            self.metallicity_fe_h = self._model_fe_h
        if self.metallicity_fe_h is None:
            self.metallicity_fe_h = metallicity.sample_fe_h(
                getattr(self.system_config, "POPULATION", None),
                scales.default_galactocentric_ly() if self.galactic_center_dist_ly is None
                else self.galactic_center_dist_ly,
                draw.Stream(f"metallicity:{self.mass!r}:{self.age!r}:{self.galactic_center_dist_ly!r}"),
                peaks=getattr(self.system_config, "FE_H_PEAKS", None))
        # GEN.786 and GEN.785: its evolutionary state and, for a white dwarf, the age since it stopped burning (no draw).
        stellar_state.set_state(self)
        # GEN.458: this star's own snow line, drawn once around the blackbody position and stored; its stream is its
        # own, like the [Fe/H] above, so no older draw moves.
        self.snow_line_au = formation.draw_snow_line_au(
            self.luminosity, draw.Stream(f"snow_line:{self.mass!r}:{self.age!r}:{self.galactic_center_dist_ly!r}"))
        # GEN.242: the disk of a young star, from its own stream like the [Fe/H] above.
        if self.yerkes_class in pre_main_sequence.STAGE_CLASSES:
            disk = pre_main_sequence.draw_disk(
                _kg_to_msun(self.mass), self.age * 1.0e9, self.luminosity / constants.SOLAR_LUMINOSITY,
                draw.Stream(f"disk:{self.mass!r}:{self.age!r}"))
            for name, value in (disk or {}).items():
                setattr(self, name, value)
        elif self.yerkes_class in debris_disk.DEBRIS_CLASSES:
            # GEN.787: the cold belt of a main-sequence star, from a stream of its own like the disk above.
            disk = debris_disk.draw_disk(_kg_to_msun(self.mass), self.luminosity / constants.SOLAR_LUMINOSITY,
                                         draw.Stream(f"debris:{self.mass!r}:{self.age!r}"))
            for name, value in (disk or {}).items():
                setattr(self, name, value)

    def _generate_from_population_model(self, initial_mass_sol=None, age_gy=None, fe_h=None):
        """
        Sets this star's properties from the population model
        (`stellarEvolution`): a Kroupa-IMF mass and a star-formation-history
        age (or the given ones), evolved to its present state. The letter,
        subclass and Yerkes class follow from the resulting temperature and
        luminosity; the radius from Stefan-Boltzmann (a white dwarf's from
        its mass-radius relation). Also fixes `age`/`lifespan` (read by
        `__init__`),
        `initial_mass_sol` (a binary secondary's mass is a fraction of it)
        and `phase_end_age_gy` (`adjust_age_for_planets` keeps the age
        inside the present phase).

        `+large_star` truncates the IMF at `LARGE_STAR_MIN_MASS_SOL`, and a
        primary that has already collapsed to a neutron star or black hole
        is redrawn (see `stellarEvolution.sample_living_star`). The config's
        `POPULATION` sets the age range and `MAX_STAR_LUMINOSITY_SOL` redraws
        any star at least that bright.
        """
        age_bias = self.system_config.AGE
        population = getattr(self.system_config, "POPULATION", None)
        fe_h_peaks = getattr(self.system_config, "FE_H_PEAKS", None)
        age_law = getattr(self.system_config, "AGE_LAW", None)
        state = None
        if initial_mass_sol is not None:
            initial_mass_sol = max(initial_mass_sol, imf.lowest_mass_sol())
            if fe_h is None:
                fe_h = birth_fe_h(draw, population, fe_h_peaks, self.galactic_center_dist_ly)
            age = age_gy if age_gy is not None else age_metallicity.tilted_age(
                fe_h, population, lambda: sample_star_age_gy(age_bias), draw)
            state = evolve_star(initial_mass_sol, age, fe_h=fe_h, population=population, fe_h_peaks=fe_h_peaks,
                                galactic_center_dist_ly=self.galactic_center_dist_ly)
            # A companion of a specified-type primary takes a fraction of
            # that star's present mass, which can be heavy enough to have
            # collapsed already at the pair's age. Lighten it until it is
            # alive rather than redraw it at another age: the pair shares
            # one age (GEN.53).
            while state is None and initial_mass_sol > imf.lowest_mass_sol():
                draw.redrawn("star.mass_lightened")
                initial_mass_sol = max(initial_mass_sol * 0.8, imf.lowest_mass_sol())
                state = evolve_star(initial_mass_sol, age, fe_h=fe_h, population=population, fe_h_peaks=fe_h_peaks,
                                    galactic_center_dist_ly=self.galactic_center_dist_ly)
        if state is None:
            # No companion given.
            habitable_host = self.system_config.HABITABLE_WORLD is True or self.system_config.INTELLIGENT_LIFE is True
            max_luminosity_sol = getattr(self.system_config, "MAX_STAR_LUMINOSITY_SOL", None)
            max_mass_sol = getattr(self.system_config, "MAX_STAR_MASS_SOL", None)
            try:
                initial_mass_sol, age, state = sample_living_star(
                    age_bias, self.system_config.LARGE_STAR, habitable_host=habitable_host,
                    population=population, max_luminosity_sol=max_luminosity_sol, max_mass_sol=max_mass_sol,
                    fe_h_peaks=fe_h_peaks, age_law=age_law, galactic_center_dist_ly=self.galactic_center_dist_ly,
                    fe_h=fe_h)
            except ValueError:
                if not (self.system_config.LARGE_STAR and (max_luminosity_sol is not None
                                                           or max_mass_sol is not None)):
                    raise
                # A large star a run's prevalence asked for (GEN.52) can be
                # brighter, or heavier, than a sector whose bright and
                # massive stars were pre-placed allows; this one stays an
                # ordinary star.
                log.choice("Large star", "dropped", "none within the sector's star limits")
                initial_mass_sol, age, state = sample_living_star(
                    age_bias, None, habitable_host=habitable_host,
                    population=population, max_luminosity_sol=max_luminosity_sol, max_mass_sol=max_mass_sol,
                    fe_h_peaks=fe_h_peaks, age_law=age_law, galactic_center_dist_ly=self.galactic_center_dist_ly,
                    fe_h=fe_h)

        state = helium_core_draw.with_helium_core(state, initial_mass_sol, age)
        self._apply_params(star_params(initial_mass_sol, age, state))
        log.choice("Star", self.type, f"population model: initial mass {initial_mass_sol:.3g} Msun, "
                                      f"age {age:.3g} Gyr")

    def _apply_params(self, params):
        """Sets this star's generated properties from a `star_params` dict
        (see `stellarEvolution.star_params`), rolling nothing."""
        self.type = params["type"]
        self.yerkes_class = params["yerkes_class"]
        self.radius = params["radius_km"]
        self.mass = params["mass_kg"]
        self.temperature = params["temperature_k"]
        self.luminosity = params["luminosity_w"]
        self.initial_mass_sol = params["initial_mass_sol"]
        self.phase_end_age_gy = params["phase_end_age_gy"]
        self.nuclear_arrival_age_gy = params.get("nuclear_arrival_age_gy")
        self._model_age_lifespan = (params["age_gy"], params["lifespan_gy"])
        self._model_fe_h = params.get("metallicity_fe_h")
        self._model_rotation_hours = params.get("rotation_period_hours")

    @classmethod
    def from_params(cls, params, system_config, name=None, **kwargs):
        """
        Builds a real `Star` from a `stellarEvolution.star_params` dict (a
        pre-placed bright star read back from storage, or a draw from
        `star_population`) without re-rolling any of its generated
        properties: only the derived ones (habitable zone, Hill sphere,
        heliosphere, galactic orbit) are computed, exactly as `__init__`
        computes them. Accepts `__init__`'s `galactic_center_dist_ly` and
        `galactic_orbital_phase_deg` keywords.
        """
        # GEN.158: a stored [Fe/H] is kept; none (or an old row) leaves the star to draw its own.
        if params.get("metallicity_fe_h") is not None:
            kwargs.setdefault("metallicity_fe_h", params["metallicity_fe_h"])
        star = cls.__new__(cls)
        StellarObject.__init__(star, system_config, name=name,
                               galactic_center_dist_ly=kwargs.get("galactic_center_dist_ly"))
        star._apply_params(params)
        star._finish_init(kwargs)
        return star

    @classmethod
    def from_record(cls, record, system_config, name=None, **kwargs):
        """
        Builds the `Star` of a member record (`member_record.MemberRecord`, GEN.1009) without evolving it: the state the tree drew the
        member from is the state the star has, so its mass, radius and temperature are the record's. Takes `from_params`'s keywords.

        Raises:
            ValueError: If the record is not a star (a brown dwarf or a compact remnant has no `star_params`).
        """
        params = record.star_params()
        if params is None:
            raise ValueError(f"a {record.kind} member has no star to build")
        return cls.from_params(params, system_config, name=name, **kwargs)

    def generate_star(self, initial_mass_sol=None, age_gy=None, fe_h=None):
        """
        Generates the physical properties of the star from the population
        model in `stellarEvolution.py`, either at random (see
        `_generate_from_population_model`) or, for a specified
        `config.STAR_TYPE`, as a star of that type the same model makes (see
        `_generate_forced_type`). `initial_mass_sol`/`age_gy` pin the random
        model's mass and age (a binary secondary: a fraction of the primary's
        initial mass, at the primary's age).
        """
        if initial_mass_sol is not None or not self.system_config.STAR_TYPE:
            self._generate_from_population_model(initial_mass_sol, age_gy, fe_h)
            return
        self._generate_forced_type(fe_h)

    def _generate_forced_type(self, fe_h=None):
        """
        Sets this star's properties from a forced `config.STAR_TYPE` (GEN.941), e.g. 'G2V': a star of that spectral
        letter, subclass and Yerkes class drawn from the population model (`forced_star.sample`), so its mass,
        luminosity, age and lifespan are the model's. A type the model never makes (a B giant) is the model's nearest
        star of that Yerkes class by temperature, at the temperature the type stands for, and says so in the log.
        `IA+` and `D` are the model's `0` and `VII` under their own names.
        """
        match = STAR_TYPE_PATTERN.fullmatch(self.system_config.STAR_TYPE.upper())
        if not match:
            raise ValueError("Invalid star type format. Expected format is e.g., G2V.")
        spectral_class, subclass_str, yerkes_class = match.groups()
        subclass = int(subclass_str)
        found = forced_star.sample(
            spectral_class, subclass, yerkes_class, age_bias=self.system_config.AGE, rng=draw,
            population=getattr(self.system_config, "POPULATION", None),
            fe_h_peaks=getattr(self.system_config, "FE_H_PEAKS", None),
            galactic_center_dist_ly=self.galactic_center_dist_ly, fe_h=fe_h)
        state = dict(found.state)
        if not found.exact:
            state["temperature_k"] = forced_star.target_temperature_k(spectral_class, subclass)
        params = star_params(found.initial_mass_sol, found.age_gy, state)
        model_class = params["yerkes_class"]
        if not found.exact or yerkes_class != model_class:
            name = YERKES_CLASS_NAMES[yerkes_class]
            color = constants.SPECTRAL_CLASS_COLORS[spectral_class]
            params["type"] = f"{spectral_class}{subclass}{yerkes_class} {color} {name} Star"
        params["yerkes_class"] = yerkes_class
        self._apply_params(params)
        if found.exact:
            log.choice("Star", self.type, f"forced type: model star of initial mass {found.initial_mass_sol:.3g} Msun, "
                                          f"age {found.age_gy:.3g} Gyr")
        else:
            log.choice("Star", self.type, f"forced type the model never makes: its nearest {model_class} star, "
                                          f"initial mass {found.initial_mass_sol:.3g} Msun, at "
                                          f"{state['temperature_k']:.0f} K")


def adjust_pair_age_for_planets(primary, secondary, primary_planets, secondary_planets, members=()):
    """
    Settles one shared age for a binary pair's two stars (GEN.53): old
    enough for every planet each star hosts (judged by that star's own
    evolutionary scales), and no older than either star's generated state
    allows. The primary's `adjust_age_for_planets` does the work, bounded
    by the companion's `age_ceiling_gy`; the companion then takes the
    result. Both stars start at the same age (the companion is born with
    the primary), so neither moves outside what it was generated as.

    Args:
        primary (Star): The heavier star.
        secondary (Star): Its companion.
        primary_planets (list): Bodies judged against the primary's scales.
        secondary_planets (list): Bodies judged against the companion's (a
            close pair passes its circumbinary planets for both).
        members (list of tuple): `(star, planets)` for each further star of the system (GEN.332): they share the age
            too, so their planets set a minimum and their states a ceiling like the companion's.
    """
    min_required = max([secondary.planet_age_requirement_gy(secondary_planets),
                        *(star.planet_age_requirement_gy(bodies) for star, bodies in members)])
    ceiling = min([secondary.age_ceiling_gy(), *(star.age_ceiling_gy() for star, _ in members)])
    primary.adjust_age_for_planets(primary_planets, min_required_age_gy=min_required, age_ceiling_gy=ceiling)
    secondary.age = primary.age
    for star, _ in members:
        star.age = primary.age


from planetgen.generation import object_timing  # noqa: E402 -- registers the generators below (PERF.92)

_kg_to_msun = units.converter("kg", "Msun")
_msun_to_kg = units.converter("Msun", "kg")
object_timing.track(Star, "star")
