# planetgen/generation/config.py

"""
Configuration settings for the generator.

This module defines the `SystemConfig` class, which encapsulates all
configuration flags and variables that control the behavior of the star system
generation process. These settings can be modified via command-line arguments
or a JSON system-specification file in `generation/run_system.py`, or directly within the
code for specific scenarios.

Most boolean-ish options are tri-state (`True`/`False`/`None`) rather than a
pair of `FORCE_X`/`NO_X` flags:
    - `True` forces the feature to be present.
    - `False` forces the feature to be absent.
    - `None` (the default) leaves it to random chance.

This tri-state shape maps directly onto the command line's `+name`/`-name`
syntax (see `run_system.process_args`): `+name` sets the option to `True`,
`-name` sets it to `False`, and omitting it leaves it `None`.
"""

from planetgen import tuning
from planetgen.util.serialization import fields_from_dict, fields_to_dict

# The SystemConfig attributes that make up a system's generation "recipe" --
# i.e. everything a `--system-file` JSON document can set (see
# `run_system.apply_system_file`) -- as opposed to `system_flavor_count` /
# `recent_flavor_texts`, which are run-time bookkeeping rather than
# configuration. `to_dict`/`from_dict` use this list so a config round-trips
# through JSON in exactly the shape `--system-file` already expects.
SERIALIZABLE_FIELDS = [
    "HABITABLE_WORLD", "ASTEROID_BELT", "COMETS", "LARGE_STAR", "MOONS",
    "MAX_PLANETS", "PLANETS", "STAR_TYPE", "NAME", "AGE", "INTELLIGENT_LIFE",
    "BINARY_SYSTEM", "WIDE_BINARY", "STARS", "MAX_STARS", "NUM_ORBITS", "SLOTS", "PREVALENCE",
]


class SystemConfig:
    FORCE_CLASS = None
    """
    str or None: A planet class letter (e.g. "O") the system must hold at least one planet of (GEN.791). A test and admin
    option for a sample that needs every class (`generation/class_sample.py`), not a recipe field: a class attribute, so a
    config holds it only when a caller sets it, it is not saved with the system and a default run never sets it. The first
    orbital slot is that class, placed in a zone that allows it (`StarSystem.calculate_distance_for_class`), and a system
    whose attempts all lost the planet reports "a planet of class X" in `StarSystem.unmet_requirements`. It cannot be
    combined with `SLOTS` (slots already name classes). Defaults to None.
    """

    def __init__(self):
        self.system_flavor_count = 0
        self.recent_flavor_texts = [] # Stores recently used flavor texts to avoid repetition
        self.TYPE_IA_SHARE = None
        """
        float or None: The share of supernova remnants that are Type Ia where the remnant is made (the sector
        profile's `type_ia_share`, GEN.319). None for the solar circle's, `tuning.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE`.
        """

        self.AMBIENT_DENSITY_CM3 = None
        """
        float or None: The density of the gas the sector's supernova remnants expand into where they sit in a
        superbubble (the cavity's `ambient_density_cm3`, GEN.298): their size law reads it. None outside a bubble:
        the remnant's own class density.
        """

        self.HABITABLE_WORLD = None
        """
        bool or None: If True, the system generation raises the odds of at
        least one habitable world (`util.request_odds`) and never forces
        one; the system may come out without and says so
        (`StarSystem.habitable_request`). If False, ensures no
        habitable world is generated. If None, left to random chance.
        Defaults to None.
        """

        self.ASTEROID_BELT = None
        """
        bool or None: If True, the system generation will attempt to force the
        creation of at least one asteroid belt. If False, ensures no asteroid
        belt is generated. If None, left to random chance. Defaults to None.
        """

        self.COMETS = None
        """
        bool or None: If True, the system generation will attempt to force the
        creation of at least one star-bound comet (see `cometData.Comet` --
        contrast the always-standalone `roguePlanetData.InterstellarComet`,
        which this flag has no effect on). If False, ensures no comet is
        generated. If None, left to random chance. Defaults to None.
        """

        self.LARGE_STAR = None
        """
        bool or None: If True, the system generation will force the creation of
        a larger, more massive star. If False, explicitly generates from the
        normal stellar population (the same as the default). If None, left to
        the default stellar population weighting. Defaults to None.
        """

        self.MOONS = None
        """
        bool or None: If True, planets generated will always have a chance at
        moons forced. If False, no planet in the system will have moons. If
        None, each planet has a 50/50 chance, unless a per-slot moon count is
        specified. Defaults to None.
        """

        self.MAX_PLANETS = None
        """
        bool or None: If True, the system generation will attempt to create the
        maximum number of planets the system can support. If False, generates
        the minimum number of objects the system's other requirements allow.
        If None, a random count is chosen. Defaults to None.
        """

        self.PLANETS = None
        """
        bool or None: If True, ensures the system has at least one planet or
        asteroid belt. If False, skips the planet generation process
        entirely, creating only a star. If None, the count (including zero)
        is left to random chance. Defaults to None.
        """

        self.STAR_TYPE = None
        """
        str or None: Specifies a precise star type to generate (e.g., 'G2V'). If None,
        a star type will be randomly generated. Defaults to None.
        """

        self.NAME = None
        """
        str or None: Specifies a name for the star system, overriding the default
        random generation. If None, a name will be randomly generated. Defaults to None.
        """

        self.AGE = None
        """
        str or None: Specifies the age of the star system, either "young" or "old".
        If None, the age will be randomly determined. Defaults to None.
        """

        self.INTELLIGENT_LIFE = None
        """
        bool or None: If True, the system generation raises the odds of at
        least one planet with intelligent life (implies HABITABLE_WORLD) and
        never forces a world, a timeline or a civilization. If False, ensures no
        planet with intelligent life is generated (also implies
        HABITABLE_WORLD). If None, left to random chance. Defaults to None.
        """

        self.BINARY_SYSTEM = None
        """
        bool or None: If True, the system generation will attempt to create a
        binary star system. If False, a single star system is generated. If
        None, the stars of the system are drawn as a tree against real
        stellar-multiplicity statistics for the primary star's own spectral
        type (the `star-count` table of `planetgen.generation.phenomena.hosts`,
        used by `StarSystem._drawn_plan`) rather than always coming out
        single. Defaults to None.
        """

        self.WIDE_BINARY = None
        """
        bool or None: Only consulted when BINARY_SYSTEM is True; ignored
        otherwise. Selects which of the two real binary-star configurations
        is generated. If True, forces an S-type (wide) binary -- two stars
        separated by tens to thousands of AU, each hosting its own
        independently-generated planets, with each star's maximum stable
        orbit limited by the companion's gravity (see
        `orbits.holman_wiegert_critical_semimajor_axis`). If False, forces a
        P-type (close/circumbinary) binary -- the original, and only,
        binary behavior this generator had before this option existed,
        where the two stars are merged into one effective star
        (`binary.BinaryStarProxy`) for the purposes of planet
        placement. If None, the pair is P-type when the tree puts its axis under
        `tuning.WIDE_BINARY_SEPARATION_MIN_AU` and S-type otherwise. Defaults to None.
        """

        self.STARS = None
        """
        int or None: Forces the number of stars in the system (`--stars N`, the system file's `stars`, GEN.267), from 1
        to MAX_STARS. 1 is a single star, 2 a pair, and more a pair with member stars (GEN.630, GEN.631). None leaves
        the count to the `star-count` table. Read by `StarSystem._drawn_plan`; carried in the config, the plan and the
        `system_configs` row so a stored config regenerates with the same count.
        `generation.run_system.stars_conflict` refuses a count that contradicts BINARY_SYSTEM or WIDE_BINARY.
        Defaults to None.
        """

        self.MAX_STARS = tuning.MAX_STARS_DEFAULT
        """
        int: The most stars the system may hold (`--max-stars N`, the system file's `max_stars`, GEN.267), the cap the
        star-count draw's tail stops at; `tuning.MAX_STARS_DEFAULT` (12) unless raised, up to
        `tuning.MAX_STARS_CEILING`. A forced STARS above it is refused. Defaults to 12.
        """

        self.NUM_ORBITS = None
        """
        int or None: Explicitly sets the number of orbital slots (planets and
        asteroid belts combined) to generate, overriding the normal
        mass-based random estimate. If None, the count is estimated/randomized
        as usual (subject to MAX_PLANETS/PLANETS). Defaults to None.
        """

        self.SLOTS = None
        """
        list or None: An optional, per-orbit list of specifications describing
        exactly what should occupy each orbital slot. Each entry is either
        `None` (meaning that slot is left to normal random generation) or a
        dict with the keys:
            - "type": "planet" or "asteroid_belt" (required).
            - "planet_class": a specific planet class letter (e.g. "M"),
              only used when type is "planet". Optional; if omitted, a class
              is chosen normally for the slot's zone.
            - "moons": an exact integer number of moons to generate for the
              planet (0 for none). Optional; if omitted, moon generation
              falls back to the MOONS setting/random chance.
        The list does not need to cover every slot; slots beyond the end of
        the list (or with a `None` entry) are generated normally. Defaults to
        None.
        """

        self.PREVALENCE = {}
        """
        dict: A sector or galaxy run's prevalences (GEN.52): feature name
        (`prevalence.FEATURES`, e.g. "comets") to a percentage deviation
        from its normal chance (+50 is 1.5 times as often, -100 never).
        Features left out keep their normal chance. The tri-state flags
        above force a feature; this only shifts how often it happens.
        "multi_star" scales the chance of five or more stars
        (`prevalence.multi_star_scale`, which `StarSystem._drawn_plan` passes as
        the `prevalence` of `multistar.draw_member_tree`).
        """

        self.POPULATION = None
        """
        str or None: The stellar population a random primary star is drawn
        from ("young", "intermediate", "old" or "bulge"; see
        `stellarEvolution.population_age_range_gy`), which sets its age
        range; a galaxy sector picks one per system from
        `galaxy.density.population_densities` at its position
        (`star_population.pick_population`). None draws from the whole
        disk's history. Not saved with the config: the star keeps its age.
        """

        self.AGE_LAW = None
        """
        galaxy.age_law.AgeLaw or None: How a random primary star's age is spread inside its population's age range
        (GEN.546): by the disk's star formation rate, from `galaxy.disk_history.population_age_laws`. None spreads it
        evenly. A galaxy sector sets it with `POPULATION`. Not saved with the config: the star keeps its age.
        """

        self.GALAXY_SHAPE = None
        """
        galaxy.density.GalaxyShape or None: The galaxy a body is made in, which gives the disc's star formation history the thorium to
        uranium ratio of its rock follows (GEN.759, `physics.rprocess_enrichment`). A galaxy sector sets it with `POPULATION`; None takes
        the default shape (`galaxy.shape_draw.default_shape`), as any work with no stored shape does. Not saved with the config.
        """

        self.FE_H_PEAKS = None
        """
        tuple or None: `((weight, mean, sigma), ...)`, the [Fe/H] a random primary star of a cluster or satellite is
        drawn from (GEN.237): the sector's own profile gives them, since those populations have no mean at a distance
        from the centre. None draws from the population's mean (`metallicity.mean_fe_h`). Not saved with the config.
        """

        self.PLACE_TERMS = None
        """
        place_terms.PlaceTerms or None: The density and velocity dispersion of the place a galaxy sector's systems are
        drawn at (GEN.103): the widest axis a wide pair keeps follows from them. None takes the Sun's neighbourhood
        (`place_terms.SOLAR_NEIGHBOURHOOD`). Not saved with the config.
        """

        self.MAX_STAR_LUMINOSITY_SOL = None
        """
        float or None: A random primary star is redrawn until it is dimmer
        than this (Lsun): a sector whose bright stars were pre-placed fills
        the rest with dimmer ones only. None for no limit. Not saved with
        the config.
        """

        self.MAX_STAR_MASS_SOL = None
        """
        float or None: A random primary star is redrawn until it was born
        lighter than this (solar masses): a sector of a galaxy whose
        scatter placed every star at or above its mass limit in advance
        draws the lighter ones only. None for no limit. Not saved with the
        config.
        """

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this config's settings, in the
        same lowercase-key shape as `generation/run_system.py`'s `--system-file` format.

        Returns:
            dict: One entry per field in `SERIALIZABLE_FIELDS`.
        """
        return fields_to_dict(self, SERIALIZABLE_FIELDS)

    @classmethod
    def from_dict(cls, data):
        """
        Builds a `SystemConfig` from a dict in the shape `to_dict()` (or a
        `--system-file` JSON document) produces. Keys not present in `data`
        are left at their normal defaults.

        Args:
            data (dict): A dict with zero or more of `SERIALIZABLE_FIELDS`'s
                         lowercase names as keys.

        Returns:
            SystemConfig: The newly built config.
        """
        config = cls()
        fields_from_dict(config, data, SERIALIZABLE_FIELDS)
        return config
