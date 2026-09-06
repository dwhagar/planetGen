# stellarObjects/config.py

"""
Configuration settings for the stellarObjects package.

This module defines the `SystemConfig` class, which encapsulates all
configuration flags and variables that control the behavior of the star system
generation process. These settings can be modified via command-line arguments
or a JSON system-specification file in `systemGen.py`, or directly within the
code for specific scenarios.

Most boolean-ish options are tri-state (`True`/`False`/`None`) rather than a
pair of `FORCE_X`/`NO_X` flags:
    - `True` forces the feature to be present.
    - `False` forces the feature to be absent.
    - `None` (the default) leaves it to random chance.

This tri-state shape maps directly onto the command line's `+name`/`-name`
syntax (see `systemGen.process_args`): `+name` sets the option to `True`,
`-name` sets it to `False`, and omitting it leaves it `None`.
"""

from .serialization import fields_from_dict, fields_to_dict

# The SystemConfig attributes that make up a system's generation "recipe" --
# i.e. everything a `--system-file` JSON document can set (see
# `systemGen.apply_system_file`) -- as opposed to `system_flavor_count` /
# `recent_flavor_texts`, which are run-time bookkeeping rather than
# configuration. `to_dict`/`from_dict` use this list so a config round-trips
# through JSON in exactly the shape `--system-file` already expects.
SERIALIZABLE_FIELDS = [
    "MARKDOWN", "HABITABLE_WORLD", "ASTEROID_BELT", "LARGE_STAR", "MOONS",
    "MAX_PLANETS", "PLANETS", "STAR_TYPE", "NAME", "AGE", "INTELLIGENT_LIFE",
    "BINARY_SYSTEM", "NUM_ORBITS", "SLOTS",
]


class SystemConfig:
    def __init__(self):
        self.system_flavor_count = 0
        self.recent_flavor_texts = [] # Stores recently used flavor texts to avoid repetition
        self.MARKDOWN = False
        """
        bool: If True, output will be formatted in Markdown. If False, output will be
        formatted in Wikitext. Defaults to False.
        """

        self.HABITABLE_WORLD = None
        """
        bool or None: If True, the system generation will attempt to force the
        creation of at least one habitable world. If False, ensures no
        habitable world is generated. If None, left to random chance.
        Defaults to None.
        """

        self.ASTEROID_BELT = None
        """
        bool or None: If True, the system generation will attempt to force the
        creation of at least one asteroid belt. If False, ensures no asteroid
        belt is generated. If None, left to random chance. Defaults to None.
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
        bool or None: If True, the system generation will attempt to force the
        creation of at least one planet with intelligent life (implies
        HABITABLE_WORLD). If False, ensures no planet with intelligent life is
        generated (also implies HABITABLE_WORLD). If None, left to random
        chance. Defaults to None.
        """

        self.BINARY_SYSTEM = None
        """
        bool or None: If True, the system generation will attempt to create a
        binary star system. If False (or None), a single star system is
        generated. Defaults to None.
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

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this config's settings, in the
        same lowercase-key shape as `systemGen.py`'s `--system-file` format.

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
