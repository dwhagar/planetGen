# stellarObjects/nebulaData.py

"""
Nebula Generation
=================

This module contains the `Nebula` class, a standalone descriptive phenomenon
representing an interstellar gas/dust cloud encountered independently of any
particular star system. Unlike `Planet`/`AsteroidBelt`, a `Nebula` never
orbits a star and is never placed into `StarSystem.planets` -- it is
generated on demand by `phenomenonGen.py`'s separate, rarer exotic-
phenomenon mode (see that module's docstring), not by
`StarSystem._generate_planets`'s normal per-slot rolls.
"""

import random

from .config import SystemConfig
from .names import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from . import program_constants
from .serialization import fields_from_dict, fields_to_dict
from .utils import generate_phoneme_salad_name, reseed_rng


class Nebula:
    """
    A basic class to store information for an interstellar nebula.

    Modeled on `AsteroidBelt`'s "simple container of already-decided
    properties" shape, but standalone rather than tied to a star system --
    see module docstring.

    Attributes:
        name (str): A generated or explicitly given name for the nebula.
        nebula_type (str): One of `program_constants.NEBULA_TYPES`'s keys
            (`"emission"`, `"reflection"`, `"planetary"`, or `"dark"`).
        radius_ly (float): The nebula's approximate radius, in light-years.
        composition (str): A descriptive composition string for this type.
        formation_cause (str): A descriptive formation-cause string for
            this type.
    """

    SERIALIZABLE_FIELDS = ["name", "nebula_type", "radius_ly", "composition", "formation_cause"]
    """Every attribute set by `__init__`, excluding `system_config` (a
    shared back-reference, threaded into `from_dict` rather than
    serialized redundantly)."""

    def __init__(self, system_config: SystemConfig, name=None):
        """
        Initializes a Nebula object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object
                (used only for its `MARKDOWN` flag, via `to_paragraph_list`).
            name (str, optional): An explicit name. Random if omitted.
        """
        reseed_rng()
        self.system_config = system_config
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)

        self.nebula_type = random.choice(list(program_constants.NEBULA_TYPES.keys()))
        type_data = program_constants.NEBULA_TYPES[self.nebula_type]
        self.radius_ly = random.uniform(*type_data["radius_range_ly"])
        self.composition = type_data["composition"]
        self.formation_cause = type_data["formation_cause"]

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this nebula's properties, per
        `SERIALIZABLE_FIELDS`.

        Returns:
            dict: One entry per field in `SERIALIZABLE_FIELDS`.
        """
        return fields_to_dict(self, self.SERIALIZABLE_FIELDS)

    @classmethod
    def from_dict(cls, data, system_config):
        """
        Reconstructs a `Nebula` from a dict in the shape `to_dict()`
        produces, without re-running generation (`__init__` is bypassed
        via `object.__new__`).

        Args:
            data (dict): A dict in the shape `to_dict()` produces.
            system_config (SystemConfig): The shared config.

        Returns:
            Nebula: The reconstructed nebula.
        """
        nebula = object.__new__(cls)
        nebula.system_config = system_config
        fields_from_dict(nebula, data, cls.SERIALIZABLE_FIELDS)
        return nebula

    def to_paragraph_list(self):
        """
        Generates a list of descriptive paragraphs for the nebula: a
        header naming it and its type, then a description of its size,
        composition, and formation.

        Returns:
            list: A list of strings, where each string is a paragraph
                  describing the nebula.
        """
        header_level = '##' if self.system_config.MARKDOWN else '=='
        header = f"{header_level} {self.name} ({self.nebula_type.capitalize()} Nebula) {header_level if not self.system_config.MARKDOWN else ''}".rstrip()

        article = "an" if self.nebula_type[0] in "aeiou" else "a"
        description = (
            f"{self.name} is {article} {self.nebula_type} nebula spanning roughly {self.radius_ly:.1f} "
            f"light-years across, composed of {self.composition}. It formed via {self.formation_cause}."
        )

        return [header, description]

    def __str__(self):
        """
        Returns a string representation of the nebula, with paragraphs
        separated by double newlines.
        """
        return "\n\n".join(self.to_paragraph_list())
