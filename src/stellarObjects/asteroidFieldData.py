# stellarObjects/asteroidFieldData.py

"""
Asteroid Field Generation
=========================

This module contains the `AsteroidField` class, a standalone descriptive
phenomenon representing a field of asteroid debris drifting in open
interstellar space, encountered independently of any particular star
system -- unlike `asteroidData.AsteroidBelt`, which always orbits a star
within a `StarSystem`. Physically these are the same kind of object
(density + mineral composition), just without a host star/orbit, so this
class reuses `asteroidData.generate_asteroid_composition`/
`format_composition_summary` rather than duplicating that logic.

Like `nebulaData.Nebula` and the other standalone exotic phenomena, an
`AsteroidField` is generated on demand by `phenomenonGen.py`'s separate,
rarer exotic-phenomenon mode, never by `StarSystem._generate_planets`'s
normal per-slot rolls.
"""

import random

from .asteroidData import format_composition_summary, generate_asteroid_composition
from .config import SystemConfig
from .names import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from . import program_constants
from .serialization import fields_from_dict, fields_to_dict
from .utils import format_galactic_orbit, generate_galactic_orbit_fields, generate_phoneme_salad_name, reseed_rng


class AsteroidField:
    """
    A basic class to store information for a standalone field of asteroid
    debris drifting in open space, bound to no star.

    Attributes:
        name (str): A generated or explicitly given name for the field.
        density (str): The density of the field ('dense', 'sparse', 'typical'
            -- the same three levels `AsteroidBelt.density` uses).
        composition (list): A list of (component, concentration) tuples,
            generated the same way `AsteroidBelt.composition` is.
        radius_ly (float): The field's approximate radius, in light-years
            (`program_constants.ASTEROID_FIELD_RADIUS_RANGE_LY`).
        galactic_orbital_speed_kms (float): Circular orbital speed around
            the galactic center, km/s -- an asteroid field is still
            gravitationally part of the galaxy even though it isn't bound
            to any star (see `utils.generate_galactic_orbit_fields`, and
            `nebulaData.Nebula`'s identical fields' docstring).
        galactic_orbital_period_gy (float): Orbital period, billions of years.
        galactic_orbital_phase_deg (float): Current angular position around
            that orbit -- advanced over time by `_db.advance_orbital_phases`.
        galactic_min_update_interval_years (float): Floating-point update
            guard for `galactic_orbital_phase_deg`.
    """

    SERIALIZABLE_FIELDS = [
        "name", "density", "radius_ly",
        "galactic_orbital_speed_kms", "galactic_orbital_period_gy",
        "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
    ]
    """Every attribute set by `__init__`, excluding `system_config` (a
    shared back-reference) and `composition` (handled separately in
    `to_dict`/`from_dict`, since it's a list of tuples -- JSON has no
    tuple type -- mirroring `AsteroidBelt.SERIALIZABLE_FIELDS`'s identical
    exclusion)."""

    def __init__(self, system_config: SystemConfig, name=None):
        """
        Initializes an AsteroidField object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object
                (used only for its `MARKDOWN` flag, via `to_paragraph_list`).
            name (str, optional): An explicit name. Random if omitted.
        """
        reseed_rng()
        self.system_config = system_config
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)

        self.density = random.choice(["dense", "sparse", "typical"])
        self.composition = generate_asteroid_composition()
        self.radius_ly = random.uniform(*program_constants.ASTEROID_FIELD_RADIUS_RANGE_LY)

        (self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy,
         self.galactic_orbital_phase_deg, self.galactic_min_update_interval_years) = \
            generate_galactic_orbit_fields()

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this field's properties, per
        `SERIALIZABLE_FIELDS`, with `composition` stored as a list of
        2-element lists (JSON has no tuple type) -- mirrors
        `AsteroidBelt.to_dict`.

        Returns:
            dict: One entry per field in `SERIALIZABLE_FIELDS`, plus
                 `composition`.
        """
        data = fields_to_dict(self, self.SERIALIZABLE_FIELDS)
        data["composition"] = [list(pair) for pair in self.composition]
        return data

    @classmethod
    def from_dict(cls, data, system_config):
        """
        Reconstructs an `AsteroidField` from a dict in the shape
        `to_dict()` produces, without re-running generation (`__init__` is
        bypassed via `object.__new__`).

        Args:
            data (dict): A dict in the shape `to_dict()` produces.
            system_config (SystemConfig): The shared config.

        Returns:
            AsteroidField: The reconstructed field.
        """
        field = object.__new__(cls)
        field.system_config = system_config
        fields_from_dict(field, data, cls.SERIALIZABLE_FIELDS)
        field.composition = [tuple(pair) for pair in data["composition"]]
        return field

    def get_composition_summary(self):
        """
        Builds the human-readable composition summary phrase, the same
        as-published-text role `AsteroidBelt.get_composition_summary`
        plays for the database persistence layer.

        Returns:
            str: The composition summary phrase.
        """
        return format_composition_summary(self.composition)

    def to_paragraph_list(self):
        """
        Generates a list of descriptive paragraphs for the field: a header
        naming it, and a description of its size, density, composition,
        and galactic motion.

        Returns:
            list: A list of strings, where each string is a paragraph
                  describing the asteroid field.
        """
        header_level = '##' if self.system_config.MARKDOWN else '=='
        header = f"{header_level} {self.name} (Asteroid Field) {header_level if not self.system_config.MARKDOWN else ''}".rstrip()

        description = (
            f"{self.name} is a {self.density} field of asteroid debris drifting in open interstellar space, "
            f"spanning roughly {self.radius_ly:.3f} light-years across, composed of {self.get_composition_summary()}. "
            f"Though bound to no star, it still orbits the galactic center at "
            f"{format_galactic_orbit(self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy)}."
        )

        return [header, description]

    def __str__(self):
        """
        Returns a string representation of the asteroid field, with
        paragraphs separated by double newlines.
        """
        return "\n\n".join(self.to_paragraph_list())
