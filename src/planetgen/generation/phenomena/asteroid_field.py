# planetgen/generation/phenomena/asteroid_field.py

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

import math
import random

from planetgen.generation.belt import format_composition_summary, generate_asteroid_composition
from planetgen.generation.config import SystemConfig
from planetgen.names.wordlists import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from planetgen import tuning
from planetgen.util import log
from planetgen.util.serialization import fields_from_dict, fields_to_dict
from planetgen.galaxy.galactic_orbit import format_galactic_orbit, generate_galactic_orbit_fields
from planetgen.names.wordsalad import generate_phoneme_salad_name
from planetgen.util.format import format_distance_ly


AU_PER_LY = 63241.077
"""float: Astronomical units in one light-year, for the size digit."""


def asteroid_field_size_digit(radius_ly):
    """The size digit of a field's class: floor(log10(radius in AU)),
    clamped to `ASTEROID_FIELD_SIZE_DIGIT_RANGE` (1-4)."""
    low, high = tuning.ASTEROID_FIELD_SIZE_DIGIT_RANGE
    return max(low, min(high, math.floor(math.log10(radius_ly * AU_PER_LY))))


def asteroid_field_class(composition_family, density, radius_ly):
    """
    A field's full class, letter plus size digit (e.g. `"C3"`).

    Args:
        composition_family (str): An `ASTEROID_FIELD_COMPOSITIONS` key.
        density (str): `"sparse"`, `"typical"` or `"dense"`.
        radius_ly (float): The field's radius, light-years.
    """
    letter = tuning.ASTEROID_FIELD_COMPOSITIONS[composition_family]["letters"][density]
    return f"{letter}{asteroid_field_size_digit(radius_ly)}"


def composition_family_for_letter(letter):
    """The composition family a class letter belongs to."""
    for family, data in tuning.ASTEROID_FIELD_COMPOSITIONS.items():
        if letter in data["letters"].values():
            return family
    raise ValueError(f"unknown asteroid field class letter {letter!r}")


def _family_composition(composition_family):
    """A (component, concentration) list like
    `generate_asteroid_composition`'s, sampled from the family's own
    components (or all asteroid components for a mixed family)."""
    components = tuning.ASTEROID_FIELD_COMPOSITIONS[composition_family]["components"]
    if components is None:
        return generate_asteroid_composition()
    concentrations = ["high", "moderate", "small", "trace"]
    chosen = random.sample(components, k=min(len(concentrations), len(components)))
    return list(zip(chosen, concentrations))


def asteroid_field_designation(field_class, sector_code, index):
    """
    An asteroid field's designation (v40): `AF <class>-<sector>-<nn>`,
    e.g. `AF C3-FE81000A2B-01`, where `<sector>` is the sector it was
    found in (its grid designation, or its name) and `<nn>` counts that
    sector's fields from 01. `AF <class>-<nn>` for one found in no
    sector.
    """
    if sector_code:
        return f"AF {field_class}-{sector_code}-{index:02d}"
    return f"AF {field_class}-{index:02d}"


class AsteroidField:
    """
    A basic class to store information for a standalone field of asteroid
    debris drifting in open space, bound to no star.

    Attributes:
        name (str): A generated or explicitly given name for the field.
        field_class (str): The class, a letter from composition and
            density plus a size digit (e.g. `"C3"`;
            `tuning.ASTEROID_FIELD_COMPOSITIONS`).
        composition_family (str): The family the letter comes from
            (`"carbonaceous"`, `"stony"`, `"metallic"`, ...).
        density (str): The density of the field ('dense', 'sparse', 'typical'
            -- the same three levels `AsteroidBelt.density` uses).
        composition (list): A list of (component, concentration) tuples,
            generated the same way `AsteroidBelt.composition` is.
        radius_ly (float): The field's approximate radius, in light-years
            (`tuning.ASTEROID_FIELD_RADIUS_RANGE_LY`).
        galactic_orbital_speed_kms (float): Circular orbital speed around
            the galactic center, km/s -- an asteroid field is still
            gravitationally part of the galaxy even though it isn't bound
            to any star (see `galactic_orbit.generate_galactic_orbit_fields`, and
            `nebulaData.Nebula`'s identical fields' docstring).
        galactic_orbital_period_gy (float): Orbital period, billions of years.
        galactic_orbital_phase_deg (float): Current angular position around
            that orbit -- advanced over time by `_db.advance_orbital_phases`.
        galactic_min_update_interval_years (float): Floating-point update
            guard for `galactic_orbital_phase_deg`.
    """

    SERIALIZABLE_FIELDS = [
        "name", "field_class", "composition_family", "density", "radius_ly",
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
        self.system_config = system_config
        # A placeholder: `_db.insert_asteroid_field` replaces it with the
        # field's designation (`asteroid_field_designation`, v40).
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)
        self.name_given = bool(name)  # a given name is kept over an object ID (GEN.64)

        self.density = random.choice(["dense", "sparse", "typical"])
        log.choice("Asteroid field density", self.density, "uniform draw among dense/sparse/typical")
        families = list(tuning.ASTEROID_FIELD_COMPOSITIONS)
        self.composition_family = random.choices(
            families,
            weights=[tuning.ASTEROID_FIELD_COMPOSITIONS[f]["frequency"] for f in families],
            k=1,
        )[0]
        log.choice("Asteroid field composition", self.composition_family,
                   "weighted draw by ASTEROID_FIELD_COMPOSITIONS frequency")
        self.composition = _family_composition(self.composition_family)
        self.radius_ly = random.uniform(*tuning.ASTEROID_FIELD_RADIUS_RANGE_LY)
        self.field_class = asteroid_field_class(self.composition_family, self.density, self.radius_ly)

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
        data = dict(data)
        if data.get("field_class") is None:
            # Saved before classes existed (schema v38): a random mineral
            # mix is honestly "mixed".
            data["composition_family"] = "mixed"
            data["field_class"] = asteroid_field_class("mixed", data["density"], data["radius_ly"])
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
        header = f"{header_level} {self.name} (Class {self.field_class} Asteroid Field) {header_level if not self.system_config.MARKDOWN else ''}".rstrip()

        family_description = tuning.ASTEROID_FIELD_COMPOSITIONS[self.composition_family]["description"]
        description = (
            f"{self.name} is a class {self.field_class} field: a {self.density} {self.composition_family} "
            f"field of asteroid debris ({family_description}) drifting in open interstellar space, "
            f"about {format_distance_ly(self.radius_ly)} in radius, composed of {self.get_composition_summary()}. "
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
