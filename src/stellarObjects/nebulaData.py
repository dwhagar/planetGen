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

import math
import random

from .config import SystemConfig
from .names import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from . import program_constants
from planetgen.util import log
from planetgen.util.serialization import fields_from_dict, fields_to_dict
from .utils import (format_distance_ly, format_galactic_orbit, format_number, generate_galactic_orbit_fields,
                    generate_phoneme_salad_name)


NEBULA_CLASS_LETTERS = tuple(
    letter for letter, data in program_constants.NEBULA_CLASSES.items()
    if data["family"] != "supernova-remnant"
)
"""tuple: The nebula classes (A-Q); R-W are supernova remnants."""

REMNANT_CLASS_LETTERS = tuple(
    letter for letter, data in program_constants.NEBULA_CLASSES.items()
    if data["family"] == "supernova-remnant"
)
"""tuple: The supernova remnant classes (R-W)."""


def _draw_in_range(low, high):
    """Log-uniform between `low` and `high`, or uniform when `low` is 0
    (a log scale can't start at nothing)."""
    if low <= 0:
        return random.uniform(low, high)
    return math.exp(random.uniform(math.log(low), math.log(high)))


def _mid_of_range(low, high):
    """The geometric middle of a range (arithmetic when `low` is 0) -- the
    typical value a migration fills in for a row saved before its class."""
    if low <= 0:
        return (low + high) / 2
    return math.sqrt(low * high)


def draw_class_contents(letter):
    """
    Draws the physical contents of one nebula or remnant class.

    Args:
        letter (str): A `program_constants.NEBULA_CLASSES` key.

    Returns:
        tuple: `(dominant_species, density_cm3, temperature_k,
            extinction_av)`.
    """
    data = program_constants.NEBULA_CLASSES[letter]
    return (
        data["species"],
        _draw_in_range(*data["density_range_cm3"]),
        _draw_in_range(*data["temperature_range_k"]),
        _draw_in_range(*data["extinction_range_av"]),
    )


def typical_class_contents(letter):
    """The same tuple as `draw_class_contents`, but each value the middle
    of its range -- what `_db`'s v38 migration fills for existing rows."""
    data = program_constants.NEBULA_CLASSES[letter]
    return (
        data["species"],
        _mid_of_range(*data["density_range_cm3"]),
        _mid_of_range(*data["temperature_range_k"]),
        _mid_of_range(*data["extinction_range_av"]),
    )


def choose_weighted_class(letters, label):
    """Picks one of `letters` by its `NEBULA_CLASSES` frequency, logging
    the choice under `label`."""
    weights = [program_constants.NEBULA_CLASSES[letter]["frequency"] for letter in letters]
    letter = random.choices(letters, weights=weights, k=1)[0]
    log.choice(label, letter, f"weighted draw among {list(letters)} by NEBULA_CLASSES frequency")
    return letter


def infer_nebula_class(nebula_type, radius_ly):
    """
    The class a nebula saved before classes existed (schema v38) most
    likely belongs to, from its family and size: the most common class of
    that family whose radius range holds `radius_ly`, else the family's
    most common class.
    """
    candidates = [letter for letter in NEBULA_CLASS_LETTERS
                  if program_constants.NEBULA_CLASSES[letter]["family"] == nebula_type]
    if not candidates:
        candidates = list(NEBULA_CLASS_LETTERS)
    by_frequency = sorted(candidates, key=lambda letter: -program_constants.NEBULA_CLASSES[letter]["frequency"])
    for letter in by_frequency:
        low, high = program_constants.NEBULA_CLASSES[letter]["radius_range_ly"]
        if low <= radius_ly <= high:
            return letter
    return by_frequency[0]


class Nebula:
    """
    A basic class to store information for an interstellar nebula.

    Modeled on `AsteroidBelt`'s "simple container of already-decided
    properties" shape, but standalone rather than tied to a star system --
    see module docstring.

    Attributes:
        name (str): A generated or explicitly given name for the nebula.
        nebula_class (str): A nebula letter class A-Q
            (`program_constants.NEBULA_CLASSES`).
        nebula_type (str): The class's family, one of
            `program_constants.NEBULA_FAMILIES`' keys (`"diffuse"`,
            `"emission"`, `"reflection"`, `"planetary"` or `"dark"`).
        radius_ly (float): The nebula's approximate radius, in light-years.
        composition (str): A descriptive composition string for the family.
        formation_cause (str): A descriptive formation-cause string for
            the family.
        dominant_species (str): What the cloud is mostly made of.
        density_cm3 (float): Particle density nH, per cm^3.
        temperature_k (float): Gas temperature, kelvin.
        extinction_av (float): Optical extinction through the cloud,
            magnitudes.
        galactic_orbital_speed_kms (float): Circular orbital speed around
            the galactic center, km/s -- a nebula is still gravitationally
            part of the galaxy even though it isn't bound to any star (see
            `utils.generate_galactic_orbit_fields`).
        galactic_orbital_period_gy (float): Orbital period, billions of
            years.
        galactic_orbital_phase_deg (float): Current angular position
            around that orbit -- the value `_db.advance_orbital_phases`
            advances over time, the same role `Star.galactic_orbital_phase_deg`
            plays.
        galactic_min_update_interval_years (float): Floating-point update
            guard for `galactic_orbital_phase_deg` (`utils.
            minimum_update_interval_years`).
    """

    SERIALIZABLE_FIELDS = [
        "name", "nebula_class", "nebula_type", "radius_ly", "composition", "formation_cause",
        "dominant_species", "density_cm3", "temperature_k", "extinction_av",
        "galactic_orbital_speed_kms", "galactic_orbital_period_gy",
        "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
    ]
    """Every attribute set by `__init__`, excluding `system_config` (a
    shared back-reference, threaded into `from_dict` rather than
    serialized redundantly)."""

    def __init__(self, system_config: SystemConfig, name=None, nebula_type=None, nebula_class=None):
        """
        Initializes a Nebula object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object
                (used only for its `MARKDOWN` flag, via `to_paragraph_list`).
            name (str, optional): An explicit name. Random if omitted.
            nebula_type (str, optional): A `NEBULA_FAMILIES` key; the
                class is then drawn within that family (sector generation
                asks for `"planetary"`, the one family with a point rate).
            nebula_class (str, optional): An explicit class letter A-Q.
                Drawn by frequency if omitted.
        """
        self.system_config = system_config
        # A draft name: `_db.insert_nebula` reserves it through the
        # system-name registry (v40), which may decorate it.
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)
        self.name_given = bool(name)  # a given name is kept over an object ID (GEN.64)

        if nebula_class is not None:
            if nebula_class not in NEBULA_CLASS_LETTERS:
                raise ValueError(f"nebula_class must be one of {list(NEBULA_CLASS_LETTERS)}, got {nebula_class!r}")
            self.nebula_class = nebula_class
        else:
            if nebula_type is not None and nebula_type not in program_constants.NEBULA_FAMILIES:
                raise ValueError(f"nebula_type must be one of {list(program_constants.NEBULA_FAMILIES)}, got {nebula_type!r}")
            letters = tuple(letter for letter in NEBULA_CLASS_LETTERS
                            if nebula_type is None
                            or program_constants.NEBULA_CLASSES[letter]["family"] == nebula_type)
            self.nebula_class = choose_weighted_class(letters, "Nebula class")
        class_data = program_constants.NEBULA_CLASSES[self.nebula_class]
        self.nebula_type = class_data["family"]
        family_data = program_constants.NEBULA_FAMILIES[self.nebula_type]
        self.radius_ly = random.uniform(*class_data["radius_range_ly"])
        self.composition = family_data["composition"]
        self.formation_cause = family_data["formation_cause"]
        (self.dominant_species, self.density_cm3, self.temperature_k,
         self.extinction_av) = draw_class_contents(self.nebula_class)

        (self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy,
         self.galactic_orbital_phase_deg, self.galactic_min_update_interval_years) = \
            generate_galactic_orbit_fields()

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
        data = dict(data)
        if data.get("nebula_class") is None:
            # Saved before classes existed (schema v38): infer one, and
            # give it its class's typical contents.
            data["nebula_class"] = infer_nebula_class(data["nebula_type"], data["radius_ly"])
            (data["dominant_species"], data["density_cm3"], data["temperature_k"],
             data["extinction_av"]) = typical_class_contents(data["nebula_class"])
        fields_from_dict(nebula, data, cls.SERIALIZABLE_FIELDS)
        return nebula

    @property
    def class_name(self):
        """The class's own name, e.g. "Classical H II region"."""
        return program_constants.NEBULA_CLASSES[self.nebula_class]["name"]

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
        header = f"{header_level} {self.name} (Class {self.nebula_class} {self.nebula_type.capitalize()} Nebula) {header_level if not self.system_config.MARKDOWN else ''}".rstrip()

        article = "an" if self.nebula_type[0] in "aeiou" else "a"
        description = (
            f"{self.name} is {article} {self.nebula_type} nebula about "
            f"{format_distance_ly(self.radius_ly)} in radius, composed of {self.composition}. It formed via {self.formation_cause}. "
            f"Though bound to no star, it still orbits the galactic center at "
            f"{format_galactic_orbit(self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy)}."
        )
        contents = (
            f"It is a class {self.nebula_class} nebula ({self.class_name.lower()}): mostly "
            f"{self.dominant_species}, at about {format_number(self.density_cm3, ',.3g')} particles per cubic centimeter "
            f"and {format_number(self.temperature_k, ',.0f')} K, dimming the stars behind it by "
            f"{self.extinction_av:.2g} magnitudes."
        )

        return [header, description, contents]

    def __str__(self):
        """
        Returns a string representation of the nebula, with paragraphs
        separated by double newlines.
        """
        return "\n\n".join(self.to_paragraph_list())
