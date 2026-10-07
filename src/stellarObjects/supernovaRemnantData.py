# stellarObjects/supernovaRemnantData.py

"""
Supernova Remnant Generation
============================

This module contains the `SupernovaRemnant` class, a standalone descriptive
phenomenon representing the expanding wreckage of a past supernova,
encountered independently of any particular star system. Like `Nebula`, it
is generated on demand by `phenomenonGen.py`'s separate exotic-phenomenon
mode, never by `StarSystem._generate_planets`'s normal per-slot rolls.

A core-collapse remnant (as opposed to a Type Ia remnant, which leaves
nothing behind -- the progenitor white dwarf is thermonuclearly disrupted
entirely) may still contain the collapsed stellar core that caused the
explosion, embedded via `compactRemnant.BlackHole`/`NeutronStar`.
"""

import math
import random

from .compactRemnant import BlackHole, NeutronStar
from .config import SystemConfig
from planetgen.names.wordlists import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from .nebulaData import REMNANT_CLASS_LETTERS, choose_weighted_class, draw_class_contents, typical_class_contents
from planetgen import tuning
from planetgen.util import log
from planetgen.util.serialization import fields_from_dict, fields_to_dict
from .utils import (format_distance_ly, format_galactic_orbit, format_number, generate_galactic_orbit_fields,
                    generate_phoneme_salad_name)


def remnant_classes_for(progenitor_type, compact_remnant_kind):
    """
    The remnant classes (R-W) a remnant with this progenitor and compact
    core can be: a Type Ia remnant is always W; a core-collapse one is any
    class whose `compact` rule allows its core (T and U need a pulsar).
    """
    if progenitor_type == "Type Ia":
        return ("W",)
    return tuple(
        letter for letter in REMNANT_CLASS_LETTERS
        if letter != "W" and compact_remnant_kind in tuning.NEBULA_CLASSES[letter]["compact"]
    )


def infer_remnant_class(morphology, progenitor_type, age_years):
    """
    The class a remnant saved before classes existed (schema v38) most
    likely belongs to, from its stored shape, progenitor and age.
    """
    if progenitor_type == "Type Ia":
        return "W"
    if morphology == "plerion":
        return "T"
    if morphology == "composite":
        return "U"
    for letter in ("R", "V"):
        low, high = tuning.NEBULA_CLASSES[letter]["age_range_years"]
        if low <= age_years <= high:
            return letter
    return "S"


class SupernovaRemnant:
    """
    A basic class to store information for a supernova remnant.

    Attributes:
        name (str): A generated or explicitly given name for the remnant.
        remnant_class (str): A remnant letter class R-W
            (`tuning.NEBULA_CLASSES`).
        morphology (str): The class's shape, one of
            `tuning.SUPERNOVA_REMNANT_MORPHOLOGIES`
            (`"shell"`, `"plerion"`, or `"composite"`).
        dominant_species (str), density_cm3 (float), temperature_k
            (float), extinction_av (float): What fills the remnant, the
            same contents `Nebula` carries.
        age_years (float): Time since the supernova, in years.
        radius_ly (float): The remnant's current radius, in light-years,
            derived from `age_years` via the Sedov-Taylor blast-wave
            relation (see `tuning.SEDOV_TAYLOR_*`).
        progenitor_type (str): `"Type Ia"` or `"core-collapse"`.
        compact_remnant (BlackHole, NeutronStar, or None): The collapsed
            stellar core left behind, if any -- always `None` for a Type
            Ia progenitor (nothing is left behind; the white dwarf is
            thermonuclearly disrupted entirely), and only sometimes
            present/detectable for a core-collapse one.
    """

    SERIALIZABLE_FIELDS = [
        "name", "remnant_class", "morphology", "age_years", "radius_ly", "progenitor_type",
        "dominant_species", "density_cm3", "temperature_k", "extinction_av", "compact_offset_ly",
        "galactic_orbital_speed_kms", "galactic_orbital_period_gy",
        "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
    ]
    """Every attribute set by `__init__` except `system_config` (a shared
    back-reference) and `compact_remnant` (a nested object, handled
    separately in `to_dict`/`from_dict` alongside its own
    `compact_remnant_kind` discriminator -- mirrors how `StarSystem.to_dict`
    nests `secondary_star`/`wide_binary` outside its own field list). See
    `nebulaData.Nebula`'s identical `galactic_orbital_*` fields' docstring
    -- a supernova remnant is likewise still gravitationally part of the
    galaxy despite being bound to no star."""

    def __init__(self, system_config: SystemConfig, name=None):
        """
        Initializes a SupernovaRemnant object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object
                (used for its `MARKDOWN` flag, and threaded into any
                embedded `compact_remnant`).
            name (str, optional): An explicit name. Random if omitted.
        """
        self.system_config = system_config
        # A draft name: `_db.insert_supernova_remnant` reserves it through
        # the system-name registry (v40), and the core follows.
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)
        self.name_given = bool(name)  # a given name is kept over an object ID (GEN.64)

        is_type_ia = random.random() < tuning.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE
        self.progenitor_type = "Type Ia" if is_type_ia else "core-collapse"
        log.choice("Supernova progenitor type", self.progenitor_type,
                   f"roll against SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE "
                   f"({tuning.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE})")

        (self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy,
         self.galactic_orbital_phase_deg, self.galactic_min_update_interval_years) = \
            generate_galactic_orbit_fields()

        # Core-collapse remnants keep a neutron star (most) or black hole,
        # offset from the center by its birth kick (compact_offset_ly
        # below); thermonuclear (Type Ia, class W) remnants have none.
        self.compact_remnant = None
        if not is_type_ia and random.random() < tuning.SUPERNOVA_CORE_COLLAPSE_REMNANT_VISIBLE_CHANCE:
            if random.random() < tuning.SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE:
                self.compact_remnant = BlackHole(system_config, name=f"{self.name} Core")
                log.choice("Core-collapse compact remnant", "black hole",
                           f"roll passed SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE "
                           f"({tuning.SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE})")
            else:
                self.compact_remnant = NeutronStar(system_config, name=f"{self.name} Core")
                log.choice("Core-collapse compact remnant", "neutron star",
                           f"roll failed SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE "
                           f"({tuning.SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE})")
        elif not is_type_ia:
            log.debug("Core-collapse compact remnant: none visible/detectable "
                       "(roll failed SUPERNOVA_CORE_COLLAPSE_REMNANT_VISIBLE_CHANCE)")

        kind = None
        if isinstance(self.compact_remnant, BlackHole):
            kind = "black_hole"
        elif isinstance(self.compact_remnant, NeutronStar):
            kind = "neutron_star"
        self.remnant_class = choose_weighted_class(
            remnant_classes_for(self.progenitor_type, kind), "Supernova remnant class")
        class_data = tuning.NEBULA_CLASSES[self.remnant_class]
        self.morphology = class_data["morphology"]
        self.age_years = random.uniform(*class_data["age_range_years"])
        self.radius_ly = (
            tuning.SEDOV_TAYLOR_RADIUS_COEFFICIENT_LY
            * (self.age_years ** tuning.SEDOV_TAYLOR_TIME_EXPONENT)
        )
        (self.dominant_species, self.density_cm3, self.temperature_k,
         self.extinction_av) = draw_class_contents(self.remnant_class)
        self.compact_offset_ly = self._draw_kick_offset_ly(kind)

    def _draw_kick_offset_ly(self, kind):
        """
        Where the compact core has drifted since the explosion, relative to
        the remnant's center, light-years: its birth kick
        (`SUPERNOVA_KICK_SPEED_RANGE_KMS`, log-uniform) times the age, in a
        random direction. `None` when there is no core.
        """
        if kind is None:
            return None
        low, high = tuning.SUPERNOVA_KICK_SPEED_RANGE_KMS[kind]
        speed_kms = math.exp(random.uniform(math.log(low), math.log(high)))
        distance_ly = speed_kms / 299792.458 * self.age_years
        z = random.uniform(-1.0, 1.0)
        phi = random.uniform(0.0, 2 * math.pi)
        r = math.sqrt(1 - z * z)
        return [distance_ly * r * math.cos(phi), distance_ly * r * math.sin(phi), distance_ly * z]

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this remnant's properties, per
        `SERIALIZABLE_FIELDS`, plus a nested `compact_remnant`
        (`BlackHole`/`NeutronStar`'s own `to_dict()`, or `None`) and its
        `compact_remnant_kind` discriminator (`"black_hole"`,
        `"neutron_star"`, or `None`).

        Returns:
            dict: One entry per field in `SERIALIZABLE_FIELDS`, plus
                 `compact_remnant_kind`/`compact_remnant`.
        """
        data = fields_to_dict(self, self.SERIALIZABLE_FIELDS)
        if isinstance(self.compact_remnant, BlackHole):
            data["compact_remnant_kind"] = "black_hole"
            data["compact_remnant"] = self.compact_remnant.to_dict()
        elif isinstance(self.compact_remnant, NeutronStar):
            data["compact_remnant_kind"] = "neutron_star"
            data["compact_remnant"] = self.compact_remnant.to_dict()
        else:
            data["compact_remnant_kind"] = None
            data["compact_remnant"] = None
        return data

    @classmethod
    def from_dict(cls, data, system_config):
        """
        Reconstructs a `SupernovaRemnant` from a dict in the shape
        `to_dict()` produces, without re-running generation (`__init__` is
        bypassed via `object.__new__`).

        Args:
            data (dict): A dict in the shape `to_dict()` produces.
            system_config (SystemConfig): The shared config.

        Returns:
            SupernovaRemnant: The reconstructed remnant.
        """
        remnant = object.__new__(cls)
        remnant.system_config = system_config
        remnant.compact_offset_ly = None
        data = dict(data)
        if data.get("remnant_class") is None:
            # Saved before classes existed (schema v38).
            data["remnant_class"] = infer_remnant_class(
                data["morphology"], data["progenitor_type"], data["age_years"])
            (data["dominant_species"], data["density_cm3"], data["temperature_k"],
             data["extinction_av"]) = typical_class_contents(data["remnant_class"])
        fields_from_dict(remnant, data, cls.SERIALIZABLE_FIELDS)

        kind = data.get("compact_remnant_kind")
        if kind == "black_hole":
            remnant.compact_remnant = BlackHole.from_dict(data["compact_remnant"], system_config)
        elif kind == "neutron_star":
            remnant.compact_remnant = NeutronStar.from_dict(data["compact_remnant"], system_config)
        else:
            remnant.compact_remnant = None

        return remnant

    def to_paragraph_list(self):
        """
        Generates a list of descriptive paragraphs for the remnant: a
        header naming it, a description of its morphology/age/size/
        progenitor, and (if present) its embedded compact remnant's own
        full paragraph list.

        Returns:
            list: A list of strings, where each string is a paragraph
                  describing the supernova remnant.
        """
        header_level = '##' if self.system_config.MARKDOWN else '=='
        header = f"{header_level} {self.name} (Class {self.remnant_class} Supernova Remnant) {header_level if not self.system_config.MARKDOWN else ''}".rstrip()

        description = (
            f"{self.name} is a {self.morphology} supernova remnant, the expanding wreckage of a "
            f"{self.progenitor_type} supernova approximately {format_number(self.age_years, ',.0f')} years ago. It now reaches "
            f"about {format_distance_ly(self.radius_ly)} from its center, and still orbits the galactic center at "
            f"{format_galactic_orbit(self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy)}."
        )

        class_name = tuning.NEBULA_CLASSES[self.remnant_class]["name"]
        paragraphs = [header, description, (
            f"It is a class {self.remnant_class} remnant ({class_name.lower()}): mostly "
            f"{self.dominant_species}, at about {format_number(self.density_cm3, ',.3g')} particles per cubic centimeter "
            f"and {format_number(self.temperature_k, ',.3g')} K."
        )]

        if self.compact_remnant is not None:
            kind_label = "black hole" if isinstance(self.compact_remnant, BlackHole) else "neutron star"
            paragraphs.append(
                f"At its center lies the collapsed core of the progenitor star: a {kind_label}, "
                f"{self.compact_remnant.name}."
            )
            paragraphs.extend(self.compact_remnant.to_paragraph_list())
        elif self.progenitor_type == "core-collapse":
            paragraphs.append(
                "Whatever compact remnant the collapse left behind, if any, has not been located within the debris."
            )

        return paragraphs

    def __str__(self):
        """
        Returns a string representation of the supernova remnant, with
        paragraphs separated by double newlines.
        """
        return "\n\n".join(self.to_paragraph_list())
