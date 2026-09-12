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

import random

from .compactRemnant import BlackHole, NeutronStar
from .config import SystemConfig
from .names import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from . import program_constants
from .serialization import fields_from_dict, fields_to_dict
from .utils import (format_galactic_orbit, generate_galactic_orbit_fields,
                    generate_phoneme_salad_name, reseed_rng)


class SupernovaRemnant:
    """
    A basic class to store information for a supernova remnant.

    Attributes:
        name (str): A generated or explicitly given name for the remnant.
        morphology (str): One of `program_constants.SUPERNOVA_REMNANT_MORPHOLOGIES`
            (`"shell"`, `"plerion"`, or `"composite"`).
        age_years (float): Time since the supernova, in years.
        radius_ly (float): The remnant's current radius, in light-years,
            derived from `age_years` via the Sedov-Taylor blast-wave
            relation (see `program_constants.SEDOV_TAYLOR_*`).
        progenitor_type (str): `"Type Ia"` or `"core-collapse"`.
        compact_remnant (BlackHole, NeutronStar, or None): The collapsed
            stellar core left behind, if any -- always `None` for a Type
            Ia progenitor (nothing is left behind; the white dwarf is
            thermonuclearly disrupted entirely), and only sometimes
            present/detectable for a core-collapse one.
    """

    SERIALIZABLE_FIELDS = [
        "name", "morphology", "age_years", "radius_ly", "progenitor_type",
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
        reseed_rng()
        self.system_config = system_config
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)

        self.morphology = random.choice(program_constants.SUPERNOVA_REMNANT_MORPHOLOGIES)
        self.age_years = random.uniform(*program_constants.SUPERNOVA_REMNANT_AGE_RANGE_YEARS)
        self.radius_ly = (
            program_constants.SEDOV_TAYLOR_RADIUS_COEFFICIENT_LY
            * (self.age_years ** program_constants.SEDOV_TAYLOR_TIME_EXPONENT)
        )

        is_type_ia = random.random() < program_constants.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE
        self.progenitor_type = "Type Ia" if is_type_ia else "core-collapse"

        (self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy,
         self.galactic_orbital_phase_deg, self.galactic_min_update_interval_years) = \
            generate_galactic_orbit_fields()

        self.compact_remnant = None
        if not is_type_ia and random.random() < program_constants.SUPERNOVA_CORE_COLLAPSE_REMNANT_VISIBLE_CHANCE:
            if random.random() < program_constants.SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE:
                self.compact_remnant = BlackHole(system_config, name=f"{self.name} Core")
            else:
                self.compact_remnant = NeutronStar(system_config, name=f"{self.name} Core")

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
        header = f"{header_level} {self.name} (Supernova Remnant) {header_level if not self.system_config.MARKDOWN else ''}".rstrip()

        description = (
            f"{self.name} is a {self.morphology} supernova remnant, the expanding wreckage of a "
            f"{self.progenitor_type} supernova approximately {self.age_years:,.0f} years ago. It now spans "
            f"roughly {self.radius_ly:.2f} light-years across, and still orbits the galactic center at "
            f"{format_galactic_orbit(self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy)}."
        )

        paragraphs = [header, description]

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
