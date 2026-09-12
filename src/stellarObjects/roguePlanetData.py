# stellarObjects/roguePlanetData.py

"""
Rogue Planets & Interstellar Comets
====================================

This module contains `RoguePlanet` and `InterstellarComet`, two standalone
descriptive phenomena representing unbound objects encountered while
surveying a sector, rather than orbiting any star. Both are generated on
demand by `phenomenonGen.py`'s separate exotic-phenomenon mode, never by
`StarSystem._generate_planets`'s normal per-slot rolls.

Kept in one module since they share a generation trigger ("an unbound
wanderer encountered in interstellar space") and a CLI sub-choice, despite
being physically distinct: a `RoguePlanet` is a planetary-mass body ejected
from (or never bound to) a star system, while an `InterstellarComet` is a
much smaller icy planetesimal passing through on a hyperbolic trajectory
(real confirmed examples: 1I/'Oumuamua, 2I/Borisov).
"""

import math
import random

from .config import SystemConfig
from .names import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from . import physical_constants, program_constants
from .serialization import fields_from_dict, fields_to_dict
from .utils import (format_galactic_orbit, generate_galactic_orbit_fields,
                    generate_phoneme_salad_name, reseed_rng)


class RoguePlanet:
    """
    A basic class to store information for a free-floating ("rogue"/nomad)
    planet with no host star.

    Attributes:
        name (str): A generated or explicitly given name.
        planet_type (str): `'t'` (terrestrial/icy) or `'g'` (gas giant),
            the same letters `Planet.body_type` uses, chosen by mass
            relative to `program_constants.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER`.
        mass_kg (float): Mass in kilograms.
        radius_km (float): Radius in kilometers.
        composition (str): A descriptive bulk-composition string.
        has_internal_heat (bool): Whether it retains detectable internal
            heat (radiogenic/primordial) worth describing.
        has_moons (bool): Whether it retains a captured companion moon.
    """

    SERIALIZABLE_FIELDS = [
        "name", "planet_type", "mass_kg", "radius_km", "composition",
        "has_internal_heat", "has_moons",
        "galactic_orbital_speed_kms", "galactic_orbital_period_gy",
        "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
    ]
    """Every attribute set by `__init__`, excluding `system_config` (a
    shared back-reference). The `galactic_orbital_*` fields (see
    `nebulaData.Nebula`'s identical fields' docstring) reflect that a
    rogue planet is unbound from any specific STAR, not from the galaxy
    itself -- it still orbits the galactic center like a lone star does."""

    def __init__(self, system_config: SystemConfig, name=None):
        """
        Initializes a RoguePlanet object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object
                (used only for its `MARKDOWN` flag, via `to_paragraph_list`).
            name (str, optional): An explicit name. Random if omitted.
        """
        reseed_rng()
        self.system_config = system_config
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)

        mass_jupiter = random.uniform(*program_constants.ROGUE_PLANET_MASS_RANGE_JUPITER)
        self.mass_kg = mass_jupiter * physical_constants.JUPITER_MASS_TO_KG

        if mass_jupiter >= program_constants.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER:
            self.planet_type = 'g'
            # Real gas giants show a near-flat mass-radius relation from
            # roughly Saturn's mass up through the deuterium-burning limit
            # -- electron-degeneracy pressure counteracts the added weight
            # of extra mass almost exactly (Chabrier & Baraffe 2000, ARA&A
            # 38:337) -- so radius is drawn as a narrow variation around
            # Jupiter's own radius rather than scaled with mass.
            self.radius_km = physical_constants.JUPITER_RADIUS_KM * random.uniform(0.8, 1.15)
            self.composition = "hydrogen and helium, similar in bulk composition to Jupiter or Saturn"
        else:
            self.planet_type = 't'
            density_range_gcm3 = physical_constants.PLANET_DENSITY["t"]
            density_kg_m3 = random.uniform(*density_range_gcm3) * 1000
            radius_m = (self.mass_kg / ((4 / 3) * math.pi * density_kg_m3)) ** (1 / 3)
            self.radius_km = radius_m / 1000
            self.composition = "rock, metal, and (at the lower end of its mass range) ice, similar in bulk composition to Earth or Mars"

        self.has_internal_heat = random.random() < program_constants.ROGUE_PLANET_INTERNAL_HEAT_CHANCE
        self.has_moons = random.random() < program_constants.ROGUE_PLANET_MOON_CHANCE

        (self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy,
         self.galactic_orbital_phase_deg, self.galactic_min_update_interval_years) = \
            generate_galactic_orbit_fields()

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this planet's properties, per
        `SERIALIZABLE_FIELDS`.

        Returns:
            dict: One entry per field in `SERIALIZABLE_FIELDS`.
        """
        return fields_to_dict(self, self.SERIALIZABLE_FIELDS)

    @classmethod
    def from_dict(cls, data, system_config):
        """
        Reconstructs a `RoguePlanet` from a dict in the shape `to_dict()`
        produces, without re-running generation (`__init__` is bypassed
        via `object.__new__`).

        Args:
            data (dict): A dict in the shape `to_dict()` produces.
            system_config (SystemConfig): The shared config.

        Returns:
            RoguePlanet: The reconstructed planet.
        """
        planet = object.__new__(cls)
        planet.system_config = system_config
        fields_from_dict(planet, data, cls.SERIALIZABLE_FIELDS)
        return planet

    def to_paragraph_list(self):
        """
        Generates a list of descriptive paragraphs for the rogue planet: a
        header naming it, and a description of its physical properties and
        drift through interstellar space.

        Returns:
            list: A list of strings, where each string is a paragraph
                  describing the rogue planet.
        """
        header_level = '##' if self.system_config.MARKDOWN else '=='
        kind_label = "Gas Giant" if self.planet_type == 'g' else "Terrestrial"
        header = f"{header_level} {self.name} (Rogue {kind_label}) {header_level if not self.system_config.MARKDOWN else ''}".rstrip()

        description = (
            f"{self.name} is a free-floating planet adrift in interstellar space, bound to no star. It has a "
            f"radius of roughly {self.radius_km:,.0f} km and is composed of {self.composition}."
        )

        sentences = [description]
        if self.has_internal_heat:
            sentences.append(
                "It retains enough residual internal heat from its formation to remain geologically active."
            )
        else:
            sentences.append("Its interior has long since cooled to the ambient temperature of deep space.")
        if self.has_moons:
            sentences.append("A smaller companion body, likely captured after ejection, still orbits it.")
        sentences.append(
            f"Though bound to no star, it still orbits the galactic center at "
            f"{format_galactic_orbit(self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy)}."
        )

        return [header, " ".join(sentences)]

    def __str__(self):
        """
        Returns a string representation of the rogue planet, with
        paragraphs separated by double newlines.
        """
        return "\n\n".join(self.to_paragraph_list())


class InterstellarComet:
    """
    A basic class to store information for a small icy body passing
    through on an unbound, hyperbolic interstellar trajectory (real
    confirmed examples: 1I/'Oumuamua, 2I/Borisov).

    Attributes:
        name (str): A generated or explicitly given name (real interstellar
            objects are cataloged as `1I/...`, `2I/...`, etc. -- this
            generator's own name is purely descriptive/narrative).
        nucleus_diameter_km (float): Nucleus diameter in kilometers.
        composition (list): A list of composition component strings,
            sampled from `program_constants.COMET_COMPOSITION`.
        velocity_kms (float): Hyperbolic excess speed, in km/s.
        is_active (bool): Whether it currently shows a coma/tail from
            sublimating ices.
    """

    SERIALIZABLE_FIELDS = [
        "name", "nucleus_diameter_km", "velocity_kms", "is_active",
        "galactic_orbital_speed_kms", "galactic_orbital_period_gy",
        "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
    ]
    """Every attribute set by `__init__`, excluding `system_config` (a
    shared back-reference) and `composition` (handled separately in
    `to_dict`/`from_dict` -- a plain list, but kept alongside the others
    for symmetry with `AsteroidBelt.composition`'s own separate handling).
    `velocity_kms` (this comet's own hyperbolic excess speed relative to
    whatever star it passes) is a separate, non-advancing descriptive
    stat from `galactic_orbital_*` (its own bulk motion around the galactic
    center, see `nebulaData.Nebula`'s identical fields' docstring) -- the
    two don't conflict, the same way a planet's `rotation_period_hours`
    doesn't conflict with its own orbital motion."""

    def __init__(self, system_config: SystemConfig, name=None):
        """
        Initializes an InterstellarComet object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object
                (used only for its `MARKDOWN` flag, via `to_paragraph_list`).
            name (str, optional): An explicit name. Random if omitted.
        """
        reseed_rng()
        self.system_config = system_config
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)

        self.nucleus_diameter_km = random.uniform(*program_constants.INTERSTELLAR_COMET_NUCLEUS_DIAMETER_RANGE_KM)
        self.velocity_kms = random.uniform(*program_constants.INTERSTELLAR_OBJECT_SPEED_KMS_RANGE)
        self.is_active = random.random() < program_constants.INTERSTELLAR_COMET_ACTIVE_CHANCE

        num_components = min(3, len(program_constants.COMET_COMPOSITION))
        self.composition = random.sample(program_constants.COMET_COMPOSITION, k=num_components)

        (self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy,
         self.galactic_orbital_phase_deg, self.galactic_min_update_interval_years) = \
            generate_galactic_orbit_fields()

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this comet's properties, per
        `SERIALIZABLE_FIELDS`, plus `composition`.

        Returns:
            dict: One entry per field in `SERIALIZABLE_FIELDS`, plus
                 `composition`.
        """
        data = fields_to_dict(self, self.SERIALIZABLE_FIELDS)
        data["composition"] = list(self.composition)
        return data

    @classmethod
    def from_dict(cls, data, system_config):
        """
        Reconstructs an `InterstellarComet` from a dict in the shape
        `to_dict()` produces, without re-running generation (`__init__` is
        bypassed via `object.__new__`).

        Args:
            data (dict): A dict in the shape `to_dict()` produces.
            system_config (SystemConfig): The shared config.

        Returns:
            InterstellarComet: The reconstructed comet.
        """
        comet = object.__new__(cls)
        comet.system_config = system_config
        fields_from_dict(comet, data, cls.SERIALIZABLE_FIELDS)
        comet.composition = list(data["composition"])
        return comet

    def get_composition_summary(self):
        """
        Builds the human-readable composition summary phrase (e.g. "water
        ice, silicate dust, and complex organic compounds"), the same
        as-published-text role `AsteroidBelt.get_composition_summary`
        plays for the database persistence layer.

        Returns:
            str: The composition summary phrase.
        """
        if not self.composition:
            return "unknown composition"
        if len(self.composition) == 1:
            return self.composition[0]
        if len(self.composition) == 2:
            return f"{self.composition[0]} and {self.composition[1]}"
        return ", ".join(self.composition[:-1]) + f", and {self.composition[-1]}"

    def to_paragraph_list(self):
        """
        Generates a list of descriptive paragraphs for the comet: a header
        naming it, and a description of its size, composition, speed, and
        current activity.

        Returns:
            list: A list of strings, where each string is a paragraph
                  describing the comet.
        """
        header_level = '##' if self.system_config.MARKDOWN else '=='
        header = f"{header_level} {self.name} (Interstellar Comet) {header_level if not self.system_config.MARKDOWN else ''}".rstrip()

        activity = (
            "It currently shows an active coma and tail as its surface ices sublimate."
            if self.is_active else
            "It shows no current activity, passing through inert."
        )
        description = (
            f"{self.name} is a small icy body on a hyperbolic, unbound trajectory through interstellar space, "
            f"with a nucleus roughly {self.nucleus_diameter_km:.2f} km across, composed of "
            f"{self.get_composition_summary()}. It is traveling at a hyperbolic excess speed of "
            f"{self.velocity_kms:.1f} km/s relative to any star it passes. {activity} Its own bulk motion "
            f"still carries it around the galactic center at "
            f"{format_galactic_orbit(self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy)}."
        )

        return [header, description]

    def __str__(self):
        """
        Returns a string representation of the comet, with paragraphs
        separated by double newlines.
        """
        return "\n\n".join(self.to_paragraph_list())
