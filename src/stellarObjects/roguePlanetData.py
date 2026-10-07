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
from .rogueSurface import ROGUE_SURFACE_FIELDS, SURFACE_REGIME_LABELS, rogue_surface_conditions
from .names import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from . import physical_constants, planetPhysics, program_constants
from planetgen.util import log
from planetgen.util.serialization import fields_from_dict, fields_to_dict
from .utils import (format_body_radius_km, format_galactic_orbit, format_number, format_speed_kms,
                    generate_galactic_orbit_fields,
                    generate_phoneme_salad_name, sample_power_law)


def format_comet_composition_summary(composition):
    """
    Builds the human-readable composition summary phrase (e.g. "water
    ice, silicate dust, and complex organic compounds") for a plain list
    of component strings, in the shape `InterstellarComet.composition`/
    `cometData.Comet.composition` both use -- extracted here, the same
    way `asteroidData.format_composition_summary` is, so `cometData.Comet`
    (a star-bound sibling of `InterstellarComet`, see that module's own
    docstring) can reuse this formatting logic instead of duplicating it.

    Args:
        composition (list): A list of component strings.

    Returns:
        str: The composition summary phrase, or "unknown composition" if
            `composition` is empty.
    """
    if not composition:
        return "unknown composition"
    if len(composition) == 1:
        return composition[0]
    if len(composition) == 2:
        return f"{composition[0]} and {composition[1]}"
    return ", ".join(composition[:-1]) + f", and {composition[-1]}"


def infer_rogue_mass_bin(mass_kg):
    """The `ROGUE_PLANET_MASS_BIN_CHOICES` bin `mass_kg` falls in, for
    rows saved before the bin was stored (schema v37)."""
    if mass_kg >= program_constants.ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER[0] * physical_constants.JUPITER_MASS_TO_KG:
        return "brown-dwarf"
    mass_earth = mass_kg / physical_constants.EARTH_MASS_TO_KG
    for name, (_low, high, _rate) in program_constants.ROGUE_PLANET_MASS_BINS.items():
        if mass_earth < high:
            return name
    return "jupiter"


def rogue_planet_classes(planet_type):
    """Every `PLANET_CLASSES` class a rogue planet of `planet_type`
    (`'t'` or `'g'`) may have: those flagged `"r": True` (GEN.8), whose
    own `"type"` matches."""
    return [code for code, data in program_constants.PLANET_CLASSES.items()
            if data.get("r") and data["type"] == planet_type]


def rogue_planet_class_candidates(planet_type, radius_km, mass_kg):
    """
    The rogue-eligible classes (`rogue_planet_classes`) a rogue planet
    fits: its radius inside the class's `radius_range` and its mass inside
    the class's mass range (`planetPhysics.planet_mass_ranges`). When none
    fits -- a rogue super-Earth past class C's 10,000 km ceiling, say --
    the eligible class whose radius range is nearest (in log radius), so
    every rogue planet still gets a class.

    Args:
        planet_type (str): `'t'` or `'g'`.
        radius_km (float): Radius in kilometers.
        mass_kg (float): Mass in kilograms.

    Returns:
        list: Class codes, never empty for `'t'` or `'g'`.
    """
    eligible = rogue_planet_classes(planet_type)
    fitting = [
        code for code in eligible
        if program_constants.PLANET_CLASSES[code]["radius_range"][0] <= radius_km
        <= program_constants.PLANET_CLASSES[code]["radius_range"][1]
        and planetPhysics.planet_mass_ranges[code][0] <= mass_kg <= planetPhysics.planet_mass_ranges[code][1]
    ]
    if fitting or not eligible:
        return fitting

    def gap(code):
        low, high = program_constants.PLANET_CLASSES[code]["radius_range"]
        return max(math.log(low / radius_km), math.log(radius_km / high), 0.0)
    return [min(eligible, key=gap)]


def choose_rogue_planet_class(planet_type, radius_km, mass_kg, mass_bin=None):
    """
    Draws a rogue planet's class (GEN.8) from `rogue_planet_class_candidates`,
    weighted by `PLANET_CLASS_PROBABILITIES` like a star's planet. A
    brown dwarf is a failed star, not a planet, so it gets none.

    Returns:
        str or None: The class code, or `None` for a brown dwarf.
    """
    if mass_bin == "brown-dwarf":
        return None
    candidates = rogue_planet_class_candidates(planet_type, radius_km, mass_kg)
    if not candidates:
        return None
    if len(candidates) == 1:
        # No draw, so a lone candidate leaves the random stream as it was.
        log.choice("Rogue planet class", candidates[0], "the only rogue class that fits its type, radius and mass")
        return candidates[0]
    return planetPhysics._choose_weighted_planet_class(candidates)


def default_rogue_planet_class(planet_type, radius_km, mass_kg, mass_bin=None):
    """
    `choose_rogue_planet_class` without the draw: the most probable
    candidate. For rows saved before rogue planets had a class (schema
    v48) and dicts written before it, so a reload gives the same answer
    every time.
    """
    if mass_bin == "brown-dwarf":
        return None
    candidates = rogue_planet_class_candidates(planet_type, radius_km, mass_kg)
    if not candidates:
        return None
    return max(candidates, key=lambda code: program_constants.PLANET_CLASS_PROBABILITIES.get(code, 0.0))


class RoguePlanet:
    """
    A basic class to store information for a free-floating ("rogue"/nomad)
    planet with no host star.

    Attributes:
        name (str): A generated or explicitly given name.
        planet_type (str): `'t'` (terrestrial/icy) or `'g'` (gas giant),
            the same letters `Planet.body_type` uses, chosen by mass
            relative to `program_constants.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER`.
        planet_class (str or None): Its `PLANET_CLASSES` letter (GEN.8),
            drawn from the classes flagged `"r"`; `None` for a brown dwarf.
        mass_kg (float): Mass in kilograms.
        radius_km (float): Radius in kilometers.
        composition (str): A descriptive bulk-composition string.
        has_internal_heat (bool): Whether it is still geologically active
            (a giant always is): its heat flow reaches
            `ROGUE_ACTIVE_HEAT_FLUX_W_M2` (from `rogueSurface`).
        has_moons (bool): Whether it retains a captured companion moon.
        age_gy, internal_heat_flux_w_m2, effective_temperature_k,
        surface_regime, surface_temperature_k, surface_pressure_pa,
        ice_shell_thickness_km, ocean_depth_km, has_liquid_water: Its
            surface conditions (`rogueSurface.rogue_surface_conditions`).
    """

    SERIALIZABLE_FIELDS = [
        "name", "planet_type", "planet_class", "mass_bin", "mass_kg", "radius_km", "composition",
        "has_internal_heat", "has_moons", *ROGUE_SURFACE_FIELDS,
        "galactic_orbital_speed_kms", "galactic_orbital_period_gy",
        "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
    ]
    """Every attribute set by `__init__`, excluding `system_config` (a
    shared back-reference). The `galactic_orbital_*` fields (see
    `nebulaData.Nebula`'s identical fields' docstring) reflect that a
    rogue planet is unbound from any specific STAR, not from the galaxy
    itself -- it still orbits the galactic center like a lone star does."""

    def __init__(self, system_config: SystemConfig, name=None, mass_bin=None):
        """
        Initializes a RoguePlanet object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object
                (used only for its `MARKDOWN` flag, via `to_paragraph_list`).
            name (str, optional): An explicit name. Random if omitted.
            mass_bin (str, optional): One of
                `program_constants.ROGUE_PLANET_MASS_BIN_CHOICES`. Omitted,
                a planet bin is drawn by its per-star rate
                (`ROGUE_PLANET_MASS_BINS`); `"brown-dwarf"` is only ever
                asked for explicitly (its own rate,
                `PHENOMENON_DENSITY_PC3["brown-dwarf"]`).
        """
        self.system_config = system_config
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)
        self.name_given = bool(name)  # a given name is kept over an object ID (GEN.64)

        if mass_bin is None:
            bins = program_constants.ROGUE_PLANET_MASS_BINS
            mass_bin = random.choices(list(bins), weights=[rate for _lo, _hi, rate in bins.values()])[0]
            log.choice("Rogue planet mass bin", mass_bin, "drawn by ROGUE_PLANET_MASS_BINS' per-star rates")
        elif mass_bin not in program_constants.ROGUE_PLANET_MASS_BIN_CHOICES:
            raise ValueError(f"mass_bin must be one of {program_constants.ROGUE_PLANET_MASS_BIN_CHOICES}, got {mass_bin!r}")
        self.mass_bin = mass_bin

        if mass_bin == "brown-dwarf":
            low, high = (m * physical_constants.JUPITER_MASS_TO_KG
                         for m in program_constants.ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER)
            self.mass_kg = math.exp(random.uniform(math.log(low), math.log(high)))
        else:
            # The same mass function the bins' rates come from (GEN.45).
            low, high, _rate = program_constants.ROGUE_PLANET_MASS_BINS[mass_bin]
            low, high = low * physical_constants.EARTH_MASS_TO_KG, high * physical_constants.EARTH_MASS_TO_KG
            self.mass_kg = sample_power_law(low, high, program_constants.ROGUE_PLANET_MASS_FUNCTION_SLOPE)
        mass_jupiter = self.mass_kg / physical_constants.JUPITER_MASS_TO_KG

        if mass_bin == "brown-dwarf":
            log.choice("Rogue planet type", "brown dwarf", "mass_bin 'brown-dwarf'")
            self.planet_type = 'g'
            # Brown dwarfs share Jupiter's near-flat mass-radius relation
            # (Chabrier & Baraffe 2000), slightly smaller when old.
            self.radius_km = physical_constants.JUPITER_RADIUS_KM * random.uniform(0.75, 1.1)
            self.composition = (
                "hydrogen and helium, a failed star that briefly fused deuterium and now glows faintly "
                "in the infrared as it cools"
            )
        elif mass_jupiter >= program_constants.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER:
            log.choice("Rogue planet type", "gas giant",
                       f"mass {mass_jupiter:.4g} Mjup >= ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER "
                       f"({program_constants.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER})")
            self.planet_type = 'g'
            # The same giant mass-radius relation as a bound giant (GEN.34,
            # GEN.60): radius grows with mass up to about Saturn's mass,
            # then stays near Jupiter's as degeneracy pressure takes over.
            self.radius_km = planetPhysics.sample_giant_radius_km(self.mass_kg)
            self.composition = "hydrogen and helium, similar in bulk composition to Jupiter or Saturn"
        else:
            log.choice("Rogue planet type", "terrestrial",
                       f"mass {mass_jupiter:.4g} Mjup < ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER "
                       f"({program_constants.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER})")
            self.planet_type = 't'
            density_range_gcm3 = physical_constants.PLANET_DENSITY["t"]
            density_kg_m3 = random.uniform(*density_range_gcm3) * 1000
            radius_m = (self.mass_kg / ((4 / 3) * math.pi * density_kg_m3)) ** (1 / 3)
            self.radius_km = radius_m / 1000
            if mass_bin == "sub-neptune":
                self.composition = (
                    "rock and ice, possibly beneath a thin hydrogen envelope, between a super-Earth and "
                    "Neptune in bulk composition"
                )
            else:
                self.composition = "rock, metal, and (at the lower end of its mass range) ice, similar in bulk composition to Earth or Mars"

        self.planet_class = choose_rogue_planet_class(self.planet_type, self.radius_km, self.mass_kg, mass_bin)

        self.has_moons = random.random() < program_constants.ROGUE_PLANET_MOON_CHANCE
        self._apply_surface_conditions(rogue_surface_conditions(
            self.mass_kg, self.radius_km, self.planet_type, self.mass_bin, self.has_moons))

        (self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy,
         self.galactic_orbital_phase_deg, self.galactic_min_update_interval_years) = \
            generate_galactic_orbit_fields()

    def _apply_surface_conditions(self, conditions):
        """Sets `ROGUE_SURFACE_FIELDS` and `has_internal_heat` from a
        `rogue_surface_conditions` result."""
        for field in ROGUE_SURFACE_FIELDS:
            setattr(self, field, conditions[field])
        self.has_internal_heat = conditions["has_internal_heat"]

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
        if getattr(planet, "mass_bin", None) is None:
            planet.mass_bin = infer_rogue_mass_bin(planet.mass_kg)
        if "planet_class" not in data:
            planet.planet_class = default_rogue_planet_class(
                planet.planet_type, planet.radius_km, planet.mass_kg, planet.mass_bin)
        if "surface_regime" not in data:
            # Written before rogues had surface conditions (schema v48):
            # computed from a generator seeded by the name, so a reload
            # always gives the same conditions.
            planet._apply_surface_conditions(rogue_surface_conditions(
                planet.mass_kg, planet.radius_km, planet.planet_type, planet.mass_bin, planet.has_moons,
                random.Random(planet.name)))
        return planet

    @property
    def kind_label(self):
        """`"Brown Dwarf"`, `"Gas Giant"` or `"Terrestrial"`."""
        if getattr(self, "mass_bin", None) == "brown-dwarf":
            return "Brown Dwarf"
        return "Gas Giant" if self.planet_type == 'g' else "Terrestrial"

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
        kind_label = self.kind_label
        if getattr(self, "planet_class", None):
            kind_label += f", Class {self.planet_class}"
        header = f"{header_level} {self.name} (Rogue {kind_label}) {header_level if not self.system_config.MARKDOWN else ''}".rstrip()

        what = "brown dwarf" if self.kind_label == "Brown Dwarf" else "planet"
        description = (
            f"{self.name} is a free-floating {what} adrift in interstellar space, bound to no star. It has a "
            f"radius of roughly {format_body_radius_km(self.system_config, self.radius_km)} and is composed of "
            f"{self.composition}."
        )

        sentences = [description]
        sentences.extend(self.surface_sentences())
        if self.has_moons:
            sentences.append("A smaller companion body, likely captured after ejection, still orbits it.")
        sentences.append(
            f"Though bound to no star, it still orbits the galactic center at "
            f"{format_galactic_orbit(self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy)}."
        )

        return [header, " ".join(sentences)]

    def surface_sentences(self):
        """
        Plain sentences on its age, own heat and surface (`rogueSurface`):
        with no star, its only warmth is what leaks out of its interior.

        Returns:
            list[str]: The sentences.
        """
        regime = self.surface_regime
        temp = self.surface_temperature_k
        sentences = [
            f"It is about {format_number(self.age_gy, ',.1f')} billion years old. With no star to warm it, its only heat is its own: "
            f"{format_number(self.internal_heat_flux_w_m2, ',.3g')} W/m\u00b2 leaking from its interior, so it glows at an effective "
            f"temperature of {format_number(self.effective_temperature_k)} K."
        ]
        if regime in ("gas-giant", "brown-dwarf"):
            sentences.append(f"It has no solid surface; at the 1 bar level the temperature is {format_number(temp)} K.")
            return sentences
        if regime == "hydrogen-envelope":
            sentences.append(
                f"A thick hydrogen envelope, opaque to heat at high pressure, blankets it, so the ground beneath "
                f"{format_number(self.surface_pressure_pa / 1e5)} bar of gas sits at {format_number(temp)} K.")
            if self.has_liquid_water and not self.ice_shell_thickness_km:
                sentences.append(f"There it holds a liquid ocean about {format_number(self.ocean_depth_km)} km deep.")
        elif regime == "bare-rock":
            sentences.append(f"Its bare rock surface has frozen to {format_number(temp)} K.")
        else:
            sentences.append(
                f"Its surface has frozen to {format_number(temp)} K, and whatever air it had lies on the ground as frost.")
        if self.ice_shell_thickness_km and self.has_liquid_water:
            sentences.append(
                f"Under an ice shell about {format_number(self.ice_shell_thickness_km)} km thick, its own heat keeps a "
                f"liquid ocean about {format_number(self.ocean_depth_km)} km deep.")
        elif self.ice_shell_thickness_km:
            sentences.append(f"Its water is frozen solid, an ice layer about "
                             f"{format_number(self.ice_shell_thickness_km)} km thick down to the rock.")
        return sentences

    @property
    def surface_regime_label(self):
        """`SURFACE_REGIME_LABELS`' text for `surface_regime`."""
        return SURFACE_REGIME_LABELS.get(self.surface_regime, self.surface_regime)

    def __str__(self):
        """
        Returns a string representation of the rogue planet, with
        paragraphs separated by double newlines.
        """
        return "\n\n".join(self.to_paragraph_list())


def interstellar_comet_designation(sector_code, index):
    """
    An interstellar comet's designation (v40): `I/<sector>-<n>`, after
    the IAU's `I/` prefix, where `<sector>` is the sector it was found in
    (its grid designation, or its name) and `<n>` counts that sector's
    interstellar comets from 1. `I/<n>` for one found in no sector.
    """
    return f"I/{sector_code}-{index}" if sector_code else f"I/{index}"


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
        self.system_config = system_config
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)
        self.name_given = bool(name)  # a given name is kept over an object ID (GEN.64)

        self.nucleus_diameter_km = random.uniform(*program_constants.INTERSTELLAR_COMET_NUCLEUS_DIAMETER_RANGE_KM)
        self.velocity_kms = random.uniform(*program_constants.INTERSTELLAR_OBJECT_SPEED_KMS_RANGE)
        self.is_active = random.random() < program_constants.INTERSTELLAR_COMET_ACTIVE_CHANCE
        log.choice("Interstellar comet activity", self.is_active,
                   f"roll against INTERSTELLAR_COMET_ACTIVE_CHANCE "
                   f"({program_constants.INTERSTELLAR_COMET_ACTIVE_CHANCE})")

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
        return format_comet_composition_summary(self.composition)

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
            f"{format_speed_kms(self.velocity_kms)} relative to any star it passes. {activity} Its own bulk motion "
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
