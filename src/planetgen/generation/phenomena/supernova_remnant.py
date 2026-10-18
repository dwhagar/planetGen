# planetgen/generation/phenomena/supernova_remnant.py

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
explosion, embedded via `phenomena.compact_remnant.BlackHole`/`NeutronStar`.
"""

import math

from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import REMNANT_CLASS_LETTERS, choose_weighted_class, draw_class_contents, typical_class_contents
from planetgen import tuning
from planetgen.util import draw
from planetgen.util import log
from planetgen.util.serialization import fields_from_dict, fields_to_dict
from planetgen.galaxy.galactic_orbit import generate_galactic_orbit_fields
from planetgen.names.wordsalad import LazySaladName


def sedov_taylor_radius_ly(age_years, density_cm3, energy_e51=None):
    """
    A supernova remnant's radius, light-years (GEN.151): the Sedov-Taylor
    blast wave `C (E51 / n)^(1/5) t^(2/5)` in the ambient density `n`
    (`density_cm3`), held at the snowplow onset radius
    `45.66 ly * E51^(2/7) * n^(-3/7)` once the remnant has cooled that far
    (`tuning.SEDOV_TAYLOR_*`, `tuning.SNOWPLOW_ONSET_*`).

    Args:
        age_years (float): Years since the explosion.
        density_cm3 (float): Ambient number density, cm^-3.
        energy_e51 (float, optional): Explosion energy in units of 1e51
            erg (`tuning.SEDOV_TAYLOR_ENERGY_E51` when omitted).
    """
    if energy_e51 is None:
        energy_e51 = tuning.SEDOV_TAYLOR_ENERGY_E51
    sedov = (tuning.SEDOV_TAYLOR_RADIUS_COEFFICIENT_LY
             * energy_e51 ** (-tuning.SEDOV_TAYLOR_DENSITY_EXPONENT)
             * density_cm3 ** tuning.SEDOV_TAYLOR_DENSITY_EXPONENT
             * age_years ** tuning.SEDOV_TAYLOR_TIME_EXPONENT)
    snowplow = (tuning.SNOWPLOW_ONSET_RADIUS_COEFFICIENT_LY
                * energy_e51 ** tuning.SNOWPLOW_ONSET_ENERGY_EXPONENT
                * density_cm3 ** tuning.SNOWPLOW_ONSET_DENSITY_EXPONENT)
    return min(sedov, snowplow)


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


class SupernovaRemnant(LazySaladName):
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
                (threaded into any
                embedded `compact_remnant`).
            name (str, optional): An explicit name. Random if omitted.
        """
        self.system_config = system_config
        # A draft name: `store.insert_supernova_remnant` reserves it through
        # the system-name registry (v40), and the core follows.
        self._init_salad_name(name)  # drawn on first read (PERF.43)
        self.name_given = bool(name)  # a given name is kept over an object ID (GEN.64)

        type_ia_chance = getattr(system_config, "TYPE_IA_SHARE", None)
        if type_ia_chance is None:
            type_ia_chance = tuning.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE
        is_type_ia = draw.random() < type_ia_chance
        self.progenitor_type = "Type Ia" if is_type_ia else "core-collapse"
        log.choice("Supernova progenitor type", self.progenitor_type,
                   f"roll against the Type Ia share of the place ({type_ia_chance:.3f}; the solar circle's is "
                   f"SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE, {tuning.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE})")

        (self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy,
         self.galactic_orbital_phase_deg, self.galactic_min_update_interval_years) = \
            generate_galactic_orbit_fields()

        # Core-collapse remnants keep a neutron star (most) or black hole,
        # offset from the center by its birth kick (compact_offset_ly
        # below); thermonuclear (Type Ia, class W) remnants have none. The
        # core is built once the remnant's age is known: it is as old.
        kind = None
        if not is_type_ia and draw.random() < tuning.SUPERNOVA_CORE_COLLAPSE_REMNANT_VISIBLE_CHANCE:
            if draw.random() < tuning.SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE:
                kind = "black_hole"
                log.choice("Core-collapse compact remnant", "black hole",
                           f"roll passed SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE "
                           f"({tuning.SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE})")
            else:
                kind = "neutron_star"
                log.choice("Core-collapse compact remnant", "neutron star",
                           f"roll failed SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE "
                           f"({tuning.SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE})")
        elif not is_type_ia:
            log.debug("Core-collapse compact remnant: none visible/detectable "
                       "(roll failed SUPERNOVA_CORE_COLLAPSE_REMNANT_VISIBLE_CHANCE)")

        self.remnant_class = choose_weighted_class(
            remnant_classes_for(self.progenitor_type, kind), "Supernova remnant class")
        class_data = tuning.NEBULA_CLASSES[self.remnant_class]
        self.morphology = class_data["morphology"]
        self.age_years = draw.uniform(*class_data["age_range_years"])
        self.compact_remnant = None
        if kind == "black_hole":
            self.compact_remnant = BlackHole(system_config, name=f"{self.name} Core")
        elif kind == "neutron_star":
            self.compact_remnant = NeutronStar(system_config, name=f"{self.name} Core",
                                               age_gy=self.age_years / 1e9)
        (self.dominant_species, self.density_cm3, self.temperature_k,
         self.extinction_av) = draw_class_contents(self.remnant_class)
        ambient_density = getattr(system_config, "AMBIENT_DENSITY_CM3", None)
        if ambient_density is None:
            ambient_density = self.density_cm3
        else:
            log.choice("Supernova remnant ambient density", f"{ambient_density:.4g} cm^-3",
                       "the remnant sits in a superbubble: its Sedov-Taylor size reads the cavity's interior density")
        self.radius_ly = sedov_taylor_radius_ly(self.age_years, ambient_density)
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
        speed_kms = math.exp(draw.uniform(math.log(low), math.log(high)))
        distance_ly = speed_kms / 299792.458 * self.age_years
        z = draw.uniform(-1.0, 1.0)
        phi = draw.uniform(0.0, 2 * math.pi)
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

from planetgen.generation import object_timing  # noqa: E402 -- registers the generators below (PERF.92)
object_timing.track(SupernovaRemnant, "supernova-remnant")
