# planetgen/generation/phenomena/rogue.py

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

from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena import hosts
from planetgen.physics.rogue_surface import ROGUE_SURFACE_FIELDS, ROGUE_SURFACE_UNITS, SURFACE_REGIME_LABELS, rogue_surface_conditions
from planetgen.physics import body_strength, brown_dwarf as brownDwarf, bulk_structure, constants, planet_classes, planets as planetPhysics, relations, rings, rogue_habitability, small_body_composition, spin
from planetgen.physics.stellar_evolution import sample_star_age_gy
from planetgen import tuning
from planetgen.util import draw
from planetgen.util import log
from planetgen.names.wordsalad import LazySaladName
from planetgen.generation.mass_object import GALACTIC_MOTION_UNITS
from planetgen.generation.planetary_object import PlanetaryObject
from planetgen.generation.small_body import SmallBody
from planetgen.util.format import format_flux_w_m2, format_number
from planetgen.util.random import sample_power_law
from planetgen.physics import unit_registry as units

_m_to_km = units.converter("m", "km")


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
    if mass_kg >= _mjup_to_kg(tuning.ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER[0]):
        return "brown-dwarf"
    mass_earth = _kg_to_mearth(mass_kg)
    for name, (_low, high, _rate) in tuning.ROGUE_PLANET_MASS_BINS.items():
        if mass_earth < high:
            return name
    return "jupiter"


ICE_KG_M3 = 917.0
OCEAN_KG_M3 = 1000.0
HIGH_PRESSURE_ICE_KG_M3 = 1300.0
"""float: The densities of the ice shell, the ocean and the high-pressure ice under it (kg/m3) that turn their thicknesses into a water mass."""


def water_mass_fraction(mass_kg, radius_km, ice_shell_km, ocean_km, hp_ice_km):
    """A rogue planet's water as a share of its mass, from the thicknesses of the layers its surface conditions hold (km, `None` for
    none): the shells' volumes at the radius times their densities, over the mass."""
    if not (ice_shell_km or ocean_km or hp_ice_km):
        return 0.0
    outer_m = radius_km * 1000.0
    kg = 0.0
    for thickness_km, density in ((ice_shell_km, ICE_KG_M3), (ocean_km, OCEAN_KG_M3), (hp_ice_km, HIGH_PRESSURE_ICE_KG_M3)):
        if thickness_km:
            inner_m = max(0.0, outer_m - thickness_km * 1000.0)
            kg += 4.0 / 3.0 * math.pi * (outer_m ** 3 - inner_m ** 3) * density
            outer_m = inner_m
    return min(1.0, kg / mass_kg)


def rocky_radius_km(mass_kg, density_kg_m3):
    """
    The radius of a rogue's rock, km (GEN.1133): the drawn density is the bulk density the body has at one Earth mass, which says its core
    mass fraction (`bulk_structure.cmf_from_mass_radius`; a rock lighter than any dry mix is the coreless one), and the radius at its own
    mass is the law of the bound planets (`bulk_structure.radius_fit_earth`, with the compression a sphere of fixed density leaves out).
    A rogue of 10 Earth masses is then 1.8 to 2.0 Earth radii, as a bound one is, and not the 2.2 to 2.4 of a body that does not compress.
    """
    earth_radius_m = constants.EARTH_RADIUS_KM * 1000.0
    radius_at_one = (_mearth_to_kg(1.0) / ((4 / 3) * math.pi * density_kg_m3)) ** (1 / 3) / earth_radius_m
    core = bulk_structure.cmf_from_mass_radius(1.0, radius_at_one)
    mass_earth = _kg_to_mearth(mass_kg)
    return bulk_structure.radius_fit_earth(mass_earth, 0.0 if core is None else core) * constants.EARTH_RADIUS_KM


def rogue_class(planet):
    """
    The class of a rogue planet (GEN.687), derived like any body's (`planet_classes.classify`) from its mass, gravity and the surface
    its own heat gives it: the temperature, the pressure of the air it holds (`surface_pressure_pa`) and the water its layers hold
    (`water_mass_fraction`). Starless, it has no flux variation, no ocean to cover its surface and no gas mixture (an air of frozen
    nitrogen or a hydrogen envelope reads as a plain air, whose shares are not stored); `None` for a brown dwarf, a failed star.
    """
    if planet.mass_bin == "brown-dwarf":
        return None
    gravity = relations.surface_gravity(planet.mass_kg, planet.radius_km * 1000.0) / constants.EARTH_GRAVITY
    return planet_classes.classify({
        "body_type": planet.planet_type, "mass_kg": planet.mass_kg, "gravity_g": gravity,
        "surface_temperature_k": planet.surface_temperature_k, "atmospheric_pressure_pa": planet.surface_pressure_pa,
        "water_mass_fraction": water_mass_fraction(planet.mass_kg, planet.radius_km, planet.ice_shell_thickness_km,
                                                   planet.ocean_depth_km, getattr(planet, "hp_ice_km", None)),
    })


class RoguePlanet(LazySaladName, PlanetaryObject):
    """
    A basic class to store information for a free-floating ("rogue"/nomad)
    planet with no host star.

    Attributes:
        name (str): A generated or explicitly given name.
        planet_type (str): `'t'` (terrestrial/icy) or `'g'` (gas giant),
            the same letters `Planet.body_type` uses, chosen by mass
            relative to `tuning.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER`.
        planet_class (str or None): Its class letter (GEN.8, GEN.687), derived from its surface (`rogue_class`);
            `None` for a brown dwarf.
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

    FIELD_UNITS = {
        "name": "text", "planet_type": "text", "planet_class": "text", "mass_bin": "text", "mass_kg": "kg",
        "radius_km": "km", "composition": "text", "has_internal_heat": "flag", "has_moons": "flag",
        "density_g_cm3": "g_cm3", "roundness_ratio": "ratio", "yield_strength_pa": "Pa", "binding_energy_j": "J",
        **rings.UNITS,
        **ROGUE_SURFACE_UNITS, **rogue_habitability.ROGUE_UNITS, **GALACTIC_MOTION_UNITS, "rotation_period_hours": "h", **spin.SPIN_UNITS,
    }
    SERIALIZABLE_FIELDS = list(FIELD_UNITS)
    """Every attribute set by `__init__`, with its native unit (GEN.570), excluding `system_config` (a
    shared back-reference). The `galactic_orbital_*` fields (see
    `nebulaData.Nebula`'s identical fields' docstring) reflect that a
    rogue planet is unbound from any specific STAR, not from the galaxy
    itself -- it still orbits the galactic center like a lone star does."""

    def __init__(self, system_config: SystemConfig, name=None, mass_bin=None):
        """
        Initializes a RoguePlanet object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object
                (its settings).
            name (str, optional): An explicit name. Random if omitted.
            mass_bin (str, optional): One of
                `tuning.ROGUE_PLANET_MASS_BIN_CHOICES`. Omitted,
                a planet bin is drawn by its per-star rate
                (`ROGUE_PLANET_MASS_BINS`); `"brown-dwarf"` is only ever
                asked for explicitly (its own rate,
                `PHENOMENON_DENSITY_PC3["brown-dwarf"]`).
        """
        self.system_config = system_config
        self._init_salad_name(name)  # drawn on first read (PERF.43)
        self.name_given = bool(name)  # a given name is kept over an object ID (GEN.64)

        if mass_bin is None:
            bins = tuning.ROGUE_PLANET_MASS_BINS
            mass_bin = draw.choices(list(bins), weights=[rate for _lo, _hi, rate in bins.values()])[0]
            log.choice("Rogue planet mass bin", mass_bin, "drawn by ROGUE_PLANET_MASS_BINS' per-star rates")
        elif mass_bin not in tuning.ROGUE_PLANET_MASS_BIN_CHOICES:
            raise ValueError(f"mass_bin must be one of {tuning.ROGUE_PLANET_MASS_BIN_CHOICES}, got {mass_bin!r}")
        self.mass_bin = mass_bin

        if mass_bin == "brown-dwarf":
            # GEN.252: the substellar mass function of the 20 pc census [R2374], not log-uniform.
            self.mass_kg = _msun_to_kg(brownDwarf.sample_mass_solar(draw))
        else:
            # The same mass function the bins' rates come from (GEN.45).
            low, high, _rate = tuning.ROGUE_PLANET_MASS_BINS[mass_bin]
            low, high = _mearth_to_kg(low), _mearth_to_kg(high)
            self.mass_kg = sample_power_law(low, high, tuning.ROGUE_PLANET_MASS_FUNCTION_SLOPE)
        mass_jupiter = _kg_to_mjup(self.mass_kg)

        if mass_bin == "brown-dwarf":
            log.choice("Rogue planet type", "brown dwarf", "mass_bin 'brown-dwarf'")
            self.planet_type = 'g'
            # The radius of the cooling track at its age (GEN.549): larger when young, nearly Jupiter's when old. The age
            # is drawn here, before the surface conditions, which are given it. It is a star's age (GEN.553): the
            # window of the population it belongs to, spread by that population's star formation rate; a brown dwarf
            # made without a place has the whole disk's window.
            brown_dwarf_age_gy = sample_star_age_gy(
                system_config.AGE, draw, getattr(system_config, "POPULATION", None), getattr(system_config, "AGE_LAW", None))
            self.radius_km = _m_to_km(
                brownDwarf.radius_sol(_kg_to_msun(self.mass_kg), brown_dwarf_age_gy) * constants.SOLAR_RADIUS_M)
            self.composition = (
                "hydrogen and helium, a failed star that briefly fused deuterium and now glows faintly "
                "in the infrared as it cools"
            )
        elif mass_jupiter >= tuning.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER:
            log.choice("Rogue planet type", "gas giant",
                       f"mass {mass_jupiter:.4g} Mjup >= ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER "
                       f"({tuning.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER})")
            self.planet_type = 'g'
            # The same giant mass-radius relation as a bound giant (GEN.34,
            # GEN.60): radius grows with mass up to about Saturn's mass,
            # then stays near Jupiter's as degeneracy pressure takes over.
            self.radius_km = planetPhysics.sample_giant_radius_km(self.mass_kg)
            if self.mass_kg < _mearth_to_kg(tuning.GIANT_RUNAWAY_MIN_MASS_EARTH):
                # GEN.1133: the core that did not run away (the bound chain's ice giant, `planet_mass.gas_outcome`)
                self.composition = ("water, ammonia and methane ices and rock under a hydrogen and helium envelope, an ice giant "
                                    "similar in bulk composition to Uranus or Neptune")
            else:
                self.composition = "hydrogen and helium, similar in bulk composition to Jupiter or Saturn"
        else:
            log.choice("Rogue planet type", "terrestrial",
                       f"mass {mass_jupiter:.4g} Mjup < ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER "
                       f"({tuning.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER})")
            self.planet_type = 't'
            density_range_gcm3 = constants.PLANET_DENSITY["t"]
            density_kg_m3 = draw.uniform(*density_range_gcm3) * 1000
            self.radius_km = rocky_radius_km(self.mass_kg, density_kg_m3)
            if mass_bin == "sub-neptune":
                self.composition = (
                    "rock and ice, possibly beneath a thin hydrogen envelope, between a super-Earth and "
                    "Neptune in bulk composition"
                )
            else:
                self.composition = "rock, metal, and (at the lower end of its mass range) ice, similar in bulk composition to Earth or Mars"

        self.has_moons = hosts.holds_something("planet-satellites", "rogue-planet")
        self._apply_surface_conditions(rogue_surface_conditions(
            self.mass_kg, self.radius_km, self.planet_type, self.mass_bin, self.has_moons,
            age_gy=brown_dwarf_age_gy if mass_bin == "brown-dwarf" else None))
        self.planet_class = rogue_class(self)    # from the surface its own heat gave it (GEN.687)
        # GEN.1162: its mean density and binding energy, and for a terrestrial one the roundness ratio and strength of its rock and the water of its shell and ocean.
        water = water_mass_fraction(self.mass_kg, self.radius_km, self.ice_shell_thickness_km, self.ocean_depth_km, getattr(self, "hp_ice_km", None)) \
            if self.planet_type == "t" else 0.0
        for name, value in body_strength.rogue_values(self, water).items():
            setattr(self, name, value)

        rings.clear(self)  # GEN.1057: a ring is what a planet's history made, and a rogue planet's moons and comets are not followed
        self.set_galactic_motion()
        # GEN.104, GEN.842: a rocky rogue planet is born spinning at the rate and angle the giant impacts that built it leave
        # (`spin.accretion_spin`, GEN.731), with no star to tide it; a giant keeps a day from the gas-giant law and a Rayleigh tilt.
        # The angle is measured from the galactic pole, the axis a rogue planet has in place of an orbit normal.
        if self.planet_type == 't':
            spin_rate, tilt = spin.accretion_spin(self.mass_kg, self.radius_km)
            self.rotation_period_hours = 2.0 * math.pi / spin_rate / 3600.0
        else:
            self.rotation_period_hours = spin.planet_day_hours(self.planet_type)
            tilt = spin.rayleigh_tilt_deg()
        spin.set_spin(self, spin.GALACTIC_POLE, tilt)

    def _apply_surface_conditions(self, conditions):
        """Sets `ROGUE_SURFACE_FIELDS`, `has_internal_heat` and the habitability score
        (`rogue_habitability.ROGUE_FIELDS`) from a `rogue_surface_conditions` result."""
        for field in ROGUE_SURFACE_FIELDS:
            setattr(self, field, conditions[field])
        self.has_internal_heat = conditions["has_internal_heat"]
        rogue_habitability.generate(self)  # the score a visitor would get on its surface; no life is drawn for a rogue

    @classmethod
    def from_dict(cls, data, system_config):
        """Reconstructs a `RoguePlanet` from a dict in the shape `to_dict()` produces, without re-running
        generation."""
        planet = cls._rebuild(data, system_config)
        if getattr(planet, "mass_bin", None) is None:
            planet.mass_bin = infer_rogue_mass_bin(planet.mass_kg)
        if "surface_regime" not in data:
            # Written before rogues had surface conditions (schema v48):
            # computed from a generator seeded by the name, so a reload
            # always gives the same conditions.
            planet._apply_surface_conditions(rogue_surface_conditions(
                planet.mass_kg, planet.radius_km, planet.planet_type, planet.mass_bin, planet.has_moons,
                draw.Stream(planet.name)))
        if "planet_class" not in data:
            planet.planet_class = rogue_class(planet)
        return planet

    @property
    def giant_kind(self):
        return self.planet_type

    @property
    def parent_kind(self):
        return None

    @property
    def kind_label(self):
        """`"Brown Dwarf"`, `"Gas Giant"` or `"Terrestrial"`."""
        if getattr(self, "mass_bin", None) == "brown-dwarf":
            return "Brown Dwarf"
        return "Gas Giant" if self.is_giant else "Terrestrial"

    @property
    def spectral_type(self):
        """A brown dwarf's L, T or Y type (`"T5"`) from its effective temperature (GEN.252), or `None` for a planet."""
        if getattr(self, "mass_bin", None) != "brown-dwarf":
            return None
        return brownDwarf.spectral_type(self.effective_temperature_k)

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
            f"{format_flux_w_m2(self.internal_heat_flux_w_m2)} leaking from its interior, so it glows at an effective "
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

def interstellar_comet_designation(sector_code, index):
    """
    An interstellar comet's designation (v40): `I/<sector>-<n>`, after
    the IAU's `I/` prefix, where `<sector>` is the sector it was found in
    (its grid designation, or its name) and `<n>` counts that sector's
    interstellar comets from 1. `I/<n>` for one found in no sector.
    """
    return f"I/{sector_code}-{index}" if sector_code else f"I/{index}"


class InterstellarComet(SmallBody):
    """
    A basic class to store information for a small icy body passing
    through on an unbound, hyperbolic interstellar trajectory (real
    confirmed examples: 1I/'Oumuamua, 2I/Borisov).

    Attributes:
        name (str): A generated or explicitly given name (real interstellar
            objects are cataloged as `1I/...`, `2I/...`, etc. -- this
            generator's own name is purely descriptive/narrative).
        nucleus_diameter_km (float): Nucleus diameter in kilometers.
        mass_kg (float): The nucleus's mass, kg: its sphere at the bulk density of its makeup (GEN.442).
        composition (list): The names of the three largest components of `makeup`.
        makeup (SmallBodyMakeup): The drawn composition (`physics/small_body_composition.py`).
        velocity_kms (float): Hyperbolic excess speed, in km/s.
        is_active (bool): Whether it currently shows a coma/tail from
            sublimating ices.
    """

    FIELD_UNITS = {
        "name": "text", "nucleus_diameter_km": "km", "mass_kg": "kg", "velocity_kms": "km_s", "is_active": "flag",
        "density_g_cm3": "g_cm3", "ice_mass_fraction": "frac", "roundness_ratio": "ratio", "yield_strength_pa": "Pa", "binding_energy_j": "J",
        **GALACTIC_MOTION_UNITS, "rotation_period_hours": "h", **spin.SPIN_UNITS,
    }
    SERIALIZABLE_FIELDS = list(FIELD_UNITS)
    """Every attribute set by `__init__`, with its native unit (GEN.570), excluding `system_config` (a
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
                (its settings).
            name (str, optional): An explicit name. Random if omitted.
        """
        super().__init__(system_config, name=name)

        self.nucleus_diameter_km = draw.uniform(*tuning.INTERSTELLAR_COMET_NUCLEUS_DIAMETER_RANGE_KM)
        self.velocity_kms = draw.uniform(*tuning.INTERSTELLAR_OBJECT_SPEED_KMS_RANGE)
        self.is_active = draw.random() < tuning.INTERSTELLAR_COMET_ACTIVE_CHANCE
        log.choice("Interstellar comet activity", self.is_active,
                   f"roll against INTERSTELLAR_COMET_ACTIVE_CHANCE "
                   f"({tuning.INTERSTELLAR_COMET_ACTIVE_CHANCE})")

        # GEN.442: an interstellar comet's ices and refractories, the wider carbon-oxide ranges of its class.
        self.makeup = small_body_composition.draw_comet(interstellar=True)
        self.composition = small_body_composition.component_names(self.makeup)
        self.mass_kg = self.makeup.mass_kg(self.nucleus_diameter_km)
        body_strength.generate_comet(self)    # GEN.1162: its density, ice share, roundness ratio, strength and binding energy

        self.set_galactic_motion()
        # GEN.104: a small body's spin, YORP-tilted from the galactic pole.
        self.rotation_period_hours = spin.small_body_period_hours()
        spin.set_spin(self, spin.GALACTIC_POLE, spin.yorp_tilt_deg())

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

from planetgen.generation import object_timing  # noqa: E402 -- registers the generators below (PERF.92)

_kg_to_mearth = units.converter("kg", "Mearth")
_kg_to_mjup = units.converter("kg", "Mjup")
_kg_to_msun = units.converter("kg", "Msun")
_mearth_to_kg = units.converter("Mearth", "kg")
_mjup_to_kg = units.converter("Mjup", "kg")
_msun_to_kg = units.converter("Msun", "kg")
object_timing.track(RoguePlanet, "rogue-planet")
object_timing.track(InterstellarComet, "interstellar-comet")
