#!/usr/bin/env python3
"""
src/tests/oddity_plausibility_cli.py -- the oddity plausibility harness (TEST.490, docs/design/systems/weird-worlds.md section 5).

Forces each oddity row on a sample of bodies (200 per row, a fixed seed) and checks that every weird outcome stays physical: the six general
invariants for every body of every oddity (energy, structure, stability, air, class derived and total, no NaN), and each row's own `limits`
(`oddities.limit_failures`). Modelled on `physical_plausibility_cli.py` and `phenomena_plausibility_cli.py`; the invariant tests are in
`test_oddity_plausibility.py`, which also holds the isolation test (the dial at 0 changes no draw but the oddity's own).

How a body gets the oddity. The table's rows are wired into the generator by GEN.1047, and each effect is implemented by its own module. Until a row
is wired, `EFFECTS` holds the effect of the row, called on a gated, finished body with the row's own stream (`oddities.stream`): the body is drawn
by the generator as it is for every body, the effect changes the named input (and the stored numbers the solvers would then make), and the
invariants read the result. A module that implements a row registers its effect in `EFFECTS` (or the harness forces it through the generator
once GEN.1047 lands: `force` is the one place that changes). The rows below `STUB_ROWS` are test rows with a stub effect: the harness runs on them
and fails on the deliberately broken ones (`BROKEN_ROWS`), which shows each invariant can fail.

Usage:
    python src/tests/oddity_plausibility_cli.py
    python src/tests/oddity_plausibility_cli.py --n 50 --rows W2 W5
    python src/tests/oddity_plausibility_cli.py --stubs        # the stub and the broken rows only
"""

import argparse
import math
import sys

from planetgen import tuning
from planetgen.generation import plausibility
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.physics import constants, mantle_stripping, oddities, planet_classes, relations, rings, surface_air
from planetgen.physics.oddities import Oddity, Rate
from planetgen.physics import bulk_draw
from planetgen.util import draw
from planetgen.util.random import log_uniform

SAMPLE = 200
"""int: Bodies forced per row (design section 5)."""
SEED = "oddity-harness"
"""str: The seed the sample is drawn from, with the row id and the body's number, so a second run draws the same bodies."""
MAX_ATTEMPTS = 40
"""int: How many gated draws a body gets before the sample records that the gate is never met."""
TEMPERATURE_SLACK = 1.01
"""float: The slack on the radiative balance of an airless body (floating point and the choice of the sub-stellar flux, which is the most it gets)."""
TIDAL_ENERGY_FRACTION = 1.0
"""float: The share of the orbital energy of the planet that tidal heating over the age of the system may spend (design 5, general invariant 1)."""
VOLATILE_ALLOWANCE = 0.09
"""float: The mass share of a body its air can hold beyond its water (the CI chondrites' carbon, nitrogen and sulfur in all, Lodders; the same value as
`tests/test_planets.py`)."""


# --- the effects (stub rows; modules register their own) --------------------------------------------------------------------------

def _spin_retrograde(body, star, stream):
    """Stub: the spin axis turned over (obliquity over 90 degrees, lever L1); nothing else changes."""
    body.axial_tilt_deg = 90.0 + 90.0 * stream.random()


def _nudge_eccentricity(body, star, stream):
    """Stub: a small bounded rise of the orbit's eccentricity (lever L2), 0.05 at most."""
    body.eccentricity = min(0.95, body.eccentricity + 0.05 * stream.random())


def _broken_nan_pressure(body, star, stream):
    body.atmospheric_pressure = float("nan")


def _broken_negative_mass(body, star, stream):
    body.mass = -abs(body.mass)


def _broken_too_dense(body, star, stream):
    body.density = 40.0   # denser than any body there is


def _broken_overheated(body, star, stream):
    body.atmospheric_pressure = 0.0
    body.atmosphere = "None"
    body.gas_partials_pa = {}
    body.surface_temperature = 5.0e4   # more than the starlight, the interior and the tides could heat an airless rock to


def _broken_tides_beyond_orbit(body, star, stream):
    body.tidal_flux_w_m2 = 1.0e9   # more energy over the age than the orbit holds


def _broken_orbit_inside_star(body, star, stream):
    body.distance = body.min_orbit_distance * 0.1


def _broken_mixture(body, star, stream):
    if body.gas_partials_pa:
        first = next(iter(body.gas_partials_pa))
        body.gas_partials_pa[first] += 0.5   # the partials no longer sum to the pressure


def _broken_class(body, star, stream):
    body.planet_class = "ZZ"


def _broken_air_beyond_inventory(body, star, stream):
    body.atmospheric_pressure = max(body.atmospheric_pressure, 1.0) * 1.0e12


def _stub_ring(body, star, stream):
    """Stub: the planet gets a ring from the ring law (lever L6) with the row's own stream."""
    cold = body.distance >= star.snow_line_au
    with draw.bound(stream):
        mass = min(log_uniform(1.0e13, 1.0e17), rings.MAX_MASS_OVER_PLANET * body.mass)   # a narrow ring of a fallen moon, the stub's one law
        ring = rings.build_ring("disrupted_moon", mass, 1.0e8, body.radius, body.mass, cold)
    rings.clear(body)
    if ring is not None:
        for name, value in ring.values().items():
            setattr(body, name, value)


def _own_ring(body, star, stream):
    _stub_ring(body, star, stream)
    if body.ring_mass_kg is None:       # the law drew none: give it the narrow ring the law would, so the broken rows have a ring to break
        body.ring_kind, body.ring_origin, body.ring_ice_fraction, body.ring_age_yr = "narrow", "disrupted_moon", 0.0, 1.0e8
        body.ring_inner_km, body.ring_outer_km = 1.5 * body.radius, 2.0 * body.radius
        body.ring_mass_kg = 1.0e13
        body.ring_optical_depth = rings.optical_depth(1.0e13, body.ring_inner_km, body.ring_outer_km, 5.0e-3, 0.0, 1.0)


def _broken_heavy_ring(body, star, stream):
    _own_ring(body, star, stream)
    body.ring_mass_kg = 1.0e-2 * body.mass


def _broken_ring_beyond_roche(body, star, stream):
    _own_ring(body, star, stream)
    body.ring_outer_km = 20.0 * body.radius


def _broken_ring_optical_depth(body, star, stream):
    _own_ring(body, star, stream)
    body.ring_optical_depth = body.ring_optical_depth * 1.0e6


def _broken_ring_near_star(body, star, stream):
    _own_ring(body, star, stream)
    body.distance = 0.05


def _broken_ring_plane(body, star, stream):
    _own_ring(body, star, stream)
    body.ring_plane_tilt_deg = 30.0


def _strip_mantle(body, star, stream):
    """
    W10 (GEN.1040): a giant impact blows the outer rock off the planet (`mantle_stripping.impact`, a head-on impact with a projectile of a third to all of its mass at 20 to 80 km/s), so
    the remnant has less mass and more iron. The density follows the iron share (the grains mixed by volume, 7.8 and 3.3 g/cm3), the radius follows the mass and the density, and the
    gravity follows both.
    """
    core = body.metal_mass_fraction if body.metal_mass_fraction is not None else mantle_stripping.EARTH_LIKE_CORE_FRACTION
    with draw.bound(stream):
        projectile = body.mass * draw.uniform(0.33, 1.0)
        speed = log_uniform(2.0e4, 8.0e4)
    remnant = mantle_stripping.impact(body.mass, core, projectile, speed)
    if remnant.mass_kg <= 0.0:       # the body is dispersed: nothing to describe, so the draw is a strip that stops short of it
        remnant = mantle_stripping.impact(body.mass, core, 0.1 * body.mass, speed)
    mix = lambda iron: 1.0 / (iron / 7.8 + (1.0 - iron) / 3.3)      # noqa: E731 -- g/cm3 of the grains mixed by volume
    body.density *= mix(remnant.core_mass_fraction) / mix(core)
    body.mass = remnant.mass_kg
    body.radius = (3.0 * body.mass / (4.0 * math.pi * body.density * 1000.0)) ** (1.0 / 3.0) / 1000.0
    body.gravity = constants.G * body.mass / (body.radius * 1000.0) ** 2 / constants.EARTH_GRAVITY
    body.metal_mass_fraction = remnant.core_mass_fraction
    if body.rock_mass_fraction is not None:
        body.rock_mass_fraction = max(0.0, 1.0 - remnant.core_mass_fraction - (body.water_mass_fraction or 0.0))


EFFECTS = {
    "strip_mantle": _strip_mantle,
    "stub_ring": _stub_ring,
    "broken_heavy_ring": _broken_heavy_ring,
    "broken_ring_beyond_roche": _broken_ring_beyond_roche,
    "broken_ring_optical_depth": _broken_ring_optical_depth,
    "broken_ring_near_star": _broken_ring_near_star,
    "broken_ring_plane": _broken_ring_plane,
    "stub_spin_retrograde": _spin_retrograde,
    "stub_nudge_eccentricity": _nudge_eccentricity,
    "broken_nan_pressure": _broken_nan_pressure,
    "broken_negative_mass": _broken_negative_mass,
    "broken_too_dense": _broken_too_dense,
    "broken_overheated": _broken_overheated,
    "broken_tides_beyond_orbit": _broken_tides_beyond_orbit,
    "broken_orbit_inside_star": _broken_orbit_inside_star,
    "broken_mixture": _broken_mixture,
    "broken_class": _broken_class,
    "broken_air_beyond_inventory": _broken_air_beyond_inventory,
}
"""dict: The effect of a row by its `apply.effect` name: `effect(body, star, stream)` changes the drawn body in place."""


def _stub(identifier, name, effect, lever, field, limits, requires=None):
    """A test row: a gated terrestrial planet, a rate of 1 (the harness forces it), the limits it must keep."""
    return Oddity(identifier, name, "planet", lever, field, effect, {},
                  tuple(requires or [{"quantity": "body_type", "equals": "t"}]), tuple(limits),
                  Rate(1.0, "test row", "RECALLED", "T9", "L4", "none: a test row"), "test row", f"test row {identifier}")


STUB_ROWS = (
    _stub("S1", "stub: spin axis turned over", "stub_spin_retrograde", "L1", "obliquity_deg", [{"quantity": "obliquity_deg", "min": 90.0}]),
    _stub("S2", "stub: a little more eccentric", "stub_nudge_eccentricity", "L2", "eccentricity", [{"quantity": "eccentricity", "lt": 0.96}]),
    _stub("S3", "stub: a ring from the ring law", "stub_ring", "L6", "ring_mass", [{"quantity": "mass_earth", "min": 0.0}],
          requires=[{"quantity": "body_type", "equals": "t"}, {"quantity": "semi_major_axis_au", "min": 0.1}]),
)
"""tuple: Test rows that must pass every invariant."""

BROKEN_ROWS = (
    _stub("B1", "broken: not a number in the air", "broken_nan_pressure", "L7", "air_branch", [{"quantity": "surface_pressure_pa", "min": 0.0}]),
    _stub("B2", "broken: a negative mass", "broken_negative_mass", "L3", "water_mass_fraction", [{"quantity": "mass_earth", "min": 0.0}]),
    _stub("B3", "broken: denser than any body", "broken_too_dense", "L3", "core_mass_fraction", [{"quantity": "density_g_cm3", "lt": 30.0}]),
    _stub("B4", "broken: hotter than the energy it has", "broken_overheated", "L4", "internal_heat_flux", [{"quantity": "surface_pressure_pa", "equals": 0.0}]),
    _stub("B5", "broken: tides beyond the orbit", "broken_tides_beyond_orbit", "L4", "tidal_power", [{"quantity": "eccentricity", "lt": 1.0}]),
    _stub("B6", "broken: an orbit inside the star", "broken_orbit_inside_star", "L2", "eccentricity", [{"quantity": "eccentricity", "lt": 1.0}]),
    _stub("B7", "broken: partial pressures that do not sum to the pressure", "broken_mixture", "L7", "air_branch", [{"quantity": "surface_pressure_pa", "min": 0.0}],
          requires=[{"quantity": "surface_pressure_pa", "min": 1.0}]),
    _stub("B8", "broken: a class that does not exist", "broken_class", "L1", "obliquity_deg", [{"quantity": "obliquity_deg", "min": 0.0}]),
    _stub("B9", "broken: more air than the body holds", "broken_air_beyond_inventory", "L7", "air_branch", [{"quantity": "surface_pressure_pa", "min": 0.0}],
          requires=[{"quantity": "surface_pressure_pa", "min": 1.0}]),
)
"""tuple: Test rows with a deliberate fault. Each must be caught by the invariant named in `BROKEN_CAUGHT_BY`."""

RING_GATE = [{"quantity": "body_type", "equals": "t"}, {"quantity": "semi_major_axis_au", "min": 0.1}]
"""list: The gate of a ring row: a terrestrial planet at 0.1 AU or more (R9: a planet inside 0.1 AU keeps no ring)."""

BROKEN_RING_ROWS = (
    _stub("R91", "broken: a ring heavier than a ten-thousandth of the planet", "broken_heavy_ring", "L6", "ring_mass", [{"quantity": "mass_earth", "min": 0.0}],
          requires=RING_GATE),
    _stub("R92", "broken: a ring far beyond the Roche limit", "broken_ring_beyond_roche", "L6", "ring_mass", [{"quantity": "mass_earth", "min": 0.0}],
          requires=RING_GATE),
    _stub("R93", "broken: an optical depth that is not the formula's", "broken_ring_optical_depth", "L6", "ring_mass", [{"quantity": "mass_earth", "min": 0.0}],
          requires=RING_GATE),
    _stub("R94", "broken: a ring inside 0.1 AU", "broken_ring_near_star", "L6", "ring_mass", [{"quantity": "mass_earth", "min": 0.0}]),
    _stub("R95", "broken: a ring plane that is not the equator", "broken_ring_plane", "L6", "ring_mass", [{"quantity": "mass_earth", "min": 0.0}],
          requires=RING_GATE),
)
"""tuple: Test rows with a deliberate ring fault (TEST.494); each must be caught by the `ring` list (`RING_CAUGHT_BY`)."""

RING_CAUGHT_BY = {"R91": "ring mass", "R92": "ring beyond the Roche limit", "R93": "ring optical depth", "R94": "ring near the star", "R95": "ring plane"}
"""dict: The ring invariant that must report each broken ring row."""

BROKEN_CAUGHT_BY = {"B1": "no NaN", "B2": "no NaN", "B3": "structure", "B4": "energy", "B5": "energy", "B6": "stability", "B7": "air", "B8": "class",
                    "B9": "air"}
"""dict: The invariant that must report each broken row."""


# --- the body ---------------------------------------------------------------------------------------------------------------------

def quantities_of(body, star):
    """The named quantities (`oddities.QUANTITIES`) a body has, for gates and limits; one the body lacks is left out (a gate on it then fails)."""
    found = {"body_type": body.body_type, "mass_earth": body.mass / constants.EARTH_MASS_TO_KG,
             "radius_earth": body.radius / 6371.0, "density_g_cm3": body.density, "semi_major_axis_au": body.distance,
             "eccentricity": body.eccentricity, "orbital_period_days": body.period * 365.25, "age_gyr": star.age,
             "host_teff_k": star.temperature, "host_mass_sun": star.mass / 1.98847e30, "host_close_binary": 0,
             "surface_pressure_pa": body.atmospheric_pressure, "obliquity_deg": body.axial_tilt_deg,
             "water_mass_fraction": body.water_mass_fraction or 0.0, "envelope_mass_fraction": body.envelope_mass_fraction or 0.0,
             "spin_period_days": abs(body.rotation_period_hours) / 24.0, "moon_mass_ratio": 0.0, "collision_count": 0}
    found["tidally_locked"] = 1 if abs(found["spin_period_days"] - found["orbital_period_days"]) <= 1e-3 * found["orbital_period_days"] else 0
    if math.isfinite(body.bond_albedo) and body.distance > 0:
        found["equilibrium_temperature_k"] = relations.equilibrium_temperature(star.luminosity, body.distance * constants.AU_M, body.bond_albedo)
    if body.metal_mass_fraction is not None:
        found["core_mass_fraction"] = body.metal_mass_fraction
    return found


def _finite_numbers(body):
    for name, value in vars(body).items():
        if isinstance(value, float) and not math.isfinite(value):
            yield name
    for name, value in (body.gas_partials_pa or {}).items():
        if not math.isfinite(value):
            yield f"gas_partials_pa[{name}]"


def general_failures(body, star):
    """
    The general invariants of design section 5 for one body: `[(invariant, text)]`, empty when the body is still physical.

    1 energy -- the surface temperature is between the cosmic background and the air solve's ceiling; an airless body is no hotter than the starlight at
    its sub-stellar point, its interior and its tides could make it; tidal power over the age of the system is no more than the orbital energy.
    2 structure -- the stored gravity and density lie inside the bounds the composition and the mass give (`plausibility`).
    3 stability -- the orbit is outside the least the star allows, inside the star's reach, and the Hill sphere is inside the orbit.
    4 air -- the pressure is not negative, the weight of the air is no more than the body holds, the partial pressures sum to the pressure.
    5 class -- the class is one of the classes and is the one the stored conditions derive.
    6 no NaN -- no stored number is not finite; mass, radius and density are positive.
    """
    failures = []
    bad = sorted(set(_finite_numbers(body)))
    if bad:
        failures.append(("no NaN", f"not finite: {bad[:5]}"))
    for name in ("mass", "radius", "density", "gravity"):
        value = getattr(body, name)
        if not (isinstance(value, (int, float)) and math.isfinite(value) and value > 0.0):
            failures.append(("no NaN", f"{name}={value!r} is not a positive finite number"))
    if any(f[0] == "no NaN" for f in failures):
        return failures      # the rest read these numbers

    temperature = body.surface_temperature
    if not (constants.COSMIC_BACKGROUND_TEMPERATURE_K * 0.999 <= temperature <= surface_air.MAX_SURFACE_TEMPERATURE_K):
        failures.append(("energy", f"surface temperature {temperature:.1f} K outside {constants.COSMIC_BACKGROUND_TEMPERATURE_K} to {surface_air.MAX_SURFACE_TEMPERATURE_K} K"))
    distance_m = body.distance * constants.AU_M
    if body.atmospheric_pressure == 0.0 and body.body_type == "t":
        starlight = star.luminosity / (4.0 * math.pi * distance_m ** 2)
        budget = starlight * (1.0 - body.bond_albedo) + (body.internal_heat_flux_w_m2 or 0.0) + (body.tidal_flux_w_m2 or 0.0)
        if constants.STEFAN_BOLTZMANN_CONSTANT * temperature ** 4 > TEMPERATURE_SLACK * budget:
            failures.append(("energy", f"an airless body at {temperature:.0f} K radiates {constants.STEFAN_BOLTZMANN_CONSTANT * temperature ** 4:.1f} W/m2 of "
                                       f"a budget of {budget:.1f} W/m2"))
    tidal = body.tidal_flux_w_m2 or 0.0
    age_s = star.age * 1.0e9 * 3.15576e7
    spent = tidal * 4.0 * math.pi * (body.radius * 1000.0) ** 2 * age_s
    orbital = TIDAL_ENERGY_FRACTION * constants.G * star.mass * body.mass / (2.0 * distance_m)
    if spent > orbital:
        failures.append(("energy", f"tidal heating over {star.age:.2f} Gyr spends {spent:.2e} J of an orbital energy of {orbital:.2e} J"))

    record = plausibility._extract_record(body, "X", is_moon=False)
    failures += [("structure", text) for text in plausibility.check_hard_invariants(record)]
    if body.body_type == "t":
        low, high = bulk_draw.DENSITY_BOUNDS_G_CM3
        if not low * 0.99 <= body.density <= high * 1.01:
            failures.append(("structure", f"density {body.density:.2f} g/cm3 outside {low} to {high}"))

    if body.distance < body.min_orbit_distance:
        failures.append(("stability", f"orbit {body.distance:.4f} AU inside the least the star allows, {body.min_orbit_distance:.4f} AU"))
    if not 0.0 < body.hill_radius < body.distance * constants.AU_TO_KM:
        failures.append(("stability", f"Hill radius {body.hill_radius:.0f} km is not inside the orbit of {body.distance * constants.AU_TO_KM:.0f} km"))
    if not 0.0 <= body.eccentricity < 1.0:
        failures.append(("stability", f"eccentricity {body.eccentricity} is not a bound orbit"))

    pressure = body.atmospheric_pressure
    if pressure < 0.0:
        failures.append(("air", f"pressure {pressure} Pa is negative"))
    if pressure > 0.0 and body.body_type == "t":
        air_kg = pressure * 4.0 * math.pi * (body.radius * 1000.0) ** 2 / (body.gravity * constants.EARTH_GRAVITY)
        if air_kg > body.mass * ((body.water_mass_fraction or 0.0) + VOLATILE_ALLOWANCE):
            failures.append(("air", f"the air weighs {air_kg:.2e} kg, more than the body holds ({body.mass * ((body.water_mass_fraction or 0.0) + VOLATILE_ALLOWANCE):.2e} kg)"))
    if body.gas_partials_pa:
        total = math.fsum(body.gas_partials_pa.values())
        if abs(total - pressure) > 1e-6 * max(pressure, 1.0) or min(body.gas_partials_pa.values()) < 0.0:
            failures.append(("air", f"the partial pressures sum to {total} Pa, not the pressure of {pressure} Pa (or a partial is negative)"))

    if body.planet_class not in tuning.PLANET_CLASSES:
        failures.append(("class", f"{body.planet_class!r} is not a class"))
    elif body.planet_class != planet_classes.classify(body):
        failures.append(("class", f"class {body.planet_class} is not the class the stored conditions derive, {planet_classes.classify(body)}"))
    return failures


RING_OUTER_OVER_ROCHE = rings.RANGES["outer_over_roche"][1]
"""float: The outer edge of a ring in Roche limits of its own material that the ring law allows (Saturn's is 1.05 with the planet's volumetric radius)."""
MIN_RING_DISTANCE_AU = 0.1
"""float: A planet nearer its star than this keeps no ring (R9, [R3957])."""
RING_RESONANCE_ROWS = ("R4",)
"""tuple: The oddity ids that name a spin-orbit resonance, which may put a ring beyond the Roche limit (R4: Quaoar, Haumea)."""


def ring_failures(body, star, markers=()):
    """
    The ring list of design section 9.3 for one body: `[(invariant, text)]`, empty when the ring (or its absence) is physical. `markers` are the ids of
    the oddities on the body.

    A ring has less than a ten-thousandth of the planet's mass (`rings.MAX_MASS_OVER_PLANET`); it ends inside the Roche limit of its own material
    (to the ring law's own margin) unless the body carries a resonance row (`RING_RESONANCE_ROWS`); its optical depth is the formula's output for the
    ring's mass, edges and kind (between the most and the least filled the law draws), never a stored input; a planet inside 0.1 AU has none; the ring
    plane is the equator (no tilt of its own is stored). A body with no ring holds no ring field.
    """
    fields = {name: getattr(body, name, None) for name in rings.BODY_FIELDS}
    if fields["ring_mass_kg"] is None:
        left = sorted(name for name, value in fields.items() if value is not None)
        return [("ring fields", f"a body with no ring mass holds {left}")] if left else []
    failures = []
    if not all(isinstance(fields[n], (int, float)) and math.isfinite(fields[n]) for n in rings.BODY_FIELDS if n not in ("ring_kind", "ring_origin")):
        return [("ring fields", f"a ring number is not finite: {fields}")]
    if fields["ring_mass_kg"] > rings.MAX_MASS_OVER_PLANET * body.mass * (1.0 + 1e-9):
        failures.append(("ring mass", f"{fields['ring_mass_kg']:.2e} kg is over {rings.MAX_MASS_OVER_PLANET:g} of the planet's {body.mass:.2e} kg"))
    inner, outer = fields["ring_inner_km"], fields["ring_outer_km"]
    if not 0.0 < inner < outer:
        failures.append(("ring edges", f"inner {inner} km, outer {outer} km"))
    if inner < body.radius:
        failures.append(("ring edges", f"the inner edge {inner:.0f} km is inside the planet's radius {body.radius:.0f} km"))
    planet_density = body.mass / (4.0 / 3.0 * math.pi * (body.radius * 1000.0) ** 3)
    roche = rings.roche_limit_km(body.radius, planet_density, rings.mixture_density_kg_m3(fields["ring_ice_fraction"]))
    if outer > RING_OUTER_OVER_ROCHE * roche * (1.0 + 1e-9) and not set(markers) & set(RING_RESONANCE_ROWS):
        failures.append(("ring beyond the Roche limit", f"outer edge {outer:.0f} km is {outer / roche:.2f} Roche limits ({roche:.0f} km) with no resonance row"))
    kind = fields["ring_kind"]
    filling = rings.RANGES["filling"][kind]
    particle = rings.RANGES["particle_radius_m"][kind]
    most = rings.optical_depth(fields["ring_mass_kg"], inner, outer, filling[0], fields["ring_ice_fraction"], particle)
    least = rings.optical_depth(fields["ring_mass_kg"], inner, outer, filling[1], fields["ring_ice_fraction"], particle)
    if not least * (1.0 - 1e-6) <= fields["ring_optical_depth"] <= most * (1.0 + 1e-6):
        failures.append(("ring optical depth", f"{fields['ring_optical_depth']:.3e} is outside the formula's {least:.3e} to {most:.3e} for this ring"))
    if body.distance < MIN_RING_DISTANCE_AU and not getattr(body, "is_moon", False):
        failures.append(("ring near the star", f"a ring at {body.distance:.3f} AU, inside {MIN_RING_DISTANCE_AU} AU"))
    if getattr(body, "ring_plane_tilt_deg", 0.0):
        failures.append(("ring plane", f"the ring plane is {body.ring_plane_tilt_deg} degrees from the equator"))
    return failures


# --- forcing a row on a sample ----------------------------------------------------------------------------------------------------

def draw_body(row, number):
    """`(body, star)` for sample number `number` of `row`: a gated body drawn by the generator, seeded from the row and the number; `None` when no
    draw in `MAX_ATTEMPTS` meets the gate."""
    for attempt in range(MAX_ATTEMPTS):
        with draw.bound(f"{SEED}:{row.id}:{number}:{attempt}"):
            star_type = draw.choice(plausibility.HOST_STAR_TYPES)
            config = SystemConfig()
            config.STAR_TYPE = star_type
            star = Star(config)
            zone = draw.choice(plausibility.ZONE_CHARS)
            body = Planet(config, star, star.habitable_zone, plausibility.distance_for_zone(star, zone), moon_count=0)
        if row in oddities.candidates(row.kind, quantities_of(body, star), rows=(row,)):
            return body, star
    return None


def force(row, number):
    """The body of sample `number` with `row` forced on it (its effect called with the row's own stream), and its star; `None` when the gate is never met.
    GEN.1047's wiring will force the row through the generator here instead."""
    drawn = draw_body(row, number)
    if drawn is None:
        return None
    body, star = drawn
    EFFECTS[row.effect](body, star, oddities.stream(row.id, f"{SEED}:{number}"))
    return body, star


def run_row(row, n=SAMPLE):
    """`{"row", "n", "gated", "general": {invariant: bodies failing}, "limits": bodies failing a limit, "examples": [text]}` for `n` forced bodies of `row`."""
    general, limits, examples, gated = {}, 0, [], 0
    for number in range(n):
        forced = force(row, number)
        if forced is None:
            continue
        gated += 1
        body, star = forced
        for invariant, text in general_failures(body, star) + ring_failures(body, star, markers=(row.id,)):
            general[invariant] = general.get(invariant, 0) + 1
            if len(examples) < 5:
                examples.append(f"{row.id} #{number} {invariant}: {text}")
        found = oddities.limit_failures(row, quantities_of(body, star))
        if found:
            limits += 1
            if len(examples) < 5:
                examples.append(f"{row.id} #{number} limit: {found}")
    return {"row": row.id, "n": n, "gated": gated, "general": general, "limits": limits, "examples": examples}


def format_report(results):
    lines = ["| row | bodies | gated | general invariants failing | limits failing |", "|---|---|---|---|---|"]
    for r in results:
        failing = ", ".join(f"{k} {v}" for k, v in sorted(r["general"].items())) or "none"
        lines.append(f"| {r['row']} | {r['n']} | {r['gated']} | {failing} | {r['limits']} |")
    for r in results:
        lines += [f"  {text}" for text in r["examples"]]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Oddity plausibility harness (TEST.490).")
    parser.add_argument("--n", type=int, default=SAMPLE, help=f"Bodies forced per row (default {SAMPLE}).")
    parser.add_argument("--rows", nargs="+", default=None, metavar="ID", help="Only these row ids (default: every shipped row).")
    parser.add_argument("--stubs", action="store_true", help="Run the stub rows and the deliberately broken rows instead of the shipped table.")
    args = parser.parse_args(argv)
    rows = STUB_ROWS + BROKEN_ROWS + BROKEN_RING_ROWS if args.stubs else tuple(oddities.ROWS)
    if args.rows:
        rows = tuple(r for r in rows if r.id in args.rows)
    unwired = [r.id for r in rows if r.effect not in EFFECTS]
    if unwired:
        print(f"No effect registered for {unwired}: the module that implements the row adds it to EFFECTS.", file=sys.stderr)
        return 2
    results = [run_row(row, args.n) for row in rows]
    print(format_report(results))
    clean = [r for r in results if r["row"] not in BROKEN_CAUGHT_BY and r["row"] not in RING_CAUGHT_BY]
    return 1 if any(r["general"] or r["limits"] for r in clean) else 0


if __name__ == "__main__":
    sys.exit(main())
