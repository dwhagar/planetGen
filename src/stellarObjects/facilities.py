# stellarObjects/facilities.py

"""
Facilities: starbases, colonies and outposts (schema v42).

A facility sits on one host. Terrestrial facilities stand on a planet or
moon, orbital ones circle a star, planet or moon inside its sphere of
influence (`orbit_limits`), asteroid ones sit in an asteroid belt (at a
random spot, circling its star from there: `belt_position`) or field,
and stand-alone ones are parked in space inside a
sector. `tuning.FACILITY_RULES` says which kinds go where.
This module holds the rules and the orbit math; `_db.add_facility` stores
one.
"""

import math
import random

from planetgen.physics import constants
from planetgen import tuning
from planetgen.physics.planets import calculate_orbital_period_years
from .utils import circular_orbital_speed_kms

PLACEMENTS = ("terrestrial", "orbital", "asteroid", "standalone")
"""tuple: Every `facilities.placement`."""

HOST_TYPES = ("star", "planet", "moon", "asteroid_belt", "asteroid_field", "space")
"""tuple: Every `facilities.host_type`."""


def allowed_kinds(placement, host_type):
    """The facility kinds `tuning.FACILITY_RULES` allows for a
    placement on a host type (empty when the pair isn't allowed at all)."""
    return tuning.FACILITY_RULES.get((placement, host_type), ())


def check_facility(kind, placement, host_type, host_body_type=None):
    """
    Checks one facility against the placement rules.

    Args:
        kind (str): A `tuning.FACILITY_KINDS` key.
        placement (str): One of `PLACEMENTS`.
        host_type (str): One of `HOST_TYPES`.
        host_body_type (str, optional): A planet's or moon's `body_type`
            (`'t'` terrestrial, `'g'` gas giant).

    Returns:
        str or None: Why it isn't allowed, or `None` when it is.
    """
    if kind not in tuning.FACILITY_KINDS:
        return f"unknown facility kind {kind!r}"
    if placement not in PLACEMENTS:
        return f"unknown placement {placement!r}"
    if host_type not in HOST_TYPES:
        return f"unknown host type {host_type!r}"
    kinds = allowed_kinds(placement, host_type)
    if not kinds:
        return f"a {placement} facility can't be placed on a {host_type.replace('_', ' ')}"
    if kind not in kinds:
        return f"a {kind} can't be a {placement} facility on a {host_type.replace('_', ' ')} (allowed: {', '.join(kinds)})"
    if placement == "terrestrial" and host_body_type != "t":
        return "a gas giant takes orbital facilities only"
    return None


def host_placements(host_type, host_body_type=None):
    """
    The placements `tuning.FACILITY_RULES` allows on a host,
    in `PLACEMENTS` order: what the facility form lists a host under. A
    gas giant takes no surface facilities (`check_facility`).
    """
    placements = [placement for placement in PLACEMENTS if allowed_kinds(placement, host_type)]
    if host_body_type != "t" and "terrestrial" in placements:
        placements.remove("terrestrial")
    return placements


def star_host(stars, star_id, binary_configuration=None, binary_separation_km=None,
              binary_heliosphere_km=None):
    """
    What an orbital facility around a star circles: `(mass_kg, radius_km,
    sphere_km)`. A close pair is orbited as one, their combined mass,
    clear of both; a single star or one star of a wide pair on its own.
    The sphere of influence is the heliosphere (the pair's shared one for
    a close pair).

    Args:
        stars (list[dict]): Every star in the system, each with `id`,
            `mass_kg`, `radius_km` and `heliosphere_radius_km`.
        star_id (int): The host star.
        binary_configuration (str, optional): `'close'`, `'wide'` or None.
        binary_separation_km (float, optional): A close pair's separation.
        binary_heliosphere_km (float, optional): A close pair's shared
            heliosphere.
    """
    star = next(star for star in stars if star["id"] == star_id)
    if binary_configuration == "close":
        mass = sum(other["mass_kg"] or 0.0 for other in stars)
        radius = (binary_separation_km or 0.0) + max(other["radius_km"] or 0.0 for other in stars)
        return mass, radius, binary_heliosphere_km or star.get("heliosphere_radius_km")
    return star["mass_kg"], star["radius_km"], star.get("heliosphere_radius_km")


def orbit_limits(host_radius_km, sphere_km=None):
    """
    The orbits a facility may take around a host, `(lowest_km,
    highest_km)`: from just above its surface
    (`tuning.FACILITY_ORBIT_FLOOR` radii) out to the edge of
    its sphere of influence (a planet's or moon's Hill sphere, a star's
    heliosphere), or `FACILITY_ORBIT_FALLBACK_RADII` radii when it has
    none stored.
    """
    lowest = host_radius_km * tuning.FACILITY_ORBIT_FLOOR
    if sphere_km is None or not math.isfinite(sphere_km) or sphere_km <= lowest:
        sphere_km = host_radius_km * tuning.FACILITY_ORBIT_FALLBACK_RADII
    return float(lowest), float(sphere_km)


def distance_from_step(step, lowest_km, highest_km):
    """
    The orbit radius at one step of the form's slider: steps 0 to
    `tuning.FACILITY_ORBIT_STEPS` run from `lowest_km` to
    `highest_km`, each step the same factor farther out (a logarithmic
    scale). Steps outside the range are clamped.
    """
    steps = tuning.FACILITY_ORBIT_STEPS
    fraction = min(max(step, 0), steps) / steps
    return lowest_km * (highest_km / lowest_km) ** fraction


def step_for_distance(distance_km, lowest_km, highest_km):
    """The slider step nearest `distance_km` (`distance_from_step`'s
    inverse), clamped to the slider."""
    steps = tuning.FACILITY_ORBIT_STEPS
    if distance_km <= lowest_km:
        return 0
    if distance_km >= highest_km:
        return steps
    return round(steps * math.log(distance_km / lowest_km) / math.log(highest_km / lowest_km))


def belt_position(lower_km, upper_km, rng=random):
    """A spot inside an asteroid belt: `(radius_km, phase_deg)`, a random
    radius between its edges and a random angle."""
    return rng.uniform(lower_km, upper_km), rng.uniform(0.0, 360.0)


def orbit_for(host_mass_kg, host_radius_km, distance_km=None, highest_km=None):
    """
    A circular orbit around a host, the way planets and moons get theirs:
    Kepler's third law for the period (`planetPhysics.
    calculate_orbital_period_years`) and `2 pi r / T` for the speed.

    Args:
        host_mass_kg (float): The host's mass.
        host_radius_km (float): The host's radius; the orbit must clear it.
        distance_km (float, optional): The orbit's radius. Defaults to
            `tuning.FACILITY_DEFAULT_ORBIT_RADII` host radii
            (or `highest_km`, if that is nearer).
        highest_km (float, optional): The edge of the host's sphere of
            influence (`orbit_limits`); the orbit must stay inside it.

    Returns:
        dict: `distance_km`, `period_years`, `orbital_speed_kms`.

    Raises:
        ValueError: If the distance isn't a finite number above the host's
            radius and inside `highest_km`, or the host has no mass.
    """
    if distance_km is None:
        distance_km = host_radius_km * tuning.FACILITY_DEFAULT_ORBIT_RADII
        if highest_km is not None:
            distance_km = min(distance_km, highest_km)
    if not (isinstance(distance_km, (int, float)) and math.isfinite(distance_km)):
        raise ValueError("orbit distance must be a finite number")
    if distance_km <= host_radius_km:
        raise ValueError(f"orbit distance {distance_km:g} km is inside the host (radius {host_radius_km:g} km)")
    if highest_km is not None and distance_km > highest_km * (1 + 1e-9):
        raise ValueError(f"orbit distance {distance_km:g} km is outside the host's sphere of influence"
                         f" ({highest_km:g} km)")
    if not host_mass_kg or host_mass_kg <= 0:
        raise ValueError("the host has no mass to orbit")
    distance_au = distance_km / constants.AU_TO_KM
    period_years = calculate_orbital_period_years(distance_au, host_mass_kg)
    return {
        "distance_km": float(distance_km),
        "period_years": period_years,
        "orbital_speed_kms": circular_orbital_speed_kms(distance_au, period_years),
    }
