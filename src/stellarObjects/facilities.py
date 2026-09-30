# stellarObjects/facilities.py

"""
Facilities: starbases, colonies and outposts (schema v42).

A facility sits on one host. Terrestrial facilities stand on a planet or
moon, orbital ones circle a star, planet or moon, asteroid ones sit in an
asteroid belt or field, and stand-alone ones are parked in space inside a
sector. `program_constants.FACILITY_RULES` says which kinds go where.
This module holds the rules and the orbit math; `_db.add_facility` stores
one.
"""

import math

from . import physical_constants, program_constants
from .planetPhysics import calculate_orbital_period_years
from .utils import circular_orbital_speed_kms

PLACEMENTS = ("terrestrial", "orbital", "asteroid", "standalone")
"""tuple: Every `facilities.placement`."""

HOST_TYPES = ("star", "planet", "moon", "asteroid_belt", "asteroid_field", "space")
"""tuple: Every `facilities.host_type`."""


def allowed_kinds(placement, host_type):
    """The facility kinds `program_constants.FACILITY_RULES` allows for a
    placement on a host type (empty when the pair isn't allowed at all)."""
    return program_constants.FACILITY_RULES.get((placement, host_type), ())


def check_facility(kind, placement, host_type, host_body_type=None):
    """
    Checks one facility against the placement rules.

    Args:
        kind (str): A `program_constants.FACILITY_KINDS` key.
        placement (str): One of `PLACEMENTS`.
        host_type (str): One of `HOST_TYPES`.
        host_body_type (str, optional): A planet's or moon's `body_type`
            (`'t'` terrestrial, `'g'` gas giant).

    Returns:
        str or None: Why it isn't allowed, or `None` when it is.
    """
    if kind not in program_constants.FACILITY_KINDS:
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


def orbit_for(host_mass_kg, host_radius_km, distance_km=None):
    """
    A circular orbit around a host, the way planets and moons get theirs:
    Kepler's third law for the period (`planetPhysics.
    calculate_orbital_period_years`) and `2 pi r / T` for the speed.

    Args:
        host_mass_kg (float): The host's mass.
        host_radius_km (float): The host's radius; the orbit must clear it.
        distance_km (float, optional): The orbit's radius. Defaults to
            `program_constants.FACILITY_DEFAULT_ORBIT_RADII` host radii.

    Returns:
        dict: `distance_km`, `period_years`, `orbital_speed_kms`.

    Raises:
        ValueError: If the distance isn't a finite number above the host's
            radius, or the host has no mass.
    """
    if distance_km is None:
        distance_km = host_radius_km * program_constants.FACILITY_DEFAULT_ORBIT_RADII
    if not (isinstance(distance_km, (int, float)) and math.isfinite(distance_km)):
        raise ValueError("orbit distance must be a finite number")
    if distance_km <= host_radius_km:
        raise ValueError(f"orbit distance {distance_km:g} km is inside the host (radius {host_radius_km:g} km)")
    if not host_mass_kg or host_mass_kg <= 0:
        raise ValueError("the host has no mass to orbit")
    distance_au = distance_km / physical_constants.AU_TO_KM
    period_years = calculate_orbital_period_years(distance_au, host_mass_kg)
    return {
        "distance_km": float(distance_km),
        "period_years": period_years,
        "orbital_speed_kms": circular_orbital_speed_kms(distance_au, period_years),
    }
