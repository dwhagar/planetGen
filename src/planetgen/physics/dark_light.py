# planetgen/physics/dark_light.py

"""
The light of the dark members of an orbit tree (GEN.1501)
=========================================================

The members of a system's tree that `Star` cannot be built as (a brown dwarf, a neutron star, a black hole, GEN.1013) hang in the tree with a mass and a radius and, until now, gave no light
to the planets of the stars beside them. They glow, and this module says how much, so that `radiation_field` can add them as sources like any other and the temperature, the sky and the
dose of a planet near one count them. Nothing is drawn here: every figure comes from the member's mass, radius and age, and for a black hole from the stars that feed it.

* **A brown dwarf** cools: its luminosity is `brown_dwarf.state(mass, age)` (GEN.252), the same law that gives a rogue one its type and its glow, at the age of the system it formed with.
* **A neutron star** cools by thermal emission from its surface: `L = 4 pi R^2 sigma T_s^4` with the surface temperature of a cooling curve of its age since the collapse (`surface_temperature_k`), the
  system's age less the life of the star it was (`NEUTRON_STAR_PROGENITOR_LIFE_GY`). The curve is two power laws through the end of the neutrino era (`NEUTRON_STAR_PHOTON_ERA_YR`,
  `NEUTRON_STAR_PHOTON_ERA_T_K`). A planet only sees it when the star is young or the planet near.
* **A black hole** shines only while it is fed. The stars of the tree feed it by their winds (`stellar_wind.wind`, the same law that sizes a star's bubble): the Bondi-Hoyle-Lyttleton capture
  `Mdot = (G M)^2 Mdot_wind / (a^2 v_w (v_w^2 + v_orb^2)^1.5)` of a star at the separation `a` from it, capped by the wind itself, turned to light at `DISK_EFFICIENCY` and capped by the
  Eddington limit. The disk's colour is that of a blackbody at the innermost stable orbit (`3 R_s`). A black hole with no wind to take gives nothing (its Hawking glow is a part in 10^20).
  Feeding by overflow of a star that fills its Roche lobe is not modelled; that needs the star's evolution, and its absence leaves the light of such a black hole at the wind's.

Stand-ins are the rows DL1 to DL7 of section 6.259 of the gap plan.
"""

import math
from dataclasses import dataclass

from planetgen.physics import brown_dwarf, constants, stellar_wind
from planetgen.physics import unit_registry as units

_au_to_m = units.converter("au", "m")
_m_to_km = units.converter("m", "km")
_msun_to_kg = units.converter("Msun", "kg")
_kg_to_msun = units.converter("kg", "Msun")

NEUTRON_STAR_PHOTON_ERA_YR = 1.0e5
"""float: The age at which a neutron star's cooling passes from the neutrino era (the core cools by neutrinos) to the photon era (the surface radiates the heat the core carries), years. [DL1, T7, L1.]"""

NEUTRON_STAR_PHOTON_ERA_T_K = 1.0e6
"""float: The surface temperature (at the star, not redshifted) at the end of the neutrino era, K. The standard and minimal cooling curves ([R3590], which tabulates them) put it at 5e5 to 2e6 K for a 1.4 solar mass star. [DL1, T7, L1.]"""

NEUTRON_STAR_PHOTON_ERA_SLOPE = -0.5
"""float: `d ln T_s / d ln t` in the photon era (the same curves fall from 1e6 K at 1e5 yr to about 1e5 K at 1e7 yr and a few 1e3 K at 1e9 yr, a slope of -0.5). [DL2, T7, L1.]"""

NEUTRON_STAR_NEUTRINO_ERA_SLOPE = -0.1
"""float: `d ln T_s / d ln t` before the photon era: the surface stays warm while the core cools by neutrinos. [DL2, T8, L0.]"""

NEUTRON_STAR_PROGENITOR_LIFE_GY = 0.03
"""float: How long the star that became the neutron star lived, Gy: a core-collapse progenitor of 8 to 25 solar masses lives 0.01 to 0.04 Gy, so the neutron star is younger than its system by this much. [DL3, T8, L0.]"""

NEUTRON_STAR_MIN_AGE_GY = 1.0e-6
"""float: The youngest age the cooling curve is read at (a star of the system's own age less its progenitor's life, never under this), Gy."""

DISK_EFFICIENCY = 0.1
"""float: The share of the rest mass of what a black hole takes in that leaves as light: 0.057 for a non-rotating hole's thin disk [R4037], up to 0.42 for a maximally rotating one. [DL4, T4, L0.]"""

ISCO_SCHWARZSCHILD_RADII = 3.0
"""float: The innermost stable circular orbit of a non-rotating black hole, in Schwarzschild radii (`6 G M / c^2`), which the disk's colour and size are taken at. [DL5, T2, L0.]"""

EDDINGTON_W_PER_SOLAR_MASS = 1.26e31
"""float: The Eddington luminosity per solar mass, W (`4 pi G M m_p c / sigma_T`). [DL6, T1, L0.]"""


@dataclass
class DarkSource:
    """
    A dark member as a source of the radiation field: what `radiation_field`'s bands read of a star (`luminosity` W, `temperature` K, `radius` km, `mass` kg), with the `kind`, the `type`
    string and the `state` the sky text names it by. It gives bolometric light only; the bands it has no value for (XUV, flares, the cumulative XUV index) are skipped, not guessed.
    """

    member_no: int
    kind: str
    mass: float
    radius: float
    luminosity: float
    temperature: float
    age: float
    type: str
    state: str
    name: str = ""

    yerkes_class = "DARK"


def surface_temperature_k(age_gy):
    """The surface temperature, K, of a neutron star `age_gy` Gy after its collapse (see the module docstring); two power laws through (`NEUTRON_STAR_PHOTON_ERA_YR`, `NEUTRON_STAR_PHOTON_ERA_T_K`)."""
    age_yr = max(age_gy, NEUTRON_STAR_MIN_AGE_GY) * 1.0e9
    slope = NEUTRON_STAR_PHOTON_ERA_SLOPE if age_yr >= NEUTRON_STAR_PHOTON_ERA_YR else NEUTRON_STAR_NEUTRINO_ERA_SLOPE
    return NEUTRON_STAR_PHOTON_ERA_T_K * (age_yr / NEUTRON_STAR_PHOTON_ERA_YR) ** slope


def neutron_star_luminosity_w(radius_km, age_gy):
    """A neutron star's thermal luminosity, W, `4 pi R^2 sigma T_s^4`, at `age_gy` Gy since its collapse."""
    radius_m = radius_km * 1.0e3
    return constants.STEFAN_BOLTZMANN_CONSTANT * 4.0 * math.pi * radius_m ** 2 * surface_temperature_k(age_gy) ** 4


def wind_capture_kg_s(black_hole_kg, separation_au, wind_kg_s, wind_m_s, star_kg):
    """
    What a black hole of `black_hole_kg` takes of the wind (`wind_kg_s`, speed `wind_m_s`) of a star of `star_kg` at `separation_au` from it, kg/s: the Bondi-Hoyle-Lyttleton rate in the star's wind
    (`rho = Mdot_w / (4 pi a^2 v_w)`) moving past the hole at `v_rel^2 = v_w^2 + v_orb^2`, never more than the wind itself.
    """
    if not (black_hole_kg > 0 and separation_au > 0 and wind_kg_s > 0 and wind_m_s > 0):
        return 0.0
    separation = _au_to_m(separation_au)
    v_orb_sq = constants.G * (star_kg + black_hole_kg) / separation
    v_rel_cubed = (wind_m_s ** 2 + v_orb_sq) ** 1.5
    captured = (constants.G * black_hole_kg) ** 2 * wind_kg_s / (separation ** 2 * wind_m_s * v_rel_cubed)
    return min(captured, wind_kg_s)


def disk_luminosity_w(black_hole_kg, captured_kg_s):
    """The light of the disk a black hole makes of `captured_kg_s` kg/s: `DISK_EFFICIENCY * Mdot * c^2`, never over the Eddington limit."""
    light = DISK_EFFICIENCY * captured_kg_s * constants.SPEED_OF_LIGHT_M_S ** 2
    return min(light, EDDINGTON_W_PER_SOLAR_MASS * _kg_to_msun(black_hole_kg))


def disk_radius_km(black_hole_kg):
    """The radius the disk's light is taken to come from, km: the innermost stable orbit, `ISCO_SCHWARZSCHILD_RADII` Schwarzschild radii."""
    return _m_to_km(ISCO_SCHWARZSCHILD_RADII * 2.0 * constants.G * black_hole_kg / constants.SPEED_OF_LIGHT_M_S ** 2)


def disk_temperature_k(luminosity_w, radius_km):
    """The blackbody temperature, K, of `luminosity_w` radiated from a disk of `radius_km` (both faces: `4 pi R^2`)."""
    area = 4.0 * math.pi * (radius_km * 1.0e3) ** 2
    return (luminosity_w / (area * constants.STEFAN_BOLTZMANN_CONSTANT)) ** 0.25


def feeders_of(stars, separations_au):
    """`[(star, separation_au)]` for the stars that may feed a black hole: those with a place in the tree (`separations_au[index]` is the separation from the hole, `None` for no place)."""
    return [(star, separations_au[index]) for index, star in enumerate(stars) if index < len(separations_au) and separations_au[index]]


def black_hole_captured_kg_s(black_hole_kg, feeders):
    """What a black hole takes in from the winds of `feeders` (`[(star, separation_au)]`), kg/s, summed."""
    total = 0.0
    for star, separation_au in feeders:
        luminosity, radius = getattr(star, "luminosity", None), getattr(star, "radius", None)
        if not (luminosity and radius and getattr(star, "mass", None)):
            continue
        try:
            wind_kg_s, wind_m_s = stellar_wind.wind(star.mass, luminosity, radius, getattr(star, "type", "") or "", getattr(star, "yerkes_class", "") or "")
        except (ValueError, ZeroDivisionError, TypeError):
            continue
        total += wind_capture_kg_s(black_hole_kg, separation_au, wind_kg_s, wind_m_s, star.mass)
    return total


def source_of(member, system_age_gy, feeders=()):
    """
    The `DarkSource` a dark member (`stored_trees.DarkMember`) is, or `None` when it gives no light (a kind this module does not light: a stripped helium star, a black hole nothing
    feeds, a member with no radius or mass).

    Args:
        member (DarkMember): The member: `member_no`, `kind`, `mass_solar`, `radius_au`.
        system_age_gy (float): The age of the system, Gy; a brown dwarf is as old as the system, a neutron star younger by its progenitor's life.
        feeders (iterable): `[(star, separation_au)]`, the stars that may feed a black hole by their winds (`feeders_of`).
    """
    if not member.mass_solar or system_age_gy is None or system_age_gy <= 0:
        return None
    mass = _msun_to_kg(member.mass_solar)
    radius_km = _m_to_km(_au_to_m(member.radius_au)) if member.radius_au else None
    if member.kind == "brown-dwarf":
        try:
            state = brown_dwarf.state(member.mass_solar, system_age_gy)
        except (ValueError, KeyError):
            return None
        radius = state["radius_sol"] * constants.SOLAR_RADIUS_M / 1.0e3
        return DarkSource(member.member_no, member.kind, mass, radius, state["luminosity_sol"] * constants.SOLAR_LUMINOSITY, state["temperature_k"], system_age_gy,
                          f"{state['spectral_type']} dwarf", "brown-dwarf", "brown dwarf")
    if member.kind == "neutron-star":
        if not radius_km:
            return None
        age = max(system_age_gy - NEUTRON_STAR_PROGENITOR_LIFE_GY, NEUTRON_STAR_MIN_AGE_GY)
        return DarkSource(member.member_no, member.kind, mass, radius_km, neutron_star_luminosity_w(radius_km, age), surface_temperature_k(age), age, "Neutron Star", "neutron-star",
                          "neutron star")
    if member.kind == "black-hole":
        captured = black_hole_captured_kg_s(mass, feeders)
        light = disk_luminosity_w(mass, captured)
        if light <= 0.0:
            return None
        disk = disk_radius_km(mass)
        return DarkSource(member.member_no, member.kind, mass, disk, light, disk_temperature_k(light, disk), system_age_gy, "Black Hole with an accretion disk", "black-hole-disk",
                          "black hole with an accretion disk")
    return None
