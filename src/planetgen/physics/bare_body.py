# planetgen/physics/bare_body.py

"""
What a body has before its air is solved: the stand-ins of the backward edges of the forward build (GEN.1548, built ahead by Lane 07)
====================================================================================================================================

The forward build makes the bare bodies first (mass, bulk composition, radius, birth orbit) and the air, the water and the surface once, on the final orbit.
Three places of the old order read the finished air of a body before there was one, and need a number of their own at the bare stage:

* `moon_orbit_bounds_km` keeps a moon 15 atmospheric scale heights clear of its planet: `scale_height_km` is that height from the bare body (its light with the
  middle of the albedo it will draw from, its gravity from mass and radius and the mean molar mass of the gas it could hold).
* The tide decay and the interior heat measure a body against its surface temperature, with `interior_heat.SURFACE_TEMPERATURE_K` (273 K) when it has none yet:
  `surface_k` is the body's own once it has a climate and else the temperature its light gives it with no air.
* `moon_history` reads the planet's surface pressure through the thermal tide. In the forward order the planet's air comes first and the moon history after
  it, so this edge needs no stand-in; `thermal_tide_is_known` says whether the planet has the pressure the tide reads, for a caller that must wait.

All three are pure functions of what a bare body stores, draw nothing and move no stream, so a body made again with the same inputs gets the same numbers.
Nothing here is called by the old order; Lane 05 switches the call sites when the build is reordered.
"""

import math

from planetgen.physics import bulk_structure, climate_draws, constants, planet_type, radiation_field

BARE_ROCKY_MOLAR_MASS_KG_MOL = 28.0e-3
"""float: The mean molar mass, kg/mol, of the gas a bare rocky body could hold (nitrogen's 28 g/mol; the Earth's air is 29, Venus's 43, Mars's 44, Titan's 28), for the
scale height of its air before the air is solved. Stand-in BB1, tier T6, tweak L0; it sets only the margin that keeps a moon clear of the air."""


def starlight_w_m2(planet, star_distance_au=None):
    """The starlight, W/m2, at `star_distance_au` (the planet's own distance when none is given; a moon passes its planet's) with the host's light and every other star's,
    as `planets.calculate_atmospheric_conditions` takes it."""
    distance = planet.distance if star_distance_au is None else star_distance_au
    return radiation_field.starlight_w_m2(planet.star, distance) * (getattr(planet, "light_factor", None) or 1.0)


def reference_albedo(planet, star_distance_au=None):
    """The middle of the Bond albedo range the body will draw from (`climate_draws.albedo_range`), which needs no air; no draw is made."""
    flux = starlight_w_m2(planet, star_distance_au)
    low, high = climate_draws.albedo_range(planet.body_type, getattr(planet, "water_mass_fraction", None), flux, getattr(planet.star, "temperature", None))
    return 0.5 * (low + high)


def equilibrium_k(planet, star_distance_au=None):
    """The temperature, K, the starlight gives the bare body with no air and the middle albedo: the stand-in for the surface temperature before there is one."""
    return climate_draws.equilibrium_temperature_k(reference_albedo(planet, star_distance_au), starlight_w_m2(planet, star_distance_au))


def surface_k(body, star_distance_au=None):
    """
    The surface temperature a body's interior is measured against: its own once its climate is known, else `equilibrium_k`. It replaces the 273 K that
    `interior_heat.surface_of` falls back to, which is wrong by hundreds of kelvin for a body far from or close to its star.
    """
    own = getattr(body, "surface_temperature", None)
    if own and own > 0.0:
        return own
    return equilibrium_k(body, star_distance_au)


def gravity_m_s2(planet):
    """The surface gravity, m/s2, from mass and radius (both are facts of the bare stage)."""
    radius_m = planet.radius * 1000.0
    return constants.G * planet.mass / (radius_m * radius_m)


def scale_height_km(planet, star_distance_au=None):
    """
    The atmospheric scale height, km, a bare body would have: `R T / (mu g)` at `equilibrium_k`, with the molar mass of a hydrogen-helium envelope for a giant or a
    body under an envelope (`bulk_structure.ENVELOPE_MEAN_MOLAR_MASS`) and `BARE_ROCKY_MOLAR_MASS_KG_MOL` otherwise. A body with no radius or mass gives 0.
    """
    if not (getattr(planet, "radius", None) and getattr(planet, "mass", None)):
        return 0.0
    envelope = not planet_type.has_ground(planet)
    molar = bulk_structure.ENVELOPE_MEAN_MOLAR_MASS if envelope else BARE_ROCKY_MOLAR_MASS_KG_MOL
    return constants.R * equilibrium_k(planet, star_distance_au) / (molar * gravity_m_s2(planet)) / 1000.0


def thermal_tide_is_known(planet):
    """Whether the planet has a surface pressure the thermal tide of `body_spin.thermal_tide_of` can read (its air is solved); False before, and the caller waits."""
    pressure = getattr(planet, "atmospheric_pressure", None)
    return bool(pressure) and pressure > 0.0
