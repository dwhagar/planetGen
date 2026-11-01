# planetgen/physics/rings.py

"""
Planetary rings (GEN.732, docs/design/composition/planetary-rings.md)
=====================================================================

A planet may hold one ring system, summed over its rings: a kind, an origin, an inner and an outer radius, a mass, a mean optical
depth and an ice share. The painter (MAP.275) reads the stored fields; nothing here draws a picture.

What the physics fixes
----------------------
* **Where it can be.** Material inside the Roche limit cannot gather into a moon, so a ring ends at it: `roche_limit_km` gives the
  fluid-body limit `2.456 R_p (rho_p / rho_m)^(1/3)` for a particle (or a clump) of density `rho_m` around a planet of radius `R_p`
  and mean density `rho_p`. The four giants of the Solar System end within a few percent of it (`tests`), with their own particle
  densities: Saturn's A ring (136 775 km) against 130 000 km for ice, Jupiter's main ring (129 000 km) against 139 000 km for rock,
  Uranus's epsilon ring (51 149 km) against 67 500 km and Neptune's Adams ring (62 933 km) against 71 300 km. A ring also ends inside the innermost
  moon that is large enough to hold it back (`MOON_CLEARANCE` of that moon's distance).
* **How thick it is.** For particles of radius `r` and density `rho_m` a ring of mass `M` over the area `A` it occupies has the mean
  optical depth `tau = 3 M / (4 rho_m r A)`; the optical depth is an output of the drawn mass, the extent and the particle size
  (`optical_depth`), not a draw of its own. With `r` = 1 m Saturn's rings (1.5e19 kg over 4.4e16 m2) come out at 0.3 and Uranus's
  rings (1.5e16 kg over 0.3 percent of their annulus) at 0.6.
* **What it is made of.** A ring beyond the snow line that has lost nothing to heat keeps its ice (Saturn's is 90 to 99 percent
  water ice); inside the snow line the ice is gone, and the ring is rock and dust. The density of the mixture is volume additive.

What is drawn (every number is a labelled stand-in; the registry rows are BC20 to BC23 of the gap plan and RO1 to RO8, RR1 to RR4)
--------------------------------------------------------------------------------------------------------------------------------
1. **Whether, and where it came from (GEN.1057).** A ring follows from history, not from a chance: `ring_origin.ring_origin` sums the rings of the moons the tide pulled inside
   the Roche limit (`moon_history`, the planet's `fallen_moons`) and of the comets that passed inside it (`ring_origin.draw_comet_passages`, drawn here from the unit's stream),
   each spread for the time since it came to be. A planet with neither has no ring, a planet inside 0.1 AU has none (R9), and a ring is the origin of the event that holds most
   of it (`disrupted_moon` or `captured_body`) and as old as that event (`ring_age_yr`). A dwarf planet with no such ring may still hold the narrow ring of its spin-orbit
   resonance (`small_body_rings`, R4).
2. **Which kind.** From the ring's mass (`kind_of_mass`): a faint one is `dusty`, a middling one `narrow`, a massive one `broad` (Jupiter, Uranus and Neptune, Saturn).
3. **Ice share, extent.** See `RANGES`.

The draw order is the order above and is part of the generator epoch.

Use:

    ring = build_ring(origin, mass_kg, age_yr, radius_km, planet_mass_kg, cold, moon_distances_km)   # a Ring or None
    apply(planet, host)                                                                              # sets the stored fields
"""

import math

from planetgen.physics import condensate_shares, constants as physical_constants, moon_composition, unit_registry as units
from planetgen.util import draw
from planetgen.util.random import log_uniform

_au_to_km = units.converter("au", "km")

KINDS = ("dusty", "narrow", "broad")
ORIGINS = ("moon_erosion", "disrupted_moon", "captured_body")

BODY_FIELDS = ("ring_kind", "ring_origin", "ring_inner_km", "ring_outer_km", "ring_mass_kg", "ring_optical_depth",
               "ring_ice_fraction", "ring_age_yr")
"""tuple: The attributes (and `planets`/`moons` columns, GEN.732) of a ring system; all `None` on a body without one."""

UNITS = {"ring_kind": "text", "ring_origin": "text", "ring_inner_km": "km", "ring_outer_km": "km", "ring_mass_kg": "kg",
         "ring_optical_depth": "ratio", "ring_ice_fraction": "frac", "ring_age_yr": "yr"}
"""dict: Each field of `BODY_FIELDS` and its native unit (GEN.570), as `RoguePlanet.FIELD_UNITS` carries it."""

ROCHE_COEFFICIENT = 2.456
"""float: The fluid Roche limit in planet radii for particles as dense as the planet (Roche 1850; [R3220]). Recalled, T2, L0."""

ICE_DENSITY_KG_M3 = 917.0
"""float: Water ice Ih, the density of a ring's ice (the repo's ice density, `constants`-free so the module stands alone)."""
ROCK_DENSITY_KG_M3 = 2500.0
"""float: The density of a ring's silicate or dark rock; a stand-in between chondritic grains (3300) and porous dust (1000 to 2000), T6, L0."""

KIND_MASS_LIMITS_KG = (1.0e12, 3.0e17)
"""tuple: A ring lighter than the first is `dusty`, one lighter than the second `narrow`, else `broad`: the middles (in log mass) of the gaps between the masses of the four Solar System
kinds (Jupiter's dust 1e7 to 1e9 kg, Uranus's rings about 1e16, Saturn's 1.5e19). Source tier T6, tweak level L3."""

RANGES = {
    "particle_radius_m": {"dusty": 1.0e-5, "narrow": 1.0, "broad": 1.0},
    "filling": {"dusty": (1.0, 1.0), "narrow": (1.0e-3, 1.0e-2), "broad": (1.0, 1.0)},
    "ice_fraction": {"dusty": (0.0, 0.0), "narrow": (0.0, 0.3), "broad": (0.90, 0.99)},
    "outer_over_roche": (0.75, 1.1),
    "inner_over_radius": (1.1, 1.7),
}
"""dict: The drawn ranges (stand-ins, T3 to T6, L0).

- `particle_radius_m`: the effective particle radius that turns a mass into an optical depth (Saturn's particles run from centimetres to
  ten metres; dust grains are about 10 micrometres).
- `filling`: the share of the annulus a ring covers, log-uniform; Uranus's rings cover about 0.3 percent of theirs.
- `ice_fraction`: the ice share by mass beyond the snow line (Saturn's rings 90 to 99 percent, [R0970], recalled; Uranus's and
  Neptune's are dark, no ice has been seen); 0 inside the snow line.
- `outer_over_roche`: the outer edge as a share of the Roche limit of the ring's own material: Saturn 1.05 (with the planet's volumetric radius), Jupiter 0.93, Neptune 0.88,
  Uranus 0.76 (COMPUTED from the edges listed above).
- `inner_over_radius`: the inner edge in planet radii: Saturn's D ring 1.11, Jupiter's halo 1.29, Uranus's 6 ring 1.64, Neptune's Galle
  ring 1.7; a drawn inner edge is held below the outer one.
"""

MOON_CLEARANCE = 0.95
"""float: A ring ends inside this share of the innermost moon's distance (a moon holds back the material beyond it; a parametrisation)."""
MAX_MASS_OVER_PLANET = 1.0e-4
"""float: A ring never holds more than this share of its planet's mass (Saturn's is 3e-8; a parametrisation)."""
MIN_WIDTH_RATIO = 1.05
"""float: A ring whose outer edge is not this many times its inner edge is not drawn (the moon clearance or a dense planet closed it)."""


class Ring:
    """
    One planet's ring system, summed.

    Attributes:
        kind (str): `KINDS`.
        origin (str): `ORIGINS`.
        inner_km, outer_km (float): The radii of its inner and outer edges, km from the planet's centre.
        mass_kg (float): Its mass.
        optical_depth (float): Its mean optical depth over the area it covers.
        ice_fraction (float): Its water-ice share by mass.
    """

    __slots__ = BODY_FIELDS

    def __init__(self, **values):
        for name in BODY_FIELDS:
            setattr(self, name, values[name])

    def values(self):
        """The attributes a body holds for this ring."""
        return {name: getattr(self, name) for name in BODY_FIELDS}


def roche_limit_km(radius_km, planet_density_kg_m3, particle_density_kg_m3):
    """The fluid Roche limit of a planet for material of `particle_density_kg_m3`, km from the planet's centre."""
    return ROCHE_COEFFICIENT * radius_km * (planet_density_kg_m3 / particle_density_kg_m3) ** (1.0 / 3.0)


def mixture_density_kg_m3(ice_fraction):
    """The density of a ring of this ice share by mass, volume additive."""
    return 1.0 / (ice_fraction / ICE_DENSITY_KG_M3 + (1.0 - ice_fraction) / ROCK_DENSITY_KG_M3)


def optical_depth(mass_kg, inner_km, outer_km, filling, ice_fraction, particle_radius_m):
    """
    The mean optical depth of a ring over the area it covers, `3 M / (4 rho r A f)` for particles of radius `particle_radius_m`
    and the mass density of the mixture, where `A` is the annulus between the edges and `f` the share of it covered.
    """
    area_m2 = math.pi * ((outer_km * 1000.0) ** 2 - (inner_km * 1000.0) ** 2) * filling
    return 3.0 * mass_kg / (4.0 * mixture_density_kg_m3(ice_fraction) * particle_radius_m * area_m2)


def kind_of_mass(mass_kg):
    """The kind of a ring of this mass: `dusty`, `narrow` or `broad` by `KIND_MASS_LIMITS_KG`."""
    dusty_below, narrow_below = KIND_MASS_LIMITS_KG
    return "dusty" if mass_kg < dusty_below else "narrow" if mass_kg < narrow_below else "broad"


def build_ring(origin, mass_kg, age_yr, radius_km, planet_mass_kg, cold, moon_distances_km=()):
    """
    The ring system of a planet whose history made a ring of `mass_kg` that came to be `age_yr` ago (`ring_origin.ring_origin`), or `None` when its edges do not fit between the planet
    and its innermost moon. The kind follows from the mass; the ice share, the outer edge and the inner edge are drawn in that order, the filling after the ice share.

    Args:
        origin (str): `ORIGINS`.
        mass_kg (float): The ring's mass now.
        age_yr (float): The years since its main event.
        radius_km, planet_mass_kg (float): The planet's radius and mass.
        cold (bool): Whether the planet orbits beyond its star's snow line, where a ring keeps its ice.
        moon_distances_km (iterable): The distances of the planet's moons from its centre.
    """
    kind = kind_of_mass(mass_kg)
    low, high = RANGES["ice_fraction"][kind]
    share = draw.uniform(low, high)
    ice_fraction = share if cold else 0.0
    filling = log_uniform(*RANGES["filling"][kind]) if RANGES["filling"][kind][0] < RANGES["filling"][kind][1] else 1.0
    planet_density = planet_mass_kg / (4.0 / 3.0 * math.pi * (radius_km * 1000.0) ** 3)
    roche = roche_limit_km(radius_km, planet_density, mixture_density_kg_m3(ice_fraction))
    outer = roche * draw.uniform(*RANGES["outer_over_roche"])
    if moon_distances_km:
        outer = min(outer, MOON_CLEARANCE * min(moon_distances_km))
    low, high = RANGES["inner_over_radius"]
    high = min(high, outer / (radius_km * MIN_WIDTH_RATIO))  # the inner edge lies inside the outer one
    if high <= low:
        return None
    inner = radius_km * draw.uniform(low, high)
    tau = optical_depth(mass_kg, inner, outer, filling, ice_fraction, RANGES["particle_radius_m"][kind])
    return Ring(ring_kind=kind, ring_origin=origin, ring_inner_km=inner, ring_outer_km=outer, ring_mass_kg=mass_kg,
                ring_optical_depth=tau, ring_ice_fraction=ice_fraction, ring_age_yr=age_yr)


def clear(body):
    """Removes a stored ring system from `body` (a body regenerated or without one)."""
    for name in BODY_FIELDS:
        setattr(body, name, None)


def apply(planet, host=None):
    """
    Stores `planet`'s ring system (`BODY_FIELDS`), replacing any it held; a belt and a body with no radius or mass have none (GEN.1057). The ring is what the planet's history made
    (`ring_origin.ring_origin`): the moons it lost to the tide (`planet.fallen_moons`, `moon_history`) and the comets that passed inside its Roche limit over the star's age, which are
    drawn here (the passages come at a rate scaled from Saturn's by the Roche sphere's cross-section, so a moon, whose Roche sphere is tiny, almost never has one). A dwarf planet that orbits a star with
    no ring from its history may have the narrow ring of its spin-orbit resonance (`small_body_rings`, R4). A moon is a planet of its own kind (GEN.836) whose star distance
    is its `host`'s, and has no fallen moons. Run once the planet's moons and its star's snow line are final.
    """
    from planetgen.physics import ring_origin, small_body_rings   # these import this module
    if planet.body_type == 'a':  # a belt is not a body
        return None
    clear(planet)
    if not planet.radius or not planet.mass:
        return None
    parent = moon_composition.parent_kind(planet.envelope_mass_fraction or 0.0, planet.mass)
    if planet.is_moon:
        cold, distances = bool(planet.born_beyond_snow_line), []
    else:
        cold, distances = planet.distance >= planet.star.snow_line_au, [_au_to_km(moon.distance) for moon in planet.moons]
    around = host if planet.is_moon and host is not None else planet
    star_mass_kg = planet.star.mass
    age_yr = planet.star.age * 1.0e9
    body = ring_origin.Planet(planet.mass, planet.radius, around.distance, star_mass_kg)
    passages = ring_origin.draw_comet_passages(body, age_yr, cold)
    history = ring_origin.ring_origin(body, age_yr, getattr(planet, "fallen_moons", ()) or (), passages)
    resonance = None
    if parent == moon_composition.DWARF_PLANET and not planet.is_moon and planet.rotation_period_hours:
        resonance = small_body_rings.resonance_ring(planet.mass, planet.radius, planet.rotation_period_hours, around.distance, star_mass_kg, cold, draw)
    ring = build_ring(history.origin, history.mass_kg, history.age_yr, planet.radius, planet.mass, cold, distances) if history.exists else None
    if ring is None and resonance is not None:
        ring = resonance.ring
    return _store(planet, ring)


def _store(body, ring):
    if ring is not None:
        for name, value in ring.values().items():
            setattr(body, name, value)
    return ring
