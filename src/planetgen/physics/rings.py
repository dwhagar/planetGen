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

What is drawn (every number is a labelled stand-in; the registry rows are BC20 to BC23 of the gap plan)
-------------------------------------------------------------------------------------------------------
1. **Whether, and which kind.** A giant always has some ring (all four have): `broad` (a massive ring of large particles, Saturn), `narrow`
   (a few dense narrow rings, Uranus and Neptune) or `dusty` (a faint dust ring from the erosion of small moons, Jupiter), in the
   ratio 1 : 2 : 1 of the Solar System's four. A terrestrial planet or a dwarf planet has a ring with probability `RING_CHANCE`
   (a transient ring of a disrupted moon, or the narrow rings of Haumea and Quaoar).
2. **Origin.** `moon_erosion` for dust, `disrupted_moon` for narrow rings, and for a broad ring `disrupted_moon` or `captured_body`
   with equal weight (a stripped differentiated moon, Canup 2010, or a captured comet, Dones 1991; no source weighs them).
3. **Ice share, mass, extent.** See `RANGES`.

The draw order is the order above and is part of the generator epoch.

Use:

    ring = draw_ring(parent_kind, radius_km, mass_kg, distance_au, snow_line_au, moon_distances_km)   # a Ring or None
    apply(planet)                                                                                       # sets the stored fields
"""

import math

from planetgen.physics import constants as physical_constants, moon_composition, unit_registry as units
from planetgen.util import draw
from planetgen.util.random import log_uniform

_au_to_km = units.converter("au", "km")

KINDS = ("dusty", "narrow", "broad")
ORIGINS = ("moon_erosion", "disrupted_moon", "captured_body")

BODY_FIELDS = ("ring_kind", "ring_origin", "ring_inner_km", "ring_outer_km", "ring_mass_kg", "ring_optical_depth",
               "ring_ice_fraction")
"""tuple: The attributes (and `planets`/`moons` columns, GEN.732) of a ring system; all `None` on a body without one."""

UNITS = {"ring_kind": "text", "ring_origin": "text", "ring_inner_km": "km", "ring_outer_km": "km", "ring_mass_kg": "kg",
         "ring_optical_depth": "ratio", "ring_ice_fraction": "frac"}
"""dict: Each field of `BODY_FIELDS` and its native unit (GEN.570), as `RoguePlanet.FIELD_UNITS` carries it."""

ROCHE_COEFFICIENT = 2.456
"""float: The fluid Roche limit in planet radii for particles as dense as the planet (Roche 1850; Murray and Dermott 1999). Recalled, T2, L0."""

ICE_DENSITY_KG_M3 = 917.0
"""float: Water ice Ih, the density of a ring's ice (the repo's ice density, `constants`-free so the module stands alone)."""
ROCK_DENSITY_KG_M3 = 2500.0
"""float: The density of a ring's silicate or dark rock; a stand-in between chondritic grains (3300) and porous dust (1000 to 2000), T6, L0."""

KIND_WEIGHTS = {"giant": (("dusty", 1.0), ("narrow", 2.0), ("broad", 1.0))}
"""dict: The kind of a giant's ring: Jupiter dusty, Uranus and Neptune narrow, Saturn broad (n = 4: T6, L0)."""
RING_CHANCE = {"terrestrial": 0.01, "dwarf-planet": 0.01}
"""dict: The chance that a terrestrial planet or a dwarf planet has a ring (none among the Solar System's terrestrial planets; Haumea, Quaoar
and Chariklo are small-body rings). A stand-in with no source, T7, L0."""
NON_GIANT_KIND = {"terrestrial": "dusty", "dwarf-planet": "narrow"}
"""dict: The kind of the ring of a body that is not a giant: dust from a disrupted moon, or narrow dense rings."""

RANGES = {
    "mass_kg": {"dusty": (1.0e7, 1.0e11), "narrow": (1.0e13, 1.0e17), "broad": (1.0e18, 1.0e21)},
    "particle_radius_m": {"dusty": 1.0e-5, "narrow": 1.0, "broad": 1.0},
    "filling": {"dusty": (1.0, 1.0), "narrow": (1.0e-3, 1.0e-2), "broad": (1.0, 1.0)},
    "ice_fraction": {"dusty": (0.0, 0.0), "narrow": (0.0, 0.3), "broad": (0.90, 0.99)},
    "outer_over_roche": (0.75, 1.1),
    "inner_over_radius": (1.1, 1.7),
}
"""dict: The drawn ranges (stand-ins, T3 to T6, L0).

- `mass_kg`: log-uniform. Jupiter's main ring and halo hold 1e7 to 1e9 kg of dust, Uranus's rings about 1e16 kg, Saturn's 1.5e19 kg
  (Iess et al. 2019, recalled); the ranges bracket each with a decade or more either side.
- `particle_radius_m`: the effective particle radius that turns a mass into an optical depth (Saturn's particles run from centimetres to
  ten metres; dust grains are about 10 micrometres).
- `filling`: the share of the annulus a ring covers, log-uniform; Uranus's rings cover about 0.3 percent of theirs.
- `ice_fraction`: the ice share by mass beyond the snow line (Saturn's rings 90 to 99 percent, Cuzzi et al. 2018, recalled; Uranus's and
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


def _draw_kind(parent):
    """The kind of a giant's ring by `KIND_WEIGHTS`: one draw."""
    kinds, weights = zip(*KIND_WEIGHTS[parent])
    pick = draw.random() * sum(weights)
    for kind, weight in zip(kinds, weights):
        pick -= weight
        if pick < 0.0:
            return kind
    return kinds[-1]


def draw_ring(parent, radius_km, mass_kg, cold, moon_distances_km=(), fallen_mass_kg=0.0):
    """
    A planet's ring system, or `None` when it has none.

    Args:
        parent (str): `moon_composition.GIANT`, `TERRESTRIAL` or `DWARF_PLANET`.
        radius_km (float): The planet's radius.
        mass_kg (float): The planet's mass.
        cold (bool): Whether the planet orbits beyond its star's snow line, where a ring keeps its ice.
        moon_distances_km (iterable): The distances of the planet's moons from its centre.
        fallen_mass_kg (float): The mass of the moons that fell inside the Roche limit (`moon_history`, GEN.959): the planet has a ring of disrupted moon,
            of that mass up to `MAX_MASS_OVER_PLANET` of the planet's (the rest rained onto the planet), whatever the chance said. The draws are the same.

    Returns:
        Ring or None: The draws are made in the order of the module docstring whether or not a ring results.
    """
    if parent == moon_composition.GIANT:
        kind = _draw_kind(parent)
    else:
        if draw.random() >= RING_CHANCE[parent] and fallen_mass_kg <= 0.0:
            return None
        kind = NON_GIANT_KIND[parent]
    origin = {"dusty": "moon_erosion", "narrow": "disrupted_moon"}.get(kind) or draw.choice(("disrupted_moon", "captured_body"))
    if fallen_mass_kg > 0.0:
        origin = "disrupted_moon"
    low, high = RANGES["ice_fraction"][kind]
    share = draw.uniform(low, high)
    ice_fraction = share if cold else 0.0
    mass = min(log_uniform(*RANGES["mass_kg"][kind]) if fallen_mass_kg <= 0.0 else fallen_mass_kg, MAX_MASS_OVER_PLANET * mass_kg)
    filling = log_uniform(*RANGES["filling"][kind]) if RANGES["filling"][kind][0] < RANGES["filling"][kind][1] else 1.0
    planet_density = mass_kg / (4.0 / 3.0 * math.pi * (radius_km * 1000.0) ** 3)
    roche = roche_limit_km(radius_km, planet_density, mixture_density_kg_m3(ice_fraction))
    outer = roche * draw.uniform(*RANGES["outer_over_roche"])
    if moon_distances_km:
        outer = min(outer, MOON_CLEARANCE * min(moon_distances_km))
    low, high = RANGES["inner_over_radius"]
    high = min(high, outer / (radius_km * MIN_WIDTH_RATIO))  # the inner edge lies inside the outer one
    if high <= low:
        return None
    inner = radius_km * draw.uniform(low, high)
    tau = optical_depth(mass, inner, outer, filling, ice_fraction, RANGES["particle_radius_m"][kind])
    return Ring(ring_kind=kind, ring_origin=origin, ring_inner_km=inner, ring_outer_km=outer, ring_mass_kg=mass,
                ring_optical_depth=tau, ring_ice_fraction=ice_fraction)


def clear(body):
    """Removes a stored ring system from `body` (a body regenerated or without one)."""
    for name in BODY_FIELDS:
        setattr(body, name, None)


def apply(planet):
    """
    Draws and stores `planet`'s ring system (`BODY_FIELDS`), replacing any it held; a belt and a body with no radius or mass have none.
    A moon draws like a planet of its own kind around no moons (GEN.836), cold when its parent is beyond the snow line. Run once the
    planet's moons and its star's snow line are final.
    """
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
    return _store(planet, draw_ring(parent, planet.radius, planet.mass, cold, distances, getattr(planet, "fallen_moon_mass_kg", 0.0) or 0.0))


def apply_rogue(rogue):
    """
    Draws and stores a rogue planet's ring system (GEN.836): the same law as a bound planet with no host and no moon distances. It is
    cold, and keeps its ice, when its effective temperature is below the snow-line temperature. A brown dwarf is no planet and has none.
    """
    clear(rogue)
    if rogue.planet_class is None:
        return None
    if rogue.planet_type == 'g':
        parent = moon_composition.GIANT
    else:
        parent = moon_composition.DWARF_PLANET if rogue.mass_kg < moon_composition.DWARF_PLANET_MAX_MASS_KG else moon_composition.TERRESTRIAL
    cold = rogue.effective_temperature_k < physical_constants.SNOW_LINE_TEMPERATURE_K
    return _store(rogue, draw_ring(parent, rogue.radius_km, rogue.mass_kg, cold))


def _store(body, ring):
    if ring is not None:
        for name, value in ring.values().items():
            setattr(body, name, value)
    return ring
