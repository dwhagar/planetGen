# planetgen/physics/giant_heat.py

"""
The internal heat of a giant, an ice giant or a brown dwarf (GEN.810, registry 6.135)
==============================================================================================

The heat a giant radiates of its own: what it leaks from the contraction and cooling of its whole body, and what the sinking of helium
adds.

**The base.** The cooling law through Jupiter and Neptune of `rogue_surface.giant_internal_flux_w_m2` (GEN.816: luminosity as the mass to the 1.65
below a Jupiter mass, falling as the age to the -1.3), at the body's mass, radius and age (`base_flux_w_m2`; its stand-ins are GEN.816's). The
giant's one-bar temperature (`giant_draw.own_flux_w_m2`) is made from the same law, before the interior exists.

**The helium rain.** A giant cold enough at one bar lets helium rain out of its molecular envelope (`giant_interior.helium_depletion`):
the helium sinks to the helium-enriched zone, the helium shell and the dilute core below, and the energy it gives up falls as heat
([R1427, R2852, R4220], who find it a third or so of Saturn's excess).
The interior stores where the helium went (`interior_layers`: every layer's mass, its helium and the radii of its top and bottom;
`helium_depletion`; `has_helium_shell`), and this module reads them back:

1. the helium of each layer against the natal share of its hydrogen and helium says which layers lost it (the envelope above) and
   which gained it (the zone, the shell and the core below);
2. the gravitational potential at the middle of each layer follows from the layers themselves, each a shell of uniform density
   (`potentials`), exactly so for a single uniform sphere;
3. the helium that moved is the sum of the gains, and the energy it gave up is that mass times the potential it fell through, from the
   losing layers' mean to the gaining layers' mean (`rain_energy_j`), weighted by the helium each moved;
4. a share of it, `RAIN_HEAT_SHARE`, leaves as heat and the rest warms the interior it falls through (virial theorem); that heat has come out
   over the body's whole age, so the power is its average over it (`RAIN_AGE_SHARE`, row GH3).

A giant with no rain (depletion 1) has no extra heat, and a giant whose helium went deeper (a shell) has more than one whose helium stayed
in a gradient. The one-bar temperature that decided the rain was made from the base alone, so the heat the rain adds is not fed back to it
(row GH4: one pass, no iteration, to keep the draws in their order).

Use:

    heat = of(body)                 # reads body.interior_layers, helium_depletion, has_helium_shell, radius, star.age
    heat.flux_w_m2                  # the heat radiated per square metre of the 1 bar surface
    heat.rain_power_w               # of it, what the helium rain gives
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, giant_interior, rogue_surface

RAIN_HEAT_SHARE = 0.5
"""float: The share of the energy the sinking helium gives up that leaves as heat; the rest warms the interior it falls through (the
virial split of gravitational energy into heat and thermal energy, row GH2)."""

RAIN_AGE_SHARE = 1.0
"""float: The share of the body's age over which the rain's heat has come out; 1 takes it as the average over the whole age, which is
the least power it can give today, since the rain began when the body cooled below the temperature where helium separates (row GH3)."""

SECONDS_PER_GYR = 1.0e9 * constants.SECONDS_PER_YEAR
"""float: Seconds in a billion years."""


@dataclass(frozen=True)
class Heat:
    """
    Attributes:
        flux_w_m2 (float): The heat radiated per square metre of the body's 1 bar surface: the base and the rain's.
        base_w_m2 (float): The base (the cooling law, `base_flux_w_m2`).
        rain_w_m2 (float): What the helium rain adds.
        rain_power_w (float): The rain's power over the whole body.
        moved_helium_earth (float): The helium that sank, Earth masses.
    """

    flux_w_m2: float
    base_w_m2: float
    rain_w_m2: float
    rain_power_w: float
    moved_helium_earth: float


def base_flux_w_m2(mass_earth, radius_m, age_gy):
    """The heat a giant of `mass_earth`, `radius_m` and `age_gy` radiates of its own before the helium rain adds to it: the cooling law of GEN.816."""
    return rogue_surface.giant_internal_flux_w_m2(mass_earth * constants.EARTH_MASS_TO_KG, radius_m / constants.KM_TO_M_FACTOR, age_gy)


def _shell_integral(mass, top, bottom):
    """The integral of dm / r over a shell of uniform density between the radii `bottom` and `top`."""
    if top <= bottom:
        return mass / top
    return 1.5 * mass * (top * top - bottom * bottom) / (top ** 3 - bottom ** 3)


def potentials(interior, radius_m):
    """
    The gravitational potential (J/kg, negative) at the middle of each layer of `interior`, outside in, treating every layer as a shell
    of uniform density between its stored radii: the pull of the mass inside the middle radius, the mass of the layer that lies deeper, and
    of the shells and the outer half of the layer above it.
    """
    kg = constants.EARTH_MASS_TO_KG
    masses = [layer.mass_earth * kg for layer in interior.layers]
    tops = [layer.top_radius_fraction * radius_m for layer in interior.layers]
    bottoms = [layer.bottom_radius_fraction * radius_m for layer in interior.layers]
    result = []
    for i, mass in enumerate(masses):
        below = math.fsum(masses[i + 1:])
        top, bottom = tops[i], bottoms[i]
        middle = ((top ** 3 + bottom ** 3) / 2.0) ** (1.0 / 3.0)
        above = math.fsum(_shell_integral(masses[j], tops[j], bottoms[j]) for j in range(i))
        outer_half = 0.5 * mass / top if top <= bottom else 1.5 * mass * (top * top - middle * middle) / (top ** 3 - bottom ** 3)
        result.append(-constants.G * ((below + mass / 2.0) / middle + above + outer_half))
    return result


def helium_moved(interior):
    """
    `[(layer, helium kg gained (+) or lost (-))]` against the natal share of each layer's hydrogen and helium: the natal share is the body's
    total helium over its total hydrogen and helium.
    """
    kg = constants.EARTH_MASS_TO_KG
    helium = interior.helium_mass_earth()
    hydrogen_helium = interior.total_mass_earth - interior.heavy_mass_earth()
    natal = helium / hydrogen_helium if hydrogen_helium > 0.0 else 0.0
    found = []
    for layer in interior.layers:
        gas = layer.mass_earth - layer.heavy_mass_earth()
        found.append(((layer.helium_mass_earth() - natal * gas) * kg))
    return found


def rain_energy_j(interior, radius_m):
    """`(helium kg that sank, energy J the sinking gave up)` from the layers of `interior`, for a body of `radius_m` (0, 0 with no rain)."""
    moves = helium_moved(interior)
    gained = [(m, i) for i, m in enumerate(moves) if m > 0.0]
    lost = [(-m, i) for i, m in enumerate(moves) if m < 0.0]
    moved = math.fsum(m for m, _ in gained)
    if not gained or not lost or moved <= 0.0:
        return 0.0, 0.0
    phi = potentials(interior, radius_m)
    phi_lost = math.fsum(m * phi[i] for m, i in lost) / math.fsum(m for m, _ in lost)
    phi_gained = math.fsum(m * phi[i] for m, i in gained) / moved
    return moved, max(0.0, moved * (phi_lost - phi_gained))


def heat(interior, radius_m, age_gy):
    """The `Heat` of a giant with the layers `interior`, radius `radius_m` (m) and age `age_gy` (Gyr)."""
    area = 4.0 * math.pi * radius_m ** 2
    base = base_flux_w_m2(interior.total_mass_earth, radius_m, age_gy)
    moved, energy = rain_energy_j(interior, radius_m)
    duration_s = max(age_gy * RAIN_AGE_SHARE, 1.0e-3) * SECONDS_PER_GYR
    power = RAIN_HEAT_SHARE * energy / duration_s
    return Heat(base + power / area, base, power / area, power, moved / constants.EARTH_MASS_TO_KG)


def of(body):
    """
    The `Heat` of the giant a body's stored columns describe (`interior_layers`, `helium_depletion`, `has_helium_shell`, its radius and its
    star's age), or `None` for a body that stores no giant interior. Without the rain (`helium_depletion` 1, or no layer that gained helium)
    it is the base alone.
    """
    interior = giant_interior.interior_of(body)
    if interior is None:
        return None
    depletion = getattr(body, "helium_depletion", None)
    raining = (depletion is not None and depletion < 1.0) or bool(getattr(body, "has_helium_shell", False))
    age_gy = getattr(getattr(body, "star", None), "age", None) or 4.6
    found = heat(interior, body.radius * constants.KM_TO_M_FACTOR, age_gy)
    if raining:
        return found
    return Heat(found.base_w_m2, found.base_w_m2, 0.0, 0.0, 0.0)
