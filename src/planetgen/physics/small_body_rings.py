# planetgen/physics/small_body_rings.py

"""
Rings of small bodies beyond the Roche limit, and the gate against rings close to a star (GEN.1056; docs/design/systems/weird-worlds.md section 9.2, rows R4 and R9)
=====================================================================================================================================================================

Two pure functions of the rings family (`rings.py`); `rings.apply` calls `resonance_ring` for a dwarf planet that orbits a star and has no ring from its history (GEN.1057). `ring_survives` is the same gate as `ring_origin.RING_MIN_DISTANCE_AU`.

R4, a resonance ring. Quaoar has a ring at 4 100 km (7.4 of its radii, well past its Roche limit), Haumea one at 2 290 km and Chariklo two at 391 and 405 km. In each the ring's
orbital period is close to three times the body's rotation period: Quaoar's ring is at the 1/3 spin-orbit resonance [R3155] and Haumea's at the 3:1 (Ortiz et al.
2017), which is the same ratio read from the other side. A resonance of that order lets a body that is not round pull on ring material at a place where the Roche limit would let it
gather into a moon, so the ring survives where an ordinary ring could not. [2ND-HAND: these are the cases of the design note; no paper was read.] The place follows from Kepler's
third law, `a = (G M (n P / 2 pi)^2)^(1/3)` for a ring that goes round `n` times slower than the body turns (`resonance_radius_km`). The ring is held to the places where it can exist:

* **Past the Roche limit** of its own material (`rings.roche_limit_km`); inside it the ordinary rings of `rings.py` are the case.
* **Inside the stable zone of the body's orbit round its star**, `HILL_FRACTION_MAX` of its Hill radius ([R1138, R1771]: 0.49 for prograde orbits).
* **Round a small body** (under `moon_composition.DWARF_PLANET_MAX_MASS_KG`); a larger body has a ring by the rule of `rings.py`.

A body that meets the conditions has the ring with probability `RESONANCE_RING_RATE` (RR1; no survey of ring searches gives a rate, three of the well observed bodies have one, so it is a
stand-in with the source's falsifier). The ring is narrow and thin: `rings.Ring` of kind `narrow`, its origin a `disrupted_moon` (the cause named for Quaoar's), its width a few percent of its radius and its
mass from the stand-in range of `RESONANCE_RING_MASS_KG` (RR4); the optical depth is an output of the mass and the width, as in `rings.py`.

R9, no ring close to a star (`ring_survives`): a planet inside `INNER_RING_LIMIT_AU` keeps no ring ([R3957], via the design note, 2ND-HAND: the rings of close planets are short lived against
the star's tides and heat). The function is a gate, the rate is zero inside the limit and the caller's own outside it.

Use:

    ring = resonance_ring(mass_kg, radius_km, spin_period_h, distance_au, host_mass_kg, cold, stream)   # a ResonanceRing or None
    ring_survives(distance_au)                                                                           # R9: False inside 0.1 AU
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants as physical_constants, moon_composition, orbits, rings
from planetgen.util.checks import finite_domain
from planetgen.util.random import log_uniform

_AU_M = 1.495978707e11
"""float: The astronomical unit in metres (IAU 2012)."""

RESONANCE_ORDER = 3
"""int: How many times slower than the body's rotation the ring goes round (Quaoar 1/3, Haumea 3:1, Chariklo near it). [Stand-in RR1, tier T4, tweak L3.]"""

RESONANCE_RING_RATE = 0.10
"""float: The chance that a small body that meets the conditions of R4 has the ring: three of the roughly twenty bodies whose shape or rotation has been mapped by occultation are known to have one,
a selection-biased count. [Stand-in RR1, tier T4, tweak L3.]"""

HILL_FRACTION_MAX = 0.5
"""float: The widest circular prograde orbit that stays bound to a body, in its Hill radii ([R1138, R1771]: 0.49). [Stand-in RR3, tier T2, tweak L0.]"""

INNER_RING_LIMIT_AU = 0.1
"""float: Inside this distance from its star a planet holds no ring (R9; [R3957] via the design note, 2ND-HAND, 'about 0.1 AU'). [Stand-in RR2, tier T3, tweak L0.]"""

RESONANCE_RING_MASS_KG = (1.0e10, 1.0e14)
"""tuple: The ring's mass, log-uniform, kg: the range that gives a mean optical depth from about 0.001 to 10 over the covered share of a thin ring at a few thousand km, the span of the dense
narrow rings measured. [Stand-in RR4, tier T6, tweak L3.]"""

RESONANCE_WIDTH_FRACTION = (0.002, 0.05)
"""tuple: The ring's width over its radius, log-uniform (Quaoar's ring is 5 to 200 km wide at 4 100 km). [Stand-in RR4, tier T6, tweak L3.]"""

NARROW = "narrow"
ORIGIN = "disrupted_moon"


@dataclass(frozen=True)
class ResonanceRing:
    """
    A small body's ring at a spin-orbit resonance.

    Attributes:
        ring (rings.Ring): The ring system (kind `narrow`, origin `disrupted_moon`).
        order (int): How many rotations of the body a ring orbit takes.
        radius_km (float): The resonance radius the ring is centred on.
        roche_km (float): The Roche limit of the ring's own material for this body; `radius_km` is past it.
    """

    ring: rings.Ring
    order: int
    radius_km: float
    roche_km: float


@finite_domain()
def resonance_radius_km(mass_kg, spin_period_h, order=RESONANCE_ORDER):
    """The radius, km from the body's centre, of the orbit whose period is `order` times the body's rotation period."""
    if not (mass_kg > 0.0 and spin_period_h > 0.0 and order > 0):
        raise ValueError("the mass, the spin period and the order must be positive")
    period_s = order * spin_period_h * 3600.0
    return (physical_constants.G * mass_kg * (period_s / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0) / 1000.0


@finite_domain()
def ring_survives(distance_au):
    """R9: whether a ring lasts at this distance from the star (False inside `INNER_RING_LIMIT_AU`)."""
    if not distance_au > 0.0:
        raise ValueError("the distance must be positive")
    return distance_au >= INNER_RING_LIMIT_AU


def resonance_ring(mass_kg, radius_km, spin_period_h, distance_au, host_mass_kg, cold, stream):
    """
    The ring of a small body at its spin-orbit resonance (R4), or `None`.

    Five draws are taken from `stream` in a fixed order whether or not a ring results (the chance, the ice share, the mass, the width, the covered share), so the stream a caller keys on the
    body is the same either way.

    Args:
        mass_kg, radius_km (float): The body's mass and mean radius.
        spin_period_h (float): Its rotation period.
        distance_au (float): Its distance from the star it orbits.
        host_mass_kg (float): The star's mass (for the Hill radius).
        cold (bool): Whether the body is beyond the star's snow line, where the ring keeps its ice.
        stream: An object with `random()` and `uniform(a, b)` (a `draw.Stream`).

    Returns:
        ResonanceRing or None
    """
    chance = stream.random()
    low, high = rings.RANGES["ice_fraction"][NARROW]
    share = stream.uniform(low, high)
    ring_mass = min(log_uniform(*RESONANCE_RING_MASS_KG, rng=stream), rings.MAX_MASS_OVER_PLANET * mass_kg)
    width = log_uniform(*RESONANCE_WIDTH_FRACTION, rng=stream)
    filling = log_uniform(*rings.RANGES["filling"][NARROW], rng=stream)
    if mass_kg >= moon_composition.DWARF_PLANET_MAX_MASS_KG or not ring_survives(distance_au) or chance >= RESONANCE_RING_RATE:
        return None
    ice_fraction = share if cold else 0.0
    centre = resonance_radius_km(mass_kg, spin_period_h)
    density = mass_kg / (4.0 / 3.0 * math.pi * (radius_km * 1000.0) ** 3)
    roche = rings.roche_limit_km(radius_km, density, rings.mixture_density_kg_m3(ice_fraction))
    hill_km = orbits.calculate_hill_sphere(distance_au * _AU_M, mass_kg, host_mass_kg) / 1000.0
    if centre <= roche or centre > HILL_FRACTION_MAX * hill_km:
        return None
    inner, outer = centre * (1.0 - width / 2.0), centre * (1.0 + width / 2.0)
    tau = rings.optical_depth(ring_mass, inner, outer, filling, ice_fraction, rings.RANGES["particle_radius_m"][NARROW])
    ring = rings.Ring(ring_kind=NARROW, ring_origin=ORIGIN, ring_inner_km=inner, ring_outer_km=outer, ring_mass_kg=ring_mass,
                      ring_optical_depth=tau, ring_ice_fraction=ice_fraction, ring_age_yr=None)
    return ResonanceRing(ring, RESONANCE_ORDER, centre, roche)
