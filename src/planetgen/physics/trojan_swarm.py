# planetgen/physics/trojan_swarm.py

"""
The two regions of a trojan swarm and where each lies at a time (GEN.1150 slice 0, stand-ins TJ9 to TJ11 of docs/design/research/gap-filling-plan.md)
===============================================================================================================================================

`trojans.py` says what may sit at L4 and L5 of a pair and how a swarm's mass splits between them. This module is the pure draw that follows it, with nothing stored and nothing wired:

* `draw_regions` takes the pair, the mass of the swarm and a stream, and gives the **leading** (L4) and **trailing** (L5) `Region`: its share of the mass, the libration amplitude of its bodies (a
  truncated Rayleigh below `trojans.TADPOLE_LIMIT_DEG`, so the region is a tadpole), the libration phase, its mean eccentricity and its mean inclination (truncated Rayleigh below `ECC_MAX` and
  `INCLINATION_MAX_DEG`). A pair that cannot hold the swarm (`trojans.stable`) gives none. The number of draws is fixed (one for the split, then four for each region, leading first), so no other
  draw moves with the mass or the pair.
* `longitude_offset_deg` and `longitude_deg` are the **libration overlay** the orbit updater adds to the secondary's mean longitude: `+-60 + A sin(2 pi t / P_lib + phase)`, `P_lib` the slow libration
  period of the pair (`trojans.libration_period`). The orbit of a trojan is the secondary's Keplerian orbit; only this term moves it off the point, so a stored epoch gives the place later.

Every constant is in `physics/data/trojans.json` with its tier and tweak level; nothing here sets how often a giant has a swarm (that is the scatter's share of the swept mass).
"""

import collections
import math

from planetgen.physics import trojans

Region = collections.namedtuple("Region", ("point", "mass", "libration_amplitude_deg", "libration_phase_deg", "eccentricity", "inclination_deg"))
"""One region of a swarm: its `point` ("L4" or "L5"), its `mass` (the unit of the masses given), the libration `amplitude` and `phase` of its bodies (degrees), and its mean `eccentricity` and `inclination_deg`
to the secondary's plane."""

LIBRATION_SIGMA_DEG = trojans.NUMBERS["libration_sigma_deg"]["value"]
"""float: The Rayleigh scale of the libration amplitude, degrees. [Stand-in TJ9, tier T4, tweak L4.]"""

ECC_SIGMA = trojans.NUMBERS["swarm_ecc_sigma"]["value"]
"""float: The Rayleigh scale of a region's mean eccentricity. [Stand-in TJ10, tier T4, tweak L4.]"""

INCLINATION_SIGMA_DEG = trojans.NUMBERS["swarm_inclination_sigma_deg"]["value"]
"""float: The Rayleigh scale of a region's mean inclination, degrees. [Stand-in TJ11, tier T4, tweak L4.]"""


def truncated_rayleigh(u, sigma, limit):
    """A Rayleigh variate of scale `sigma` cut at `limit`, from the uniform `u` in 0..1 by its inverse distribution: one draw, no rejection."""
    top = 1.0 - math.exp(-0.5 * (limit / sigma) ** 2)
    return sigma * math.sqrt(-2.0 * math.log(1.0 - u * top))


def draw_regions(host_mass, secondary_mass, trojan_mass, rng):
    """
    The regions of a swarm of `trojan_mass` at L4 and L5 of a pair (masses in one unit), or `()` when the pair cannot hold it (`trojans.stable`) or the swarm has no mass.

    Args:
        host_mass (float): The host's mass.
        secondary_mass (float): The secondary's mass (a giant, a brown dwarf, a large moon).
        trojan_mass (float): The swarm's mass (the scatter's trojan share of the swept mass, less what the neighbours take).
        rng: Anything with `random()`; the `draw` module or a `Stream`.

    Returns:
        tuple: The leading and then the trailing `Region`; draws exactly nine numbers from `rng`, or none for a swarm it cannot hold.
    """
    if trojan_mass <= 0.0 or not trojans.stable(host_mass, secondary_mass, trojan_mass):
        return ()
    leading_mass, trailing_mass = trojans.swarm_split(trojan_mass, trojans.leading_share(rng.random()))
    found = []
    for point, mass in (("L4", leading_mass), ("L5", trailing_mass)):
        amplitude = truncated_rayleigh(rng.random(), LIBRATION_SIGMA_DEG, trojans.TADPOLE_LIMIT_DEG)
        phase = 360.0 * rng.random()
        eccentricity = truncated_rayleigh(rng.random(), ECC_SIGMA, trojans.ECC_MAX)
        inclination = truncated_rayleigh(rng.random(), INCLINATION_SIGMA_DEG, trojans.INCLINATION_MAX_DEG)
        found.append(Region(point, mass, amplitude, phase, eccentricity, inclination))
    return tuple(found)


def longitude_offset_deg(region, years, orbital_period_years, mu):
    """
    The mean longitude of a region's bodies from the secondary's, degrees, `years` after the epoch: the point's `+-60` plus the libration `A sin(2 pi t / P_lib + phase)`, with `P_lib` the slow libration
    period of a pair of mass parameter `mu` whose secondary takes `orbital_period_years` to go round (`trojans.libration_period`).

    Raises:
        trojans.TrojanError: For a pair that holds nothing (no libration period).
    """
    period = trojans.libration_period(orbital_period_years, mu)
    if period is None:
        raise trojans.TrojanError(f"a pair of mass parameter {mu!r} holds no tadpole")
    swing = region.libration_amplitude_deg * math.sin(2.0 * math.pi * years / period + math.radians(region.libration_phase_deg))
    return trojans.longitude_offset_deg(region.point) + swing


def longitude_deg(region, secondary_longitude_deg, years, orbital_period_years, mu):
    """The mean longitude of a region's bodies, degrees in 0..360: the secondary's plus `longitude_offset_deg`."""
    return (secondary_longitude_deg + longitude_offset_deg(region, years, orbital_period_years, mu)) % 360.0
