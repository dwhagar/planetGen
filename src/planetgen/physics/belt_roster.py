# planetgen/physics/belt_roster.py

"""
Where each stored member of a belt is, and the numbers that pick what it is made of (GEN.1070, docs/design/systems/belt-evolution.md sections 3 and 5)
==============================================================================================================================================

The belt ledger (`belt_ledger`, GEN.1061) says how heavy each stored body is, and `belt_member` (GEN.1067) says what a body at a distance is made of. This module places each stored
body: its semi-major axis, its eccentricity and inclination, and the three angles that put it on its orbit, with the two uniform numbers that pick its taxon and, if light, its porosity.
Everything comes from one stream of the belt's own (`draw_roster`): the same number of draws whatever the belt, so the member of rank `n` is the same body whatever the belt's mass
or the ledger's length, and the rest of the run's draws are not touched.

1. **The semi-major axis** is drawn where the belt's solids are: the cumulative mass of the disk's solids across the stretch (`belt_mass.primordial_mass_earth`, `CDF_POINTS` log-spaced
   points) is inverted by one uniform number. A belt thick near its inner edge puts more of its bodies there.
2. **The eccentricity** is a Rayleigh draw whose rms is the width rule of `belt_mass` (the particles' eccentricity from the belt's width, `e = (upper - lower) / (upper + lower)`, the extent
   of the orbits of a belt that stirs itself), cut at `E_MAX`. **The inclination** is a Rayleigh draw with the same rms in radians (`belt_mass`'s `e/I = 1`), cut at `I_MAX_DEG`.
3. **The angles** (ascending node, argument of periapsis, mean anomaly) are uniform over the circle.

Limits. A belt stirred by a giant has a larger `e` than its width gives (`belt_mass`'s own limit); a body in a resonance or a gap keeps no trace of it here (the gaps the giants
cleared, `belt_origin`, are not cut out of the solids profile); the bodies do not know the planets near them. Stand-ins MR1 to MR3 of `gap-filling-plan.md` section 6.167.

Use:

    entries = draw_roster(count, lower_au, upper_au, star_mass_kg, snow_line_au, metal_ratio, stream)    # one Entry per stored body, heaviest first
"""

import math
from dataclasses import dataclass

from planetgen.physics import belt_ledger, belt_mass
from planetgen.physics.belt_page import Placement
from planetgen.util.checks import finite_domain

CDF_POINTS = 48
"""int: The log-spaced points across a belt at which the cumulative mass of the disk's solids is tabulated (numerical). [Stand-in MR1, tier T9, tweak L0.]"""

CDF_STEPS = 8
"""int: The log-spaced steps `belt_mass.primordial_mass_earth` takes between the belt's inner edge and each tabulated point (numerical). [Stand-in MR1, tier T9, tweak L0.]"""

E_MAX = 0.5
"""float: The largest eccentricity a stored member is given (a body above it would cross the belt's whole width and the next body's orbit). [Stand-in MR3, tier T9, tweak L0.]"""

I_MAX_DEG = 60.0
"""float: The largest inclination a stored member is given, degrees (the Rayleigh tail above it is cut). [Stand-in MR3, tier T9, tweak L0.]"""

E_OVER_I = 1.0
"""float: The ratio of the particles' rms eccentricity to their rms inclination in radians, 1 as `belt_mass` takes it (Wyatt et al.'s `e/I = 1`). [Stand-in MR2, tier T1, tweak L0.]"""

DRAWS_PER_BODY = 8
"""int: The uniform numbers the stream gives each body: taxon, porosity, semi-major axis, eccentricity, inclination, node, periapsis, anomaly."""


class BeltRosterError(ValueError):
    """A belt or a count the roster cannot be drawn for."""


@dataclass(frozen=True)
class Entry:
    """
    One stored body's place.

    Attributes:
        placement (belt_page.Placement): The semi-major axis in AU and the two numbers that pick the taxon and the porosity.
        eccentricity (float): 0 to `E_MAX`.
        inclination_deg (float): 0 to `I_MAX_DEG`.
        node_deg, periapsis_deg, anomaly_deg (float): Longitude of the ascending node, argument of periapsis and mean anomaly, degrees, 0 to 360.
    """

    placement: Placement
    eccentricity: float
    inclination_deg: float
    node_deg: float
    periapsis_deg: float
    anomaly_deg: float


def rms_eccentricity(lower_au, upper_au):
    """The rms eccentricity of the particles of a belt of these limits, from its width (`belt_mass`'s rule)."""
    return (upper_au - lower_au) / (upper_au + lower_au)


def rayleigh(u, rms):
    """A Rayleigh draw of rms `rms` from the uniform number `u` (the scale is `rms / sqrt(2)`)."""
    return rms / math.sqrt(2.0) * math.sqrt(-2.0 * math.log(1.0 - u))


def cumulative_table(lower_au, upper_au, star_mass_kg, snow_line_au, metal_ratio):
    """`(edges, masses)`: `CDF_POINTS + 1` log-spaced distances across the stretch and the cumulative mass of the disk's solids inside each, Earth masses (0 at the inner edge)."""
    step = (math.log(upper_au) - math.log(lower_au)) / CDF_POINTS
    edges = [math.exp(math.log(lower_au) + k * step) for k in range(CDF_POINTS + 1)]
    edges[-1] = upper_au
    masses = [0.0] + [belt_mass.primordial_mass_earth(lower_au, edge, star_mass_kg, snow_line_au, metal_ratio, steps=CDF_STEPS) for edge in edges[1:]]
    return edges, masses


def semi_major_axis_au(table, u):
    """The distance below which a share `u` (0 to 1) of the disk's solids lies: the table inverted, log-linear between points."""
    edges, masses = table
    target = u * masses[-1]
    lo, hi = 0, len(masses) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if masses[mid] <= target:
            lo = mid
        else:
            hi = mid
    span = masses[hi] - masses[lo]
    fraction = 0.0 if span <= 0.0 else (target - masses[lo]) / span
    return math.exp(math.log(edges[lo]) + fraction * (math.log(edges[hi]) - math.log(edges[lo])))


@finite_domain()
def draw_roster(count, lower_au, upper_au, star_mass_kg, snow_line_au, metal_ratio, stream):
    """
    The `Entry` of each of `count` stored bodies (rank 1 first) of a belt between `lower_au` and `upper_au`.

    Args:
        count (int): How many stored bodies the belt's ledger holds (at most `belt_ledger.LEDGER_MAX_BODIES`).
        stream: An object with `random()` (a `draw.Stream` of the belt's own); `DRAWS_PER_BODY * belt_ledger.LEDGER_MAX_BODIES` draws are taken, whatever `count` is.
    """
    if not 0.0 < lower_au < upper_au:
        raise BeltRosterError("the belt needs 0 < lower < upper")
    if not 0 <= count <= belt_ledger.LEDGER_MAX_BODIES:
        raise BeltRosterError(f"a belt stores 0 to {belt_ledger.LEDGER_MAX_BODIES} bodies, not {count}")
    numbers = [[stream.random() for _ in range(DRAWS_PER_BODY)] for _ in range(belt_ledger.LEDGER_MAX_BODIES)]
    table = cumulative_table(lower_au, upper_au, star_mass_kg, snow_line_au, metal_ratio)
    e_rms = rms_eccentricity(lower_au, upper_au)
    i_rms = e_rms / E_OVER_I
    entries = []
    for u_taxon, u_porosity, u_a, u_e, u_i, u_node, u_peri, u_anomaly in numbers[:count]:
        inclination = min(math.degrees(rayleigh(u_i, i_rms)), I_MAX_DEG)
        entries.append(Entry(Placement(semi_major_axis_au(table, u_a), u_taxon, u_porosity), min(rayleigh(u_e, e_rms), E_MAX), inclination,
                             360.0 * u_node, 360.0 * u_peri, 360.0 * u_anomaly))
    return entries
