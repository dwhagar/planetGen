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
4. **The population** (GEN.1092). The members of an outer disk that a giant acted on belong to the scatter step's populations (`belt_populations`): each draws one and its semi-major axis,
   eccentricity and inclination follow from it (cold ones stirred little, resonant ones at the giant's resonances, scattered, crossing and detached ones with their perihelion set), and the
   tag is stored with the member. A belt without those shares draws as 1 and 2 say and has no tag.

Limits. A belt stirred by a giant has a larger `e` than its width gives (`belt_mass`'s own limit); a body in a resonance or a gap keeps no trace of it here (the gaps the giants
cleared, `belt_origin`, are not cut out of the solids profile); the bodies do not know the planets near them. Stand-ins MR1 to MR3 of `gap-filling-plan.md` section 6.182.

Use:

    entries = draw_roster(count, lower_au, upper_au, star_mass_kg, snow_line_au, metal_ratio, stream)    # one Entry per stored body, heaviest first
"""

import math
from dataclasses import dataclass

from planetgen.physics import belt_ledger, belt_mass, belt_populations
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

DRAWS_PER_BODY = 11
"""int: The uniform numbers the stream gives each body: taxon, porosity, semi-major axis, eccentricity, inclination, node, periapsis, anomaly, then the population and its two choices."""


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
        tag (str or None): The scatter population the body belongs to (`belt_populations`), None where the belt has none.
    """

    placement: Placement
    eccentricity: float
    inclination_deg: float
    node_deg: float
    periapsis_deg: float
    anomaly_deg: float
    tag: str | None = None


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


def mass_below(table, a_au):
    """The solids of the disk inside `a_au` (Earth masses), the table read forward: log-linear between points."""
    edges, masses = table
    if a_au <= edges[0]:
        return 0.0
    if a_au >= edges[-1]:
        return masses[-1]
    lo, hi = 0, len(edges) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if edges[mid] <= a_au:
            lo = mid
        else:
            hi = mid
    fraction = (math.log(a_au) - math.log(edges[lo])) / (math.log(edges[hi]) - math.log(edges[lo]))
    return masses[lo] + fraction * (masses[hi] - masses[lo])


def _distance_at(table, target):
    """The distance below which `target` Earth masses of the disk's solids lie: the table inverted, log-linear between points."""
    edges, masses = table
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


def semi_major_axis_au(table, u, cleared=()):
    """
    The distance below which a share `u` (0 to 1) of the disk's solids lies: the table inverted, log-linear between points.

    `cleared` is a list of `(low, high)` stretches the giants emptied (GEN.1092, `belt_origin.BeltOrigin.cleared_spans`): the share `u` is then of the solids outside them, so no body is
    placed in one.
    """
    edges, masses = table
    spans = sorted((max(a, edges[0]), min(b, edges[-1])) for a, b in cleared if min(b, edges[-1]) > max(a, edges[0]))
    if not spans:
        return _distance_at(table, u * masses[-1])
    removed = [(mass_below(table, a), mass_below(table, b)) for a, b in spans]
    target = u * (masses[-1] - sum(high - low for low, high in removed))
    for low, high in removed:                # the solids outside the stretches, counted from the star: each stretch passed adds its own mass to the table's
        if target + 0.0 < low:
            break
        target += high - low
    return _distance_at(table, target)


@finite_domain()
def draw_roster(count, lower_au, upper_au, star_mass_kg, snow_line_au, metal_ratio, stream, populations=None, giant_au=None, cleared=()):
    """
    The `Entry` of each of `count` stored bodies (rank 1 first) of a belt between `lower_au` and `upper_au`.

    Args:
        count (int): How many stored bodies the belt's ledger holds (at most `belt_ledger.LEDGER_MAX_BODIES`).
        stream: An object with `random()` (a `draw.Stream` of the belt's own); `DRAWS_PER_BODY * belt_ledger.LEDGER_MAX_BODIES` draws are taken, whatever `count` is.
        populations (tuple or None): `(tag, weight)` of the populations that hold members (`belt_populations.member_weights`); None for a belt the scatter step did not act on, whose bodies have no tag.
        giant_au (float or None): The distance of the giant that swept the disk, for the populations that depend on it (None: all the members are cold).
        cleared (sequence): The `(low, high)` stretches the giants cleared of everything (GEN.1092); the members of the cold population are not placed in them.
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
    for u_taxon, u_porosity, u_a, u_e, u_i, u_node, u_peri, u_anomaly, u_pop, u_pa, u_pe in numbers[:count]:
        a_au, e_scale, e_fixed, tag = semi_major_axis_au(table, u_a, cleared), 1.0, None, None
        if populations is not None:
            tag, a_au, e_scale, e_fixed = belt_populations.orbit(belt_populations.pick(populations, u_pop), a_au, lower_au, upper_au, giant_au, u_pa, u_pe)
        inclination = min(math.degrees(rayleigh(u_i, i_rms * (e_scale if e_scale is not None else 1.0))), I_MAX_DEG)
        eccentricity = e_fixed if e_fixed is not None else min(rayleigh(u_e, e_rms * e_scale), E_MAX)
        entries.append(Entry(Placement(a_au, u_taxon, u_porosity), eccentricity, inclination, 360.0 * u_node, 360.0 * u_peri, 360.0 * u_anomaly, tag))
    return entries
