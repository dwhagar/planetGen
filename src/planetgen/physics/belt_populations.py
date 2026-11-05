# planetgen/physics/belt_populations.py

"""
The population each stored member of an outer disk belongs to, and the orbit that gives it (GEN.1092)
======================================================================================================

The scatter step (`scatter.scatter_shares`, GEN.1089) says what share of an outer disk the giants left on cold orbits, excited, held in resonances, scattered, left crossing the giants or
lifted out of reach. A stored member of that belt (`belt_roster`, GEN.1070) is one body of the disk, so it belongs to one of those populations with the chance of the population's share of the mass
(the biggest bodies are drawn from the same size law whichever population holds them), and its orbit follows from which one it is:

* `cold_classical`: stirred only as much as the belt's width allows, `COLD_STIRRING` times the rms eccentricity and inclination the whole belt has.
* `hot_classical`: the belt's own rms eccentricity and inclination.
* `resonant`: at the giant's exterior 3:2 or 2:1 resonance (`scatter.RESONANCES`, `scatter.resonance_au`), eccentricity drawn from `RESONANT_E_RANGE`; a belt that does not reach the resonance
  keeps the body cold.
* `scattered_disk`: perihelion between the giant's orbit and the belt's inner edge (the zone the giant scatters), so the eccentricity follows from the semi-major axis.
* `centaurs`: perihelion `CENTAUR_PERIHELION_OVER_GIANT` times the giant's distance, inside it: the orbit crosses the giant's, and the semi-major axis stays in the belt's first stretch.
* `detached`: perihelion `DETACHED_PERIHELION_OVER_INNER_EDGE` times the belt's inner edge, beyond the giant's reach, with the semi-major axis the belt's.

The populations that left the disk (`ejected`, `trojans`, `implanted_inner_belt`, `captured_moons`) hold no member. A belt without a giant has all its members cold. All the numbers are stand-ins
(registry 6.185).
"""

import math

from planetgen.physics import scatter

IN_DISK = ("cold_classical", "hot_classical", "resonant", "scattered_disk", "centaurs", "detached")
"""tuple: The scatter populations whose bodies stay in the disk and so can be stored members."""

COLD_STIRRING = 0.25
"""float: The cold classical members' eccentricity and inclination over the belt's rms (a cold disk stirred only by its own bodies). [Stand-in BPO1, tier Confident educated guess from inferred physics, tweak L4.]"""

RESONANT_E_RANGE = (0.05, 0.35)
"""tuple: The eccentricity a body held in a resonance is drawn from, uniform. [Stand-in BPO2, tier Marginal educated guess from inferred physics, tweak L4.]"""

CENTAUR_PERIHELION_OVER_GIANT = (0.5, 1.0)
"""tuple: The perihelion of a body crossing the giant's orbit, over the giant's distance, uniform. [Stand-in BPO3, tier Marginal educated guess from inferred physics, tweak L4.]"""

DETACHED_PERIHELION_OVER_INNER_EDGE = (1.0, 1.6)
"""tuple: The perihelion of a body lifted out of the giant's reach, over the belt's inner edge, uniform. [Stand-in BPO4, tier Marginal educated guess from inferred physics, tweak L4.]"""

DETACHED_MIN_A_OVER_Q = 1.05
"""float: A detached body's semi-major axis is at least this many times its perihelion (a detached orbit is an eccentric one, not a circle). [Stand-in BPO5, tier T9, tweak L0.]"""

E_MAX_BOUND = 0.9
"""float: The largest eccentricity of a scattered, crossing or detached body (a larger one has its perihelion inside the star's room). [Stand-in BPO6, tier T9, tweak L0.]"""


def member_weights(shares):
    """`(tag, weight)` of the populations that hold stored members, the weights the shares of the scatter step (none: all `cold_classical`)."""
    if not shares:
        return (("cold_classical", 1.0),)
    weights = tuple((tag, max(float(shares.get(tag, 0.0)), 0.0)) for tag in IN_DISK)
    return weights if sum(weight for _, weight in weights) > 0.0 else (("cold_classical", 1.0),)


def pick(weights, u):
    """The population a uniform number `u` (0 to 1) falls in, by the weights."""
    total = sum(weight for _, weight in weights)
    target = u * total
    running = 0.0
    for tag, weight in weights:
        running += weight
        if target < running:
            return tag
    return weights[-1][0]


def _between(low_high, u):
    return low_high[0] + u * (low_high[1] - low_high[0])


def orbit(tag, base_au, lower_au, upper_au, giant_au, u_a, u_e):
    """
    `(tag, a_au, e_scale, e_fixed)` for a body of population `tag` whose draw from the belt's solids gave `base_au`: the semi-major axis the population puts it at and either
    `e_scale` (the multiple of the belt's rms eccentricity to draw it from, `e_fixed` None) or `e_fixed` (the eccentricity itself). A population the belt cannot hold the body in
    (no giant, no resonance inside the belt, no room for the perihelion) is cold: the tag returned is the one that applies.

    Args:
        base_au (float): The body's semi-major axis from the belt's solids.
        lower_au, upper_au (float): The belt's edges.
        giant_au (float or None): The distance of the giant the disk was swept by, if there was one.
        u_a, u_e (float): Two uniform numbers, 0 to 1, for the choices the population makes.
    """
    cold = ("cold_classical", base_au, COLD_STIRRING, None)
    if tag == "cold_classical" or giant_au is None:
        return cold
    if tag == "hot_classical":
        return ("hot_classical", base_au, 1.0, None)
    if tag == "resonant":
        places = [scatter.resonance_au(giant_au, p, q) for p, q in scatter.RESONANCES]
        places = [place for place in places if lower_au <= place <= upper_au]
        if not places:
            return cold
        return ("resonant", places[min(int(u_a * len(places)), len(places) - 1)], None, _between(RESONANT_E_RANGE, u_e))
    if tag == "scattered_disk":
        perihelion = giant_au + u_e * (lower_au - giant_au)
        return ("scattered_disk", base_au, None, _bound(perihelion, base_au))
    if tag == "centaurs":
        a = lower_au * (1.0 + 0.5 * u_a)
        a = min(a, upper_au)
        perihelion = giant_au * _between(CENTAUR_PERIHELION_OVER_GIANT, u_e)
        return ("centaurs", a, None, _bound(perihelion, a))
    if tag == "detached":
        perihelion = lower_au * _between(DETACHED_PERIHELION_OVER_INNER_EDGE, u_e)
        if perihelion * DETACHED_MIN_A_OVER_Q > upper_au:
            return cold
        a = min(max(base_au, perihelion * DETACHED_MIN_A_OVER_Q), upper_au)
        return ("detached", a, None, _bound(perihelion, a))
    raise ValueError(f"unknown population {tag!r}")


def _bound(perihelion_au, a_au):
    """The eccentricity that puts the perihelion at `perihelion_au` for semi-major axis `a_au`, at most `E_MAX_BOUND`, never negative."""
    return min(max(1.0 - perihelion_au / a_au, 0.0), E_MAX_BOUND)
