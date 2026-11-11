# planetgen/physics/moon_resonances.py

"""
Resonances between a planet's moons: capture on the way out, and the eccentricity they force (GEN.1449 part 2)
==============================================================================================================

`moon_orbits.evolve` moves every moon under its planet's tide and leaves the moons blind to each other. The tide drives an inner moon outward faster than an outer one
(`da/dt` goes as `a^(-11/2)`), so two moons that start apart approach in period ratio, and a pair that reaches a first-order resonance of 2:1, 3:2 or 4:3 with little eccentricity
is **captured**: the inner moon pushes the outer along as it drifts, the ratio stays where it was caught and the two keep a forced eccentricity (Io and Europa, Europa and Ganymede, the
Saturnian Enceladus and Dione). Without it the model lets two moons close to a period ratio of 1.1, which no moon system shows.

**Capture** (`capture`). A pair has *swept* a resonance `p:q` when the ratio of its birth periods lay above `p/q` (with the offset `CAPTURE_OFFSET`) and the ratio the integrator ends
with lies below it; the first resonance met on the way in is the largest `p/q` below the birth ratio. A swept pair is captured and sits at `p/q (1 + CAPTURE_OFFSET)`. The pair
(and a chain of pairs, Io, Europa and Ganymede) keeps the total orbital angular momentum the integrator gave it, `sum m sqrt(a)` (the planet's spin is unchanged), shared by the
resonance condition `sqrt(a_k) = K_k sqrt(a_0)` with `K` the product of the ratios' cube roots up the chain, so `sqrt(a_0) = sum m sqrt(a) / sum m K`: a light outer moon is pushed out almost
for free and a heavy one holds the inner moon back. If that would hold the inner moon inside the place it was born (the outer moon is too heavy to be pushed) the pair is not captured.

**The forced eccentricity** (`give_forced_eccentricities`). Every pair of moons within `resonant_chains.RESONANCE_TOLERANCE` of a first-order ratio, captured or just close, is given the
eccentricity `resonant_chains.pair_forced_eccentricities` finds (the Lithwick, Xie and Wu forced eccentricity, two bodies being enough). A moon already more eccentric keeps its own. A moon
flexed by it has its tide made again (a rocky one's air and heat, `refinish_rocky`), so Io's heat follows from Europa.

**What it does not do.** Second-order resonances (5:3) capture nothing and force nothing; a pair is captured with certainty (the adiabatic limit for a first-order resonance at a low
eccentricity; there is no escape by eccentricity or inclination, no libration amplitude and no three-body resonance); the angular momentum shared inside a captured pair is the orbital
`m sqrt(a)` of the integrator's finals, not the integral of the shared migration; and the pulls of the moons on each other outside a resonance stay out of the tide's history. Stand-ins
are TS9 to TS12 of section 6.260 of the gap plan.
"""

import math
from dataclasses import dataclass

from planetgen.physics import resonant_chains

CAPTURE_OFFSET = 0.002
"""float: How far wide of exact resonance a captured pair sits, as a share of `p/q` (Io and Europa are within 0.001, Enceladus and Dione within 0.003 of 2:1; the Kepler multis
sit at 0.01 to 0.02 because disks and stellar tides add detuning). Source tier T6, tweak level L1."""

FIRST_ORDER = tuple((p, q) for p, q in sorted(resonant_chains.TARGET_RATIOS, key=lambda pq: -pq[0] / pq[1]) if p - q == 1)
"""tuple: The resonances that capture, `p:q` by falling ratio (2:1, 3:2, 4:3)."""


@dataclass(frozen=True)
class Capture:
    """A captured pair: `inner` and `outer` are their places in the list given, in the resonance `p:q`."""

    inner: int
    outer: int
    p: int
    q: int


def swept_resonance(birth_ratio, final_ratio):
    """The `(p, q)` a pair swept, going from the period ratio `birth_ratio` to `final_ratio`, or `None`: the first first-order resonance the falling ratio meets."""
    for p, q in FIRST_ORDER:
        edge = p / q * (1.0 + CAPTURE_OFFSET)
        if birth_ratio > edge > final_ratio:
            return p, q
    return None


def capture(masses_kg, birth_distances, distances):
    """
    Where moons end once the pairs that swept a first-order resonance are held in it.

    Args:
        masses_kg (sequence): The moons' masses, in order of the final distance.
        birth_distances (sequence): Their semi-major axes when born (any length unit).
        distances (sequence): The semi-major axes the integrator ended with, rising (the same unit).

    Returns:
        tuple: `(distances, captures)`: the new distances in the same unit and order, and the `Capture`s made.
    """
    count = len(distances)
    swept = []
    for i in range(count - 1):
        birth_ratio = (birth_distances[i + 1] / birth_distances[i]) ** 1.5
        final_ratio = (distances[i + 1] / distances[i]) ** 1.5
        swept.append(swept_resonance(birth_ratio, final_ratio))
    result, made = list(distances), []
    i = 0
    while i < count - 1:
        if swept[i] is None:
            i += 1
            continue
        j = i
        while j < count - 1 and swept[j] is not None:
            j += 1
        members = list(range(i, j + 1))
        factors = [1.0]
        for k in range(i, j):
            p, q = swept[k]
            factors.append(factors[-1] * (p / q * (1.0 + CAPTURE_OFFSET)) ** (1.0 / 3.0))
        momentum = sum(masses_kg[m] * math.sqrt(distances[m]) for m in members)
        inner_root = momentum / sum(masses_kg[m] * f for m, f in zip(members, factors))
        if inner_root * inner_root >= birth_distances[i] * (1.0 - 1.0e-9):
            for m, f in zip(members, factors):
                result[m] = (inner_root * f) ** 2
            made += [Capture(k, k + 1, *swept[k]) for k in range(i, j)]
        i = j
    return result, tuple(made)


def forced_eccentricities(planet):
    """The eccentricity the neighbouring moons of `planet` force on each, in the order of `planet.moons` (0.0 for a moon near no first-order ratio), before the moon keeps its own if larger."""
    moons = list(planet.moons)
    if len(moons) < 2 or not planet.mass > 0.0:
        return [0.0] * len(moons)
    usable = [m for m in moons if getattr(m, "period", None) and m.mass and m.mass > 0.0]
    if len(usable) < 2:
        return [0.0] * len(moons)
    found = resonant_chains.pair_forced_eccentricities([m.period * 365.25 for m in usable], [m.mass for m in usable], planet.mass)
    by_id = {id(m): e for m, e in zip(usable, found)}
    return [by_id.get(id(m), 0.0) for m in moons]


def give_forced_eccentricities(planet):
    """Raises the eccentricity of each moon of `planet` to the one its resonant neighbours force on it (a moon already above keeps its own) and puts the moon on its new ellipse; returns the moons that
    moved by more than `planets.ECCENTRICITY_SETTLE_TOLERANCE`, which have their tide made again by the final-age refresh."""
    from planetgen.physics import planets
    moved = []
    for moon, forced in zip(list(planet.moons), forced_eccentricities(planet)):
        if forced - (moon.eccentricity or 0.0) <= planets.ECCENTRICITY_SETTLE_TOLERANCE:
            continue
        moon.eccentricity = forced
        planets.update_orbital_position(moon)
        if getattr(moon, "_age_drawn", None) is not None:
            moon._age_drawn = 0.0      # the air, class and heat of a rocky moon are made again once, on the orbit it ends with (`planets.stale_air`)
        moon._resonance_e = forced
        moved.append(moon)
    return moved
