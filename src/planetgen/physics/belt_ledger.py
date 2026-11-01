# planetgen/physics/belt_ledger.py

"""
The ledger of a belt: which bodies are its own objects, and how much mass is left over for the rest (GEN.1061; docs/design/systems/belt-evolution.md sections 3.1 and 11, belt-sources.md)
====================================================================================================================================================================================

A planetesimal belt or outer disc holds far more bodies than anyone could store. Boss's rule (4:27 PM EDT 2026-10-30) is a **mass cut**: store every body above `L1_MASS_KG`; if those
are fewer than the bodies needed to hold `LEDGER_COVERAGE` of the belt's mass, add the next heaviest until they do, or until `LEDGER_MAX_BODIES` are stored. The belt's mass is a budget that is spent, never exceeded:
the stored bodies plus the `remainder_kg` (everything below the cut, the L2 and L3 tiers) add up to the belt's mass.

1. **The size law.** The cumulative number of bodies above a mass `m` is a piecewise power law, `N(>m) proportional to m^-b` with `b = q / 3` for a cumulative diameter exponent `q` (mass goes as diameter cubed), one exponent
   for each stretch between the breaks. The lowest stretch must have `q < 3` and the highest `q > 3`, so the mass is finite at both ends, and the belt's mass fixes the number at the top break: the mass of each stretch is
   `m dN` integrated (`_table`). Two laws are held: the icy one from the hot Kuiper belt (`ICY_LAW`, [R1439], READ) and a rocky main-belt one (`ROCKY_LAW`, [R2108], READ, and a stand-in top).
   Ivezic et al. give the *differential* index of `dN/dD`: 2.3 below 5 km and 4 from 5 to 40 km; the cumulative exponents used here are one less, 1.3 and 3.0 (the size-law note had them as cumulative, corrected in GEN.1061).
2. **The top end is an order statistic.** The bodies of a belt are a Poisson process in mass, so the heaviest is `m` with `N(>m) = G1`, the next `N(>m) = G2`, and so on, with `Gk` the sum of `k` unit exponential draws. The draw is
   inverse of the law, `mass_at_count`, so rank `k` always uses the `k`-th exponential: **ranks are stable**. The same stream with a lighter belt (an older one) gives the same bodies in the same order, each lighter by the same factor
   (`M^(3 / q_large)` above the top break), which is what lets the evolved belt keep an object's ID.
3. **Mass is conserved.** The draw can by chance put more mass in the stored bodies than the belt has. If the stored bodies come to more than `LEDGER_MAX_SHARE` of the belt, all of them are scaled down together (keeping the order),
   so `remainder_kg` is always at least `1 - LEDGER_MAX_SHARE` of the belt. The scale is rare for a belt with many bodies and it is reported (`Ledger.scaled`).
4. **The largest body and the belt's mass are settled together.** The collisional mass of `belt_mass.belt_mass_earth` depends on `Dc`, the largest body, and the largest body depends on the mass the ledger is drawn for.
   `solve_belt` finds the `Dc` that is the largest body of its own ledger: the ledger of the mass left for `Dc` has a heaviest body of `Dc`. The largest body of the primordial belt is an upper bound (the mass left is never
   more) and the left mass rises with `Dc` more slowly than the body does, so the difference changes sign once and the root is found by bisection on `log Dc` to `SOLVE_TOLERANCE`. (Plain repetition of "mass from `Dc`, `Dc` from mass"
   converges too, at about a tenth of the gap per pass; the design note said "iterate once".) A young belt keeps the primordial mass whatever `Dc` is and is settled in one pass.

The Solar System's main belt (about 10 to 30 bodies above 1e19 kg, Ceres 9.4e20 kg) and the Kuiper belt (the 31 largest bodies about 40 percent of 1.97e-2 Earth masses, [R3539]) are test targets, never inputs.

Limits. One law for the whole belt: a belt of mixed populations (hot and cold classical, scattered) is a sum of these draws, which is the work of the outer-disc items (GEN.1071). The top of the rocky law (the break above 40 km and the slope past it) is a stand-in tuned
to the count of the largest asteroids (BR20, L3). The top body is not limited by what a young disc can grow (Kenyon and Bromley's 1 750 km, BR10, is a check, not a cap). Nothing calls this yet: GEN.1069 stores the ledger and
GEN.1071 and GEN.1092 wire the outer disc and the scatter; they ride merge trains.

Use:

    ledger = draw_ledger(ROCKY_LAW, belt_mass_kg, stream)                          # a Ledger for a given mass
    result = solve_belt(primordial_earth, lower_au, upper_au, star_mass_solar, age_myr, ICY_LAW, stream)   # a BeltLedger: Dc and the mass settled together
"""

import bisect
import math
from dataclasses import dataclass
from functools import lru_cache

from planetgen.physics import belt_mass, constants
from planetgen.util.checks import finite_domain

L1_MASS_KG = 1.0e19
"""float: Every body above this mass is its own stored object (about 200 km across for rock). [Stand-in BR17, tier T9, tweak L0.]"""

LEDGER_COVERAGE = 0.9
"""float: The share of the belt's mass the stored bodies must hold; the heaviest bodies below the cut are added until they do. [Stand-in BR18, tier T9, tweak L0.]"""

LEDGER_MAX_BODIES = 300
"""int: The most bodies stored for one belt. [Stand-in BR18, tier T9, tweak L0.]"""

LEDGER_MAX_SHARE = 0.95
"""float: The most of the belt's mass the stored bodies may hold; if a draw puts more in them they are scaled down together, so the rest of the belt keeps at least `1 - this`. [Stand-in BR21, tier T9, tweak L0.]"""

SOLVE_TOLERANCE = 1.0e-3
"""float: The relative width of the bracket on `Dc` at which `solve_belt` stops. [Stand-in BR21, tier T9, tweak L0.]"""

SOLVE_MIN_KM = 1.0e-3
"""float: The lower end of the bracket on `Dc` of `solve_belt`, far below any belt's largest body. [Stand-in BR21, tier T9, tweak L0.]"""


@dataclass(frozen=True)
class SizeLaw:
    """
    A piecewise cumulative size law of a belt.

    Attributes:
        breaks_km (tuple): The diameters where the exponent changes, ascending.
        exponents (tuple): The cumulative diameter exponent `q` of each stretch, lowest diameters first (one more than `breaks_km`); the first must be under 3 and the last over 3 so the mass is finite.
        density_kg_m3 (float): The bulk density that turns the breaks' diameters and the bodies' masses into one another.
    """

    breaks_km: tuple
    exponents: tuple
    density_kg_m3: float

    def __post_init__(self):
        breaks, exponents = self.breaks_km, self.exponents
        if not (len(exponents) == len(breaks) + 1 and breaks and all(b > 0.0 for b in breaks) and all(x < y for x, y in zip(breaks, breaks[1:]))
                and exponents[0] < 3.0 < exponents[-1] and all(q > 0.0 for q in exponents) and self.density_kg_m3 > 0.0):
            raise ValueError("a size law needs ascending positive breaks, one more positive exponent than breaks, the first under 3 and the last over 3, and a positive density")


ICY_LAW = SizeLaw(breaks_km=(121.0,), exponents=(1.0, 4.35), density_kg_m3=1000.0)
"""SizeLaw: The hot Kuiper belt's: bright-end slope alpha1 = 0.87 (q = 5 alpha = 4.35) and faint-end slope near 0.2 (q = 1.0), break H_B = 7.7, which is 121 km across for an albedo of 0.10 and 1 g/cm3
([R1439], READ; the albedo is an assumption, 99 to 192 km across 0.04 to 0.15). [Stand-in BR19, tier T1, tweak L1.]"""

ROCKY_LAW = SizeLaw(breaks_km=(5.0, 200.0), exponents=(1.3, 3.0, 4.0), density_kg_m3=belt_mass.GRAIN_DENSITY_KG_M3)
"""SizeLaw: The rocky belt's: cumulative q = 1.3 below 5 km and 3.0 from 5 km up ([R2108], astro-ph/0105511, READ: differential indices 2.3 and 4, measured over 0.4 to 5 km and 5 to 40 km), and a
steeper q = 4 above a break at 200 km, a stand-in (a stretch of q = 3 holds equal mass in every doubling of size, so the law must steepen somewhere; the break is chosen so that a belt of the main belt's mass has
some tens of bodies above 1e19 kg). [Stand-in BR20, tier T9, tweak L3.]"""


@lru_cache(maxsize=None)
def _table(law):
    """
    The law's bookkeeping, relative to one body above the top break: `(masses, exponents, relative, mass)` with `masses` the breaks as masses, `exponents` the `b = q / 3` of each stretch, `relative[j]` the
    number above the j-th break over the number above the top break, and `mass` the whole law's mass in the same unit (so a belt's top-break count is its mass over `mass`).
    """
    masses = tuple(math.pi / 6.0 * (d * 1000.0) ** 3 * law.density_kg_m3 for d in law.breaks_km)
    exponents = tuple(q / 3.0 for q in law.exponents)
    top = len(masses) - 1
    relative = [1.0] * len(masses)
    for j in range(top - 1, -1, -1):
        relative[j] = relative[j + 1] * (masses[j] / masses[j + 1]) ** -exponents[j + 1]
    total = 0.0
    for j in range(len(exponents)):
        low = masses[j - 1] if j > 0 else 0.0
        high = masses[j] if j < len(masses) else math.inf
        b = exponents[j]
        anchor = relative[j] * masses[j] ** b if j < len(masses) else masses[top] ** b
        if j == 0:
            total += anchor * b / (1.0 - b) * high ** (1.0 - b)
        elif j == len(masses):
            total += anchor * b / (b - 1.0) * low ** (1.0 - b)
        elif abs(b - 1.0) < 1e-12:
            total += anchor * b * math.log(high / low)
        else:
            total += anchor * b / (1.0 - b) * (high ** (1.0 - b) - low ** (1.0 - b))
    return masses, exponents, tuple(relative), total


@finite_domain()
def top_count(law, belt_mass_kg):
    """The number of bodies above the top break, so that the law's mass adds up to the belt's mass."""
    if not belt_mass_kg > 0.0:
        raise ValueError("the belt's mass must be positive")
    return belt_mass_kg / _table(law)[3]


@finite_domain()
def count_above(law, belt_mass_kg, mass_kg):
    """`N(>m)`: the expected number of bodies heavier than `mass_kg` in a belt of this mass."""
    if not mass_kg > 0.0:
        raise ValueError("the mass must be positive")
    masses, exponents, relative, _ = _table(law)
    j = bisect.bisect_right(masses, mass_kg)           # the stretch the mass is in: j breaks lie below it
    if j == len(masses):
        share = (mass_kg / masses[-1]) ** -exponents[j]
    else:
        share = relative[j] * (mass_kg / masses[j]) ** -exponents[j]
    return top_count(law, belt_mass_kg) * share


@finite_domain()
def mass_at_count(law, belt_mass_kg, count):
    """The inverse of `count_above`: the mass with `count` bodies expected above it (`count = 1` is the typical largest body)."""
    if not count > 0.0:
        raise ValueError("the count must be positive")
    masses, exponents, relative, _ = _table(law)
    share = count / top_count(law, belt_mass_kg)
    if share <= 1.0:
        return masses[-1] * share ** (-1.0 / exponents[-1])
    j = len(masses) - 1
    while j > 0 and share > relative[j - 1]:
        j -= 1
    return masses[j] * (share / relative[j]) ** (-1.0 / exponents[j])


@finite_domain()
def diameter_km(mass_kg, density_kg_m3):
    """The diameter of a sphere of this mass and density."""
    if not (mass_kg > 0.0 and density_kg_m3 > 0.0):
        raise ValueError("the mass and the density must be positive")
    return (6.0 * mass_kg / (math.pi * density_kg_m3)) ** (1.0 / 3.0) / 1000.0


@dataclass(frozen=True)
class Ledger:
    """
    The stored bodies of a belt.

    Attributes:
        masses_kg (tuple): The stored bodies' masses, heaviest first; the index is the rank, stable for a stream whatever the belt's mass.
        belt_mass_kg (float): The belt's mass, the budget.
        remainder_kg (float): The mass of everything not stored (belt mass minus the bodies); at least `1 - LEDGER_MAX_SHARE` of the belt.
        above_cut (int): How many bodies are above `L1_MASS_KG` before any scaling.
        scaled (bool): Whether the bodies had to be scaled down to stay inside `LEDGER_MAX_SHARE`.
        law (SizeLaw): The size law the draw followed.
    """

    masses_kg: tuple
    belt_mass_kg: float
    remainder_kg: float
    above_cut: int
    scaled: bool
    law: SizeLaw

    @property
    def covered(self):
        """The share of the belt's mass in the stored bodies."""
        return 1.0 - self.remainder_kg / self.belt_mass_kg

    @property
    def largest_km(self):
        """The diameter of the heaviest body."""
        return diameter_km(self.masses_kg[0], self.law.density_kg_m3)


def draw_exponentials(stream):
    """The belt's `LEDGER_MAX_BODIES` unit exponential draws, the same number from `stream` whatever the belt (so a run's stream does not depend on the belt's mass)."""
    return [-math.log(1.0 - stream.random()) for _ in range(LEDGER_MAX_BODIES)]


@finite_domain()
def ledger_from(law, belt_mass_kg, exponentials):
    """The ledger of a belt of this mass from the `draw_exponentials` of a stream (the heaviest body is rank 0)."""
    raw, gap = [], 0.0
    for e in exponentials:
        gap += e
        raw.append(mass_at_count(law, belt_mass_kg, gap))
    above = sum(1 for m in raw if m >= L1_MASS_KG)
    n = min(above, LEDGER_MAX_BODIES)
    held = sum(raw[:n])
    while n < LEDGER_MAX_BODIES and held < LEDGER_COVERAGE * belt_mass_kg:
        held += raw[n]
        n += 1
    scale = min(1.0, LEDGER_MAX_SHARE * belt_mass_kg / held) if held > 0.0 else 1.0
    masses = tuple(m * scale for m in raw[:n])
    return Ledger(masses, belt_mass_kg, belt_mass_kg - sum(masses), above, scale < 1.0, law)


def draw_ledger(law, belt_mass_kg, stream):
    """The ledger of a belt of `belt_mass_kg` kg whose bodies follow `law`, drawn from `stream` (a `draw.Stream`; `LEDGER_MAX_BODIES` draws of `random()`)."""
    return ledger_from(law, belt_mass_kg, draw_exponentials(stream))


@dataclass(frozen=True)
class BeltLedger:
    """
    A belt's ledger settled with its collisional mass.

    Attributes:
        ledger (Ledger): The stored bodies for the settled mass.
        mass_earth (float): The belt's mass, Earth masses.
        largest_km (float): `Dc`, the diameter of the heaviest body, which the mass was worked out for (to `SOLVE_TOLERANCE`).
        passes (int): How many ledgers were drawn on the way.
    """

    ledger: Ledger
    mass_earth: float
    largest_km: float
    passes: int


@finite_domain()
def solve_belt(primordial_earth, lower_au, upper_au, star_mass_solar, age_myr, law, stream):
    """
    The belt's ledger and mass, with `Dc` and `belt_mass.belt_mass_earth` settled together (step 4 above).

    Args:
        primordial_earth (float): What formed between the limits, Earth masses (`belt_mass.primordial_mass_earth`).
        lower_au, upper_au (float): The belt's limits.
        star_mass_solar (float): The star's mass.
        age_myr (float): The belt's age.
        law (SizeLaw): `ICY_LAW` or `ROCKY_LAW`.
        stream: An object with `random()` (a `draw.Stream`); `LEDGER_MAX_BODIES` draws are taken, whatever the belt.

    Returns:
        BeltLedger
    """
    exponentials = draw_exponentials(stream)
    kg_per_earth = constants.EARTH_MASS_TO_KG
    passes = 0

    def left_kg(largest_km):
        return belt_mass.belt_mass_earth(primordial_earth, lower_au, upper_au, star_mass_solar, age_myr, largest_km) * kg_per_earth

    def ledger_for(largest_km):
        nonlocal passes
        passes += 1
        return ledger_from(law, left_kg(largest_km), exponentials)

    high = ledger_from(law, primordial_earth * kg_per_earth, exponentials).largest_km
    passes += 1
    young = ledger_for(high)
    if abs(young.belt_mass_kg - primordial_earth * kg_per_earth) <= SOLVE_TOLERANCE * young.belt_mass_kg:
        return BeltLedger(young, young.belt_mass_kg / kg_per_earth, young.largest_km, passes)
    low = SOLVE_MIN_KM
    while high / low > 1.0 + SOLVE_TOLERANCE:
        middle = math.sqrt(low * high)
        if ledger_for(middle).largest_km > middle:
            low = middle
        else:
            high = middle
    ledger = ledger_for(math.sqrt(low * high))
    return BeltLedger(ledger, ledger.belt_mass_kg / kg_per_earth, ledger.largest_km, passes)
