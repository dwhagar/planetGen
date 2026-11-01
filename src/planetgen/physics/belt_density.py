# planetgen/physics/belt_density.py

"""
How crowded a belt is: objects per square AU above a stated size, and mass per square AU (GEN.1093)
=====================================================================================================

The asteroid belt information page (`belt-evolution.md` section 12) gives each belt a **density factor**: how many objects above 1, 10 and 100 km a square astronomical
unit of the belt holds, and how much mass. This module is that calculation as pure functions of a belt's stored quantities. Nothing calls it yet (UX.240, the belt page
data builder, reads it; GEN.1061's ledger will hand it the individually kept bodies).

**The size law.** The population below the individually kept bodies follows a broken cumulative power law, `N(>D)` proportional to `D^-q` in each stretch between breaks,
the two stretches joined where they meet (`SizeLaw`: the floor, the breaks and the exponents). Three laws ship as stand-ins: the main belt ([R2108], differential exponents 2.3 from
0.4 to 5 km and 4 above, so `q` 1.3 and 3.0), and the hot and cold classical Kuiper belt ([R1439], `q = 5 alpha`: 4.35 or 7.5 above the break and 1.0 below, the breaks 121 and 175 km at
albedo 0.10). The law gives the **shape**; the number of objects comes from the mass.

**The ledger.** A belt's stored mass is the budget (`belt-evolution.md` section 3): the bodies kept one by one (every one above the cut rule's mass) are spent first, what is
left is the tail, and the law is normalised so that the tail, from the law's floor up to the smallest kept body, holds exactly that mass. Then the count above a
diameter is the bodies kept above it (exact, they are listed) plus the tail's count above it, and the mass of the listed bodies and the tail adds to the belt's mass with
nothing over and nothing short. Without a ledger the tail runs up to `largest_km` (the law's top end).

**What the density is.** The belt's mean over its annulus, `pi (outer^2 - inner^2)` square AU: the generator gives a belt no radial profile, so a position in the belt has the
mean. A diameter under the law's floor is answered with the count at the floor (and flagged), since the law says nothing about dust.

Units: objects per square AU is the registry unit `au-2` and mass per square AU `kg_au2` (`unit_registry`), so no number shows without its unit.
Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.99.
"""

import math
from dataclasses import dataclass

from planetgen.util.checks import finite_domain

OBJECTS_UNIT = "au-2"
"""str: The `unit_registry` key of an object count per square AU."""
MASS_UNIT = "kg_au2"
"""str: The `unit_registry` key of a mass per square AU."""

DEFAULT_DIAMETERS_KM = (1.0, 10.0, 100.0)
"""tuple: The diameters the page states the density above (Boss, 4:37 PM EDT 2026-10-30). Source tier T9, tweak level L0."""
BULK_DENSITY_KG_M3 = 2500.0
"""float: The rock density that turns a mass into a diameter where the caller gives none (the design note's, `belt-evolution.md` section 11.1). Source tier T6, tweak level L1."""
_Q_THREE_TOLERANCE = 1e-9


class BeltDensityError(ValueError):
    """An input the model cannot use."""


@dataclass(frozen=True)
class SizeLaw:
    """
    A broken cumulative power law, `N(>D)` proportional to `D^-q` between breaks.

    Attributes:
        name (str): What it is the law of.
        breaks_km (tuple): Ascending diameters: the floor (below it the law says nothing), the breaks, and the top.
        slopes (tuple): The `q` of each stretch between neighbouring breaks, one fewer than the breaks.
    """

    name: str
    breaks_km: tuple
    slopes: tuple

    def __post_init__(self):
        if len(self.breaks_km) < 2 or len(self.slopes) != len(self.breaks_km) - 1:
            raise BeltDensityError("a size law has a floor, a top and one exponent for each stretch between its breaks")
        if any(not math.isfinite(b) or b <= 0.0 for b in self.breaks_km) or any(b >= c for b, c in zip(self.breaks_km, self.breaks_km[1:])):
            raise BeltDensityError("the breaks of a size law must be positive and ascending")
        if any(not math.isfinite(q) or q <= 0.0 for q in self.slopes):
            raise BeltDensityError("the exponents of a size law must be positive")

    @property
    def floor_km(self):
        return self.breaks_km[0]

    @property
    def top_km(self):
        return self.breaks_km[-1]


MAIN_BELT = SizeLaw("main belt", (0.4, 5.0, 1000.0), (1.3, 3.0))
"""SizeLaw: [R2108] (abstract READ): `D^-2.3` for 0.4 to 5 km and `D^-4` for 5 to 40 km, with 530,000 objects above 1 km. The abstract does not say whether the exponents
are of the cumulative or the differential distribution; the paper's section 6 defines the index for `dN/dD` (a summary of the PDF, to be checked verbatim), so here they are
differential and the cumulative `q` is one less: 1.3 and 3.0 (the repo's size-law table took them as cumulative). Read as cumulative the law leaves 3 asteroids above 40 km
where the belt has some hundreds; read as differential, about 130. The 3.0 is carried to the top, where the ledger takes over. Source tier T1 below 40 km, T6 above, tweak level L0."""
KUIPER_HOT = SizeLaw("hot Kuiper belt", (1.0, 121.0, 2500.0), (1.0, 4.35))
"""SizeLaw: [R1439] (READ), hot population: `alpha` 0.87 bright and about 0.2 faint, break `H` 7.7 (121 km at albedo 0.10), cumulative `q = 5 alpha`; the faint side is
set to 1.0 from the table's note that the faint exponents are about 1.0 ([R4107] give 0.75 below 1 to 2 km). Source tier T1 bright side, T6 faint side, tweak level L0."""
KUIPER_COLD = SizeLaw("cold Kuiper belt", (1.0, 175.0, 1500.0), (1.0, 7.5))
"""SizeLaw: [R1439] (READ), cold classical population: `alpha` 1.5, break `H` 6.9 (175 km at albedo 0.10), `q = 7.5` bright side. Source tier T1 bright side, T6 faint side, tweak level L0."""


@dataclass(frozen=True)
class Count:
    """
    The objects above one diameter.

    Attributes:
        diameter_km (float): The diameter asked.
        objects (float): The number above it in the whole belt.
        per_au2 (float): The same per square AU of the annulus (`OBJECTS_UNIT`).
        listed (int): Of `objects`, the bodies kept one by one.
        below_floor (bool): The diameter was under the law's floor and the count is the one at the floor.
    """

    diameter_km: float
    objects: float
    per_au2: float
    listed: int
    below_floor: bool


@dataclass(frozen=True)
class BeltDensity:
    """
    The density factor of a belt.

    Attributes:
        area_au2 (float): The annulus.
        mass_per_au2_kg (float): The belt's mass over it (`MASS_UNIT`).
        listed_mass_kg, tail_mass_kg (float): The bodies kept one by one and the tail; they add to the belt's mass (the tail is zero if the listed bodies are all of it).
        tail_top_km (float): The diameter the tail runs up to.
        counts (tuple): A `Count` for each diameter asked.
    """

    area_au2: float
    mass_per_au2_kg: float
    listed_mass_kg: float
    tail_mass_kg: float
    tail_top_km: float
    counts: tuple


class _Tail:
    """The law between its floor and `top_km`, normalised to `mass_kg`. `shape(D)` is the law's cumulative count above `D` with 1 at the top, so the tail holds
    `scale * (shape(D) - 1)` objects above `D`."""

    def __init__(self, law, top_km, mass_kg, density_kg_m3):
        self.floor_km = law.floor_km
        self.top_km = top_km
        self.pieces = []                                         # (lower, upper, q, shape at upper), from the top down
        shape_upper = 1.0
        for lower, upper, q in reversed(list(zip(law.breaks_km, law.breaks_km[1:], law.slopes))):
            if lower >= top_km:
                continue
            upper = min(upper, top_km)
            self.pieces.append((lower, upper, q, shape_upper))
            shape_upper *= (upper / lower) ** q
        self.scale = 0.0
        if self.pieces and mass_kg > 0.0:
            self.scale = mass_kg / (math.pi / 6.0 * density_kg_m3 * 1e9 * self.mass_shape(self.floor_km, top_km))

    def shape(self, diameter_km):
        for lower, upper, q, shape_upper in self.pieces:
            if diameter_km >= lower:
                return shape_upper * (min(diameter_km, upper) / upper) ** -q
        return self.pieces[-1][3] * (self.floor_km / self.pieces[-1][1]) ** -self.pieces[-1][2] if self.pieces else 1.0

    def objects_above(self, diameter_km):
        """The tail's objects above a diameter (clamped to the tail's range)."""
        if not self.scale:
            return 0.0
        return self.scale * (self.shape(min(max(diameter_km, self.floor_km), self.top_km)) - 1.0)

    def mass_shape(self, low_km, high_km):
        """The integral of `D^3 dN` between two diameters, km^3, in the units of `scale`."""
        total = 0.0
        for lower, upper, q, shape_upper in self.pieces:
            a, b = max(lower, low_km), min(upper, high_km)
            if a >= b:
                continue
            at_b = shape_upper * (b / upper) ** -q
            exponent = 3.0 - q
            ratio = math.log(b / a) if abs(exponent) < _Q_THREE_TOLERANCE else -math.expm1(-exponent * math.log(b / a)) / exponent
            total += q * at_b * b ** 3 * ratio
        return total

    def mass_between(self, low_km, high_km, density_kg_m3):
        """The tail's mass between two diameters, kg."""
        return self.scale * math.pi / 6.0 * density_kg_m3 * 1e9 * self.mass_shape(low_km, high_km)


def _diameter_km(mass_kg, density_kg_m3):
    return (6.0 * mass_kg / (math.pi * density_kg_m3)) ** (1.0 / 3.0) / 1e3


def _tail(law, belt_mass_kg, density_kg_m3, listed_diameters_km, listed_masses_kg, largest_km):
    if len(listed_diameters_km) != len(listed_masses_kg):
        raise BeltDensityError("a mass for each listed diameter")
    if any(not d > 0.0 or not math.isfinite(d) for d in listed_diameters_km) or any(not m > 0.0 or not math.isfinite(m) for m in listed_masses_kg):
        raise BeltDensityError("listed diameters and masses must be positive")
    listed_mass = math.fsum(listed_masses_kg)
    top = min(listed_diameters_km) if listed_diameters_km else (largest_km if largest_km is not None else law.top_km)
    if not top > 0.0:
        raise BeltDensityError("the largest body must have a positive size")
    return listed_mass, _Tail(law, min(top, law.top_km) if not listed_diameters_km else top, max(belt_mass_kg - listed_mass, 0.0), density_kg_m3)


@finite_domain()
def belt_density(inner_au, outer_au, belt_mass_kg, law=MAIN_BELT, density_kg_m3=BULK_DENSITY_KG_M3, listed_diameters_km=(), listed_masses_kg=(),
                 diameters_km=DEFAULT_DIAMETERS_KM, largest_km=None):
    """
    The density factor of a belt: the objects above each diameter and the mass, per square AU of its annulus.

    Args:
        inner_au, outer_au (float): The belt's limits.
        belt_mass_kg (float): Its stored mass (the listed bodies and the tail together).
        law (SizeLaw): The size law of the tail.
        density_kg_m3 (float): The bulk density that gives the mass of an object of a diameter.
        listed_diameters_km, listed_masses_kg (sequence): The bodies kept one by one (the ledger's), by diameter and mass.
        diameters_km (sequence): The diameters to count above.
        largest_km (float, optional): With no listed bodies, the diameter the tail runs up to (default the law's top).

    Returns:
        BeltDensity: The density, with a `Count` for each diameter asked.
    """
    if not 0.0 <= inner_au < outer_au:
        raise BeltDensityError("a belt runs from its inner to its outer limit")
    if belt_mass_kg < 0.0 or not density_kg_m3 > 0.0:
        raise BeltDensityError("the mass must not be negative and the density must be positive")
    area = math.pi * (outer_au ** 2 - inner_au ** 2)
    listed_mass, tail = _tail(law, belt_mass_kg, density_kg_m3, listed_diameters_km, listed_masses_kg, largest_km)
    tail_mass = tail.mass_between(tail.floor_km, tail.top_km, density_kg_m3) if tail.scale else 0.0
    counts = []
    for diameter in diameters_km:
        if not diameter > 0.0:
            raise BeltDensityError("a diameter must be positive")
        listed = sum(1 for d in listed_diameters_km if d > diameter)
        objects = listed + tail.objects_above(diameter)
        counts.append(Count(float(diameter), objects, objects / area, listed, diameter < law.floor_km))
    return BeltDensity(area, belt_mass_kg / area, listed_mass, tail_mass, tail.top_km, tuple(counts))


@finite_domain()
def tail_bin(low_km, high_km, belt_mass_kg, law=MAIN_BELT, density_kg_m3=BULK_DENSITY_KG_M3, listed_diameters_km=(), listed_masses_kg=(), largest_km=None):
    """
    The tail between two diameters, for the page's size table: `(objects, mass_kg)` of the objects that are not kept one by one. The bins of a table that cover the law
    from its floor to the tail's top add up to the tail, and with the listed bodies to the belt's mass.
    """
    if not 0.0 < low_km < high_km:
        raise BeltDensityError("a bin runs from a smaller to a larger diameter")
    _, tail = _tail(law, belt_mass_kg, density_kg_m3, listed_diameters_km, listed_masses_kg, largest_km)
    if not tail.scale:
        return 0.0, 0.0
    return tail.objects_above(low_km) - tail.objects_above(high_km), tail.mass_between(low_km, high_km, density_kg_m3)


def diameter_of(mass_kg, density_kg_m3=BULK_DENSITY_KG_M3):
    """The diameter in km of a sphere of the mass and density."""
    if not mass_kg > 0.0 or not density_kg_m3 > 0.0:
        raise BeltDensityError("mass and density must be positive")
    return _diameter_km(mass_kg, density_kg_m3)
