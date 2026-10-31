# planetgen/physics/resonance_gaps.py

"""
The map a belt keeps of its giants: resonance zones, chaotic zones, the size below which Yarkovsky drift has emptied it, and how stirred it is (GEN.1062; docs/design/systems/belt-evolution.md sections 2.3, 2.4 and 4 steps 2 and 5, belt-sources.md BR1 to BR3)
==========================================================================================================================================================================================================================================================

A belt that has giants beside it is not the annulus it formed as. Where a giant's orbit period is a small-number ratio of a body's, the body is shaken out; within a distance of the giant's own orbit nothing lasts; and the bodies small
enough to be moved by sunlight drift into those places and go. This is a pure function of the giants' orbits and the belt's limits and age. Nothing calls it yet (the belt wiring, GEN.1069 and later, reads the map; the clearance ledger GEN.1014 can use the same zones).

1. **Resonance zones.** A body whose period is a small-number ratio of a giant's sits at a resonance. The label `p:q` has the longer period first: `3:1` is the giant's period over the body's for a body inside the giant's orbit
   (`a = a_giant (q / p)^(2/3)`; Jupiter's gives 2.50 AU) and the body's over the giant's for a body outside it (`a = a_giant (p / q)^(2/3)`; Neptune's `3:2` gives 39.4 AU). `order` is `p - q`. Only resonances of `MAX_ORDER` or less are listed.
2. **Width.** A first-order resonance has the width `(da / a)_res = 3.2 mu^(2/3) j^(1/3)` for particles of low eccentricity and `7.78 (e mu j)^(1/2)` above (Mustill and Wyatt 2012, arXiv:1110.1282, READ, their Eqs. 3 and 9; `j` is the larger number
   of the ratio, `mu` the giant over the star); the larger of the two is taken. A resonance of higher order is narrower by a factor `e^((order - 1) / 2)`, the order dependence of a resonance's libration width (`RGP2`, a stand-in: no paper was opened for the coefficient).
3. **Chaotic zones.** A giant clears the zone round its orbit of half-width `da / a = 1.3 mu^(2/7)` for particles of low eccentricity (Wisdom 1980, as Mustill and Wyatt confirm by numerical integration: 1.35 mu^0.278) and `1.8 (e mu)^(1/5)` above the
   critical eccentricity `0.21 mu^(3/7)` (their Eq. 10 and 11); above `e_max = 2.1 mu^(1/4)` the zone is where orbits cross, `da / a = e`. The belt's `survivors` are its intervals outside every giant's zone.
4. **Clearing time.** A resonance zone is empty once the belt is older than `RESONANCE_CLEARING_MYR` (the design's "1 to 10 Myr", Minton and Malhotra 2010 for the later, logarithmic tail, `belt-sources.md` BR3); a younger belt keeps them.
5. **The Yarkovsky cut.** Sunlight moves a body by `da/dt = 9.3e-4 au km / My / D x (1.1 au / a)^2` (Vokrouhlicky et al. 2015, READ: Bennu's measured drift scaled with 1/D and the flux, `belt-sources.md` BR1). Bodies that drift as far as half the spacing between
   the belt's barriers (its edges after the chaotic zones, and the resonance zones) within the belt's age are fed into them and lost: `cut_diameter_km` is the diameter at which that distance is just reached, so it grows with the age (`D_cut` proportional to `t`)
   and a body above it stays. The Solar System's belt gives about 5 km; bodies above 40 km do not move in 4.5 Gyr.
6. **Stirring.** A giant on an eccentric orbit forces eccentricity on the belt: `e_forced = (5/4) (a / a_giant) e_giant / (1 - e_giant^2)` (the leading secular term; Wyatt et al. 1999, recalled, not opened, `RGP5`), the strongest giant's at the belt's middle (`min(a / a_giant, a_giant / a)` stands for the ratio on either side). The belt's eccentricity
   is the larger of that and the width of its own annulus `(upper - lower) / (upper + lower)` (what `belt_mass` already uses); the zones and resonance widths are worked out for the giants' forcing alone (at least `ZONE_ECCENTRICITY_FLOOR`), not for the annulus's spread, which for a wide belt is not a free eccentricity.

The Solar System is a test target: Jupiter's 3:1, 5:2, 7:3 and 2:1 at 2.50, 2.82, 2.96 and 3.27 AU, and Neptune's chaotic zone (27.7 to 32.4 AU) and 3:2 and 2:1 at 39.4 and 47.7 AU, are checked in `tests/test_resonance_gaps.py`, never fed in.

Limits. The widths are the resonance's libration widths (the region, not the part cleared of bodies, which is smaller: the 2:1 comes out about 0.4 AU wide against a Kirkwood gap of 0.1 to 0.2); resonances of higher order than `MAX_ORDER`, the
secular resonances (the nu6 that empties the inner belt), the sweeping of resonances while the giants migrate, three-body resonances and the resonances of two giants together are not here. A resonance that holds bodies stably (the 3:2 Hildas) is listed as a resonance;
whether it clears or traps depends on the eccentricity and is not decided. Bodies inside a resonance are not removed from the survivors: the zones are listed, the survivors are cut only at the chaotic zones.

Use:

    giants = [Giant(a_au, mass_kg, eccentricity), ...]
    belt = belt_map(giants, star_mass_kg, lower_au, upper_au, age_myr)      # a BeltMap
    belt.resonances, belt.zones, belt.survivors, belt.cut_diameter_km, belt.eccentricity
"""

import math
from dataclasses import dataclass

from planetgen.util.checks import finite_domain

CHAOTIC_ZONE_COEFFICIENT = 1.3
"""float: The half-width of a giant's chaotic zone over its distance is this times `mu^(2/7)` for particles of low eccentricity (Wisdom 1980; Mustill and Wyatt 2012: 1.3, numerically 1.35). [Stand-in RGP1, tier T1, tweak L0.]"""

EXTENDED_ZONE_COEFFICIENT = 1.8
"""float: The same above the critical eccentricity: `1.8 (e mu)^(1/5)` (Mustill and Wyatt 2012, Eq. 10). [Stand-in RGP1, tier T1, tweak L0.]"""

CRITICAL_ECCENTRICITY_COEFFICIENT = 0.21
"""float: The critical eccentricity is this times `mu^(3/7)` (Mustill and Wyatt 2012, Eq. 11). [Stand-in RGP1, tier T1, tweak L0.]"""

MAX_ECCENTRICITY_COEFFICIENT = 2.1
"""float: Above this times `mu^(1/4)` particles are removed by close encounters (Mustill and Wyatt 2012, Eq. 12). [Stand-in RGP1, tier T1, tweak L0.]"""

LOW_ECCENTRICITY_WIDTH = 3.2
"""float: A first-order resonance's width over its distance is this times `mu^(2/3) j^(1/3)` for low eccentricity (Wisdom 1980; Mustill and Wyatt 2012, Eq. 3). [Stand-in RGP1, tier T1, tweak L0.]"""

HIGH_ECCENTRICITY_WIDTH = 7.78
"""float: The same above the critical eccentricity: this times `(e mu j)^(1/2)` (Mustill and Wyatt 2012, Eq. 9). [Stand-in RGP1, tier T1, tweak L0.]"""

ORDER_WIDTH_EXPONENT = 0.5
"""float: A resonance of order `k` is narrower than a first-order one by `e^(EXPONENT (k - 1))`: the libration width of a resonance goes as `e^(k/2)`. [Stand-in RGP2, tier T6, tweak L1.]"""

ZONE_ECCENTRICITY_FLOOR = 0.005
"""float: The least free eccentricity the zones and widths are worked out for, a cold belt's (the giants' forcing replaces it when larger). [Stand-in RGP7, tier T6, tweak L1.]"""

MAX_ORDER = 4
"""int: Resonances of this order (`p - q`) and less are listed. [Numerical, stand-in RGP4, tier T9, tweak L0.]"""

MAX_P = 12
"""int: The largest `p` of a listed resonance; a bound on the list for a belt close to a giant, where the chaotic zone is the better description. [Numerical, stand-in RGP4, tier T9, tweak L0.]"""

RESONANCE_CLEARING_MYR = 10.0
"""float: A belt older than this has emptied its resonance zones (the design note's "1 to 10 Myr"; Minton and Malhotra 2010 for the tail after it). [Stand-in RGP3, tier T2, tweak L1.]"""

YARKOVSKY_COEFFICIENT_AU_KM_PER_MYR = 9.3e-4
"""float: `da/dt` times the diameter, au km/Myr, at `YARKOVSKY_REFERENCE_AU` (Bennu: 493 m across, 18.95e-4 au/Myr at 1.10 au; Vokrouhlicky et al. 2015, arXiv:1502.01249, READ). [Stand-in BR1 of 6.95, tier T2, tweak L1; row RGP6.]"""

YARKOVSKY_REFERENCE_AU = 1.1
"""float: The distance the drift coefficient is measured at; the drift is taken to go as the flux, `1/a^2`. [Stand-in RGP6, tier T2, tweak L1.]"""

FORCED_ECCENTRICITY_COEFFICIENT = 1.25
"""float: The eccentricity a giant forces on a belt: this times `(a / a_giant) e_giant / (1 - e_giant^2)` (the leading secular term; Wyatt et al. 1999, recalled, not opened). [Stand-in RGP5, tier T2, tweak L0.]"""

INTERIOR = "interior"
EXTERIOR = "exterior"


@dataclass(frozen=True)
class Giant:
    """
    A giant planet that acts on the belt.

    Attributes:
        a_au (float): Its semi-major axis.
        mass_kg (float): Its mass.
        eccentricity (float): Its eccentricity.
    """

    a_au: float
    mass_kg: float
    eccentricity: float = 0.0

    def __post_init__(self):
        if not (self.a_au > 0.0 and self.mass_kg > 0.0 and 0.0 <= self.eccentricity < 1.0):
            raise ValueError("a giant needs a positive distance and mass and an eccentricity in [0, 1)")


@dataclass(frozen=True)
class Resonance:
    """One mean-motion resonance zone of a giant.

    Attributes:
        giant (int): The index of the giant.
        label (str): `p:q`, the longer period first.
        p, q (int): The ratio.
        order (int): `p - q`.
        side (str): `interior` (the body is inside the giant's orbit) or `exterior`.
        a_au (float): The resonance's distance from the star.
        width_au (float): The width of the zone.
    """

    giant: int
    label: str
    p: int
    q: int
    order: int
    side: str
    a_au: float
    width_au: float

    @property
    def inner_au(self):
        return self.a_au - 0.5 * self.width_au

    @property
    def outer_au(self):
        return self.a_au + 0.5 * self.width_au


@dataclass(frozen=True)
class ChaoticZone:
    """The zone round a giant's orbit that nothing survives in.

    Attributes:
        giant (int): The index of the giant.
        inner_au, outer_au (float): Its edges.
        relative_half_width (float): `da / a`.
        regime (str): `classical` (low eccentricity), `extended` or `crossing`.
    """

    giant: int
    inner_au: float
    outer_au: float
    relative_half_width: float
    regime: str


@dataclass(frozen=True)
class BeltMap:
    """
    What the giants leave of a belt.

    Attributes:
        resonances (tuple): The `Resonance` zones inside the belt, by distance.
        zones (tuple): The `ChaoticZone`s that overlap the belt.
        survivors (tuple): The `(low, high)` intervals of the belt outside every chaotic zone.
        gaps_open (bool): Whether the belt is old enough that its resonance zones are empty.
        spacing_au (float): The median distance between the belt's barriers.
        cut_diameter_km (float): The Yarkovsky cut: smaller bodies have drifted into a barrier and are gone.
        forced_eccentricity (float): The eccentricity the giants force at the belt's middle.
        eccentricity (float): The belt's eccentricity: the larger of that and its annulus's.
    """

    resonances: tuple
    zones: tuple
    survivors: tuple
    gaps_open: bool
    spacing_au: float
    cut_diameter_km: float
    forced_eccentricity: float
    eccentricity: float


def _mu(giant, star_mass_kg):
    return giant.mass_kg / star_mass_kg


@finite_domain()
def chaotic_zone(mu, eccentricity):
    """The half-width of a giant's chaotic zone over its distance, and the regime, for a mass ratio `mu` and particle eccentricity: `(da / a, regime)`."""
    if not (mu > 0.0 and eccentricity >= 0.0):
        raise ValueError("the mass ratio must be positive and the eccentricity not negative")
    classical = CHAOTIC_ZONE_COEFFICIENT * mu ** (2.0 / 7.0)
    critical = CRITICAL_ECCENTRICITY_COEFFICIENT * mu ** (3.0 / 7.0)
    if eccentricity <= critical:
        return classical, "classical"
    if eccentricity > MAX_ECCENTRICITY_COEFFICIENT * mu ** 0.25:
        return max(eccentricity, EXTENDED_ZONE_COEFFICIENT * (eccentricity * mu) ** 0.2), "crossing"
    return EXTENDED_ZONE_COEFFICIENT * (eccentricity * mu) ** 0.2, "extended"


@finite_domain()
def resonance_width(mu, p, q, eccentricity):
    """The width of the `p:q` resonance over its distance for a mass ratio `mu` and particle eccentricity (first-order law, narrowed by `ORDER_WIDTH_EXPONENT` for higher orders)."""
    if not (mu > 0.0 and p > q > 0 and eccentricity >= 0.0):
        raise ValueError("the mass ratio must be positive, p over q above 1 and the eccentricity not negative")
    low = LOW_ECCENTRICITY_WIDTH * mu ** (2.0 / 3.0) * p ** (1.0 / 3.0)
    high = HIGH_ECCENTRICITY_WIDTH * math.sqrt(eccentricity * mu * p)
    width = max(low, high)
    return width * eccentricity ** (ORDER_WIDTH_EXPONENT * (p - q - 1))


@finite_domain()
def drift_au_per_myr(diameter_km, distance_au):
    """The Yarkovsky drift of a body of this diameter at this distance, au per Myr (unsigned)."""
    if not (diameter_km > 0.0 and distance_au > 0.0):
        raise ValueError("the diameter and the distance must be positive")
    return YARKOVSKY_COEFFICIENT_AU_KM_PER_MYR / diameter_km * (YARKOVSKY_REFERENCE_AU / distance_au) ** 2


@finite_domain()
def cut_diameter_km(distance_au, age_myr, reach_au):
    """The diameter whose Yarkovsky drift over `age_myr` is `reach_au` at `distance_au`: smaller bodies reach it."""
    if not (distance_au > 0.0 and age_myr >= 0.0 and reach_au > 0.0):
        raise ValueError("the distance and reach must be positive and the age not negative")
    return YARKOVSKY_COEFFICIENT_AU_KM_PER_MYR * (YARKOVSKY_REFERENCE_AU / distance_au) ** 2 * age_myr / reach_au


def _ratios():
    out = []
    for p in range(2, MAX_P + 1):
        for q in range(1, p):
            if p - q <= MAX_ORDER and math.gcd(p, q) == 1:
                out.append((p, q))
    return out


@finite_domain()
def forced_eccentricity(giants, distance_au):
    """The largest eccentricity any giant forces at this distance (the leading secular term)."""
    if not distance_au > 0.0:
        raise ValueError("the distance must be positive")
    values = [FORCED_ECCENTRICITY_COEFFICIENT * min(distance_au / g.a_au, g.a_au / distance_au) * g.eccentricity / (1.0 - g.eccentricity ** 2) for g in giants]
    return max(values, default=0.0)


def belt_map(giants, star_mass_kg, lower_au, upper_au, age_myr):
    """
    The resonance zones, chaotic zones, survivors, Yarkovsky cut and eccentricity of a belt between two distances with these giants.

    Args:
        giants (sequence of Giant): The planets that act on the belt (typically the giants; every body in the sequence is treated as a perturber).
        star_mass_kg (float): The star's mass.
        lower_au, upper_au (float): The belt's limits, as formed.
        age_myr (float): The belt's age.

    Returns:
        BeltMap
    """
    if not (0.0 < lower_au < upper_au and star_mass_kg > 0.0 and age_myr >= 0.0):
        raise ValueError("the belt needs 0 < lower < upper, a star with mass and an age not negative")
    middle = 0.5 * (lower_au + upper_au)
    annulus = (upper_au - lower_au) / (upper_au + lower_au)
    forced = forced_eccentricity(giants, middle)
    eccentricity = max(annulus, forced)
    particles = max(forced, ZONE_ECCENTRICITY_FLOOR)       # the zones and widths are for the particles' free eccentricity, not the annulus's spread
    zones, resonances = [], []
    for index, giant in enumerate(giants):
        mu = _mu(giant, star_mass_kg)
        half, regime = chaotic_zone(mu, particles)
        inner, outer = giant.a_au * (1.0 - half), giant.a_au * (1.0 + half)
        if outer > lower_au and inner < upper_au:
            zones.append(ChaoticZone(index, inner, outer, half, regime))
        for p, q in _ratios():
            for side, a in ((INTERIOR, giant.a_au * (q / p) ** (2.0 / 3.0)), (EXTERIOR, giant.a_au * (p / q) ** (2.0 / 3.0))):
                if lower_au <= a <= upper_au and not (inner <= a <= outer):
                    width = resonance_width(mu, p, q, particles) * a
                    resonances.append(Resonance(index, f"{p}:{q}", p, q, p - q, side, a, width))
    resonances.sort(key=lambda r: (r.a_au, r.giant))
    survivors, start = [], lower_au
    for zone in sorted(zones, key=lambda z: z.inner_au):
        if zone.inner_au > start:
            survivors.append((start, min(zone.inner_au, upper_au)))
        start = max(start, zone.outer_au)
    if start < upper_au:
        survivors.append((start, upper_au))
    barriers = sorted({edge for low, high in survivors for edge in (low, high)} | {r.a_au for r in resonances})
    spacings = sorted(b - a for a, b in zip(barriers, barriers[1:]))
    spacing = spacings[len(spacings) // 2] if spacings else upper_au - lower_au
    reach = 0.5 * spacing
    cut = cut_diameter_km(middle, age_myr, reach) if reach > 0.0 else 0.0
    return BeltMap(tuple(resonances), tuple(zones), tuple(survivors), age_myr >= RESONANCE_CLEARING_MYR, spacing, cut, forced, eccentricity)
