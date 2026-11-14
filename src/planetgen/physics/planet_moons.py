# planetgen/physics/planet_moons.py

"""
Moons of planet size around the giants (GEN.1157, stand-ins PSM1 to PSM4 of docs/design/research/gap-filling-plan.md section 6.276)
=========================================================================================================================================

How often a giant planet holds a moon as big as a planet (`PLANET_MOON_MIN_MASS_KG`, a tenth of an Earth mass or more), and how big, from the two ways a giant can get one. A pure module (slice 0), read by `planets.generate_moons` and
`planet_kind.moon_mass_ceiling_kg` since the wiring slice (an epoch). Nothing is stored: the captured moon is an ordinary moon row with origin `exchange`. Design: docs/design/systems/planet-size-moons.md. Report: reports/gen1157-planet-size-moons-2026-11-12.md.

* **The circumplanetary disk** makes regular moons whose total is about 1e-4 to 6e-4 of the planet (`moon_composition.total_moon_mass_ratio` at the ends of `ALPHA_OVER_F`, [R0642], [R0643:
  Eq. 4.15]). No disk moon can pass that share of its parent, so `disk_moon_ceiling_kg` is the ceiling the generator holds a giant's moons to (`planet_kind.moon_mass_ceiling_kg`); a planet-size moon of the disk needs a parent of
  a few Jupiter masses (the whole disk of a Jupiter holds 0.2 Earth masses), the super-Jovians of [R1899]. A disk moon is drawn with the disk, not here.
* **Capture** (a giant that meets a binary, keeps one and flings the other away, [R0037, R5229]; or a planet pulled in during a migration instability, [R5230]; the candidate
  Kepler-1625 b-i of [R4322] is the case studied, and [R1900] read its signal another way) is rare and not tied to the disk: `exchange_chance` of a giant has one, `exchange_ceiling_kg` is the most it can hold (a share of the parent, `EXCHANGE_MAX_RATIO`) and
  `draw_exchange_moon` draws its mass from the generator's own moon-mass law (`tuning.MOON_MASS_FUNCTION_SLOPE`) between `PLANET_MOON_MIN_MASS_KG` and that ceiling.
* **Giant impact** is the rocky parents' route (`tuning.ROCKY_IMPACT_MOON_CHANCE`, GEN.962) and is not repeated here.

The moon itself is then made by the ordinary body generator (no special case, same classes and life rules as a planet); this module gives only its mass and origin.

Use:

    disk_moon_ceiling_kg(10 * JUPITER_MASS_KG)       # 1.2e25 kg, about two Earth masses
    exchange_chance(5 * JUPITER_MASS_KG, hill_radius_m, orbital_speed_m_s)         # 0.014, or 0 for a parent no partner can leave a moon to
    draw_exchange_moon(JUPITER_MASS_KG, hill_radius_m, orbital_speed_m_s, limit_mass_kg, rng)     # PlanetMoon("exchange", mass) or None
"""

import collections
import math

from planetgen import tuning
from planetgen.physics import moon_composition
from planetgen.physics.constants import EARTH_MASS_TO_KG, G

JUPITER_MASS_KG = moon_composition.JUPITER_MASS_KG

PLANET_MOON_MIN_MASS_KG = 0.1 * EARTH_MASS_TO_KG
"""float: The lightest moon that counts as a moon of planet size: a tenth of an Earth mass (Mars is 0.107). [Stand-in PSM1, tier "from needed output", tweak L0.]"""

DISK_TOTAL_RATIO_MAX = moon_composition.total_moon_mass_ratio(moon_composition.ALPHA_OVER_F[1])
"""float: The most the moons of one circumplanetary disk can hold of their parent: Canup and Ward's total at the top of the alpha/f range, 6.4e-4 ([R0643: Eq. 4.15]). [Stand-in PSM2, tier T1, tweak L0.]"""

CAPTURE_CHANCE = 1.0 / 70.0
"""float: The share of giants that hold a moon of planet size from exchange capture: one candidate, Kepler-1708 b-i, in the 70 cool giants of the survey [R2365]. A candidate, not a confirmed moon, and
the survey found no more, so this is a loose upper figure. [Stand-in PSM3, tier "marginal unverified", tweak L2.]"""

EXCHANGE_ENCOUNTER_HILL_FRACTION = 19.1 / 261.0
"""float: The closest approach of the binary to the parent (`b` of Eq. 1 to 3 of [R5229]) in the parent's Hill radius: the 19.1 Jupiter radii of the worked case of Kepler-1625 b (10 Jupiter masses at 0.87 au round 1.079 solar masses, a Hill
radius of 261 Jupiter radii). One encounter distance for every parent, since the paper gives one case. [Stand-in PSM4, tier "educated guess from physics", tweak L1, replaces the L3 cap `EXCHANGE_MAX_RATIO` of GEN.1157.]"""

EXCHANGE_V_INFINITY_FRACTION = 1.0
"""float: The speed of the binary relative to the parent at infinity (`v_inf` of Eq. 2 of [R5229]) in the parent's speed round its star: [R5229] tests 1 to 30 km/s for Kepler-1625 b, 30 km/s being its Keplerian speed, the hardest case
(a faster encounter needs a heavier partner), so the bound it gives is the careful end. [Stand-in PSM5, tier "educated guess from physics", tweak L0.]"""

APOAPSIS_HILL_FRACTION = 0.5
"""float: The widest apoapsis a captured moon keeps, in the parent's Hill radius: [R5229] (Eq. 3) takes a stable orbit as one whose apoapsis is under half the Hill radius, the same share `moon_orbit_bounds_km` gives the moons. [Derived, tier T3.]"""

PlanetMoon = collections.namedtuple("PlanetMoon", ("origin", "mass_kg"))
"""A moon of planet size: its `origin` ("exchange") and its mass in kg."""


def disk_moon_ceiling_kg(parent_mass_kg):
    """The most mass the moons of a circumplanetary disk can hold between them, kg: `DISK_TOTAL_RATIO_MAX` of the parent. No single disk moon is heavier.

    Raises:
        ValueError: The mass is not positive.
    """
    if not parent_mass_kg > 0.0:
        raise ValueError(f"the parent mass {parent_mass_kg!r} is not positive")
    return DISK_TOTAL_RATIO_MAX * parent_mass_kg


def encounter_limit_kg(parent_mass_kg, partner_mass_kg, hill_radius_m, orbital_speed_m_s):
    """
    Eq. 1 of [R5229] (the framework of Williams 2013), the heaviest mass a parent captures from a binary encounter whose other member leaves with `partner_mass_kg`, kg:

        M_cap < 3 M_p (pi G M_esc / (2 b v_enc dv))^(3/2) - M_esc,

    with the encounter distance `b` a fixed share of the Hill radius (`EXCHANGE_ENCOUNTER_HILL_FRACTION`), the speed at infinity `v_inf` the parent's own speed round the star (`EXCHANGE_V_INFINITY_FRACTION`), the speed at
    the encounter `v_enc = sqrt(v_esc^2 + v_inf^2)` over the parent's escape speed at `b`, and `dv` the least kick that leaves the moon bound on an orbit whose semimajor axis `a` is half of (half the Hill radius plus `b`)
    (Eq. 2 and 3 of [R5229]). 0 for a parent with no Hill radius or a partner too light. Checked against [R5229]: a 17 Earth-mass moon of Kepler-1625 b needs a 5 to 50 Earth-mass partner at 1 to 30 km/s.
    """
    if not (parent_mass_kg > 0.0 and partner_mass_kg > 0.0 and hill_radius_m > 0.0):
        return 0.0
    b = EXCHANGE_ENCOUNTER_HILL_FRACTION * hill_radius_m
    v_inf = EXCHANGE_V_INFINITY_FRACTION * orbital_speed_m_s
    v_esc_squared = 2.0 * G * parent_mass_kg / b
    v_enc = math.sqrt(v_esc_squared + v_inf * v_inf)
    a = 0.5 * (APOAPSIS_HILL_FRACTION * hill_radius_m + b)
    kick = v_enc - math.sqrt(G * parent_mass_kg * (2.0 / b - 1.0 / a))
    if kick <= 0.0:       # no kick is needed: the encounter itself leaves the moon bound, and the equation has no limit
        return math.inf
    bound = 3.0 * parent_mass_kg * (math.pi * G * partner_mass_kg / (2.0 * b * v_enc * kick)) ** 1.5 - partner_mass_kg
    return max(bound, 0.0)


def exchange_bound_kg(parent_mass_kg, partner_mass_kg, hill_radius_m, orbital_speed_m_s):
    """
    The heaviest moon a parent keeps from a binary encounter, kg: `encounter_limit_kg` (Eq. 1 of [R5229]) of the parent and the partner that leaves, and never more than the partner: the lighter member of a binary is the one
    kept, the heavier carries off the energy ([R5229] needs a partner of 5 to 50 Earth masses to leave a Neptune-mass moon). 0 for a parent with no Hill radius.
    """
    return min(encounter_limit_kg(parent_mass_kg, partner_mass_kg, hill_radius_m, orbital_speed_m_s), max(partner_mass_kg, 0.0))


def minimum_partner_kg(parent_mass_kg, hill_radius_m, orbital_speed_m_s):
    """
    The lightest partner that lets a parent keep a moon of planet size (`PLANET_MOON_MIN_MASS_KG`), kg, or `None` where even a partner as heavy as the parent cannot (a fast encounter, a tiny Hill sphere).
    The bound rises with the partner, so a bisection in the logarithm finds it.
    """
    if not (parent_mass_kg > 0.0 and hill_radius_m > 0.0):
        return None
    low, high = PLANET_MOON_MIN_MASS_KG, parent_mass_kg
    if high <= low or exchange_bound_kg(parent_mass_kg, high, hill_radius_m, orbital_speed_m_s) < PLANET_MOON_MIN_MASS_KG:
        return None
    log_low, log_high = math.log(low), math.log(high)
    for _ in range(60):
        middle = 0.5 * (log_low + log_high)
        if exchange_bound_kg(parent_mass_kg, math.exp(middle), hill_radius_m, orbital_speed_m_s) < PLANET_MOON_MIN_MASS_KG:
            log_low = middle
        else:
            log_high = middle
    return math.exp(log_high)


def exchange_chance(parent_mass_kg, hill_radius_m, orbital_speed_m_s):
    """The chance that a giant of this mass holds a moon of planet size from exchange capture: `CAPTURE_CHANCE`, or 0 where no partner no heavier than the parent can leave one (`minimum_partner_kg`)."""
    if not parent_mass_kg > 0.0:
        raise ValueError(f"the parent mass {parent_mass_kg!r} is not positive")
    return CAPTURE_CHANCE if minimum_partner_kg(parent_mass_kg, hill_radius_m, orbital_speed_m_s) is not None else 0.0


def _law_mass_kg(low_kg, high_kg, u):
    """A mass between `low_kg` and `high_kg` from the uniform `u`: `dN/dlogM ~ M^-s` (`tuning.MOON_MASS_FUNCTION_SLOPE`, the generator's own moon-mass law), by its inverse distribution."""
    if not 0.0 <= u <= 1.0:
        raise ValueError(f"the uniform {u!r} is outside 0..1")
    if high_kg <= low_kg:
        return low_kg
    slope = tuning.MOON_MASS_FUNCTION_SLOPE
    low_term, high_term = low_kg ** -slope, high_kg ** -slope
    return (low_term + u * (high_term - low_term)) ** (-1.0 / slope)


def draw_exchange_moon(parent_mass_kg, hill_radius_m, orbital_speed_m_s, limit_mass_kg, rng):
    """
    The moon of planet size a giant of `parent_mass_kg` holds by exchange capture, or `None`. Always three numbers from `rng` (anything with `random()`): the roll against `exchange_chance`, then the partner of the binary
    (its mass from the generator's moon-mass law between the lightest partner that leaves a moon of planet size and the parent), then the moon (from the same law between `PLANET_MOON_MIN_MASS_KG` and the bound of the
    encounter, `exchange_bound_kg`, held to `limit_mass_kg`, the heaviest moon the generator lets the parent hold), so the stream never moves with the parent.

    Returns:
        PlanetMoon or None: The moon, or `None` for no capture (including a parent no partner can leave a moon to).
    """
    roll, partner_roll, mass_roll = rng.random(), rng.random(), rng.random()
    lightest = minimum_partner_kg(parent_mass_kg, hill_radius_m, orbital_speed_m_s)
    if lightest is None or roll >= CAPTURE_CHANCE:
        return None
    partner = _law_mass_kg(lightest, parent_mass_kg, partner_roll)
    ceiling = min(exchange_bound_kg(parent_mass_kg, partner, hill_radius_m, orbital_speed_m_s), limit_mass_kg)
    return PlanetMoon("exchange", _law_mass_kg(PLANET_MOON_MIN_MASS_KG, ceiling, mass_roll))
