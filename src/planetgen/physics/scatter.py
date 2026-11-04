# planetgen/physics/scatter.py

"""
The scatter step: where the solids of an outer disc end up once the giants have acted on it (GEN.1089, docs/design/systems/belt-origins.md section 6, registry 6.83)
=====================================================================================================================================================================

A cold outer disc beyond the last giant is one flat ring until a giant acts on it. This module is the rule that says what the giant did, as shares of the disc's mass,
from the giant's mass, radius and orbit, how far it moved outward and whether the system had an instability. It is a rule over the architecture's events with stored
mass shares and no N-body run (Boss, 2026-10-30). The Solar System is a test target only.

**The rule.**

1. *Swept band.* A giant of mass ratio `q = m / M*` at `a` has a chaotic zone of half-width `w = ZONE_COEFFICIENT a q^(2/7)` ([R4727], [R] recalled) in which bodies are
   scattered in a few orbits. An outward-migrating giant sweeps the band it crossed, from its starting `a - migration` to `a`, and the zone with it: the disc inside
   `[a - migration - w, a + w]` is the swept mass. A giant that never moved sweeps only its zone.
2. *Resonances.* The giant's outward migration carries its exterior resonances `P:1` outward through the disc (Malhotra): the 3:2 at `a (3/2)^(2/3)` = 1.31 `a` and the 2:1 at 1.59 `a`.
   `RESONANCE_CAPTURE` of the mass in each swept span (beyond the zone) is captured, the resonant share.
3. *What is left* beyond the band and the resonances stays on cold, low-eccentricity orbits (`cold_classical`); after an instability a share (`HOT_SHARE` times the instability,
   0 to 1) is excited, the `hot_classical` population.
4. *The swept mass* is ejected in the proportion `1 - exp(-Theta / THETA_EJECT)`, `Theta = q a / R` the Safronov number of the giant (the escape speed squared over twice the orbital
   speed squared at its orbit: where it is large a scattered body is thrown out, where it is near 1 it stays bound), and the rest is split between the
   `scattered_disc` (perihelion near the giant), the `centaurs` (crossing the giants), the `detached` (perihelion lifted out of reach), the `trojans` (1:1), the share implanted in the
   inner belt, and the share captured as irregular moons (`SWEPT_SPLIT`).

A giant-free disc keeps all its mass cold. The shares sum to 1. Defaults pending Boss and the registry rows (6.83, SX1 to SX5) say what is stand-in.

5. *Losses of the architecture* (GEN.1119, registry 6.98). Every body the architecture loses, the chain's sub-minimum bodies, the planets a giant sweeps or throws out and the mass a
   merge does not keep, is shared among `ejected`, `collided` and `bound` by `loss_fates` from the Safronov number at the place of the encounter, `Theta = q d / R` with the
   perturber's mass ratio `q`, its radius `R` and the body's own distance `d` (the escape speed over the orbital speed there: Theta far below 1 a scattered body mostly collides or
   stays, far above it is thrown out; [R3867], recalled). The ratio of ejection to collision goes as `Theta^EJECTION_THETA_EXPONENT`, and a body the perturber does not reach (outside
   `reach_band_au`) stays bound. No draw is made: the shares are expected values and the rogue builder (GEN.1111) draws the discrete fate.

Use:

    shares = scatter_shares(star_mass_solar, giant_mass_earth, giant_au, giant_radius_earth, migration_au, disc_inner_au, disc_edge_au, instability=0.0)
    history_line(shares, giant_au, migration_au)
"""

import math

from planetgen.physics import constants, unit_registry as units

_mearth_to_kg = units.converter("Mearth", "kg")
_msun_to_kg = units.converter("Msun", "kg")

ZONE_COEFFICIENT = 1.3
"""float: The half-width of the chaotic zone over `a q^(2/7)` ([R4727], recalled from the design's table, row SX1)."""

DISC_SLOPE = 1.5
"""float: The surface density of the disc goes as `a^-DISC_SLOPE` (the minimum-mass nebula, [R1854], recalled; row SX1)."""

RESONANCES = ((3, 2), (2, 1))
"""tuple: The exterior resonances `(p, q)` that capture (row SX2)."""

RESONANCE_CAPTURE = 0.25
"""float: The share of the mass in a swept resonance span that is captured (row SX2, a stand-in)."""

THETA_EJECT = 10.0
"""float: The Safronov number at which `1 - 1/e` of the swept mass is ejected (row SX3, a stand-in)."""

HOT_SHARE = 0.4
"""float: The share of the unswept cold disc that is excited, hot classical, after a full instability (row SX4, a stand-in)."""

SWEPT_SPLIT = (("scattered_disc", 0.62), ("centaurs", 0.25), ("detached", 0.10), ("trojans", 0.002), ("implanted_inner_belt", 0.02), ("captured_moons", 0.008))
"""tuple: How the swept mass that is not ejected is shared (sums to 1, row SX5, a stand-in)."""

EJECTION_THETA_EXPONENT = 2.0
"""float: Ejection to collision goes as Theta to this power ([R3867], recalled; the cross-section for a hit grows as 1 + Theta, a throw-out needs the escape speed, row EJ1)."""

LOSS_BOUND_WEIGHT = 1.0
"""float: The weight of a scattered body that is neither hit nor thrown out, against a collision weight of 1 (row EJ1, a stand-in)."""

DISC_LIFETIME_YR = 3.0e6
"""float: The age at which the gas disc is gone and the giants have finished moving, the time of an architecture loss ([R1750], recalled; row EJ2)."""

POPULATIONS = ("cold_classical", "hot_classical", "resonant", "scattered_disc", "centaurs", "detached", "trojans", "ejected", "implanted_inner_belt", "captured_moons")
"""tuple: The keys of `scatter_shares`."""

_AU = constants.AU_TO_M
_R_EARTH_M = constants.EARTH_RADIUS_KM * 1.0e3


def mass_ratio(giant_mass_earth, star_mass_solar):
    """`q = m / M*`."""
    return _mearth_to_kg(giant_mass_earth) / (_msun_to_kg(star_mass_solar))


def zone_half_width_au(giant_au, giant_mass_earth, star_mass_solar):
    """`w = 1.3 a q^(2/7)`, AU: the half-width of the giant's chaotic zone."""
    return ZONE_COEFFICIENT * giant_au * mass_ratio(giant_mass_earth, star_mass_solar) ** (2.0 / 7.0)


def safronov_number(giant_mass_earth, giant_radius_earth, giant_au, star_mass_solar):
    """`Theta = q a / R`, the Safronov number of the giant at its orbit."""
    return mass_ratio(giant_mass_earth, star_mass_solar) * giant_au * _AU / (giant_radius_earth * _R_EARTH_M)


def resonance_au(giant_au, p, q):
    """The exterior `p:q` resonance of a giant, AU: `a (p/q)^(2/3)`."""
    return giant_au * (p / q) ** (2.0 / 3.0)


def disc_fraction(lower_au, upper_au, disc_inner_au, disc_edge_au, slope=DISC_SLOPE):
    """The share of the disc's mass between `lower_au` and `upper_au` (clipped to the disc) for a surface density `a^-slope`."""
    low, high = max(lower_au, disc_inner_au), min(upper_au, disc_edge_au)
    if high <= low:
        return 0.0

    def integral(x):
        power = 2.0 - slope
        return math.log(x) if abs(power) < 1e-9 else x ** power / power

    return (integral(high) - integral(low)) / (integral(disc_edge_au) - integral(disc_inner_au))


def scatter_shares(star_mass_solar, giant_mass_earth, giant_au, giant_radius_earth, migration_au, disc_inner_au, disc_edge_au, instability=0.0):
    """
    The shares of an outer disc's mass in each of `POPULATIONS`, summing to 1.

    Args:
        star_mass_solar (float): The star's mass.
        giant_mass_earth, giant_au, giant_radius_earth (float): The outermost giant's mass, final orbit and radius; a mass of 0 or less means no giant.
        migration_au (float): How far outward the giant moved (0 or negative: it did not move outward).
        disc_inner_au, disc_edge_au (float): The disc's extent.
        instability (float): 0 to 1, how much the architecture's giants scattered each other (a Nice-type event is 1).
    """
    shares = dict.fromkeys(POPULATIONS, 0.0)
    if giant_mass_earth <= 0.0:
        shares["cold_classical"] = 1.0
        return shares
    width = zone_half_width_au(giant_au, giant_mass_earth, star_mass_solar)
    start = giant_au - max(migration_au, 0.0)
    swept = disc_fraction(start - width, giant_au + width, disc_inner_au, disc_edge_au)
    resonant = 0.0
    if migration_au > 0.0:
        for p, q in RESONANCES:
            low, high = resonance_au(start, p, q), resonance_au(giant_au, p, q)
            resonant += RESONANCE_CAPTURE * disc_fraction(max(low, giant_au + width), high, disc_inner_au, disc_edge_au)
    resonant = min(resonant, 1.0 - swept)
    left = 1.0 - swept - resonant
    hot = left * HOT_SHARE * max(0.0, min(instability, 1.0))
    ejected = swept * (1.0 - math.exp(-safronov_number(giant_mass_earth, giant_radius_earth, giant_au, star_mass_solar) / THETA_EJECT))
    bound = swept - ejected
    shares.update(cold_classical=left - hot, hot_classical=hot, resonant=resonant, ejected=ejected)
    for key, share in SWEPT_SPLIT:
        shares[key] = bound * share
    return shares


def scattered_share(shares):
    """The share of the disc that ends scattered, ejected or implanted: everything but the cold, hot classical and resonant populations."""
    return 1.0 - shares["cold_classical"] - shares["hot_classical"] - shares["resonant"]


def history_line(shares, giant_au, migration_au):
    """One sentence for the belt page."""
    percent = round(100.0 * scattered_share(shares))
    if migration_au > 0.0:
        return f"when the outer giant moved from {giant_au - migration_au:.0f} to {giant_au:.0f} AU it scattered {percent} percent of the disc"
    if percent > 0:
        return f"the outer giant at {giant_au:.0f} AU never moved but its chaotic zone scattered {percent} percent of the disc"
    return "no giant reached the disc, which kept all its solids on cold orbits"


def reach_band_au(giant_au, formation_au, giant_mass_earth, star_mass_solar):
    """
    `(low, high)` AU, the stretch a giant reaches: its chaotic zone about its final orbit, widened to the whole path it moved along when it formed
    at `formation_au` (it scattered what it passed, as in `scatter_shares`).
    """
    width = zone_half_width_au(giant_au, giant_mass_earth, star_mass_solar)
    return min(giant_au, formation_au) - width, max(giant_au, formation_au) + width


def loss_fates(theta):
    """`(ejected, collided, bound)`, shares summing to 1, of a body a perturber of Safronov number `theta` scatters: weights `theta^EJECTION_THETA_EXPONENT`, 1 and `LOSS_BOUND_WEIGHT`."""
    ejected = max(theta, 0.0) ** EJECTION_THETA_EXPONENT
    total = ejected + 1.0 + LOSS_BOUND_WEIGHT
    return ejected / total, 1.0 / total, LOSS_BOUND_WEIGHT / total


def hill_radius_au(mass_earth, distance_au, star_mass_solar):
    """The Hill radius of a body, AU: `d (m / 3 M*)^(1/3)`; the moons that can be kept through an encounter sit inside a fraction of it."""
    return distance_au * (mass_ratio(mass_earth, star_mass_solar) / 3.0) ** (1.0 / 3.0)
