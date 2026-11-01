# planetgen/physics/rogue_collision.py

"""
The collision of two rogue planets (GEN.1099; docs/design/systems/interstellar-asteroid-fields.md 2 and 3, docs/design/research/interstellar-field-sources.md)
==================================================================================================================================================

Pure functions, no draws. Two rogues almost never meet (about one collision in the life of the Galaxy), and when they do, what is left depends on the speed:
this module gives the rate, the speed at which the pair breaks up, what is left and how fast the pieces fly. The collision law is [R2633], read in full for GEN.1098 (stand-in rows IF1 to IF8 of `gap-filling-plan.md` section 6.103):

* **Rate.** `Gamma = n sigma v` per rogue, with `sigma = pi (R1 + R2)^2 (1 + v_esc^2 / v^2)` (`cross_section_m2`, gravitational focusing with the escape speed of the pair
  at touching); for the Galaxy `0.5 sigma v (rogues per star)^2 integral(n_star^2 dV)` (`galaxy_rate_per_year`). A relative speed of Maxwellian dispersion `s` weights
  each speed by `sigma v` (`speed_fractions`, closed form), which gives the shares of the collisions below, in and above a speed window.
* **Disruption.** The catastrophic disruption energy of strengthless planets is `Q*_RD = c* (4/5) pi rho1 G RC1^2` for equal masses, with `c* = 1.9`, `rho1 = 1000 kg/m3` and
  `RC1` the radius of the combined mass at `rho1`; other mass ratios scale as `V* = ((gamma + 1)^2 / (4 gamma))^(1/(3 mu_bar)) V*_(gamma=1)` (`disruption_speed_ms`). Two
  Earths: `V*` is 22.6 km/s. An oblique impact needs more: only the part of the projectile that overlaps the target counts (`interacting_share`,
  `oblique_disruption_energy_j_kg`), and the angle comes from `sin theta cos theta` (`impact_angle_cdf`).
* **What is left.** The universal law `Mlr / Mtot = 0.5 - 0.5 (Q_R / Q*_RD - 1)` between `Q_R = 0` and `1.8 Q*_RD`, a power law of index -1.5 beyond (`largest_remnant_fraction`);
  the second-largest remnant and a cumulative tail of index 2.85 (`second_remnant_fraction`); `regime` names the outcome.
* **How fast.** The mass of the pieces per unit speed (in escape speeds of the combined mass) falls as `10^(A - S v)` with `A = 0.3 - 0.3 Mlr / Mtot`; `unbound_fraction` is the part
  that leaves the pair's escape speed and `unbound_rms_speed_ms` the root mean square of the speed it keeps at infinity, the dispersion of the field it makes.
* **Shock.** `peak_shock_pressure_pa`: the planar impact of two equal rocks, particle speed half the impact speed, the linear Hugoniot of forsterite.

What the law does not have: melt and vapour (it is a gravity-only N-body fit at speeds of up to a few km/s, used here far outside that range: IF14), and the transition from merging
to hit-and-run of a grazing impact (a grazing impact below the disruption energy is taken as hit-and-run, IF14). `rubble_speed_window_ms` is the window of a large remnant plus rubble, from `V*` to the speed at which the
largest remnant falls to a tenth of the mass. Nothing reads the module yet: the field host rule and the columns (GEN.1103 and GEN.1104) wire it.
"""

import math

from planetgen.physics import constants

RHO1_KG_M3 = 1000.0
"""float: The reference density of the disruption law: `RC1` is the radius of the combined mass at this density ([R2633]; IF3, tier T1, L0)."""
C_STAR_PLANET = 1.9
"""float: The disruption offset of strengthless planet-size bodies, in units of the binding energy (1.9 +/- 0.3; IF3)."""
MU_BAR_PLANET = 0.36
"""float: The coupling exponent of planet-size bodies (0.36 +/- 0.01; IF3). Between 1/3 (momentum) and 2/3 (energy)."""
C_STAR_RUBBLE = 5.0
"""float: The disruption offset of small bodies with strength or rubble (5 +/- 2; IF3). A rogue with a thick ice shell sits between the two."""
MU_BAR_RUBBLE = 0.37
"""float: The coupling exponent of small bodies (IF3)."""
SUPER_CATASTROPHIC_RATIO = 1.8
"""float: `Q_R / Q*_RD` at which the largest remnant falls to a tenth of the mass and the universal law gives way to a power law (IF1, IF2)."""
SUPER_CATASTROPHIC_REMNANT = 0.1
"""float: The largest-remnant fraction at `SUPER_CATASTROPHIC_RATIO`."""
SUPER_CATASTROPHIC_INDEX = -1.5
"""float: The index of the largest remnant against `Q_R / Q*_RD` beyond it (-1.2 to -1.5 in the laboratory; IF2)."""
FRAGMENT_TAIL_INDEX = 2.85
"""float: The cumulative index of the fragments below the largest remnant (the analytic model needs it under 3; IF5)."""
MERGER_REMNANT = 0.95
"""float: A non-grazing impact that leaves at least this share of the mass in the largest remnant is a merger. STAND-IN IF14 (tier T9, L1): the law has no sharp edge."""
FORSTERITE_DENSITY_KG_M3 = 3220.0
"""float: Zero-pressure density of forsterite, the rock of the shock estimate (IF13, tier T3, L0: recalled)."""
FORSTERITE_BULK_SOUND_SPEED_MS = 5500.0
"""float: The intercept of the linear shock-speed law `U_s = C + s u_p` of forsterite (IF13)."""
FORSTERITE_SLOPE = 1.3
"""float: The slope `s` of that law (IF13)."""

REGIMES = ("merger", "hit_and_run", "partial_erosion", "catastrophic", "super_catastrophic")
"""tuple: The outcomes in order of the energy that makes them (a grazing impact skips to hit-and-run below the disruption energy)."""

PC_PER_YEAR_AT_ONE_MS = constants.SECONDS_PER_YEAR / constants.PARSEC_M
"""float: Parsecs travelled in a year at 1 m/s."""


class RogueCollisionError(ValueError):
    """A masses, radii or speeds argument is not positive."""


def _positive(**values):
    for name, value in values.items():
        if not value > 0.0:
            raise RogueCollisionError(f"{name} {value!r} is not positive")


# --- rate ----------------------------------------------------------------------------------------------------------------

def escape_speed_touching_ms(mass1_kg, mass2_kg, radius1_m, radius2_m):
    """The escape speed of a pair at touching, `(2 G (M1 + M2) / (R1 + R2))^(1/2)`, m/s (11.2 km/s for two Earths)."""
    _positive(mass1_kg=mass1_kg, mass2_kg=mass2_kg, radius1_m=radius1_m, radius2_m=radius2_m)
    return math.sqrt(2.0 * constants.G * (mass1_kg + mass2_kg) / (radius1_m + radius2_m))


def cross_section_m2(mass1_kg, mass2_kg, radius1_m, radius2_m, speed_ms):
    """The collision cross-section with focusing, `pi (R1 + R2)^2 (1 + v_esc^2 / v^2)`, m2, for a relative speed at infinity of `speed_ms`."""
    _positive(speed_ms=speed_ms)
    v_esc = escape_speed_touching_ms(mass1_kg, mass2_kg, radius1_m, radius2_m)
    return math.pi * (radius1_m + radius2_m) ** 2 * (1.0 + (v_esc / speed_ms) ** 2)


def rate_per_rogue_per_year(number_density_pc3, mass1_kg, mass2_kg, radius1_m, radius2_m, speed_ms):
    """`Gamma = n sigma v`, collisions per rogue per year, for a number density of partners in rogues per pc3 (2.3e-23 to 3.5e-23 at 0.78 per pc3 for two Earths at 20 to 80 km/s)."""
    sigma_pc2 = cross_section_m2(mass1_kg, mass2_kg, radius1_m, radius2_m, speed_ms) / constants.PARSEC_M ** 2
    return number_density_pc3 * sigma_pc2 * speed_ms * PC_PER_YEAR_AT_ONE_MS


def galaxy_rate_per_year(rogues_per_star, star_pair_integral_pc3, mass1_kg, mass2_kg, radius1_m, radius2_m, speed_ms):
    """
    Collisions per year in one component of the Galaxy, `0.5 sigma v (rogues per star)^2 I`, with `I = integral(n_star^2 dV)` in stars squared per pc3
    (3.4e10 for the exponential disc and 5.1e10 for the bulge of the design). Two Earth-size pairs, the disc at 40 km/s and the bulge at 120 km/s: 8.9e-11 a year in all, about one in
    the life of the Galaxy.
    """
    _positive(rogues_per_star=rogues_per_star, star_pair_integral_pc3=star_pair_integral_pc3)
    sigma_pc2 = cross_section_m2(mass1_kg, mass2_kg, radius1_m, radius2_m, speed_ms) / constants.PARSEC_M ** 2
    return 0.5 * sigma_pc2 * speed_ms * PC_PER_YEAR_AT_ONE_MS * rogues_per_star ** 2 * star_pair_integral_pc3


def speed_fractions(dispersion_ms, escape_speed_ms, edges_ms):
    """
    The share of the collisions in each speed bin, for a Maxwellian relative speed of one-dimensional dispersion `dispersion_ms`, each speed weighted by `sigma v` with the focusing
    of a pair of escape speed `escape_speed_ms`: density `v^2 exp(-v^2 / 2 s^2) (v + v_esc^2 / v)`. `edges_ms` ascending give `len(edges) + 1` shares that add up to 1
    (below the first edge, between each pair, above the last). With `s` = 25 km/s and edges 15, 30 and 60 km/s: 3, 17, 60 and 20 percent; with 50 km/s: 0.2, 1.6, 15 and 83.
    """
    _positive(dispersion_ms=dispersion_ms)
    if escape_speed_ms < 0.0 or list(edges_ms) != sorted(edges_ms) or any(edge <= 0.0 for edge in edges_ms):
        raise RogueCollisionError("the edges must be positive and ascending and the escape speed not negative")
    two_s2 = 2.0 * dispersion_ms ** 2
    total = 2.0 * dispersion_ms ** 4 + escape_speed_ms ** 2 * dispersion_ms ** 2

    def below(edge):
        u = edge ** 2 / two_s2
        return (2.0 * dispersion_ms ** 4 * (1.0 - math.exp(-u) * (1.0 + u)) + escape_speed_ms ** 2 * dispersion_ms ** 2 * (1.0 - math.exp(-u))) / total

    cumulative = [0.0] + [below(edge) for edge in edges_ms] + [1.0]
    return tuple(upper - lower for lower, upper in zip(cumulative, cumulative[1:]))


# --- the impact ----------------------------------------------------------------------------------------------------------

def combined_radius_m(total_mass_kg):
    """`RC1`, the radius of the combined mass at `RHO1_KG_M3`, m (14 180 km for two Earths)."""
    _positive(total_mass_kg=total_mass_kg)
    return (3.0 * total_mass_kg / (4.0 * math.pi * RHO1_KG_M3)) ** (1.0 / 3.0)


def combined_escape_speed_ms(total_mass_kg):
    """`V_esc = (2 G Mtot / RC1)^(1/2)`, the unit of the fragment speeds, m/s (10.6 km/s for two Earths)."""
    return math.sqrt(2.0 * constants.G * total_mass_kg / combined_radius_m(total_mass_kg))


def specific_impact_energy_j_kg(mass1_kg, mass2_kg, speed_ms):
    """`Q_R = 0.5 mu V^2 / Mtot`, J/kg (equation 1)."""
    _positive(mass1_kg=mass1_kg, mass2_kg=mass2_kg, speed_ms=speed_ms)
    return 0.5 * mass1_kg * mass2_kg / (mass1_kg + mass2_kg) * speed_ms ** 2 / (mass1_kg + mass2_kg)


def disruption_speed_ms(target_mass_kg, projectile_mass_kg, c_star=C_STAR_PLANET, mu_bar=MU_BAR_PLANET):
    """
    The head-on impact speed at which the largest remnant is half the mass, `V*`, m/s (equations 22 and 30): `(32 pi c* / 5)^(1/2) (rho1 G)^(1/2) RC1` for equal masses, times
    `((gamma + 1)^2 / (4 gamma))^(1/(3 mu_bar))` for a projectile of mass ratio `gamma` (at most 1). 22.6 km/s for two Earths, 15.2 for two 0.3 Earths, 10.5 for two 0.1 Earths.
    """
    _positive(target_mass_kg=target_mass_kg, projectile_mass_kg=projectile_mass_kg)
    big, small = max(target_mass_kg, projectile_mass_kg), min(target_mass_kg, projectile_mass_kg)
    gamma = small / big
    equal = math.sqrt(32.0 * math.pi * c_star / 5.0) * math.sqrt(RHO1_KG_M3 * constants.G) * combined_radius_m(big + small)
    return ((gamma + 1.0) ** 2 / (4.0 * gamma)) ** (1.0 / (3.0 * mu_bar)) * equal


def disruption_energy_j_kg(target_mass_kg, projectile_mass_kg, c_star=C_STAR_PLANET, mu_bar=MU_BAR_PLANET):
    """`Q*_RD` for a head-on impact, J/kg: the energy of an impact at `disruption_speed_ms` (6.4e7 J/kg for two Earths)."""
    return specific_impact_energy_j_kg(target_mass_kg, projectile_mass_kg, disruption_speed_ms(target_mass_kg, projectile_mass_kg, c_star, mu_bar))


def critical_impact_parameter(target_radius_m, projectile_radius_m):
    """`b_crit = R / (R + r)`: an impact with a larger `b = sin theta` is grazing."""
    _positive(target_radius_m=target_radius_m, projectile_radius_m=projectile_radius_m)
    return target_radius_m / (target_radius_m + projectile_radius_m)


def impact_angle_cdf(impact_parameter):
    """The share of impacts with an impact parameter `b = sin theta` or less, `b^2` (the angle follows `sin theta cos theta d theta`); half have `b` under 0.707."""
    if not 0.0 <= impact_parameter <= 1.0:
        raise RogueCollisionError(f"the impact parameter {impact_parameter!r} is outside 0 to 1")
    return impact_parameter ** 2


def interacting_share(target_radius_m, projectile_radius_m, impact_parameter):
    """`alpha`, the share of the projectile's mass that overlaps the target (equations 7 to 11): 1 unless the impact is oblique enough to shear part of it off."""
    _positive(target_radius_m=target_radius_m, projectile_radius_m=projectile_radius_m)
    if not 0.0 <= impact_parameter <= 1.0:
        raise RogueCollisionError(f"the impact parameter {impact_parameter!r} is outside 0 to 1")
    r = projectile_radius_m
    length = (target_radius_m + r) * (1.0 - impact_parameter)
    if length >= 2.0 * r:
        return 1.0
    return max((3.0 * r * length ** 2 - length ** 3) / (4.0 * r ** 3), 0.0)


def oblique_disruption_energy_j_kg(target_mass_kg, projectile_mass_kg, target_radius_m, projectile_radius_m, impact_parameter, c_star=C_STAR_PLANET, mu_bar=MU_BAR_PLANET):
    """`Q'*_RD = (mu / mu_alpha)^(2 - 3 mu_bar / 2) Q*_RD`, J/kg (equation 15): the disruption energy of an impact at `impact_parameter`; at least `Q*_RD`."""
    alpha = interacting_share(target_radius_m, projectile_radius_m, impact_parameter)
    if alpha <= 0.0:
        return math.inf
    mu = target_mass_kg * projectile_mass_kg / (target_mass_kg + projectile_mass_kg)
    mu_alpha = alpha * projectile_mass_kg * target_mass_kg / (alpha * projectile_mass_kg + target_mass_kg)
    return (mu / mu_alpha) ** (2.0 - 1.5 * mu_bar) * disruption_energy_j_kg(target_mass_kg, projectile_mass_kg, c_star, mu_bar)


# --- what is left --------------------------------------------------------------------------------------------------------

def largest_remnant_fraction(energy_ratio):
    """
    `Mlr / Mtot` for `energy_ratio = Q_R / Q*_RD`: the universal law `0.5 - 0.5 (ratio - 1)` up to 1.8 (1 for no energy, 0.5 at the threshold, 0.1 at 1.8), the power law
    `0.1 (ratio / 1.8)^-1.5` beyond (equations 5 and 44).
    """
    if energy_ratio < 0.0:
        raise RogueCollisionError(f"the energy ratio {energy_ratio!r} is negative")
    if energy_ratio <= SUPER_CATASTROPHIC_RATIO:
        return min(1.0, 0.5 - 0.5 * (energy_ratio - 1.0))
    return SUPER_CATASTROPHIC_REMNANT * (energy_ratio / SUPER_CATASTROPHIC_RATIO) ** SUPER_CATASTROPHIC_INDEX


def second_remnant_fraction(largest_fraction, beta=FRAGMENT_TAIL_INDEX, largest_count=1, second_count=2):
    """
    `Mslr / Mtot` from equation 37: `(3 - beta)(1 - Nlr Mlr / Mtot) / (Nslr beta)` (0.013 of the mass at the disruption threshold). `largest_count` and `second_count` are 2 and 4 for the
    hit-and-run of equal masses.
    """
    if not 0.0 <= largest_fraction <= 1.0 or not 0.0 < beta < 3.0:
        raise RogueCollisionError("the remnant fraction must lie in 0 to 1 and beta in 0 to 3")
    return max((3.0 - beta) * (1.0 - largest_count * largest_fraction) / (second_count * beta), 0.0)


def rubble_speed_window_ms(target_mass_kg, projectile_mass_kg, c_star=C_STAR_PLANET, mu_bar=MU_BAR_PLANET):
    """
    `(V*, V* sqrt(1.8))`, m/s: the head-on speeds that leave a largest remnant of a half down to a tenth of the mass and rubble round it (22.6 to 30.3 km/s for two Earths). Slower,
    the pair merges or hits and runs; faster, the largest remnant is under a tenth.
    """
    v_star = disruption_speed_ms(target_mass_kg, projectile_mass_kg, c_star, mu_bar)
    return v_star, v_star * math.sqrt(SUPER_CATASTROPHIC_RATIO)


def regime(energy_ratio, grazing=False):
    """
    The outcome of an impact of `energy_ratio = Q_R / Q'*_RD` (`grazing`: `b` over `critical_impact_parameter`), one of `REGIMES`: a grazing impact under the disruption energy is
    `hit_and_run`; otherwise `merger` while the largest remnant keeps `MERGER_REMNANT` of the mass, `partial_erosion` to the threshold, `catastrophic` to 1.8 and
    `super_catastrophic` beyond.
    """
    if energy_ratio < 0.0:
        raise RogueCollisionError(f"the energy ratio {energy_ratio!r} is negative")
    if energy_ratio >= SUPER_CATASTROPHIC_RATIO:
        return "super_catastrophic"
    if energy_ratio >= 1.0:
        return "catastrophic"
    if grazing:
        return "hit_and_run"
    return "merger" if largest_remnant_fraction(energy_ratio) >= MERGER_REMNANT else "partial_erosion"


# --- how fast ------------------------------------------------------------------------------------------------------------

def _decay_rate(largest_fraction):
    """`lambda = S ln 10 = 10^A / (Mrem / Mtot)` of the fragment speed law `10^(A - S v)` per unit of escape speed (equations 38 to 43, read with dm/dv per unit speed: IF6)."""
    if not 0.0 <= largest_fraction < 1.0:
        raise RogueCollisionError("a remnant fraction of 1 or more leaves no fragments")
    return 10.0 ** (0.3 - 0.3 * largest_fraction) / (1.0 - largest_fraction)


def unbound_fraction(largest_fraction, threshold_escape_speeds=1.0):
    """
    The share of the total mass whose speed relative to the largest remnant is over `threshold_escape_speeds` of the escape speed of the combined mass: `(1 - Mlr / Mtot) exp(-lambda x)`.
    0.11 at the start of the super-catastrophic regime and 0.136 at most (IF6; GEN.1099 reading of equation 43).
    """
    return (1.0 - largest_fraction) * math.exp(-_decay_rate(largest_fraction) * threshold_escape_speeds)


def unbound_rms_speed_ms(largest_fraction, total_mass_kg):
    """
    The root mean square of the speed the unbound pieces keep at infinity, `V_esc sqrt(2 / lambda + 2 / lambda^2)`, m/s: the dispersion of the field they make (the
    one-dimensional dispersion is this over the square root of 3). About 13 km/s for two Earths with a tenth of the mass in the largest remnant.
    """
    lam = _decay_rate(largest_fraction)
    return combined_escape_speed_ms(total_mass_kg) * math.sqrt(2.0 / lam + 2.0 / lam ** 2)


def peak_shock_pressure_pa(speed_ms):
    """
    The peak shock pressure of a planar impact of two equal rocks, Pa: `rho0 (C + s u_p) u_p` with the particle speed `u_p` half the impact speed, forsterite's `C` and `s`
    (IF13). 70 GPa at 5 km/s, 600 GPa at 20 km/s, 4.3 TPa at 60 km/s.
    """
    _positive(speed_ms=speed_ms)
    particle = 0.5 * speed_ms
    return FORSTERITE_DENSITY_KG_M3 * (FORSTERITE_BULK_SOUND_SPEED_MS + FORSTERITE_SLOPE * particle) * particle
