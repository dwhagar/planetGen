# planetgen/physics/star_pair_tides.py

"""
Mutual tides, spin and apsidal precession between the two stars of a pair (GEN.1467)
====================================================================================

The stars of a binary raise tides on each other. Three things follow, none of which a star drawn alone can have:

* **Spin-orbit locking.** The equilibrium tide drags a star's spin toward the orbit's. A star's spin approaches the pseudo-synchronous rate `Omega_ps` (Hut 1981) with a
  synchronisation time `t_sync`: `Omega(t) = Omega_ps + (Omega_0 - Omega_ps) exp(-t / t_sync)`. The time is Zahn's (1977, 1989) for a star with a convective envelope, in the
  form of Rasio et al. (1996) and Hurley et al. (2002) equation 21 on, and Zahn's (1975) dynamical tide for a star with a radiative envelope (Hurley et al. 2002 equations 41 and
  42). A giant has a convective envelope; a white dwarf, a neutron star and a black hole keep the spin they were drawn with.
* **Spin-orbit alignment.** The same torque turns the spin axis onto the orbit normal; the angle between them falls as `exp(-t / t_align)` with `t_align = ALIGN_FACTOR t_sync`.
* **Apsidal precession.** The orbit's periapsis turns because of the stars' tidal bulges, their rotational bulges and general relativity (Sterne 1939; Hilditch 2001): with `n`
  the mean motion, `omega_dot = n (c_rot + c_tide) + omega_dot_GR`, `c_rot = k2 (Omega / n)^2 (1 + m_other / m_this) (R / a)^5 g(e)` and `c_tide = 15 k2 (m_other / m_this)
  (R / a)^5 f(e)` for each star, `g = (1 - e^2)^-2` and `f = (1 + 3 e^2 / 2 + e^4 / 8) (1 - e^2)^-5`. Relativity is `stability.gr_precession_period_yr`.

Everything here is pure: it reads a star's mass, radius, temperature, class and spin and the pair's axis and eccentricity, and returns numbers. `generation.system` applies them
to the stars of a pair (a hierarchy is taken a pair node at a time, each star against the mass of the other side) and notes the precession in the effect ledger (effect 6).
The orbit is held as it is today (the period and eccentricity cuts of GEN.1011), so no stored orbit moves; the circular cut below 7 days stays (`multistar.CIRCULAR_BELOW_DAYS`).

Stand-ins are KT1 to KT8 of `docs/design/research/gap-filling-plan.md` section 6.273: the radius of gyration (KT1), the convective envelope share by mass (KT2), the convective and
radiative choice (KT3), the first-order relaxation over the whole age with the orbit held and the breakup cap (KT4), the white dwarf, neutron star and black hole exclusions (KT5), the
alignment time (KT6), the apsidal constants by kind (KT7), and the pair-by-pair reading of a hierarchy with no third-body term (KT8).
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, magnetic_geometry, spin, stability

GYRATION_SQUARED = 0.1
"""float: `I / (M R^2)`, the star's moment of inertia over `M R^2` (Hurley et al. 2002 take 0.1 for every star; a stand-in, row KT1, tier Confident Unverified, tweak L0)."""

ALIGN_FACTOR = 1.0
"""float: The alignment time over the synchronisation time (Hut 1981 gives the same order for the two; a stand-in, row KT6, tier Marginal educated guess from physics, tweak L0)."""

ZAHN_RADIATIVE_E2_COEFFICIENT = 1.592e-9
ZAHN_RADIATIVE_E2_EXPONENT = 2.84
"""float: The tidal coefficient of a radiative envelope, `E2 = 1.592e-9 (M / Msun)^2.84` (Zahn 1975 as fitted by Hurley et al. 2002 equation 42)."""

ENVELOPE_MASS_SHARE = ((0.35, 1.0), (0.6, 0.25), (1.0, 0.03), (1.3, 0.003))
"""tuple: `(mass in solar masses, convective envelope share of the mass)` points of a dwarf, interpolated in the log of the share: fully convective at 0.35 and below, the Sun's 2 to 3
percent, none worth counting at the Kraft break (a stand-in, row KT2, tier Marginal educated guess from physics, tweak L1)."""

GIANT_ENVELOPE_MASS_SHARE = 0.5
"""float: A giant's convective envelope share of its mass (a stand-in, row KT2)."""

RADIATIVE_MASS_SOL = 1.3
"""float: A dwarf above this mass, or hotter than the Kraft break's middle, has a radiative envelope (a stand-in, row KT3, tier Confident educated guess from physics, tweak L0)."""

APSIDAL_CONSTANT_FULLY_CONVECTIVE = 0.14
APSIDAL_CONSTANT_SOLAR = 0.03
APSIDAL_CONSTANT_MASSIVE = 0.004
APSIDAL_CONSTANT_GIANT = 0.005
"""float: The apsidal motion constant `k2` of a star: 0.14 for a fully convective one (a polytrope of index 1.5), 0.03 at a solar mass, 0.004 at 5 solar masses and over, 0.005 for a
giant (centrally condensed); the log is interpolated in the mass between (observed `log k2` of main-sequence stars runs from about -1.5 at a solar mass to about -2.5 at 10;
Claret and Gimenez 2010; stand-ins, row KT7, tier Marginal unverified, tweak L2)."""

KT_MASS_POINTS = ((0.35, APSIDAL_CONSTANT_FULLY_CONVECTIVE), (1.0, APSIDAL_CONSTANT_SOLAR), (5.0, APSIDAL_CONSTANT_MASSIVE))

COMPACT_CLASSES_PREFIXES = ("D", "N", "X")
"""tuple: Yerkes classes that begin a white dwarf (D), a neutron star (N) or a black hole (X): they keep the spin they were drawn with."""

NEVER = math.inf


def _mass_sol(mass_kg):
    return mass_kg / constants.SOLAR_MASS_TO_KG


def _yr(seconds):
    return seconds / constants.SECONDS_PER_YEAR


def _interpolate_log(points, x):
    """`y` at `x` from `(x, y)` points, linear in `x` and in `log y`, held at the ends."""
    if x <= points[0][0]:
        return points[0][1]
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x <= x1:
            return math.exp(math.log(y0) + (x - x0) / (x1 - x0) * (math.log(y1) - math.log(y0)))
    return points[-1][1]


def is_compact(yerkes_class):
    """Whether a star is a white dwarf, a neutron star or a black hole by its class string."""
    return bool(yerkes_class) and yerkes_class.strip().startswith(COMPACT_CLASSES_PREFIXES)


def has_convective_envelope(star):
    """Whether the star's tide is the convective one: a giant, or a dwarf of the mass and temperature below the radiative limit."""
    cls = (star.yerkes_class or "").strip()
    if cls in spin.GIANT_CLASSES:
        return True
    return _mass_sol(star.mass) <= RADIATIVE_MASS_SOL and spin.kraft_weight(star.temperature) < 0.5


def envelope_mass_share(star):
    """The convective envelope's share of the star's mass."""
    if (star.yerkes_class or "").strip() in spin.GIANT_CLASSES:
        return GIANT_ENVELOPE_MASS_SHARE
    return _interpolate_log(ENVELOPE_MASS_SHARE, _mass_sol(star.mass))


def apsidal_constant(star):
    """The star's apsidal motion constant `k2` (see `APSIDAL_CONSTANT_*`)."""
    if (star.yerkes_class or "").strip() in spin.GIANT_CLASSES:
        return APSIDAL_CONSTANT_GIANT
    mass = _mass_sol(star.mass)
    if mass >= KT_MASS_POINTS[-1][0]:
        return APSIDAL_CONSTANT_MASSIVE
    return _interpolate_log(KT_MASS_POINTS, mass)


def pseudo_synchronous_ratio(eccentricity):
    """`Omega_ps / n`: the spin rate an eccentric orbit's tide settles a star at, over the mean motion (Hut 1981 equation 42)."""
    e2 = eccentricity * eccentricity
    return (1.0 + 7.5 * e2 + 45.0 / 8.0 * e2 * e2 + 5.0 / 16.0 * e2 ** 3) / ((1.0 + 3.0 * e2 + 3.0 / 8.0 * e2 * e2) * (1.0 - e2) ** 1.5)


def sync_time_yr(star, partner_mass_kg, axis_au, spin_period_hours=None):
    """
    The time (years) in which the tide of a partner of mass `partner_mass_kg` at a mean separation `axis_au` takes a star's spin a factor `e` of the way to synchronous.

    A white dwarf, a neutron star and a black hole are `NEVER` (inf). The convective form is `1 / t = 3 (k / T) q^2 (R / a)^6 / GYRATION_SQUARED` with `k / T = (2 / 21) (f_conv /
    tau_conv) (M_env / M)`, `tau_conv` the convective turnover time of `magnetic_geometry.convective_turnover_days`, `f_conv = min(1, (P_tid / 2 tau_conv)^2)` and `P_tid` the
    period of the tide in the star's frame; the radiative form is `1 / t = 5 2^(5/3) sqrt(G M / R^3) q^2 (1 + q)^(5/6) E2 (R / a)^(17/2) / GYRATION_SQUARED`.
    """
    if is_compact(star.yerkes_class) or not star.mass or not star.radius or not axis_au or not partner_mass_kg:
        return NEVER
    q = partner_mass_kg / star.mass
    radius_m = star.radius * 1000.0
    a_m = axis_au * constants.AU_TO_M
    ratio = radius_m / a_m
    if ratio >= 1.0:
        return NEVER
    if has_convective_envelope(star):
        tau_yr = magnetic_geometry.convective_turnover_days(_mass_sol(star.mass)) / 365.25
        orbital_period_yr = stability.orbital_period_yr(axis_au, _mass_sol(star.mass + partner_mass_kg))
        spin_hours = spin_period_hours or star.rotation_period_hours
        spin_rate = 1.0 / (spin_hours / 24.0 / 365.25) if spin_hours else 0.0
        tide_rate = abs(1.0 / orbital_period_yr - spin_rate)
        tide_period_yr = 1.0 / tide_rate if tide_rate > 0.0 else NEVER
        f_conv = min(1.0, (tide_period_yr / (2.0 * tau_yr)) ** 2)
        k_over_t = (2.0 / 21.0) * (f_conv / tau_yr) * envelope_mass_share(star)       # per year
        rate_per_yr = 3.0 * k_over_t * q * q * ratio ** 6 / GYRATION_SQUARED
    else:
        e2 = ZAHN_RADIATIVE_E2_COEFFICIENT * _mass_sol(star.mass) ** ZAHN_RADIATIVE_E2_EXPONENT
        dynamical = math.sqrt(constants.G * star.mass / radius_m ** 3) * constants.SECONDS_PER_YEAR      # per year
        rate_per_yr = 5.0 * 2.0 ** (5.0 / 3.0) * dynamical * q * q * (1.0 + q) ** (5.0 / 6.0) * e2 * ratio ** 8.5 / GYRATION_SQUARED
    return 1.0 / rate_per_yr if rate_per_yr > 0.0 else NEVER


@dataclass(frozen=True)
class SpinResult:
    """What a pair's tide leaves of a star's spin."""
    sync_time_yr: float            # inf when the star cannot be tidally spun
    rotation_period_hours: float   # the spin period at the star's age (the drawn one when untouched)
    obliquity_factor: float        # the share of the drawn angle between spin and orbit normal that is left
    locked: bool                   # the spin is within 1 percent of the pseudo-synchronous rate
    changed: bool


def spin_state(star, partner_mass_kg, axis_au, eccentricity, age_yr):
    """
    The star's spin at `age_yr` under the tide of its partner: the relaxation toward `Omega_ps` over `sync_time_yr` (the orbit held as it is today), never faster than
    `spin.BREAKUP_FRACTION` of breakup. A star the tide does not reach (`sync_time_yr` over 100 ages) comes back unchanged.
    """
    drawn_hours = star.rotation_period_hours
    t_sync = sync_time_yr(star, partner_mass_kg, axis_au, drawn_hours)
    if not math.isfinite(t_sync) or t_sync > 100.0 * age_yr or not drawn_hours:
        return SpinResult(t_sync, drawn_hours, 1.0, False, False)
    orbital_period_hours = stability.orbital_period_yr(axis_au, _mass_sol(star.mass + partner_mass_kg)) * 365.25 * 24.0
    omega_orbit = 2.0 * math.pi / orbital_period_hours
    omega_ps = pseudo_synchronous_ratio(eccentricity) * omega_orbit
    omega_0 = 2.0 * math.pi / drawn_hours
    decay = math.exp(-age_yr / t_sync)
    omega = omega_ps + (omega_0 - omega_ps) * decay
    hours = 2.0 * math.pi / omega
    cap_hours = spin.period_hours_from_speed(star.radius, spin.BREAKUP_FRACTION * spin.breakup_speed_kms(star.mass, star.radius))
    hours = max(hours, cap_hours)
    obliquity = math.exp(-age_yr / (ALIGN_FACTOR * t_sync))
    locked = abs(omega - omega_ps) <= 0.01 * omega_ps
    return SpinResult(t_sync, hours, obliquity, locked, True)


def aligned_axis(axis, normal, factor):
    """
    The unit vector `axis` turned toward the unit vector `normal`, in the plane of the two, until the angle between them is `factor` times what it was (`factor` 1 leaves it, 0 puts it
    on the normal): the spin axis of a star after the pair's tide has taken part of its misalignment away (`SpinResult.obliquity_factor`).
    """
    dot = max(-1.0, min(1.0, sum(a * n for a, n in zip(axis, normal))))
    angle = math.acos(dot)
    if angle < 1e-12 or factor >= 1.0:
        return tuple(axis)
    across = tuple(a - dot * n for a, n in zip(axis, normal))
    length = math.sqrt(sum(c * c for c in across))
    if length < 1e-12:
        return tuple(axis)
    new_angle = angle * max(factor, 0.0)
    return tuple(math.cos(new_angle) * n + math.sin(new_angle) * c / length for n, c in zip(normal, across))


@dataclass(frozen=True)
class Precession:
    """The apsidal precession rates (radians a year) of a pair's orbit, by cause, and their sum."""
    relativity: float
    tide: float
    rotation: float

    @property
    def total(self):
        return self.relativity + self.tide + self.rotation

    @property
    def period_yr(self):
        return 2.0 * math.pi / self.total if self.total > 0.0 else NEVER


def apsidal_precession(star_a, star_b, axis_au, eccentricity):
    """
    The apsidal precession of the orbit of two stars about each other, `Precession` rates in radians a year: the tidal and rotational bulges of both (`k2` of each, its spin
    `rotation_period_hours`) and general relativity. A compact star has no bulge worth counting (its `R / a` is tiny) but still carries its mass for relativity.
    """
    total_mass = star_a.mass + star_b.mass
    period_yr = stability.orbital_period_yr(axis_au, _mass_sol(total_mass))
    mean_motion = 2.0 * math.pi / period_yr
    e2 = eccentricity * eccentricity
    g_e = (1.0 - e2) ** -2
    f_e = (1.0 + 1.5 * e2 + e2 * e2 / 8.0) * (1.0 - e2) ** -5
    a_km = axis_au * constants.AU_TO_M / 1000.0
    tide = rotation = 0.0
    for this, other in ((star_a, star_b), (star_b, star_a)):
        if not this.radius or is_compact(this.yerkes_class):
            continue
        ratio5 = (this.radius / a_km) ** 5
        k2 = apsidal_constant(this)
        mass_ratio = other.mass / this.mass
        spin_hours = this.rotation_period_hours
        omega = 2.0 * math.pi / (spin_hours / 24.0 / 365.25) if spin_hours else 0.0
        rotation += k2 * (omega / mean_motion) ** 2 * (1.0 + mass_ratio) * ratio5 * g_e * mean_motion
        tide += 15.0 * k2 * mass_ratio * ratio5 * f_e * mean_motion
    gr_period = stability.gr_precession_period_yr(axis_au, eccentricity, _mass_sol(total_mass), period_yr)
    return Precession(2.0 * math.pi / gr_period, tide, rotation)
