"""
Stellar activity (GEN.86): each star's X-ray and XUV output, whether its
corona is still saturated, how often it flares, and the XUV it has put out
over its life. The design is docs/design/activity-magnetism-radiation-hydrosphere.md
sections 2.1 to 2.4; every rate here is a design default good to about
+/-0.5 dex.

1. **Coronal X-rays.** A young star's `L_X/L_bol` is saturated at
   `SATURATED_LX_LBOL` (7.4e-4 for G, K and M; falling through F to the
   1e-7 of A, B and O stars, whose X-rays come from wind shocks). After the
   saturation age `t_sat` it decays as `(t / t_sat)^-b`. Both `t_sat` and
   `b` depend on mass (`ACTIVITY_TABLE`): fully convective M dwarfs stay
   saturated for billions of years and then decline slowly. Each star draws
   a log-normal scatter on `t_sat` (`T_SAT_SIGMA_DEX`) and on its decayed
   level (`DECAY_SIGMA_DEX`).
2. **XUV.** The EUV (100 to 920 A) follows the X-rays,
   `log L_EUV = 4.80 + 0.860 log L_X` (erg/s; Sanz-Forcada 2011), and
   `L_XUV = L_X + L_EUV`. A hot photosphere adds its own blackbody output
   above 13.6 eV (`blackbody_fraction_above`), which dominates for O and B
   stars and hot white dwarfs.
3. **Flares.** `log10 N(>1e33 erg)` per year by class and age
   (`FLARE_TABLE`), scattered +/-0.5 dex; the frequency slope `alpha` is
   uniform in 1.8 to 2.2. A, F, B and O stars flare at 1e-4 per year;
   giants, white dwarfs and compact remnants not at all.
4. **Compact remnants.** A white dwarf's XUV is its photosphere's; a
   neutron star's is its thermal X-rays or a thousandth of its spin-down
   power (magnetic dipole), whichever is larger; a black hole's is its
   accretion disk's luminosity (or its Hawking glow without one).
5. **Fluence.** `xuv_fluence_j` integrates `L_XUV` over the star's age; a
   planet's exposure index divides it by its distance squared and by a
   1 solar mass star's at 1 AU over 5 Gyr (`REFERENCE_FLUENCE_J`).
"""

import math

from planetgen.physics import constants
from planetgen.util import draw

STAR_ACTIVITY_FIELDS = (
    "log_lx_lbol", "l_xuv_w", "xuv_saturated", "flare_n33_per_yr", "flare_alpha", "xuv_fluence_j",
    "lethal_event_rate_per_gyr",
)
"""tuple: What `generate_activity` sets on a star, as stored in `stars`."""

ERG_PER_J = 1e7

SN_RATE_AT_SUN_PER_GYR = 1.5
SN_SCALE_LENGTH_PC = 1450.0
SUN_GALACTIC_RADIUS_PC = 8000.0
SN_MAX_BOOST = 250.0

ACTIVITY_TABLE = (
    # mass (Msun), t_sat (Gyr), decay exponent b
    (0.1, 4.0, 1.0),
    (0.2, 3.0, 1.0),
    (0.3, 2.0, 1.0),
    (0.45, 1.0, 1.75),
    (0.6, 0.4, 2.0),
    (0.8, 0.15, 2.0),
    (1.0, 0.10, 2.0),
)
"""tuple: Section 2.1's saturation age and decay exponent by mass, linear
in between and held at the ends."""

SATURATED_LX_LBOL = ((1.0, 7.4e-4), (1.4, 5e-5), (2.0, 1e-7))
"""tuple: `(mass, saturated L_X/L_bol)`: 7.4e-4 up to 1 solar mass (Wright
2011, Jackson 2012), log-linear to 5e-5 at 1.4 and to the 1e-7 of hot stars'
wind shocks at 2, held beyond."""

T_SAT_SIGMA_DEX = 0.3
"""float: Per-star scatter of the saturation age (Tu et al. 2015: slow and
fast starters differ tenfold)."""

DECAY_SIGMA_DEX = 0.4
"""float: Per-star scatter of `L_X` after saturation (section 2.1)."""

GIANT_LX_LBOL = 1e-7
"""float: A giant's or supergiant's `L_X/L_bol` (past the coronal dividing
line, no saturated dynamo)."""

FLARE_TABLE = (
    # lowest mass of the class, then log10 N33 at <0.1, 0.1-0.6, 0.6-2, 2-6, >6 Gyr
    (0.85, (1.0, 0.0, -1.3, -2.5, -2.8)),    # G (0.85 to 1.1)
    (0.6, (1.3, 0.5, -0.8, -2.0, -2.3)),     # K
    (0.35, (2.0, 1.6, 0.8, -0.3, -0.7)),     # early M
    (0.15, (2.3, 2.2, 1.8, 0.9, 0.5)),       # mid M
    (0.0, (2.5, 2.5, 2.2, 0.8, 0.5)),        # late M
)
"""tuple: Section 2.4's `log10 N(>1e33 erg)` per year, heaviest class
first."""

FLARE_AGE_EDGES_GY = (0.1, 0.6, 2.0, 6.0)
FLARE_MAX_MASS = 1.1
"""float: Above this the star flares like an A or F star (`HOT_STAR_LOG_N33`)."""
HOT_STAR_LOG_N33 = -4.0
FLARE_SCATTER_DEX = 0.5
FLARE_ALPHA_RANGE = (1.8, 2.2)

XUV_EDGE_EV = 13.6
"""float: The XUV band starts at the hydrogen ionization edge (912 A)."""

REFERENCE_AGE_GY = 5.0
FLUENCE_STEPS = 64

DIPOLE_EDOT_TO_LX = 1e-3
"""float: A rotation-powered pulsar's X-ray efficiency (Becker and
Truemper 1997)."""


def _interp(table, mass):
    """Linear interpolation of `(t_sat, b)` in `ACTIVITY_TABLE`, held at
    the ends."""
    if mass <= table[0][0]:
        return table[0][1:]
    for (m0, *v0), (m1, *v1) in zip(table, table[1:]):
        if mass <= m1:
            f = (mass - m0) / (m1 - m0)
            return tuple(a + f * (b - a) for a, b in zip(v0, v1))
    return table[-1][1:]


def saturated_lx_lbol(mass_sol):
    """Section 2.1's saturated `L_X/L_bol` at `mass_sol`."""
    if mass_sol <= SATURATED_LX_LBOL[0][0]:
        return SATURATED_LX_LBOL[0][1]
    for (m0, r0), (m1, r1) in zip(SATURATED_LX_LBOL, SATURATED_LX_LBOL[1:]):
        if mass_sol <= m1:
            f = (mass_sol - m0) / (m1 - m0)
            return 10 ** (math.log10(r0) + f * (math.log10(r1) - math.log10(r0)))
    return SATURATED_LX_LBOL[-1][1]


def lx_lbol_at(age_gy, saturated, t_sat_gy, decay, offset_dex=0.0):
    """`L_X/L_bol` at `age_gy`: `saturated` until `t_sat_gy`, then decaying
    as `(age / t_sat)^-decay`, times `10^offset_dex`."""
    if age_gy <= t_sat_gy:
        return saturated
    return saturated * (age_gy / t_sat_gy) ** -decay * 10 ** offset_dex


def coronal_xuv_w(l_x_w):
    """`L_X + L_EUV`, W, with `log L_EUV = 4.80 + 0.860 log L_X` (erg/s)."""
    if l_x_w <= 0.0:
        return 0.0
    l_euv_erg = 10 ** (4.80 + 0.860 * math.log10(l_x_w * ERG_PER_J))
    return l_x_w + l_euv_erg / ERG_PER_J


def blackbody_fraction_above(energy_ev, temperature_k):
    """The share of a blackbody's output above `energy_ev` (the series of
    the Planck integral)."""
    if not temperature_k or temperature_k <= 0.0:
        return 0.0
    x = energy_ev * constants.ELECTRON_VOLT_J / (constants.BOLTZMANN * temperature_k)
    if x > 700.0:
        return 0.0
    total = 0.0
    for n in range(1, 60):
        term = math.exp(-n * x) * (x ** 3 / n + 3 * x ** 2 / n ** 2 + 6 * x / n ** 3 + 6 / n ** 4)
        total += term
        if term < 1e-12 * total:
            break
    return min(1.0, total * 15.0 / math.pi ** 4)


def flare_log_n33(mass_sol, age_gy):
    """Section 2.4's table value of `log10 N(>1e33 erg)` per year."""
    if mass_sol > FLARE_MAX_MASS:
        return HOT_STAR_LOG_N33
    column = sum(1 for edge in FLARE_AGE_EDGES_GY if age_gy >= edge)
    for lowest, row in FLARE_TABLE:
        if mass_sol >= lowest:
            return row[column]
    return FLARE_TABLE[-1][1][column]


def _fluence_j(luminosity_w, temperature_k, mass_sol, age_gy, t_sat_gy, decay, offset_dex):
    """`L_XUV` integrated from 0 to `age_gy`, J: the saturated stretch in one
    piece, the decay on log-spaced steps (midpoint rule)."""
    photosphere = luminosity_w * blackbody_fraction_above(XUV_EDGE_EV, temperature_k)
    saturated = saturated_lx_lbol(mass_sol)
    gy_s = 1e9 * constants.SECONDS_PER_YEAR
    sat_end = min(age_gy, t_sat_gy)
    total = (coronal_xuv_w(saturated * luminosity_w) + photosphere) * sat_end * gy_s
    if age_gy > t_sat_gy:
        ratio = (age_gy / t_sat_gy) ** (1.0 / FLUENCE_STEPS)
        t0 = t_sat_gy
        for _ in range(FLUENCE_STEPS):
            t1 = t0 * ratio
            mid = math.sqrt(t0 * t1)
            l_x = lx_lbol_at(mid, saturated, t_sat_gy, decay, offset_dex) * luminosity_w
            total += (coronal_xuv_w(l_x) + photosphere) * (t1 - t0) * gy_s
            t0 = t1
    return total


def reference_fluence_j():
    """A 1 solar mass, solar-luminosity star's XUV over `REFERENCE_AGE_GY`,
    with no scatter: the exposure index's unit (section 2.3)."""
    t_sat, decay = _interp(ACTIVITY_TABLE, 1.0)
    return _fluence_j(constants.SOLAR_LUMINOSITY, 5772.0, 1.0, REFERENCE_AGE_GY, t_sat, decay, 0.0)


REFERENCE_FLUENCE_J = reference_fluence_j()


def lethal_event_rate_per_gyr(galactic_distance_ly):
    """GEN.87 rule 7: supernovae within 8 pc per Gyr at this distance from
    the galactic centre (1.5 at the Sun's 8 kpc, up to 250 times that),
    or `None` for a star with no galaxy placement."""
    if galactic_distance_ly is None:
        return None
    radius_pc = galactic_distance_ly * constants.LY_TO_M / constants.PARSEC_M
    boost = min(SN_MAX_BOOST, math.exp((SUN_GALACTIC_RADIUS_PC - radius_pc) / SN_SCALE_LENGTH_PC))
    return SN_RATE_AT_SUN_PER_GYR * boost


def set_galactic_hazard(star):
    """Sets `lethal_event_rate_per_gyr` from the star's galaxy placement."""
    star.lethal_event_rate_per_gyr = lethal_event_rate_per_gyr(getattr(star, "galactic_center_dist_ly", None))


def generate_activity(star):
    """Sets `STAR_ACTIVITY_FIELDS` on an ordinary star (main sequence,
    subgiant, subdwarf, giant or white dwarf)."""
    set_galactic_hazard(star)
    mass_sol = star.mass / constants.SOLAR_MASS_TO_KG
    luminosity = star.luminosity or 0.0
    age = max(star.age or 0.0, 1e-4)
    photosphere = luminosity * blackbody_fraction_above(XUV_EDGE_EV, star.temperature)
    yerkes = star.yerkes_class or ""
    if yerkes == "VII":
        _set_remnant(star, photosphere, age)
        return
    t_sat, decay = _interp(ACTIVITY_TABLE, mass_sol)
    t_sat *= 10 ** draw.gauss(0.0, T_SAT_SIGMA_DEX)
    offset = draw.gauss(0.0, DECAY_SIGMA_DEX)
    log_n33 = flare_log_n33(mass_sol, age) + draw.uniform(-FLARE_SCATTER_DEX, FLARE_SCATTER_DEX)
    star.flare_alpha = draw.uniform(*FLARE_ALPHA_RANGE)
    if yerkes in ("IV", "V", "VI"):
        ratio = lx_lbol_at(age, saturated_lx_lbol(mass_sol), t_sat, decay, offset)
        star.xuv_saturated = age <= t_sat
        star.flare_n33_per_yr = 10 ** log_n33
        star.xuv_fluence_j = _fluence_j(luminosity, star.temperature, mass_sol, age, t_sat, decay, offset)
    else:
        # A giant has left the saturated dynamo behind; its early fluence
        # is the main-sequence one, which this does not track.
        ratio = GIANT_LX_LBOL
        star.xuv_saturated = False
        star.flare_n33_per_yr = 0.0
        star.xuv_fluence_j = (coronal_xuv_w(ratio * luminosity) + photosphere) * age * 1e9 * constants.SECONDS_PER_YEAR
    star.log_lx_lbol = math.log10(ratio)
    star.l_xuv_w = coronal_xuv_w(ratio * luminosity) + photosphere


def _set_remnant(star, l_xuv_w, age_gy):
    """A remnant's activity: no corona or flares, `l_xuv_w` held over its
    age (its cooling is not integrated)."""
    star.log_lx_lbol = None
    star.xuv_saturated = False
    star.flare_n33_per_yr = 0.0
    star.flare_alpha = None
    star.l_xuv_w = l_xuv_w
    star.xuv_fluence_j = l_xuv_w * age_gy * 1e9 * constants.SECONDS_PER_YEAR


def neutron_star_xuv_w(neutron_star):
    """A neutron star's XUV, W: its thermal surface output or
    `DIPOLE_EDOT_TO_LX` of its magnetic-dipole spin-down power
    (`B^2 R^6 Omega^4 / 6 c^3`, cgs), whichever is larger."""
    radius_cm = neutron_star.radius * 1e5
    omega = 2 * math.pi / (neutron_star.spin_period_ms / 1000.0)
    c_cm = constants.SPEED_OF_LIGHT_M_S * 100.0
    edot_erg = neutron_star.magnetic_field_gauss ** 2 * radius_cm ** 6 * omega ** 4 / (6 * c_cm ** 3)
    return max(neutron_star.luminosity, DIPOLE_EDOT_TO_LX * edot_erg / ERG_PER_J)


def generate_remnant_activity(remnant, l_xuv_w):
    """`STAR_ACTIVITY_FIELDS` for a neutron star or black hole (no draws)."""
    set_galactic_hazard(remnant)
    _set_remnant(remnant, l_xuv_w, max(remnant.age or 0.0, 1e-4))


def xuv_flux_earth(star, distance_au):
    """XUV flux at `distance_au` in units of Earth's today
    (`constants.EARTH_XUV_FLUX_W_M2`)."""
    l_xuv = getattr(star, "l_xuv_w", None)
    if l_xuv is None or not distance_au:
        return None
    d_m = distance_au * constants.AU_TO_M
    return l_xuv / (4 * math.pi * d_m ** 2) / constants.EARTH_XUV_FLUX_W_M2


def exposure_index(star, distance_au):
    """Cumulative XUV at `distance_au` against a 1 solar mass star's at
    1 AU over 5 Gyr (section 2.3: 2 or less Earth-like, 2 to 8 elevated, 8
    to 30 high, over 30 extreme)."""
    fluence = getattr(star, "xuv_fluence_j", None)
    if fluence is None or not distance_au:
        return None
    return fluence / distance_au ** 2 / REFERENCE_FLUENCE_J


def flare_irradiation_index(star, distance_au):
    """`N33 / d_AU^2` (section 2.3: Earth 0.003, Proxima b about 2,100)."""
    rate = getattr(star, "flare_n33_per_yr", None)
    if rate is None or not distance_au:
        return None
    return rate / distance_au ** 2
