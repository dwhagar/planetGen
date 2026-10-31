"""
Stellar activity (GEN.86): each star's X-ray and XUV output, whether its
corona is still saturated, how often it flares, and the XUV it has put out
over its life. The design is docs/design/habitability/activity-magnetism-radiation-hydrosphere.md
sections 2.1 to 2.4; every rate here is a design default good to about
+/-0.5 dex.

1. **Coronal X-rays.** A young star's `L_X/L_bol` is saturated at
   `SATURATED_LX_LBOL` (7.4e-4 for G, K and M; falling through F to the
   1e-7 of A, B and O stars, whose X-rays come from wind shocks). After the
   saturation age `t_sat` it decays as `(t / t_sat)^-b`. `t_sat` comes from
   `saturation_time` (GEN.493, Johnstone et al. 2021) by the star's mass and
   its own draw of the initial-rotation percentile: fully convective M dwarfs
   stay saturated for billions of years. `b` depends on mass
   (`DECAY_EXPONENTS`), and each star draws a log-normal scatter on its
   decayed level (`DECAY_SIGMA_DEX`).
2. **XUV.** The EUV (100 to 920 A) follows the X-rays,
   `log L_EUV = 4.80 + 0.860 log L_X` (erg/s), and
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
5. **Rotation (GEN.807).** The activity-rotation relation (Wright et al. 2011): X-rays saturate while the Rossby number `Ro = P_rot / tau_conv` is
   under 0.13 and fall as `Ro^-2.7` after, so a star that spins faster than a star of its mass and age spins stays saturated longer by the factor
   `(P_nominal / P)^(2.7 / b)` of its saturation age (`rotation_saturation_scale`), and flares correspondingly more often. The saturation
   age of 1 stays the mean (it is the observed one of M dwarfs, whose spin-down the gyrochronology law does not follow); a star's own
   rotation (`spin.star_rotation_period_hours`, drawn before) moves it about that mean. A star the law does not cover (a hot star, a giant, a
   white dwarf) keeps its saturation age.
6. **Fluence.** `xuv_fluence_j` integrates `L_XUV` over the star's age; a
   planet's exposure index divides it by its distance squared and by a
   1 solar mass star's at 1 AU over 5 Gyr (`REFERENCE_FLUENCE_J`).
"""

import math

from planetgen.galaxy import scales
from planetgen.physics import constants, saturation_time, spin, stellar_library, stellar_spectrum, unit_registry as units
from planetgen.util import draw

_kg_to_msun = units.converter("kg", "Msun")

_au_to_m = units.converter("au", "m")
_ly_to_m = units.converter("ly", "m")
_m_to_pc = units.converter("m", "pc")
_ms_to_s = units.converter("ms", "s")
_gyr_to_yr = units.converter("Gyr", "yr")
_yr_to_s = units.converter("yr", "s")
_gyr_to_s = units.converter("Gyr", "s")

STAR_ACTIVITY_UNITS = {
    "log_lx_lbol": "dex", "l_xuv_w": "W", "xuv_saturated": "flag", "flare_n33_per_yr": "1/yr", "flare_alpha": "ratio",
    "xuv_fluence_j": "J", "lethal_event_rate_per_gyr": "1/Gy", "uv_c_share": "ratio",
}
"""dict: What `generate_activity` sets on a star, as stored in `stars`, with each attribute's native unit."""

STAR_ACTIVITY_FIELDS = tuple(STAR_ACTIVITY_UNITS)
"""tuple: The names of `STAR_ACTIVITY_UNITS`."""

ERG_PER_J = 1e7

SN_RATE_AT_SUN_PER_GYR = 1.5
SN_SCALE_LENGTH_PC = 1450.0
SUN_GALACTIC_RADIUS_PC = 8000.0
SN_MAX_BOOST = 250.0

# Gamma-ray bursts (GEN.189): a beamed long burst within about 2 kpc whose
# fluence at the planet is 30 to 100 kJ/m2 strips about a third of the ozone
# (Thomas et al. 2005, ApJ 634:509), about one per 0.5 Gyr at the Sun
# (Melott et al. 2004, Int. J. Astrobiol. 3:55). Long bursts need low-metal
# massive stars and follow star formation outward, so the rate falls
# off beyond the Sun with the supernova scale length but is not raised
# inside it, where metal-rich stars rarely make one.
GRB_LETHAL_RATE_AT_SUN_PER_GYR = 2.0

# Active galactic nucleus (GEN.189): there is no catalogued rate, so this is a
# design estimate, not a measurement. Sgr A* is quiet now, but the Fermi
# bubbles and the X-ray echoes of the molecular clouds (Su et al. 2010;
# Ponti et al. 2013) show Seyfert-like episodes in the last few Myr. Take
# 1e44 erg/s for 1 Myr every 100 Myr: the fluence 1e44 * 3.15e13 s / (4 pi r^2)
# is 30 kJ/m2 at 2.9 kpc, so a planet that close to the centre loses its
# ozone once per episode.
AGN_EPISODES_PER_GYR = 10.0
AGN_LETHAL_RADIUS_PC = 2900.0

DECAY_EXPONENTS = (
    # mass (Msun), decay exponent b
    (0.1, 1.0),
    (0.2, 1.0),
    (0.3, 1.0),
    (0.45, 1.35),
    (0.6, 1.8),
    (0.8, 1.45),
    (1.0, 1.45),
)
"""tuple: Section 2.1's decay exponent after saturation by mass, linear in
between and held at the ends. Fully convective stars (0.3 and below) decline
slowly (Wright 2018). From 0.45 up the exponent is the one that holds each
mass's `L_X/L_bol` at `REFERENCE_AGE_GY` where the earlier table put it, now
that `t_sat` is the median of Johnstone et al. 2021 and not the table's own:
1.45 at 1 solar mass is the Skumanich band (1.35 to 1.5) and gives the Sun
1.3e27 erg/s today. [Stand-in AB1, tier T9, tweak L3.]"""

PERCENTILE_RANGE = (2.0, 98.0)
"""tuple: The initial-rotation percentiles a star draws uniformly between: the range
the tracks of Johnstone et al. 2021 offer, uniform because the distribution is the
observed one (row p6)."""

MEDIAN_PERCENTILE = 50.0

SATURATED_LX_LBOL = ((1.0, 7.4e-4), (1.4, 5e-5), (2.0, 1e-7))
"""tuple: `(mass, saturated L_X/L_bol)`: 7.4e-4 up to 1 solar mass (Wright
2011, Jackson 2012), log-linear to 5e-5 at 1.4 and to the 1e-7 of hot stars'
wind shocks at 2, held beyond."""

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

ROTATION_SLOPE = 2.7
"""float: The exponent of the Rossby number in the unsaturated `L_X/L_bol` (Wright et al. 2011: -2.70 +/- 0.13), which sets how far a star's
own spin moves its saturation age (GEN.807)."""

XUV_EDGE_EV = 13.6
"""float: The XUV band starts at the hydrogen ionization edge (912 A)."""

REFERENCE_AGE_GY = 5.0
FLUENCE_STEPS = 64

DIPOLE_EDOT_TO_LX = 1e-3
"""float: A rotation-powered pulsar's X-ray efficiency."""


def decay_exponent(mass_sol):
    """Linear interpolation of the decay exponent in `DECAY_EXPONENTS`, held at the ends."""
    table = DECAY_EXPONENTS
    if mass_sol <= table[0][0]:
        return table[0][1]
    for (m0, b0), (m1, b1) in zip(table, table[1:]):
        if mass_sol <= m1:
            return b0 + (mass_sol - m0) / (m1 - m0) * (b1 - b0)
    return table[-1][1]


def saturation_age_gy(mass_sol, percentile=MEDIAN_PERCENTILE):
    """The age, Gyr, a star of `mass_sol` and initial-rotation `percentile` leaves X-ray saturation (GEN.493; the median unless a percentile is given)."""
    return saturation_time.saturation_time_myr(mass_sol, percentile) / 1000.0


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


def coronal_xuv_at(mass_sol, age_gy, luminosity_w):
    """`L_XUV` of the corona of a `mass_sol` star at `age_gy`, W, at the nominal (unscattered) saturation age, the median of
    `saturation_time`, and the decay of `DECAY_EXPONENTS`; the photosphere's share is not added (the history of an escape model needs the shape of the decay)."""
    t_sat, decay = saturation_age_gy(mass_sol), decay_exponent(mass_sol)
    return coronal_xuv_w(lx_lbol_at(age_gy, saturated_lx_lbol(mass_sol), t_sat, decay) * luminosity_w)


def effective_saturation_age_gy(mass_sol, age_gy, log_lx_lbol=None):
    """
    The saturation age that reproduces a star's own `L_X/L_bol` at `age_gy` in the decay of its mass (GEN.804): a star that spun slowly or fast, or that
    drew a high or low decay offset, ran along the decay of its mass shifted in time, and its stored `log_lx_lbol` fixes the shift. The nominal
    age without the stored value, and the older of the two for a star still saturated (its saturation has not ended, so its end is not known).
    """
    nominal, decay = saturation_age_gy(mass_sol), decay_exponent(mass_sol)
    if log_lx_lbol is None or decay <= 0.0 or age_gy <= 0.0:
        return nominal
    ratio = 10 ** log_lx_lbol / saturated_lx_lbol(mass_sol)
    if ratio >= 1.0:
        return max(age_gy, nominal)
    return age_gy * ratio ** (1.0 / decay)


def coronal_xuv_with(mass_sol, at_gy, luminosity_w, t_sat_gy):
    """`coronal_xuv_at` along the decay with a given saturation age (`effective_saturation_age_gy` for a star's own)."""
    decay = decay_exponent(mass_sol)
    return coronal_xuv_w(lx_lbol_at(at_gy, saturated_lx_lbol(mass_sol), t_sat_gy, decay) * luminosity_w)


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
    gy_s = _gyr_to_s(1.0)
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
    t_sat, decay = saturation_age_gy(1.0), decay_exponent(1.0)
    return _fluence_j(constants.SOLAR_LUMINOSITY, 5772.0, 1.0, REFERENCE_AGE_GY, t_sat, decay, 0.0)


REFERENCE_FLUENCE_J = reference_fluence_j()


def _galactic_radius_pc(galactic_distance_ly):
    return _m_to_pc(_ly_to_m(galactic_distance_ly))


def supernova_rate_per_gyr(galactic_distance_ly):
    """GEN.87 rule 7: supernovae within 8 pc per Gyr at this distance from
    the galactic centre (1.5 at the Sun's 8 kpc, up to 250 times that)."""
    boost = min(SN_MAX_BOOST, math.exp(
        (SUN_GALACTIC_RADIUS_PC - _galactic_radius_pc(galactic_distance_ly)) / SN_SCALE_LENGTH_PC))
    return SN_RATE_AT_SUN_PER_GYR * boost


def grb_rate_per_gyr(galactic_distance_ly):
    """GEN.189: gamma-ray bursts per Gyr that reach this distance from the
    centre with an ozone-stripping fluence: `GRB_LETHAL_RATE_AT_SUN_PER_GYR`
    at the Sun's distance and inward, falling off outward."""
    boost = min(1.0, math.exp(
        (SUN_GALACTIC_RADIUS_PC - _galactic_radius_pc(galactic_distance_ly)) / SN_SCALE_LENGTH_PC))
    return GRB_LETHAL_RATE_AT_SUN_PER_GYR * boost


def agn_rate_per_gyr(galactic_distance_ly):
    """GEN.189: nuclear-activity episodes per Gyr that strip the ozone of a
    planet this far from the centre: `AGN_EPISODES_PER_GYR` inside
    `AGN_LETHAL_RADIUS_PC`, none beyond."""
    return AGN_EPISODES_PER_GYR if _galactic_radius_pc(galactic_distance_ly) <= AGN_LETHAL_RADIUS_PC else 0.0


def lethal_event_rate_per_gyr(galactic_distance_ly):
    """Ozone-stripping galactic events per Gyr at this distance from the
    galactic centre: nearby supernovae (GEN.87), gamma-ray bursts and the
    nucleus (GEN.189), or `None` for a star with no galaxy placement."""
    if galactic_distance_ly is None:
        return None
    return (supernova_rate_per_gyr(galactic_distance_ly) + grb_rate_per_gyr(galactic_distance_ly)
            + agn_rate_per_gyr(galactic_distance_ly))


def set_galactic_hazard(star):
    """Sets `lethal_event_rate_per_gyr` from the star's galaxy placement; a star with none is rated where a system with no place in the galaxy
    is taken to be (`scales.default_galactocentric_ly`, where its metallicity and orbit are also taken), so the column is filled for every star."""
    distance = getattr(star, "galactic_center_dist_ly", None)
    star.lethal_event_rate_per_gyr = lethal_event_rate_per_gyr(scales.default_galactocentric_ly() if distance is None else distance)


def rotation_saturation_scale(nominal_period_hours, period_hours, decay):
    """
    The factor on a star's saturation age from its own rotation (GEN.807): the unsaturated `L_X/L_bol` falls as `Ro^-2.7`, the Rossby number
    follows the period, and the period grows as `t^0.55` (`spin.GYRO_N`), so a star whose period is `nominal / x` is as active as a nominal star
    `x^(2.7 / b)` times younger in the decay `b`. 1 for a star with no period or no nominal one.
    """
    if not nominal_period_hours or not period_hours or period_hours <= 0.0 or decay <= 0.0:
        return 1.0
    return (nominal_period_hours / period_hours) ** (ROTATION_SLOPE / decay)


def generate_activity(star):
    """Sets `STAR_ACTIVITY_FIELDS` on an ordinary star (main sequence,
    subgiant, subdwarf, giant or white dwarf)."""
    set_galactic_hazard(star)
    mass_sol = _kg_to_msun(star.mass)
    luminosity = star.luminosity or 0.0
    age = max(star.age or 0.0, 1e-4)
    photosphere = luminosity * blackbody_fraction_above(XUV_EDGE_EV, star.temperature)
    yerkes = star.yerkes_class or ""
    if yerkes == "VII":
        _set_remnant(star, photosphere, age)
        return
    t_sat, decay = saturation_age_gy(mass_sol, draw.uniform(*PERCENTILE_RANGE)), decay_exponent(mass_sol)
    offset = draw.gauss(0.0, DECAY_SIGMA_DEX)
    log_n33 = flare_log_n33(mass_sol, age) + draw.uniform(-FLARE_SCATTER_DEX, FLARE_SCATTER_DEX)
    star.flare_alpha = draw.uniform(*FLARE_ALPHA_RANGE)
    saturated_ratio = saturated_lx_lbol(mass_sol)
    spin_scale = rotation_saturation_scale(
        spin.nominal_rotation_period_hours(star.mass, star.temperature, age) if star.temperature else None,
        getattr(star, "rotation_period_hours", None), decay)
    if spin_scale != 1.0:     # GEN.807: the star's own spin moves its activity, and its flares, about the table's mean
        log_n33 += math.log10(lx_lbol_at(age, saturated_ratio, t_sat * spin_scale, decay) / lx_lbol_at(age, saturated_ratio, t_sat, decay))
        t_sat *= spin_scale
    if yerkes in ("IV", "V", "VI", "PMS", "protostar"):  # a young star is saturated (GEN.242)
        ratio = lx_lbol_at(age, saturated_lx_lbol(mass_sol), t_sat, decay, offset)
        star.xuv_saturated = bool(age <= t_sat)
        star.flare_n33_per_yr = 10 ** log_n33
        star.xuv_fluence_j = _fluence_j(luminosity, star.temperature, mass_sol, age, t_sat, decay, offset)
    else:
        # A giant has left the saturated dynamo behind; its early fluence
        # is the main-sequence one, which this does not track.
        ratio = GIANT_LX_LBOL
        star.xuv_saturated = False
        star.flare_n33_per_yr = 0.0
        star.xuv_fluence_j = _yr_to_s(_gyr_to_yr((coronal_xuv_w(ratio * luminosity) + photosphere) * age))
    star.log_lx_lbol = math.log10(ratio)
    star.l_xuv_w = coronal_xuv_w(ratio * luminosity) + photosphere
    star.uv_c_share = _uv_c_share(star)


def uv_shares_of(star):
    """
    The UV-C, UV-B and UV-A power of an ordinary star as multiples of its blackbody's (GEN.971, `stellar_library.uv_shares`): the
    model atmosphere of its effective temperature and surface gravity plus the far ultraviolet of its corona, which is its coronal
    XUV (`log_lx_lbol`; the photosphere's own XUV is in the model). None when the star lacks one of those quantities or the library
    is missing and the run disregards that.
    """
    mass, radius = getattr(star, "mass", None), getattr(star, "radius", None)
    temperature, luminosity = getattr(star, "temperature", None), getattr(star, "luminosity", None)
    if not (mass and radius and temperature and luminosity):
        return None
    log_ratio = getattr(star, "log_lx_lbol", None)
    corona = 0.0 if log_ratio is None else coronal_xuv_w(10 ** log_ratio * luminosity) / luminosity
    return stellar_library.uv_shares(mass, radius, temperature, corona)


SUN_AGE_GY = 4.57
"""float: The age of the Sun the far-ultraviolet flux of `far_uv_flux_ratio` is measured against."""
SUN_TEMPERATURE_K = 5772.0
SUN_LOGG = 4.438
"""float: The effective temperature (IAU 2015) and the surface gravity of the Sun, the model atmosphere that stands for the reference star."""


def _sun_far_uv_flux_w_m2():
    """The 91.2 to 200 nm flux the Sun gives at 1 AU in the same model as any other star (`far_uv_flux_ratio` is 1 there by construction)."""
    luminosity = constants.SOLAR_LUMINOSITY
    corona = coronal_xuv_at(1.0, SUN_AGE_GY, luminosity) / luminosity
    return stellar_library.far_uv_power_w(SUN_TEMPERATURE_K, SUN_LOGG, luminosity, corona) / (4.0 * math.pi * constants.AU_M ** 2)


def far_uv_flux_ratio(star, distance_au):
    """
    The 91.2 to 200 nm flux of `star` at `distance_au`, relative to the quiet Sun at 1 AU (GEN.982): the number `physics.photolysis.rate_per_s`
    takes as `flux_ratio`. The star's model-atmosphere share plus its corona's far ultraviolet (`stellar_library.far_uv_power_w`) over the
    area of the sphere at that distance, over the Sun's own flux in the same model. `None` when the star lacks a mass, radius, temperature
    or luminosity, or the distance is not above zero.
    """
    temperature, luminosity = getattr(star, "temperature", None), getattr(star, "luminosity", None)
    if not (getattr(star, "mass", None) and getattr(star, "radius", None) and temperature and luminosity and distance_au and distance_au > 0.0):
        return None
    log_ratio = getattr(star, "log_lx_lbol", None)
    corona = 0.0 if log_ratio is None else coronal_xuv_w(10 ** log_ratio * luminosity) / luminosity
    power = stellar_library.far_uv_power_w(temperature, stellar_library.surface_logg(star.mass, star.radius), luminosity, corona)
    return power / (4.0 * math.pi * (distance_au * constants.AU_M) ** 2) / _sun_far_uv_flux_w_m2()


def _uv_c_share(star):
    """The star's `uv_c_share`: the library's (GEN.971); the seeded GEN.443 share, from the star's own stream so no draw of the run's
    moves, for a run that disregards a missing library."""
    shares = uv_shares_of(star)
    if shares is not None:
        return shares["uv_c"]
    return stellar_spectrum.draw_uv_c_share(draw.Stream(f"uv_c_share:{star.mass!r}:{star.age!r}:{star.temperature!r}"))


def _set_remnant(star, l_xuv_w, age_gy):
    """A remnant's activity: no corona or flares, `l_xuv_w` held over its
    age (its cooling is not integrated)."""
    star.log_lx_lbol = None
    star.xuv_saturated = False
    star.flare_n33_per_yr = 0.0
    star.flare_alpha = None
    star.uv_c_share = None
    star.l_xuv_w = l_xuv_w
    star.xuv_fluence_j = _yr_to_s(_gyr_to_yr(l_xuv_w * age_gy))


def neutron_star_spin_down_w(neutron_star):
    """A neutron star's magnetic-dipole spin-down power, W (`B^2 R^6 Omega^4 / 6 c^3`, cgs)."""
    radius_cm = neutron_star.radius * 1e5
    omega = 2 * math.pi / _ms_to_s(neutron_star.spin_period_ms)
    c_cm = constants.SPEED_OF_LIGHT_M_S * 100.0
    edot_erg = neutron_star.magnetic_field_gauss ** 2 * radius_cm ** 6 * omega ** 4 / (6 * c_cm ** 3)
    return edot_erg / ERG_PER_J


def neutron_star_xuv_w(neutron_star):
    """A neutron star's XUV, W: its thermal surface output or `DIPOLE_EDOT_TO_LX` of its spin-down power
    (`neutron_star_spin_down_w`), whichever is larger."""
    return max(neutron_star.luminosity, DIPOLE_EDOT_TO_LX * neutron_star_spin_down_w(neutron_star))


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
    d_m = _au_to_m(distance_au)
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
