# planetgen/physics/radiation.py

"""
Surface radiation dose, UV and galactic hazards (GEN.87).

The rules are docs/design/activity-magnetism-radiation-hydrosphere.md
section 4 (design defaults). A rocky body's yearly dose, mSv/yr, is

    surface = gcr * helio_mult + sep + ground

1. **Cosmic rays (`dose_gcr_msv_yr`).** The dose under a column of air `X =
   P0 / g` (g/cm2) is log-log interpolated through the Moon (X 0, 500
   mSv/yr), Mars (16.4, 237) and a 100 g/cm2 point (110), then falls with an
   attenuation length of 164 g/cm2 (sea level, 1033, gives 0.37 against
   Earth's 0.39 observed), floored at 1e-3 of sea level. A dipole cuts it
   to 0.6 at most: `1 - 0.4 min(1, log10(Rc0 / 1 GV) / log10(15))`, the
   cutoff rigidity `Rc0 = 14.9 GV (M_dip f_dip / 7.768e22) (R_E / R_p)^3`.
2. **Compressed heliosphere (`dose_helio_mult`).** `1 + 1.5 exp(-X / 100)
   c`, at most 2.5, `c` running from 0 in open space to 1 when a cloud
   presses the star's heliosphere to 1 AU or less (`heliosphere_compression`).
   Open space at generation; `store.refresh_containment` updates it when
   the system is found inside a nebula or remnant.
3. **Stellar particles (`dose_sep_msv_yr`).** A free-space `20 mSv/yr
   sqrt(N33 / N33_sun) / d_AU^2`, times `exp(-p_cut(X) / 0.2 GV)` (the
   rigidity of the proton whose range is X: `R = 0.0022 E^1.77`) and the
   polar-cap share `1 - sqrt(1 - sqrt(R_eff / Rc0))` with `R_eff = max(0.5
   p_cut, 0.2 GV)` (1 without a field).
4. **Ground (`dose_ground_msv_yr`).** `0.48 mSv/yr a_rad H(age) / H(4.5
   Gyr)`, `a_rad` log-uniform 0.5 to 2 (`tuning.ROGUE_RADIOGENIC_ABUNDANCE_RANGE`),
   plus radon, 0.5 mSv/yr times `a_rad`, with air (1 kPa or more) and land.
5. **UV (`uv_surface_index`).** Earth's DNA-weighted surface UV is 1. A
   star's 200 to 300 nm blackbody power against the Sun's (floored at 0.1:
   cool dwarfs have chromospheric UV) times the bolometric flux against
   Earth's, times `(O3 / O3_Earth)^-1.6` (at most 1e3, bare ground), with
   `O3 / O3_Earth = sqrt(pO2 / 21.2 kPa)` (planetGen's Chapman-like
   default; none below 0.02 kPa).
6. **Ozone loss (`ozone_loss_flag`).** True when there is ozone and either
   the system's lethal galactic events exceed one per 100 Myr or particle
   events deliver 100 mSv/yr at the ozone layer (X = 10 g/cm2).
7. **Galactic hazard (`lethal_event_rate_per_gyr`, on the star).** A
   supernova within 8 pc, 1.5 per Gyr at 8 kpc from the centre, times
   `exp((8 kpc - R) / 1.45 kpc)` (250x at the centre).

A planet of a pulsar or a black hole is rated by GEN.89, not here.
"""

import math

from planetgen import tuning
from planetgen.physics import constants, magnetism, rogue_surface
from planetgen.util.random import log_uniform

PLANCK_J_S = 2 * math.pi * constants.REDUCED_PLANCK

BODY_FIELDS = (
    "surface_dose_msv_yr", "dose_gcr_msv_yr", "dose_sep_msv_yr", "dose_ground_msv_yr", "dose_helio_mult",
    "uv_surface_index", "ozone_loss_flag",
)
"""tuple: What `generate_dose` sets on a planet or moon, as stored in
`planets` and `moons`."""

GCR_POINTS = ((0.0, 500.0), (16.4, 237.0), (100.0, 110.0))
"""tuple: `(X g/cm2, mSv/yr)`: the Moon, Mars and a high-altitude point
(Chang'E-4, MSL RAD, a design value)."""
GCR_ATTENUATION_G_CM2 = 164.0
GCR_FLOOR_SHARE = 1e-3
SEA_LEVEL_G_CM2 = 1033.0

EARTH_DIPOLE_A_M2 = 7.768e22
EARTH_CUTOFF_GV = 14.9
MAGNETIC_CUT = 0.4
MAGNETIC_FULL_GV = 15.0

HELIO_MAX_EXTRA = 1.5
HELIO_FADE_G_CM2 = 100.0

SEP_FREE_MSV_YR = 20.0
SEP_R0_GV = 0.2
SEP_MIN_RIGIDITY_GV = 0.2
N33_SUN_PER_YR = 10 ** -2.5
PROTON_REST_MEV = 938.272
RANGE_COEFFICIENT = 0.0022
RANGE_EXPONENT = 1.77

GROUND_MSV_YR = 0.48
RADON_MSV_YR = 0.5
RADON_MIN_PRESSURE_PA = 1000.0
EARTH_AGE_GY = 4.5

UV_FLOOR = 0.1
UV_EXPONENT = -1.6
UV_MAX = 1e3
O3_MIN_PO2_KPA = 0.02
EARTH_PO2_KPA = 21.2
OZONE_LAYER_G_CM2 = 10.0
OZONE_SEP_MSV_YR = 100.0
SUN_TEMPERATURE_K = 5772.0

LETHAL_EVENT_LIMIT_PER_GYR = 10.0
"""float: One lethal event per 100 Myr: the rate above which galactic hazards
flag ozone loss (Boss's default, GEN.87)."""


def _interpolate_gcr(column_g_cm2):
    """Rule 1's dose for an unshielded sky, mSv/yr, at column mass X."""
    x = column_g_cm2
    points = GCR_POINTS
    if x >= points[-1][0]:
        dose = points[-1][1] * math.exp(-(x - points[-1][0]) / GCR_ATTENUATION_G_CM2)
        floor = GCR_FLOOR_SHARE * points[-1][1] * math.exp(-(SEA_LEVEL_G_CM2 - points[-1][0]) / GCR_ATTENUATION_G_CM2)
        return max(dose, floor)
    for (x0, d0), (x1, d1) in zip(points, points[1:]):
        if x <= x1:
            f = (math.log1p(x) - math.log1p(x0)) / (math.log1p(x1) - math.log1p(x0))
            return math.exp(math.log(d0) + f * (math.log(d1) - math.log(d0)))
    return points[-1][1]


def gcr_dose_msv_yr(column_g_cm2, cutoff_gv=0.0):
    """The cosmic-ray dose under `column_g_cm2`, with a dipole of
    equatorial cutoff `cutoff_gv` (0 for none)."""
    return _interpolate_gcr(column_g_cm2) * magnetic_factor(cutoff_gv)


def magnetic_factor(cutoff_gv):
    """Rule 1's global-average cut of the cosmic-ray dose by a dipole."""
    if cutoff_gv <= 1.0:
        return 1.0
    return 1.0 - MAGNETIC_CUT * min(1.0, math.log10(cutoff_gv) / math.log10(MAGNETIC_FULL_GV))


def cutoff_rigidity_gv(dipole_moment_a_m2, radius_km):
    """`Rc0` for a dipole moment (A m2, the dipole part) on a body of
    `radius_km`."""
    if not dipole_moment_a_m2 or radius_km <= 0:
        return 0.0
    return EARTH_CUTOFF_GV * dipole_moment_a_m2 / EARTH_DIPOLE_A_M2 * (constants.EARTH_RADIUS_KM / radius_km) ** 3


def heliosphere_compression(open_radius_au, compressed_radius_au):
    """0 in open space to 1 when the heliosphere is pressed to 1 AU or
    less (log scale)."""
    if not compressed_radius_au or compressed_radius_au >= open_radius_au or open_radius_au <= 1.0:
        return 0.0
    return min(1.0, math.log(open_radius_au / compressed_radius_au) / math.log(open_radius_au))


def helio_multiplier(column_g_cm2, compression):
    """Rule 2's cosmic-ray multiplier (1 to 2.5)."""
    return 1.0 + HELIO_MAX_EXTRA * math.exp(-column_g_cm2 / HELIO_FADE_G_CM2) * min(max(compression, 0.0), 1.0)


def proton_cutoff_gv(column_g_cm2):
    """The rigidity, GV, of the proton whose range is `column_g_cm2`."""
    if column_g_cm2 <= 0:
        return 0.0
    kinetic = (column_g_cm2 / RANGE_COEFFICIENT) ** (1.0 / RANGE_EXPONENT)
    return math.sqrt(kinetic ** 2 + 2.0 * kinetic * PROTON_REST_MEV) / 1000.0


def sep_dose_msv_yr(free_msv_yr, column_g_cm2, cutoff_gv):
    """Rule 3's stellar-particle dose under `column_g_cm2` for a free-space
    dose and a dipole's `cutoff_gv`."""
    p_cut = proton_cutoff_gv(column_g_cm2)
    r_eff = max(0.5 * p_cut, SEP_MIN_RIGIDITY_GV)
    polar = 1.0 if cutoff_gv <= r_eff else 1.0 - math.sqrt(1.0 - math.sqrt(r_eff / cutoff_gv))
    return free_msv_yr * math.exp(-p_cut / SEP_R0_GV) * polar


def free_sep_msv_yr(flare_n33_per_yr, distance_au):
    """Rule 3's free-space dose at `distance_au`."""
    if not flare_n33_per_yr or not distance_au:
        return 0.0
    return SEP_FREE_MSV_YR * math.sqrt(flare_n33_per_yr / N33_SUN_PER_YR) / distance_au ** 2


def ground_dose_msv_yr(age_gy, mass_kg, radius_km, abundance, has_radon):
    """Rule 4's crust and radon dose."""
    now = rogue_surface.radiogenic_flux_w_m2(mass_kg, radius_km, age_gy)
    then = rogue_surface.radiogenic_flux_w_m2(mass_kg, radius_km, EARTH_AGE_GY)
    heat = now / then if then > 0 else 1.0
    return abundance * (GROUND_MSV_YR * heat + (RADON_MSV_YR if has_radon else 0.0))


def _band_fraction(temperature_k):
    """The share of a blackbody's power between 200 and 300 nm."""
    if temperature_k <= 0:
        return 0.0
    c1 = 2 * PLANCK_J_S * constants.SPEED_OF_LIGHT_M_S ** 2
    c2 = PLANCK_J_S * constants.SPEED_OF_LIGHT_M_S / (constants.BOLTZMANN * temperature_k)
    steps = 40
    low, high = 200e-9, 300e-9
    width = (high - low) / steps
    total = 0.0
    for i in range(steps + 1):
        wavelength = low + i * width
        x = c2 / wavelength
        radiance = 0.0 if x > 700 else c1 / (wavelength ** 5 * (math.exp(x) - 1.0))
        total += radiance * (0.5 if i in (0, steps) else 1.0) * width
    return total * math.pi / (constants.STEFAN_BOLTZMANN_CONSTANT * temperature_k ** 4)


_SUN_BAND = _band_fraction(SUN_TEMPERATURE_K)


def uv_ratio(temperature_k):
    """A star's 200 to 300 nm power per unit luminosity against the Sun's,
    floored at `UV_FLOOR`."""
    if not temperature_k:
        return UV_FLOOR
    return max(_band_fraction(temperature_k) / _SUN_BAND, UV_FLOOR)


def ozone_ratio(p_o2_kpa):
    """O3 against Earth's, from the oxygen partial pressure."""
    if not p_o2_kpa or p_o2_kpa < O3_MIN_PO2_KPA:
        return 0.0
    return min(math.sqrt(p_o2_kpa / EARTH_PO2_KPA), 1.2)


def uv_index(star, distance_au, p_o2_kpa):
    """Rule 5, or `None` without a star temperature."""
    temperature = getattr(star, "temperature", None)
    luminosity = getattr(star, "luminosity", None)
    if not temperature or not luminosity or not distance_au:
        return None
    flux = luminosity / constants.SOLAR_LUMINOSITY / distance_au ** 2
    ozone = ozone_ratio(p_o2_kpa)
    shield = UV_MAX if ozone <= 0 else min(UV_MAX, ozone ** UV_EXPONENT)
    return uv_ratio(temperature) * flux * shield


def _clear(body):
    for field in BODY_FIELDS:
        setattr(body, field, None)


def generate_dose(body):
    """Sets `BODY_FIELDS` on a planet or moon (call after its atmosphere,
    water and magnetic field and once `magnetism.update_exposure` has given
    its distance). A gas giant has no surface (every field `None`). Draws
    the crust's `a_rad` once."""
    _clear(body)
    if body.body_type == "g":
        return
    star = body.star
    distance = getattr(body, "_star_distance_au", None)
    pressure_pa = getattr(body, "atmospheric_pressure", None) or 0.0
    column = column_g_cm2(pressure_pa, getattr(body, "gravity", None))
    moment = (getattr(body, "magnetic_moment_a_m2", None) or 0.0) * magnetism.body_dipole_share(body)
    cutoff = cutoff_rigidity_gv(moment, body.radius)
    gcr = gcr_dose_msv_yr(column, cutoff)
    free = free_sep_msv_yr(getattr(star, "flare_n33_per_yr", None), distance)
    sep = sep_dose_msv_yr(free, column, cutoff)
    abundance = log_uniform(*tuning.ROGUE_RADIOGENIC_ABUNDANCE_RANGE)
    radon = pressure_pa >= RADON_MIN_PRESSURE_PA and (getattr(body, "land_fraction", None) or 0.0) > 0.0
    ground = ground_dose_msv_yr(getattr(star, "age", None) or EARTH_AGE_GY, body.mass, body.radius, abundance, radon)
    p_o2 = getattr(body, "p_o2_kpa", None) or 0.0
    body.dose_gcr_msv_yr = gcr
    body.dose_sep_msv_yr = sep
    body.dose_ground_msv_yr = ground
    body.dose_helio_mult = 1.0
    body.surface_dose_msv_yr = gcr + sep + ground
    body.uv_surface_index = uv_index(star, distance, p_o2)
    hazard = getattr(star, "lethal_event_rate_per_gyr", None)
    ozone_sep = sep_dose_msv_yr(free, OZONE_LAYER_G_CM2, cutoff)
    body.ozone_loss_flag = bool(ozone_ratio(p_o2) > 0.0 and (
        (hazard is not None and hazard > LETHAL_EVENT_LIMIT_PER_GYR) or ozone_sep >= OZONE_SEP_MSV_YR))


def column_g_cm2(pressure_pa, gravity_g):
    """X = P0 / g, g/cm2, from a stored pressure (Pa) and gravity (Earth
    g)."""
    gravity = (gravity_g or 0.0) * constants.EARTH_GRAVITY
    return pressure_pa / gravity * 0.1 if gravity > 0 and pressure_pa else 0.0


def total_dose_msv_yr(gcr, helio_mult, sep, ground):
    """`surface_dose_msv_yr` from its parts."""
    return gcr * helio_mult + sep + ground
