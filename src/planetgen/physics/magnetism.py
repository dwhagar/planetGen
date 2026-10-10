"""
Planetary magnetic fields and the star's radiation at the planet (GEN.86).

The rules are docs/design/activity-magnetism-radiation-hydrosphere.md
section 3.5 (design defaults):

1. **Rocky bodies.** A dynamo needs mass of at least 0.05 Earth and a core
   mass fraction of at least 0.15. The core fraction comes from the bulk
   density (`core_mass_fraction`): the density is uncompressed with
   `1 + 0.25 sqrt(M)`, mixed from rock (3.3 g/cm3) and iron (7.9), and
   scaled so Earth gives 0.325. The dynamo lasts `4 Gyr * M^0.5 *
   U(0.5, 1.5)` (half as long under a core fraction of 0.2); the moment is
   `7.8e22 A m2 * M * (c / 0.325)`, log-normal 0.4 dex, fading as
   `(1 - age / tau)^0.3`. A small, iron-rich body (Mercury) gets 3e19,
   log-normal 0.5 dex. Whether the lid is stagnant is not modelled.
2. **Rotation.** Slow rotation makes the field multipolar: `Ro_l = 0.09 (P /
   24 h)^1.23` and `f_dip = 0.15 + 0.85 / (1 + (Ro_l / 0.12)^4)`. A tidally
   locked M-dwarf planet keeps about 0.15 of its moment in the dipole.
3. **Giants.** Always magnetized: log-log between Neptune (3e24 at 17 Earth
   masses), Saturn (5e25 at 95) and Jupiter (1.5e27 at 318), linear in mass
   above, log-normal 0.4 dex.
4. **Moons.** A rocky moon of at least 0.02 Earth masses is Ganymede-like
   with a chance of 0.1 (1.3e20, log-normal 0.7 dex); any other has none.
5. **Magnetopause.** `R_mp / R_p = 9.75 (M_dip f_dip / M_earth)^(1/3)
   (p_sw / 2.24 nPa)^(-1/6)`, the wind pressure scaling with the star's
   surface X-ray flux to the 1.3 (capped at 10x solar, 1x under 0.4 solar
   masses) over distance squared. `None` without a field.

`dipole_class`: none, weak (under 0.01 Earth's), Earth-like, strong (over 5
Earth's) or multipolar (`f_dip` under 0.5).
"""

import math

from planetgen.physics import activity, constants
from planetgen.util import draw

BODY_FIELDS = (
    "magnetic_moment_a_m2", "dipole_class", "magnetopause_rp",
    "xuv_flux_earth", "xuv_exposure_index", "flare_irradiation_index",
)
"""tuple: What this module sets on a planet or moon, as stored in `planets`
and `moons`."""

DIPOLE_CLASSES = ("none", "weak", "earth-like", "strong", "multipolar")

EARTH_MOMENT_A_M2 = 7.8e22
EARTH_CORE_FRACTION = 0.325
ROCK_DENSITY = 3.3
IRON_DENSITY = 7.9
EARTH_DENSITY = 5.51
MIN_DYNAMO_MASS_EARTH = 0.05
MIN_CORE_FRACTION = 0.15
DYNAMO_LIFETIME_GY = 4.0
MERCURY_MOMENT_A_M2 = 3e19
GANYMEDE_MOMENT_A_M2 = 1.3e20
GANYMEDE_MIN_MASS_EARTH = 0.02
GANYMEDE_CHANCE = 0.1
GIANT_ANCHORS = ((17.0, 3e24), (95.0, 5e25), (318.0, 1.5e27))
"""tuple: `(mass in Earth masses, dipole moment A m2)` for Neptune, Saturn
and Jupiter."""

EARTH_STANDOFF_RP = 9.75
EARTH_WIND_PA = 2.24e-9
SOLAR_LX_W = 1.4e20
"""float: The Sun's X-ray luminosity today (1.4e27 erg/s, what
`physics.activity` gives at 4.57 Gyr): the wind pressure's unit."""
WIND_CAP_SOLAR = 10.0
M_DWARF_WIND_CAP = 1.0
M_DWARF_WIND_MASS = 0.4


def _log_normal(sigma_dex):
    return 10 ** draw.gauss(0.0, sigma_dex)


def core_mass_fraction(mass_earth, density_g_cm3):
    """The iron core's share of a rocky body's mass, from its bulk density
    (scaled so Earth gives 0.325)."""
    def raw(mass, density):
        uncompressed = density / (1.0 + 0.25 * math.sqrt(max(mass, 0.0)))
        share = (1 / ROCK_DENSITY - 1 / uncompressed) / (1 / ROCK_DENSITY - 1 / IRON_DENSITY)
        return min(max(share, 0.0), 1.0)
    return min(1.0, raw(mass_earth, density_g_cm3) * EARTH_CORE_FRACTION / raw(1.0, EARTH_DENSITY))


def dipole_share(rotation_period_hours):
    """`f_dip`: the share of the moment in the dipole, falling for slow
    rotators (rule 2)."""
    if not rotation_period_hours:
        return 1.0
    rossby = 0.09 * (rotation_period_hours / 24.0) ** 1.23
    return 0.15 + 0.85 / (1.0 + (rossby / 0.12) ** 4)


def _share(body):
    """`dipole_share` for a rocky planet; giants and moons keep a dipole."""
    if body.body_type == "g" or body.is_moon:
        return 1.0
    return dipole_share(body.rotation_period_hours)


def _giant_moment(mass_earth):
    anchors = GIANT_ANCHORS
    if mass_earth <= anchors[0][0]:
        return anchors[0][1] * mass_earth / anchors[0][0]
    if mass_earth >= anchors[-1][0]:
        return anchors[-1][1] * mass_earth / anchors[-1][0]
    for (m0, b0), (m1, b1) in zip(anchors, anchors[1:]):
        if mass_earth <= m1:
            f = (math.log(mass_earth) - math.log(m0)) / (math.log(m1) - math.log(m0))
            return math.exp(math.log(b0) + f * (math.log(b1) - math.log(b0)))
    return anchors[-1][1]


def _rocky_moment(body, mass_earth, age_gy):
    """Rule 1's moment for a rocky planet, or 0."""
    core = core_mass_fraction(mass_earth, getattr(body, "density", None) or EARTH_DENSITY)
    lifetime = DYNAMO_LIFETIME_GY * math.sqrt(mass_earth) * draw.uniform(0.5, 1.5)
    strength = _log_normal(0.4)
    if mass_earth < MIN_DYNAMO_MASS_EARTH or core < MIN_CORE_FRACTION:
        return 0.0
    if core < 0.2:
        lifetime *= 0.5
    if age_gy >= lifetime:
        return 0.0
    if mass_earth < 0.1 and core > 0.4:
        return MERCURY_MOMENT_A_M2 * _log_normal(0.5)
    return (EARTH_MOMENT_A_M2 * mass_earth * (core / EARTH_CORE_FRACTION) * strength
            * (1.0 - age_gy / lifetime) ** 0.3)


def generate_field(body):
    """Draws `magnetic_moment_a_m2` and `dipole_class` for a planet or moon
    (call once its rotation period is known), and its magnetopause when
    `update_exposure` has given its distance from the star."""
    mass_earth = body.mass / constants.EARTH_MASS_TO_KG
    age_gy = getattr(body.star, "age", None) or 0.0
    if body.body_type == "g":
        moment = _giant_moment(mass_earth) * _log_normal(0.4)
    elif body.is_moon:
        ganymede = mass_earth >= GANYMEDE_MIN_MASS_EARTH and draw.random() < GANYMEDE_CHANCE
        moment = GANYMEDE_MOMENT_A_M2 * _log_normal(0.7) if ganymede else 0.0
    else:
        moment = _rocky_moment(body, mass_earth, age_gy)
    share = _share(body)
    body.magnetic_moment_a_m2 = moment
    if moment <= 0.0:
        body.dipole_class = "none"
    elif share < 0.5:
        body.dipole_class = "multipolar"
    elif moment < 0.01 * EARTH_MOMENT_A_M2:
        body.dipole_class = "weak"
    elif moment > 5 * EARTH_MOMENT_A_M2:
        body.dipole_class = "strong"
    else:
        body.dipole_class = "earth-like"
    distance = getattr(body, "_star_distance_au", None)
    body.magnetopause_rp = magnetopause_rp(body, distance) if distance else None


def wind_pressure_pa(star, distance_au):
    """The stellar wind's pressure at `distance_au` (rule 5)."""
    if (getattr(star, "radius", None) is None or getattr(star, "luminosity", None) is None
            or getattr(star, "log_lx_lbol", None) is None):
        return EARTH_WIND_PA / distance_au ** 2
    l_x = 10 ** star.log_lx_lbol * star.luminosity
    radius_ratio = constants.SOLAR_RADIUS_M / (star.radius * 1000.0)
    surface_flux = (l_x / SOLAR_LX_W) * radius_ratio ** 2
    mass_sol = star.mass / constants.SOLAR_MASS_TO_KG
    cap = M_DWARF_WIND_CAP if mass_sol < M_DWARF_WIND_MASS else WIND_CAP_SOLAR
    return EARTH_WIND_PA * min(surface_flux ** 1.3, cap) / distance_au ** 2


def magnetopause_rp(body, star_distance_au):
    """Rule 5's standoff, planet radii, or `None` with no field."""
    moment = getattr(body, "magnetic_moment_a_m2", None)
    if not moment:
        return None
    share = _share(body)
    pressure = wind_pressure_pa(body.star, star_distance_au)
    return (EARTH_STANDOFF_RP * (moment * share / EARTH_MOMENT_A_M2) ** (1 / 3)
            * (pressure / EARTH_WIND_PA) ** (-1 / 6))


def update_exposure(body, star_distance_au):
    """The star-distance values (no draws): XUV flux, exposure index, flare
    irradiation index and, once the field is drawn, the magnetopause.
    `star_distance_au` is the body's distance from its star (a moon's
    planet's)."""
    body._star_distance_au = star_distance_au
    star = body.star
    body.xuv_flux_earth = activity.xuv_flux_earth(star, star_distance_au)
    body.xuv_exposure_index = activity.exposure_index(star, star_distance_au)
    body.flare_irradiation_index = activity.flare_irradiation_index(star, star_distance_au)
    body.magnetopause_rp = magnetopause_rp(body, star_distance_au) if star_distance_au else None
