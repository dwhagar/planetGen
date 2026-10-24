"""
Planetary magnetic fields and the star's radiation at the planet (GEN.86).

The rules are docs/design/habitability/activity-magnetism-radiation-hydrosphere.md
section 3.5 (design defaults):

1. **Rocky bodies.** A dynamo needs mass of at least 0.05 Earth and a core
   mass fraction of at least 0.15. The core fraction comes from the bulk
   density (`core_mass_fraction`): the density is uncompressed with
   `1 + 0.25 sqrt(M)`, mixed from rock (3.3 g/cm3) and iron (7.9), and
   scaled so Earth gives 0.325. The dynamo lasts `4 Gyr * M^0.5 *
   U(0.5, 1.5)` (half as long under a core fraction of 0.2); the moment is
   `7.8e22 A m2 * M * (c / 0.325)`, log-normal 0.4 dex, fading as
   `(1 - age / tau)^0.3`. A small, iron-rich body (Mercury) gets 3e19,
   log-normal 0.5 dex. A planet with a **stagnant lid** (one rigid plate, no
   plate tectonics: Venus, Mars) cools its core slowly, so its dynamo lasts
   0.4 as long (GEN.177). Plate tectonics (a mobile lid) needs a liquid surface
   ocean to lubricate the plates: a body of 0.5 to 5 Earth masses with one has
   a 0.6 chance of a mobile lid, a heavier one 0.3 (the literature is
   unsettled about super-Earths, so the odds fall rather than rise), and any
   other rocky planet is stagnant (`mobile_lid_chance`).
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
6. **Induced field of an ocean moon** (GEN.178). A moon with a conducting
   liquid ocean answers its parent's changing field with an induced dipole
   (Europa and Callisto). The parent's field at the moon's orbit,
   `mu0 M_p f_dip / (4 pi a^3)`, is seen to oscillate at the synodic period
   `1 / |1/P_rot - 1/P_orb|` with `EXCITATION_SHARE` (0.4) of its strength,
   which gives Europa's ~200 nT amplitude for Jupiter's 1.56e27 A m2. The
   ocean answers as a thin conducting shell, `|A| = x / sqrt(1 + x^2)`,
   `x = omega mu0 Sigma R / 2` for the shell's conductance `Sigma = sigma h`
   at the shell's radius R (Europa: A = 0.97 [R]), and the induced dipole is
   `2 pi R_moon^3 A B_exc / mu0`, stored as `induced_moment_a_m2` beside the
   intrinsic `magnetic_moment_a_m2`. `sigma` is a planetGen default per ocean
   class (`OCEAN_CONDUCTIVITY_S_M`; Europa's magnesium-sulfate ocean 2.75 S/m
   [Hand and Chyba 2007]). Nothing is drawn.

`dipole_class`: none, weak (under 0.01 Earth's), Earth-like, strong (over 5
Earth's) or multipolar (`f_dip` under 0.5).
"""

import math

from planetgen import tuning
from planetgen.physics import activity, constants
from planetgen.util import draw
from planetgen.physics import unit_registry as units

_h_to_s = units.converter("h", "s")
_yr_to_s = units.converter("yr", "s")
_km_to_m = units.converter("km", "m")
_au_to_m = units.converter("au", "m")

BODY_FIELDS = (
    "magnetic_moment_a_m2", "dipole_class", "induced_moment_a_m2", "magnetopause_rp",
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
STAGNANT_LID_LIFETIME_FACTOR = 0.4
MOBILE_LID_MIN_MASS_EARTH = 0.5
MOBILE_LID_MAX_MASS_EARTH = 5.0
MOBILE_LID_CHANCE = 0.6
MOBILE_LID_HEAVY_CHANCE = 0.3
"""float: Rule 1's lid regime (GEN.177; design note section 3.5 rule 2): the dynamo
lifetime factor of a stagnant lid, and the chance that a planet with a liquid
surface ocean has plate tectonics, by mass (planetGen defaults)."""
MERCURY_MOMENT_A_M2 = 3e19
GANYMEDE_MOMENT_A_M2 = 1.3e20
GANYMEDE_MIN_MASS_EARTH = 0.02
GANYMEDE_CHANCE = 0.1
GIANT_ANCHORS = ((17.0, 3e24), (95.0, 5e25), (318.0, 1.5e27))
"""tuple: `(mass in Earth masses, dipole moment A m2)` for Neptune, Saturn
and Jupiter."""

VACUUM_PERMEABILITY = 4e-7 * math.pi
"""float: mu0, T m/A."""
EXCITATION_SHARE = 0.4
"""float: The share of the parent's equatorial dipole field that reaches a moon
as a changing field (tilt, plasma sheet): Jupiter's 520 nT at Europa gives the
~200 nT amplitude Galileo saw. A planetGen default."""
OCEAN_CONDUCTIVITY_S_M = {
    "ice-sealed": 0.5, "chloride brine": 5.0, "acid sulfate": 10.0, "soda": 3.0, "neutral": 2.75, None: 2.75,
}
"""dict: Electrical conductivity of an ocean by `ocean_class`, S/m: a cold
dilute ocean 0.5, Europa's MgSO4 2.75 (Hand and Chyba 2007), sodium carbonate
3, saturated chloride 5, sulfuric acid 10 (planetGen defaults)."""
INDUCTION_HYDROSPHERES = ("ice-covered ocean", "surface ocean", "hycean")
"""tuple: The states with a liquid ocean for an induced field."""

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


def body_dipole_share(body):
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


def mobile_lid_chance(body, mass_earth):
    """The chance that a rocky planet has plate tectonics (a mobile lid): only
    a body with a liquid surface ocean can, and only from `MOBILE_LID_MIN_MASS_EARTH`
    (smaller worlds lose their heat too fast)."""
    if getattr(body, "hydrosphere", None) != "surface ocean" or mass_earth < MOBILE_LID_MIN_MASS_EARTH:
        return 0.0
    return MOBILE_LID_CHANCE if mass_earth <= MOBILE_LID_MAX_MASS_EARTH else MOBILE_LID_HEAVY_CHANCE


def _rocky_moment(body, mass_earth, age_gy):
    """Rule 1's moment for a rocky planet, or 0."""
    core = core_mass_fraction(mass_earth, getattr(body, "density", None) or EARTH_DENSITY)
    lifetime = DYNAMO_LIFETIME_GY * math.sqrt(mass_earth) * draw.uniform(0.5, 1.5)
    strength = _log_normal(0.4)
    mobile_lid = draw.random() < mobile_lid_chance(body, mass_earth)
    if mass_earth < MIN_DYNAMO_MASS_EARTH or core < MIN_CORE_FRACTION:
        return 0.0
    if mass_earth < 0.1 and core > 0.4:
        return MERCURY_MOMENT_A_M2 * _log_normal(0.5)  # before the lifetime test: Mercury's field is alive at 4.5 Gyr
    if not mobile_lid:
        lifetime *= STAGNANT_LID_LIFETIME_FACTOR
    if core < 0.2:
        lifetime *= 0.5
    if age_gy >= lifetime:
        return 0.0
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
    share = body_dipole_share(body)
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


def synodic_period_s(parent_rotation_hours, orbital_period_years):
    """The period at which a moon sees the parent's field turn, s; `None` when
    the parent turns with the moon (tidally locked: nothing oscillates)."""
    if not parent_rotation_hours or not orbital_period_years:
        return None
    rate = abs(1.0 / _h_to_s(parent_rotation_hours)
               - 1.0 / _yr_to_s(orbital_period_years))
    return 1.0 / rate if rate > 1e-12 else None


def induced_moment(parent_moment_a_m2, orbit_radius_m, synodic_s, moon_radius_m, shell_conductance_s, shell_radius_m):
    """
    The induced dipole moment of a moon's ocean, A m2 (rule 6), and the
    response `|A|` it gave, as `(moment, response)`.

    Args:
        parent_moment_a_m2: The parent's dipole moment (times its dipole share).
        orbit_radius_m: The moon's distance from the parent.
        synodic_s: The excitation period (`synodic_period_s`), `None` for none.
        moon_radius_m: The moon's radius.
        shell_conductance_s: The ocean's conductance `sigma h`, siemens.
        shell_radius_m: The ocean shell's radius (the moon's less its ice).
    """
    if not parent_moment_a_m2 or not synodic_s or shell_conductance_s <= 0.0 or orbit_radius_m <= 0.0:
        return 0.0, 0.0
    excitation_t = EXCITATION_SHARE * VACUUM_PERMEABILITY * parent_moment_a_m2 / (4.0 * math.pi * orbit_radius_m ** 3)
    omega = 2.0 * math.pi / synodic_s
    x = omega * VACUUM_PERMEABILITY * shell_conductance_s * shell_radius_m / 2.0
    response = x / math.sqrt(1.0 + x * x)
    return 2.0 * math.pi * moon_radius_m ** 3 * response * excitation_t / VACUUM_PERMEABILITY, response


def generate_induced_field(moon, parent):
    """Sets `induced_moment_a_m2` on a moon from its parent planet's field
    and its own ocean (rule 6; no draws): 0 for a moon with no liquid ocean,
    no parent field or a parent that turns with it. Call once the moon's
    water and orbit are set and the parent's field exists."""
    moon.induced_moment_a_m2 = 0.0
    if getattr(moon, "hydrosphere", None) not in INDUCTION_HYDROSPHERES:
        return
    depth_m = _km_to_m(getattr(moon, "ocean_depth_km", None) or 0.0)
    if depth_m <= 0.0:
        return
    parent_moment = (getattr(parent, "magnetic_moment_a_m2", None) or 0.0) * body_dipole_share(parent)
    synodic = synodic_period_s(getattr(parent, "rotation_period_hours", None), moon.period)
    radius_m = _km_to_m(moon.radius)
    shell_radius_m = max(radius_m - _km_to_m(getattr(moon, "ice_shell_km", None) or 0.0), 0.0)
    conductivity = OCEAN_CONDUCTIVITY_S_M.get(getattr(moon, "ocean_class", None), OCEAN_CONDUCTIVITY_S_M[None])
    moon.induced_moment_a_m2 = induced_moment(
        parent_moment, _au_to_m(moon.distance), synodic, radius_m, conductivity * depth_m, shell_radius_m)[0]


PULSING_TYPES = ("young", "millisecond", "magnetar")


def pulsar_wind_pressure_pa(edot_w, distance_au):
    """The momentum flux of a pulsar wind that carries a spin-down power `edot_w`, Pa, at `distance_au`:
    `xi Edot / (4 pi r^2 c)` (GEN.317)."""
    distance_m = _au_to_m(distance_au)
    return tuning.PULSAR_WIND_XI * edot_w / (4 * math.pi * distance_m ** 2 * constants.SPEED_OF_LIGHT_M_S)


def wind_pressure_pa(star, distance_au):
    """The stellar wind's pressure at `distance_au` (rule 5); a pulsing neutron star's is its pulsar wind's
    (`pulsar_wind_pressure_pa`)."""
    if getattr(star, "pulsar_type", None) in PULSING_TYPES and distance_au:
        return pulsar_wind_pressure_pa(activity.neutron_star_spin_down_w(star), distance_au)
    if (getattr(star, "radius", None) is None or getattr(star, "luminosity", None) is None
            or getattr(star, "log_lx_lbol", None) is None):
        return EARTH_WIND_PA / distance_au ** 2
    l_x = 10 ** star.log_lx_lbol * star.luminosity
    radius_ratio = constants.SOLAR_RADIUS_M / _km_to_m(star.radius)
    surface_flux = (l_x / SOLAR_LX_W) * radius_ratio ** 2
    mass_sol = star.mass / constants.SOLAR_MASS_TO_KG
    cap = M_DWARF_WIND_CAP if mass_sol < M_DWARF_WIND_MASS else WIND_CAP_SOLAR
    return EARTH_WIND_PA * min(surface_flux ** 1.3, cap) / distance_au ** 2


def magnetopause_rp(body, star_distance_au):
    """Rule 5's standoff, planet radii, or `None` with no field."""
    moment = getattr(body, "magnetic_moment_a_m2", None)
    if not moment:
        return None
    share = body_dipole_share(body)
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
