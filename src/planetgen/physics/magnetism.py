"""
Planetary magnetic fields and the star's radiation at the planet (GEN.86).

The rules are docs/design/habitability/activity-magnetism-radiation-hydrosphere.md
section 3.5 (design defaults):

1. **Rocky bodies and moons.** The dynamo is the core's (GEN.838, `core_heat`): the interior heat model keeps a core temperature beside the mantle's,
   the core loses heat across the core-mantle boundary, and the dynamo runs while that heat and the gravitational energy of the inner core's
   growth leave power over the core's own adiabatic conduction (and the liquid shell is thick enough). The moment is the dipole of the
   buoyancy-flux scaling of the liquid shell (`dynamo.field`: the power, the shell's radii and density, the spin), already reduced by the
   dipole's share, with a log-normal scatter of `FIELD_SCATTER_DEX`. A dead core has no moment; a moon's core is run the same way (GEN.967 replaces
   the Ganymede-like chance). A small, iron-rich body (Mercury) whose bulk density says so (mass 0.05 to 0.1 Earth, core fraction over 0.4) gets 3e19,
   log-normal 0.5 dex, whatever the model says: its iron-sulfide core melts at lower temperatures than the pure iron the model has, and Mercury's
   field is alive at 4.5 Gyr. Plate tectonics (a mobile lid) needs a liquid surface ocean to lubricate the plates: a body of 0.5 to 5 Earth masses
   with one has a 0.6 chance of a mobile lid, a heavier one 0.3 (the literature is unsettled about super-Earths, so the odds fall rather than rise),
   and any other rocky planet is stagnant (`mobile_lid_chance`). The regime is drawn once, with the body's interior (`interior_heat.generate_lid`,
   GEN.721), and stored as `lid_regime`; the mantle's heat flow (and so the core's) and the ground dose read it.
2. **Rotation.** Slow rotation makes the field multipolar: the local Rossby number of the liquid shell (`dynamo.local_rossby`) gates the dipole,
   `f_dip = 0.15 + 0.85 / (1 + (Ro_l / 0.12)^4)` (`dynamo.dipole_fraction`), and the stored moment is the dipole's, so it already holds the share. A tidally
   locked M-dwarf planet keeps about 0.15 of its field in the dipole.
3. **Giants.** The field is made in the convecting layer of the stored interior (the metallic envelope of a gas giant or brown dwarf, the ionic water
   of an ice giant's mantle), driven by a share of the heat the giant radiates (`giant_dynamo.giant_field`, `giant_heat.of`). A giant with no placed
   layer has no field.
4. **Moons.** See rule 1: a moon's core is run through the same model as a planet's.
5. **Magnetopause.** `R_mp / R_p = 9.75 (M_dip f_dip / M_earth)^(1/3)
   (p_sw / 2.24 nPa)^(-1/6)`, the wind pressure scaling with the star's
   surface X-ray flux to the 1.3 (capped at 10x solar, 1x under 0.4 solar
   masses) over distance squared. `None` without a field.
6. **Induced field of an ocean moon** (GEN.178). A moon with a conducting
   liquid ocean answers its parent's changing field with an induced dipole
   (Europa and Callisto). The parent's field at the moon's orbit,
   `mu0 M_p / (4 pi a^3)` (the stored moment is the dipole's), is seen to oscillate at the synodic period
   `1 / |1/P_rot - 1/P_orb|` with `EXCITATION_SHARE` (0.4) of its strength,
   which gives Europa's ~200 nT amplitude for Jupiter's 1.56e27 A m2. The
   ocean answers as a thin conducting shell, `|A| = x / sqrt(1 + x^2)`,
   `x = omega mu0 Sigma R / 2` for the shell's conductance `Sigma = sigma h`
   at the shell's radius R (Europa: A = 0.97 [R]), and the induced dipole is
   `2 pi R_moon^3 A B_exc / mu0`, stored as `induced_moment_a_m2` beside the
   intrinsic `magnetic_moment_a_m2`. `sigma` is a planetGen default per ocean
   class (`OCEAN_CONDUCTIVITY_S_M`; Europa's magnesium-sulfate ocean 2.75 S/m
   [[R1785]]). Nothing is drawn.
7. **Induced field of a metal core and of a planet in its star's field** (GEN.981). The iron core of a body (`interior_heat.core_of`, whatever its dynamo
   does: solid iron conducts as well as liquid) answers a changing field as a conducting sphere (`sphere_response`: `|F| = |1 - 3/z^2 + 3 cot(z)/z|`,
   `z = (1 - i) a / delta`, skin depth `delta = sqrt(2 / (omega mu0 sigma))`, `sigma` `CORE_CONDUCTIVITY_S_M`), which is 1 for a core many skin depths
   across (an Earth core at a day's period is 20 000 of them) and falls as `omega mu0 sigma a^2 / 15` for a thin one. A moon's core feels the
   field its ocean lets through (`sqrt(1 - A^2)` of it) and its moment adds to the ocean's. A rocky planet is excited by the star's field at its orbit,
   `IMF_EARTH_T` (5 nT) times the square root of the wind's pressure over Earth's (the field's energy a fixed share of the wind's), `STELLAR_EXCITATION_SHARE`
   of it perpendicular to the spin, turning at the planet's synodic rotation (a tidally locked planet at the star's synodic rotation instead); its ocean
   (a surface ocean's conductance scaled by the share of the surface it covers) and its core answer as a moon's do, and the moment is stored as `induced_moment_a_m2` beside the intrinsic one. A giant has none.

`dipole_class`: none, weak (under 0.01 Earth's), Earth-like, strong (over 5
Earth's) or multipolar (the dipole's share under 0.5).
"""

import cmath
import math

from planetgen import tuning
from planetgen.physics import activity, constants, core_heat, dynamo, giant_dynamo, giant_interior
from planetgen.util import draw
from planetgen.physics import unit_registry as units

_kg_to_mearth = units.converter("kg", "Mearth")
_kg_to_msun = units.converter("kg", "Msun")

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
LID_MOBILE = "mobile"
"""str: What a body's stored `lid_regime` is when it has plate tectonics (`interior_heat.generate_lid`, GEN.721)."""
MOBILE_LID_MIN_MASS_EARTH = 0.5
MOBILE_LID_MAX_MASS_EARTH = 5.0
MOBILE_LID_CHANCE = 0.6
MOBILE_LID_HEAVY_CHANCE = 0.3
"""float: Rule 1's lid regime (GEN.177; design note section 3.5 rule 2): the chance
that a planet with a liquid surface ocean has plate tectonics, by mass (planetGen defaults)."""
MERCURY_MOMENT_A_M2 = 3e19
CORE_POWER_SCALE = dynamo.EARTH_DYNAMO_POWER_W / core_heat.EARTH_DYNAMO_W
"""float: The core model's power as the dynamo scaling takes it: the model's Earth drives its dynamo with 2.4 TW (`core_heat.EARTH_DYNAMO_W`) and the scaling's Earth
with 1.15 TW (`dynamo.EARTH_DYNAMO_POWER_W`), so the ratio puts the model's Earth on the scaling's. Registry row DY6."""
MERCURY_SCATTER_DEX = 0.5
"""float: The log-normal scatter of the stand-in moment of a small, iron-rich body (GEN.86's, dex)."""
FIELD_SCATTER_DEX = 0.07
"""float: The log-normal scatter of the dynamo scaling, dex: the 145 numerical dynamos fit the buoyancy-flux law within about 15 percent ([R3355],
second-hand, `dynamo.DIPOLE_COEFFICIENT`). Registry row DY5."""
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
dilute ocean 0.5, Europa's MgSO4 2.75 [R1785], sodium carbonate
3, saturated chloride 5, sulfuric acid 10 (planetGen defaults)."""
CORE_CONDUCTIVITY_S_M = 1.0e6
"""float: The electrical conductivity of an iron-alloy core, S/m: Earth's outer core at its pressure and temperature, 1.0 to 1.5 million ([R1628, R3597] find about half), the same for a solid core. Registry row IN1."""
IMF_EARTH_T = 5e-9
"""float: The field a star's wind carries at Earth's distance, T: 5 nT at 1 AU in the Sun's (the interplanetary field of the Parker spiral). Registry row IN2."""
STELLAR_EXCITATION_SHARE = 0.5
"""float: The share of the wind's field that a planet sees as a changing field: the part perpendicular to its spin axis, averaged over the tilt of the field
(the mean of sin of an angle on a sphere is 0.785; the Parker spiral and the sector structure take it lower). Registry row IN3."""
SPHERE_SKIN_DEPTHS_MAX = 1e3
"""float: A core this many skin depths across answers with exactly 1 (the series of the exact response is no longer needed)."""
INDUCTION_HYDROSPHERES = ("ice-covered ocean", "surface ocean", "hycean")
"""tuple: The states with a liquid ocean for an induced field (a `brine seep` is a puddle, not a conducting ocean layer; GEN.908)."""

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


def mobile_lid_chance(body, mass_earth):
    """The chance that a rocky planet has plate tectonics (a mobile lid): only
    a body with a liquid surface ocean can, and only from `MOBILE_LID_MIN_MASS_EARTH`
    (smaller worlds lose their heat too fast)."""
    if getattr(body, "hydrosphere", None) != "surface ocean" or mass_earth < MOBILE_LID_MIN_MASS_EARTH:
        return 0.0
    return MOBILE_LID_CHANCE if mass_earth <= MOBILE_LID_MAX_MASS_EARTH else MOBILE_LID_HEAVY_CHANCE


def _core_field(body, core):
    """The `dynamo.Field` of a rocky body's core (`core_heat.Core`, `interior_heat.core_of`), or `None` for a body with no metal or a dead core."""
    if core is None or not core.running:
        return None
    return dynamo.field(core.dynamo_w * CORE_POWER_SCALE, core.radius_m, core.inner_fraction * core.radius_m, core.density_kg_m3, giant_dynamo.spin_rad_s(body.rotation_period_hours))


def _giant_field(body):
    """The `dynamo.Field` of a giant from the layers it stores and the heat it radiates (`giant_dynamo.giant_field`), or `None` when it has no placed
    region that convects."""
    interior = giant_interior.interior_of(body)
    from planetgen.physics import giant_heat  # late: giant_heat reaches interior_heat, which reads this module's LID_MOBILE
    heat = giant_heat.of(body)
    if interior is None or heat is None or not body.rotation_period_hours:
        return None
    radius_m = body.radius * constants.KM_TO_M_FACTOR
    found = giant_dynamo.giant_field(interior, radius_m, heat.flux_w_m2 * 4.0 * math.pi * radius_m ** 2, giant_dynamo.spin_rad_s(body.rotation_period_hours))
    return None if found is None else found.field


def _moment(body, mass_earth, core):
    """The dipole moment and the dipole's share of the dynamo that runs in `body`: `(0, 1)` without one. The rocky core and the giant's
    metallic or ionic layer are the only sources (`dynamo.field`); Mercury's iron-rich small body keeps the one stand-in value."""
    if body.body_type == "g":
        found = _giant_field(body)
    else:
        density = getattr(body, "density", None) or EARTH_DENSITY
        if not body.is_moon and MIN_DYNAMO_MASS_EARTH <= mass_earth < 0.1 and core_mass_fraction(mass_earth, density) > 0.4:
            return MERCURY_MOMENT_A_M2 * _log_normal(MERCURY_SCATTER_DEX), 1.0
        found = _core_field(body, core)
    scatter = _log_normal(FIELD_SCATTER_DEX)
    if found is None or found.moment_a_m2 <= 0.0:
        return 0.0, 1.0
    return found.moment_a_m2 * scatter, found.dipole_fraction


def generate_field(body, core=None):
    """Draws `magnetic_moment_a_m2` (the dipole's moment, already reduced by the dipole's share) and `dipole_class` for a planet or moon (call once its
    rotation period is known), and its magnetopause when `update_exposure` has given its distance from the star. `core` is a rocky planet's or moon's
    `core_heat.Core` (`interior_heat.core_of`, GEN.838); a giant's field comes from its stored layers and radiated heat (GEN.967)."""
    mass_earth = _kg_to_mearth(body.mass)
    moment, share = _moment(body, mass_earth, core)
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


def sphere_response(omega, conductivity_s_m, radius_m):
    """
    `|F|`, the size of the answer of a homogeneous conducting sphere to a uniform field that oscillates at `omega` rad/s (0 to 1, rule 7): the induced
    dipole is `2 pi R^3 |F| B / mu0`. `F = 1 - 3/z^2 + 3 cot(z)/z` for `z = (1 - i) R / delta`, the skin depth `delta = sqrt(2 / (omega mu0 sigma))`:
    `omega mu0 sigma R^2 / 15` for a sphere much thinner than a skin depth, 1 for one much thicker. 0 for no frequency, conductor or size.
    """
    if not omega or omega <= 0.0 or conductivity_s_m <= 0.0 or radius_m <= 0.0:
        return 0.0
    skin_depth = math.sqrt(2.0 / (omega * VACUUM_PERMEABILITY * conductivity_s_m))
    if radius_m / skin_depth > SPHERE_SKIN_DEPTHS_MAX:
        return 1.0
    z = complex(1.0, -1.0) * radius_m / skin_depth
    return min(abs(1.0 - 3.0 / z ** 2 + 3.0 / (cmath.tan(z) * z)), 1.0)


def core_induced_moment(core, excitation_t, synodic_s, transmission=1.0):
    """The induced dipole moment of a metal core (`core_heat.Core`) in a field of amplitude `excitation_t` T that turns at the period `synodic_s`, A m2 (rule 7);
    `transmission` is the share of the field that reaches it past an ocean. 0 for no core or excitation."""
    if core is None or not excitation_t or not synodic_s:
        return 0.0
    response = sphere_response(2.0 * math.pi / synodic_s, CORE_CONDUCTIVITY_S_M, core.radius_m)
    return 2.0 * math.pi * core.radius_m ** 3 * response * transmission * excitation_t / VACUUM_PERMEABILITY


def _ocean_and_core_moment(body, excitation_t, synodic_s):
    """The induced moment of a body's ocean and its core in a changing field (rules 6 and 7): the ocean answers as a thin shell, the core as a sphere in the field the
    ocean lets through."""
    from planetgen.physics import interior_heat  # late: interior_heat reads this module's LID_MOBILE
    radius_m = _km_to_m(body.radius)
    moment, response = 0.0, 0.0
    depth_m = _km_to_m(getattr(body, "ocean_depth_km", None) or 0.0)
    if getattr(body, "hydrosphere", None) in INDUCTION_HYDROSPHERES and depth_m > 0.0 and excitation_t and synodic_s:
        shell_radius_m = max(radius_m - _km_to_m(getattr(body, "ice_shell_km", None) or 0.0), 0.0)
        conductivity = OCEAN_CONDUCTIVITY_S_M.get(getattr(body, "ocean_class", None), OCEAN_CONDUCTIVITY_S_M[None])
        omega = 2.0 * math.pi / synodic_s
        coverage = getattr(body, "ocean_fraction", None) if body.hydrosphere == "surface ocean" else None      # an ocean under ice or a hycean one is global
        coverage = 1.0 if coverage is None else min(max(coverage, 0.0), 1.0)
        x = omega * VACUUM_PERMEABILITY * conductivity * depth_m * coverage * shell_radius_m / 2.0
        response = x / math.sqrt(1.0 + x * x)
        moment = 2.0 * math.pi * radius_m ** 3 * response * excitation_t / VACUUM_PERMEABILITY
    core = interior_heat.core_of(body) if getattr(body, "body_type", None) != "g" else None
    return moment + core_induced_moment(core, excitation_t, synodic_s, math.sqrt(max(1.0 - response ** 2, 0.0)))


def generate_induced_field(moon, parent):
    """Sets `induced_moment_a_m2` on a moon from its parent planet's field,
    its own ocean and its metal core (rules 6 and 7; no draws): 0 for a moon with neither,
    no parent field or a parent that turns with it. Call once the moon's
    water, orbit and core are set and the parent's field exists."""
    moon.induced_moment_a_m2 = 0.0
    parent_moment = (getattr(parent, "magnetic_moment_a_m2", None) or 0.0)
    synodic = synodic_period_s(getattr(parent, "rotation_period_hours", None), moon.period)
    orbit_radius_m = _au_to_m(moon.distance)
    if not parent_moment or not synodic or orbit_radius_m <= 0.0:
        return
    excitation_t = EXCITATION_SHARE * VACUUM_PERMEABILITY * parent_moment / (4.0 * math.pi * orbit_radius_m ** 3)
    moon.induced_moment_a_m2 = _ocean_and_core_moment(moon, excitation_t, synodic)


def generate_stellar_induction(planet):
    """Sets `induced_moment_a_m2` on a rocky planet from its star's field at its orbit (rule 7; no draws): the wind's field `IMF_EARTH_T` scaled by the wind's
    pressure, seen at the planet's synodic rotation (the star's, for a planet that turns with its orbit), answered by its ocean and its metal core. 0 for a body
    with neither or with no orbit; a giant is left as it is. Call once the water, the core and the distance from the star are set."""
    if planet.body_type == "g":
        return
    planet.induced_moment_a_m2 = 0.0
    distance = getattr(planet, "_star_distance_au", None)
    period = getattr(planet, "period", None)
    if not distance or not period:
        return
    synodic = synodic_period_s(getattr(planet, "rotation_period_hours", None), period)
    if synodic is None:
        synodic = synodic_period_s(getattr(planet.star, "rotation_period_hours", None), period)
    if synodic is None:
        return
    excitation_t = STELLAR_EXCITATION_SHARE * IMF_EARTH_T * math.sqrt(wind_pressure_pa(planet.star, distance) / EARTH_WIND_PA)
    planet.induced_moment_a_m2 = _ocean_and_core_moment(planet, excitation_t, synodic)


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
    mass_sol = _kg_to_msun(star.mass)
    cap = M_DWARF_WIND_CAP if mass_sol < M_DWARF_WIND_MASS else WIND_CAP_SOLAR
    return EARTH_WIND_PA * min(surface_flux ** 1.3, cap) / distance_au ** 2


def magnetopause_rp(body, star_distance_au):
    """Rule 5's standoff in radii of `body` (GEN.843: Earth's 9.75 Earth radii scaled to the field and the wind, then divided by the
    body's own radius), or `None` with no field."""
    moment = getattr(body, "magnetic_moment_a_m2", None)
    if not moment:
        return None
    pressure = wind_pressure_pa(body.star, star_distance_au)
    standoff_earth_radii = (EARTH_STANDOFF_RP * (moment / EARTH_MOMENT_A_M2) ** (1 / 3)
                            * (pressure / EARTH_WIND_PA) ** (-1 / 6))
    return standoff_earth_radii * constants.EARTH_RADIUS_KM / body.radius


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
