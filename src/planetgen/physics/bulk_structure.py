# planetgen/physics/bulk_structure.py

"""
Radius, density and gravity from the bulk fractions (GEN.429, docs/design/bulk-composition.md sections 4 and 7, the
contracts of docs/design/composition-spine-contracts.md)
=============================================================================================================

Pure functions, no draws and no state. A body's mass and its `BulkFractions` (`physics/bulk_composition.py`) give its
radius at the 1 bar level (Decided, Boss, 1:24 AM EDT 2026-10-17), and from the radius its density, surface gravity and
escape speed. Mass, radius and gravity are never drawn apart.

One reference solver and three fast forms:

- **The reference solver** `layered_structure` is part B's two-layer integrator (docs/design/body-composition-parts/B-rocky.md
  section 3.8, verified there against a scipy reference built from the full mantle tables: 0.20 percent in radius, 0.28 in
  the moment of inertia, 0.93 in the core-mantle pressure over 90 cases) extended with a water layer on top whose
  density is the 300 K isotherm of liquid water (IAPWS-95 through the NIST WebBook, to 950 MPa), the ice VI, VII and X
  reference points of the materials table (GEN.427) and [R3992]'s polytrope above 1 TPa (the way part C built its
  three-layer reference, within 1.5 percent of [R4881] for water-bearing cases from 0.5 to 50 Earth masses).
  Hydrostatic balance in s = m^(1/3), fourth-order Runge-Kutta, shooting on the central pressure to 1 bar at the surface.
  Validated 0.1 to 20 Earth masses; it is also used out of the closed forms' ranges (a moon, a dwarf planet) with the
  same equations, and then the caller reads `Structure.method`.
- **Dry rock and iron**: `radius_fit_earth`, B's fit to the solver, 0.1 to 20 Earth masses, core mass fraction 0 to 0.8,
  mantle Fe/(Fe+Mg) 0 to 0.3, core density factor 0.8 to 1.0; maximum 0.40 percent, rms 0.10 percent against the solver.
- **Water-bearing**: `water_world_radius_earth`, C's fit to its 175-point three-layer grid, 0.3 to 20 Earth masses, water
  fraction 0 to 1, core fraction of the dry part 0 to 0.5; rms 1.2 and at most 3.6 percent.
- **Small differentiated bodies**: `small_body_radius_earth`, a fit to the solver for 1.5e-4 to 0.3 Earth masses (rms 0.57,
  at most 2.0 percent), so a Mars, a Mercury or a Ganymede has a radius without running the solver.
- **An envelope on a core**: `envelope_radius_earth`, [R2741] Eqs. 1 and 3 (1 to 20 Earth masses, envelope
  0.01 to 20 percent, flux 0.1 to 1000 Earth's, age 0.1 to 10 Gyr; rms 6.3 percent at 1 Gyr against their tables), at the
  20 mbar level; the 1 bar radius adds ln(50) scale heights.

Below `SMALL_BODY_MASS_KG` a body is a heap of grains: `mixing_density_kg_m3` is the zero-pressure mixing rule with the
porosity the caller gives (the Bierson and Nimmo law and its taper are the small-body lane's, GEN.442).

A giant (envelope more than `GIANT_ENVELOPE_FRACTION` of the mass, or over 20 Earth masses of core) is outside this module:
`structure` raises `OutsideRange` and the giant's radius comes from the giant law (GEN.460).

Units: SI inside, Earth masses and radii in the names that say so; `Structure` carries km, g/cm3 and g.
"""

import math
from bisect import bisect_right
from dataclasses import dataclass

from planetgen.physics import constants, materials

G = constants.G
M_EARTH = constants.EARTH_MASS_TO_KG
R_EARTH = constants.EARTH_RADIUS_KM * 1000.0

P_SURFACE_PA = 1.0e5
"""float: The 1 bar reference level of the radius."""

SOLVER_MASS_EARTH = (1.0e-4, 50.0)
"""tuple: The masses (Earth masses) `layered_structure` takes; its validation is 0.1 to 20 (part B 3.8) and 0.5 to 50 for water
(part C 2.4), and below 0.1 it has no check beyond the equations (bulk-composition.md section 12)."""

CLOSED_FORM_RANGES = {
    "dry": {"mass_earth": (0.1, 20.0), "cmf": (0.0, 0.8), "x_fe": (0.0, 0.30), "core_scale": (0.8, 1.0)},
    "water": {"mass_earth": (0.3, 20.0), "water": (0.0, 1.0), "cmf": (0.0, 0.5)},
    "envelope": {"mass_earth": (1.0, 20.0), "envelope": (0.0001, 0.20), "flux_earth": (0.1, 1000.0), "age_gy": (0.1, 10.0)},
}
"""dict: The range each closed form was fitted over; outside it `structure` uses the solver (or refuses)."""

SMALL_BODY_MASS_KG = 1.0e21
"""float: Below this a body is a heap of grains with a porosity (bulk-composition.md decision 25; the endpoint is a stand-in
until viscous compaction replaces the taper, small-bodies-composition.md S9)."""

WATER_NEGLIGIBLE = 1.0e-3
"""float: Water and ice together under this share of the mass leave the radius where the dry body has it (the closed forms for water
start at 0.3 Earth masses and the dry ones at 0.1), so `structure` treats the body as dry."""

GIANT_ENVELOPE_FRACTION = 0.5
"""float: An envelope share above this makes the body a giant."""

# ---------------------------------------------------------------------------------------------------------------------
# Constants (sources in the module text and in part B 3.8, part C 2.3 to 2.5)
# ---------------------------------------------------------------------------------------------------------------------

CORE = (8264.2, 157.07e9, 5.6167)
"""tuple: The core's Vinet fit to hcp iron, (rho0 kg/m3, K0 Pa, K0'), [R1096] below 234.4 [R1752] above,
0.3 percent in density from 1 to 1000 GPa (part B 3.8)."""

_MANTLE_LOG_P = (5.5, 8.0, 9.0, 9.4, 9.7, 9.95, 10.1, 10.25, 10.35, 10.45, 10.6, 10.8, 11.0, 11.15, 11.3,
                 11.6, 11.9, 12.2, 12.45)
_MANTLE_RHO = (3055, 3060, 3097, 3154, 3234, 3361, 3457, 3880, 4049, 4351, 4494, 4748, 5091, 5534, 5951,
               7042, 8560, 10631, 12902)
"""The mantle's density along a 1600 K adiabat at Fe/(Fe+Mg) = 0 and Mg/Si = 1 by log10 P (Pa), from a Gibbs minimisation with
the [R4233] end-members (part B 3.8)."""
_MANTLE_LN_RHO = tuple(math.log(rho) for rho in _MANTLE_RHO)
_MANTLE_TOP_SLOPE = (_MANTLE_LN_RHO[-1] - _MANTLE_LN_RHO[-2]) / (_MANTLE_LOG_P[-1] - _MANTLE_LOG_P[-2])
MANTLE_IRON_FACTOR = 0.264
"""float: The mantle's density rises by 1 + this times Fe/(Fe+Mg)."""

# The water layer, the 300 K isotherm: (P Pa, rho kg/m3). Liquid water from IAPWS-95 (NIST WebBook fluid properties, read
# 2026-10-23), ice VI and VII and X from the materials table's reference points ([R2204] splines), the
# Seager polytrope above 1 TPa scaled to join (part C 2.4).
_WATER_POINTS = (
    (1.0e5, 996.56), (5.01e7, 1017.9), (1.001e8, 1037.2), (1.501e8, 1054.9), (2.001e8, 1071.1), (2.501e8, 1086.2),
    (3.001e8, 1100.2), (3.501e8, 1113.3), (4.001e8, 1125.7), (4.501e8, 1137.4), (5.001e8, 1148.5), (5.501e8, 1159.1),
    (6.001e8, 1169.2), (6.501e8, 1178.9), (7.001e8, 1188.2), (7.501e8, 1197.2), (8.001e8, 1205.8), (8.501e8, 1214.1),
    (9.001e8, 1222.2), (9.501e8, 1230.0),
)
_SEAGER_WATER = (1460.0, 0.00311, 0.513)
"""tuple: rho0 (kg/m3), c, n of [R3992: Table 3] for H2O, rho = rho0 + c P^n (P in Pa)."""

ICE_OTHER_DENSITY_KG_M3 = 1000.0
"""float: The grain density of the other ices (N2 950, NH3 850, CO 1000, CO2 1000 to 1600, CH4 500; bulk-composition.md section 6:
close to water ice)."""
ENVELOPE_MEAN_MOLAR_MASS = 2.35e-3
"""float: kg/mol of a solar H/He envelope (mean molecular weight 2.35), for the scale height at the 1 bar level."""


class OutsideRange(ValueError):
    """The body is outside what this module can give a radius for."""


@dataclass(frozen=True)
class Structure:
    """
    What a body's mass and fractions give.

    Attributes:
        radius_km (float): At the 1 bar level.
        density_g_cm3 (float): Mean density.
        gravity_g (float): Surface gravity in Earth's.
        escape_kms (float): Escape speed.
        moment_of_inertia (float or None): C/(M R^2), where the solver ran.
        cmb_pressure_gpa (float or None): At the core-mantle boundary, where the solver ran and the body has a core and a mantle.
        method (str): `closed-dry`, `closed-water`, `closed-small`, `closed-envelope`, `solver` or `mixing`.
    """

    radius_km: float
    density_g_cm3: float
    gravity_g: float
    escape_kms: float
    moment_of_inertia: float | None
    cmb_pressure_gpa: float | None
    method: str


# ---------------------------------------------------------------------------------------------------------------------
# Equations of state
# ---------------------------------------------------------------------------------------------------------------------

def core_density(pressure_pa, scale=1.0):
    """The iron core's density at `pressure_pa` from the Vinet fit (inverted by bisection on x = (rho0/rho)^(1/3)), times the
    light-element density factor `scale`."""
    rho0, k0, k0_prime = CORE
    if pressure_pa <= 0.0:
        return rho0 * scale
    low, high = 0.05, 1.0
    for _ in range(60):
        x = 0.5 * (low + high)
        p = 3.0 * k0 * (1.0 - x) / (x * x) * math.exp(1.5 * (k0_prime - 1.0) * (1.0 - x))
        if p > pressure_pa:
            low = x
        else:
            high = x
    return rho0 / (0.5 * (low + high)) ** 3 * scale


def mantle_density(pressure_pa, x_fe=0.0):
    """The silicate mantle's density at `pressure_pa`, with mantle Fe/(Fe+Mg) `x_fe`."""
    log_p = math.log10(max(pressure_pa, 1.0))
    if log_p <= _MANTLE_LOG_P[0]:
        ln_rho = _MANTLE_LN_RHO[0]
    elif log_p >= _MANTLE_LOG_P[-1]:
        ln_rho = _MANTLE_LN_RHO[-1] + _MANTLE_TOP_SLOPE * (log_p - _MANTLE_LOG_P[-1])
    else:
        i = bisect_right(_MANTLE_LOG_P, log_p) - 1
        t = (log_p - _MANTLE_LOG_P[i]) / (_MANTLE_LOG_P[i + 1] - _MANTLE_LOG_P[i])
        ln_rho = _MANTLE_LN_RHO[i] + t * (_MANTLE_LN_RHO[i + 1] - _MANTLE_LN_RHO[i])
    return math.exp(ln_rho) * (1.0 + MANTLE_IRON_FACTOR * x_fe)


def _water_nodes():
    """The water isotherm as (ln P, ln rho) nodes: liquid, then the reference points of ice VI, VII and X (300 K), then the
    polytrope scaled to join at 1 TPa."""
    points = list(_WATER_POINTS)
    for key in ("ice_vi", "ice_vii", "ice_x"):
        for pressure_mpa, temperature_k, density, *_ in materials.material(key).eos_of("reference-points").params["points"]:
            if abs(temperature_k - 300.0) < 0.5 and pressure_mpa * 1.0e6 > points[-1][0]:
                points.append((pressure_mpa * 1.0e6, density))
    last_pressure, last_density = points[-1]
    rho0, c, n = _SEAGER_WATER
    scale = last_density / (rho0 + c * last_pressure ** n)
    for exponent in (12.1, 12.5, 13.0, 13.5, 14.0, 15.0, 16.0):
        pressure = 10.0 ** exponent
        points.append((pressure, scale * (rho0 + c * pressure ** n)))
    return tuple((math.log(p), math.log(rho)) for p, rho in points)


_WATER_NODES = None


def water_density(pressure_pa):
    """The water layer's density at `pressure_pa` on the 300 K isotherm (log-log interpolation between the nodes)."""
    global _WATER_NODES
    if _WATER_NODES is None:
        _WATER_NODES = _water_nodes()
    nodes = _WATER_NODES
    ln_p = math.log(max(pressure_pa, 1.0e5))
    if ln_p >= nodes[-1][0]:
        (p0, r0), (p1, r1) = nodes[-2], nodes[-1]
    else:
        i = max(bisect_right([node[0] for node in nodes], ln_p) - 1, 0)
        (p0, r0), (p1, r1) = nodes[i], nodes[i + 1]
    return math.exp(r0 + (r1 - r0) * (ln_p - p0) / (p1 - p0))


# ---------------------------------------------------------------------------------------------------------------------
# The reference solver
# ---------------------------------------------------------------------------------------------------------------------

def _rk4_layer(state, s_from, s_to, steps, density_of):
    """Integrates [radius, pressure, moment] from s_from to s_to in s = m^(1/3), so dm = 3 s^2 ds."""
    h = (s_to - s_from) / steps

    def derivatives(s, r, p):
        rho = density_of(max(p, 1.0))
        return (3.0 * s * s / (4.0 * math.pi * r * r * rho),
                -3.0 * G * s ** 5 / (4.0 * math.pi * r ** 4),
                2.0 * r * r * s * s)

    s = s_from
    r, p, c = state
    for _ in range(steps):
        k1 = derivatives(s, r, p)
        k2 = derivatives(s + h / 2, r + h / 2 * k1[0], p + h / 2 * k1[1])
        k3 = derivatives(s + h / 2, r + h / 2 * k2[0], p + h / 2 * k2[1])
        k4 = derivatives(s + h, r + h * k3[0], p + h * k3[1])
        r += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        p += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        c += h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        s += h
    return [r, p, c]


def _integrate(central_pressure, mass_kg, layers, steps):
    """Integrates outward from the centre through `layers`, a list of (mass fraction, density function), innermost first.
    Returns the end state and the pressure at the top of the first layer when there is more than one."""
    first_density = layers[0][1](central_pressure)
    s0 = (1e-6 * mass_kg) ** (1.0 / 3.0)
    r0 = s0 * (3.0 / (4.0 * math.pi * first_density)) ** (1.0 / 3.0)
    p0 = central_pressure - (3.0 * G / (8.0 * math.pi)) * (4.0 * math.pi * first_density / 3.0) ** (4.0 / 3.0) * s0 * s0
    state = [r0, p0, 0.0]
    s = s0
    covered = 0.0
    boundary_pressure = None
    for index, (fraction, density_of) in enumerate(layers):
        covered += fraction
        s_to = (covered * mass_kg) ** (1.0 / 3.0)
        state = _rk4_layer(state, s, s_to, steps[index], density_of)
        s = s_to
        if index == 0 and len(layers) > 1:
            boundary_pressure = state[1]
    return state, boundary_pressure


def layered_structure(mass_earth, metal, rock, water, x_fe=0.11, core_scale=1.0, surface_pressure_pa=P_SURFACE_PA,
                      core_steps=60, mantle_steps=160, water_steps=160):
    """
    Radius, central and core-mantle pressure, gravity, escape speed and C/(M R^2) of a body with an iron core, a silicate
    mantle and a water layer, by shooting on the central pressure (Illinois regula falsi).

    Args:
        mass_earth (float): Earth masses, `SOLVER_MASS_EARTH`.
        metal, rock, water (float): Mass fractions of the whole body; they are normalised to add up to 1 (any other ice is
            water-like and goes in `water`).
        x_fe (float): Mantle Fe/(Fe+Mg), molar.
        core_scale (float): The core's density factor from its light elements and its temperature.
        surface_pressure_pa (float): The level the radius is read at.

    Returns:
        dict: `radius` (Earth radii), `central_pressure`, `cmb_pressure` (the top of the core, `None` without a core and a
        mantle) in Pa, `gravity` (m/s2), `escape` (m/s), `moment_of_inertia`.

    Raises:
        OutsideRange: The mass is outside `SOLVER_MASS_EARTH`, or the fractions are all zero.
    """
    if not SOLVER_MASS_EARTH[0] <= mass_earth <= SOLVER_MASS_EARTH[1]:
        raise OutsideRange(f"{mass_earth!r} Earth masses is outside the solver's {SOLVER_MASS_EARTH}")
    total = metal + rock + water
    if total <= 0.0:
        raise OutsideRange("the layer fractions are all zero")
    mass_kg = mass_earth * M_EARTH
    parts = [(metal / total, lambda p: core_density(p, core_scale), core_steps),
             (rock / total, lambda p: mantle_density(p, x_fe), mantle_steps),
             (water / total, water_density, water_steps)]
    parts = [part for part in parts if part[0] > 1e-9]
    layers = [(fraction, density_of) for fraction, density_of, _ in parts]
    steps = [count for _, _, count in parts]

    def residual(log_pc):
        return _integrate(math.exp(log_pc), mass_kg, layers, steps)[0][1] - surface_pressure_pa

    a, b = math.log(1.0e5), math.log(2.0e15)
    fa, fb = residual(a), residual(b)
    side = 0
    c = b
    for _ in range(80):
        c = (a * fb - b * fa) / (fb - fa)
        fc = residual(c)
        if abs(fc) < 2.0e3 or abs(b - a) < 1.0e-9:
            break
        if fc * fb > 0:
            b, fb = c, fc
            if side == -1:
                fa *= 0.5
            side = -1
        else:
            a, fa = c, fc
            if side == 1:
                fb *= 0.5
            side = 1
    (radius, _, moment), boundary = _integrate(math.exp(c), mass_kg, layers, steps)
    return {
        "radius": radius / R_EARTH,
        "central_pressure": math.exp(c),
        "cmb_pressure": boundary if len(layers) > 1 else None,
        "gravity": G * mass_kg / radius ** 2,
        "escape": math.sqrt(2.0 * G * mass_kg / radius),
        "moment_of_inertia": moment / (mass_kg * radius * radius),
    }


# ---------------------------------------------------------------------------------------------------------------------
# The fast forms
# ---------------------------------------------------------------------------------------------------------------------

_A = ((0.075176, -0.23494, -0.033795, -0.044414),
      (0.28771, -0.0099019, 0.0021907, -0.0049603),
      (-0.0052534, -0.0041485, 0.0081412, -0.009653),
      (-0.0011712, 0.0006952, -0.0026943, 0.0033964))
_B = (-0.11559, 0.086052, -0.0067087, 0.0040513)
_E = (0.25003, 0.18349, 0.022991, 0.013879)


def radius_fit_earth(mass_earth, cmf, x_fe=0.11, core_scale=0.9):
    """
    Part B's fast radius of a dry rocky planet in Earth radii: ln(R) = sum a_ij u^i c^j + x_fe (b0 + b1 c + b2 u + b3 u c) +
    (1 - core_scale) c (e0 + e1 c + e2 u + e3 u c), u = ln M, c = `cmf`, fitted to `layered_structure` on 1500 random draws.
    Valid 0.1 to 20 Earth masses, cmf 0 to 0.8, x_fe 0 to 0.3, core_scale 0.8 to 1 (`CLOSED_FORM_RANGES`).
    """
    u = math.log(mass_earth)
    s = sum(_A[i][j] * u ** i * cmf ** j for i in range(4) for j in range(4))
    s += x_fe * (_B[0] + _B[1] * cmf + _B[2] * u + _B[3] * u * cmf)
    s += (1.0 - core_scale) * cmf * (_E[0] + _E[1] * cmf + _E[2] * u + _E[3] * u * cmf)
    return math.exp(s)


def cmf_from_mass_radius(mass_earth, radius_earth, x_fe=0.11, core_scale=0.9, low=0.0, high=0.9):
    """
    The core mass fraction (of the dry part) that gives `radius_earth` at `mass_earth`, by bisection on `radius_fit_earth`
    (the radius falls as the core grows). `None` when the radius is outside the coreless and `high` limits.
    """
    def excess(cmf):
        return radius_fit_earth(mass_earth, cmf, x_fe, core_scale) - radius_earth

    if excess(low) < 0.0 or excess(high) > 0.0:
        return None
    for _ in range(50):
        mid = 0.5 * (low + high)
        if excess(mid) > 0.0:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


WATER_FIT = (1.0688, 0.2317, 0.4492, -0.1597)
"""tuple: Part C's closed form R = (a - b CMF (1 - x)) M^(1/3.7) (1 + c x + d x^2), a fit to its 175-point grid (fit_params2.json)."""
WATER_FIT_EXPONENT = 1.0 / 3.7


def water_world_radius_earth(mass_earth, cmf, water):
    """
    Part C's closed form of a body with water mass fraction `water` of the whole body and iron share `cmf` of the dry part,
    in Earth radii, cold, 1 bar. Valid 0.3 to 20 Earth masses; cold: a planet above the runaway-greenhouse limit holds its water
    as a steam atmosphere and is larger ([R4448], not coded).
    """
    a, b, c, d = WATER_FIT
    return (a - b * cmf * (1.0 - water)) * mass_earth ** WATER_FIT_EXPONENT * (1.0 + c * water + d * water * water)




SMALL_FIT_RANGE = {"mass_earth": (1.5e-4, 0.3), "cmf": (0.0, 0.6), "water": (0.0, 1.0)}
"""dict: What `SMALL_FIT` was fitted over (mantle Fe/(Fe+Mg) 0.11 and core density factor 0.9 are fixed in it)."""

SMALL_FIT = (
    (0, 0, 0, 0.049425130042842665),
    (0, 0, 1, 0.31363754824882184),
    (0, 0, 2, -0.01906252315657364),
    (0, 1, 0, -0.14752683174150732),
    (0, 1, 1, 0.3260477559995153),
    (0, 1, 2, 0.10002166443621706),
    (0, 2, 0, -0.06722043936874862),
    (0, 2, 1, 0.8051686034778492),
    (0, 2, 2, -1.6624387836019827),
    (1, 0, 0, 0.2686285620165081),
    (1, 0, 1, -0.13783424891100762),
    (1, 0, 2, 0.17418538396869027),
    (1, 1, 0, 0.019123959930845494),
    (1, 1, 1, 0.2732751713095798),
    (1, 1, 2, -0.18985619096388456),
    (1, 2, 0, 0.011947360531166129),
    (1, 2, 1, 0.08018234756937437),
    (1, 2, 2, -0.5981870670191283),
    (2, 0, 0, -0.010273379662690707),
    (2, 0, 1, -0.02358510505287554),
    (2, 0, 2, 0.03869118380625455),
    (2, 1, 0, 0.00303876227938904),
    (2, 1, 1, 0.08290739618144893),
    (2, 1, 2, -0.08190719662812826),
    (2, 2, 0, 0.00275121718531722),
    (2, 2, 1, -0.00895794989297094),
    (2, 2, 2, -0.07111833254864358),
    (3, 0, 0, -0.0005373091562476395),
    (3, 0, 1, -0.0014064529197319597),
    (3, 0, 2, 0.002612483666679266),
    (3, 1, 0, 0.00014152854430338202),
    (3, 1, 1, 0.006482801312019726),
    (3, 1, 2, -0.007428094208973751),
    (3, 2, 0, 0.00020551908567645327),
    (3, 2, 1, -0.0022821619375062163),
    (3, 2, 2, -0.0010074383086345895),
)
"""tuple: (i, j, k, a_ijk) of ln R/R_earth = sum a_ijk u^i c^j w^k, u = ln M/M_earth, c = metal share of the dry part, w = water
share: a fit to `layered_structure` on 600 random bodies, rms 0.57 and at most 2.0 percent on 150 fresh ones
(`scripts/fit_small_body_radius.py`). It gives a Mercury, a Mars or a Ganymede a radius without the 0.09 s solver."""


def small_body_radius_earth(mass_earth, cmf, water=0.0):
    """The fast radius, Earth radii, of a differentiated body of 1.5e-4 to 0.3 Earth masses (`SMALL_FIT_RANGE`)."""
    u = math.log(mass_earth)
    return math.exp(sum(a * u ** i * cmf ** j * water ** k for i, j, k, a in SMALL_FIT))


def core_radius_earth(core_mass_earth):
    """[R2741: Eq. 1]: the radius of a rock and iron core of 2:1 rock to iron, within about 2 percent."""
    return core_mass_earth ** 0.25


def envelope_thickness_earth(mass_earth, envelope, flux_earth, age_gy):
    """[R2741: Eq. 3]: the envelope's thickness in Earth radii, to the 20 mbar level."""
    return 2.06 * mass_earth ** -0.21 * (envelope / 0.05) ** 0.59 * flux_earth ** 0.044 * (age_gy / 5.0) ** -0.18


def envelope_radius_earth(mass_earth, envelope, flux_earth=1.0, age_gy=4.5, temperature_k=None, gravity_m_s2=None):
    """
    A rock and iron core of mass M(1 - envelope) under an H/He envelope of mass fraction `envelope`: R_core + R_env ([R2741] Eqs. 1 and 3), in Earth radii, at the 20 mbar level, or at the 1 bar level when `temperature_k` and
    `gravity_m_s2` are given (the radius moves by ln(50) scale heights H = R T / (mu g), part C and bulk-composition.md 4.1).
    """
    radius = core_radius_earth(mass_earth * (1.0 - envelope)) + envelope_thickness_earth(mass_earth, envelope, flux_earth, age_gy)
    if temperature_k is not None and gravity_m_s2 is not None:
        scale_height = constants.R * temperature_k / (ENVELOPE_MEAN_MOLAR_MASS * gravity_m_s2)
        radius -= math.log(50.0) * scale_height / R_EARTH
    return radius


def mixing_density_kg_m3(metal, rock, water, ice, porosity=0.0):
    """
    The zero-pressure mixing rule of a heap of grains: 1 / sum(w_i / rho_i) times (1 - porosity), with the grain densities
    of the materials table (iron, enstatite and water ice) and `ICE_OTHER_DENSITY_KG_M3` for the other ices.
    """
    grains = ((metal, materials.material("iron").prop("rho_grain_g_cm3").value * 1000.0),
              (rock, materials.material("enstatite").prop("rho_grain_g_cm3").value * 1000.0),
              (water, materials.material("water_ice").prop("rho_grain_g_cm3").value * 1000.0),
              (ice, ICE_OTHER_DENSITY_KG_M3))
    volume = sum(fraction / density for fraction, density in grains if fraction > 0.0)
    return (1.0 - (porosity or 0.0)) / volume


# ---------------------------------------------------------------------------------------------------------------------
# From fractions and mass to a structure
# ---------------------------------------------------------------------------------------------------------------------

def _inside(ranges, **values):
    return all(low <= values[name] <= high for name, (low, high) in ranges.items() if name in values)


def _finish(radius_km, mass_kg, method, moment_of_inertia=None, cmb_pressure_gpa=None):
    gravity = G * mass_kg / (radius_km * 1000.0) ** 2
    return Structure(
        radius_km=radius_km,
        density_g_cm3=mass_kg / (4.0 / 3.0 * math.pi * (radius_km * 1000.0) ** 3) / 1000.0,
        gravity_g=gravity / constants.EARTH_GRAVITY,
        escape_kms=math.sqrt(2.0 * G * mass_kg / (radius_km * 1000.0)) / 1000.0,
        moment_of_inertia=moment_of_inertia,
        cmb_pressure_gpa=cmb_pressure_gpa,
        method=method,
    )


def structure(bulk, mass_kg, *, flux_earth=1.0, age_gy=4.5, temperature_k=None, x_fe=None, core_scale=1.0,
              prefer_solver=False):
    """
    The radius, density, gravity and escape speed of a body of mass `mass_kg` and bulk composition `bulk`, at the 1 bar level.

    Args:
        bulk (BulkFractions): The composition. The mantle's Fe/(Fe+Mg) is `bulk.rock_feo` converted to moles unless `x_fe` is given;
            the core's density factor is `core_scale` (light elements and temperature, part B 2.4 and 3.3 step 6).
        mass_kg (float): The mass.
        flux_earth, age_gy (float): For an envelope (Lopez and Fortney).
        temperature_k (float or None): For the 1 bar level of an envelope; default is the equilibrium value at `flux_earth`.
        prefer_solver (bool): Use `layered_structure` even where a closed form applies (tests, edit views).

    Raises:
        OutsideRange: A giant, or a mass the solver cannot take.
    """
    if bulk.envelope > GIANT_ENVELOPE_FRACTION:
        raise OutsideRange(f"an envelope of {bulk.envelope!r} of the mass is a giant, outside this module (`physics/planets.py` has the giant radius law)")
    if mass_kg < SMALL_BODY_MASS_KG:
        density = mixing_density_kg_m3(bulk.metal, bulk.rock, bulk.water, bulk.ice, bulk.porosity)
        radius_km = (3.0 * mass_kg / (4.0 * math.pi * density)) ** (1.0 / 3.0) / 1000.0
        return _finish(radius_km, mass_kg, "mixing")

    mass_earth = mass_kg / M_EARTH
    dry = bulk.metal + bulk.rock
    cmf = bulk.metal / dry if dry > 0.0 else 0.0
    water_like = bulk.water + bulk.ice
    if water_like < WATER_NEGLIGIBLE:
        water_like = 0.0
    if x_fe is None:
        x_fe = x_fe_from_feo_mass_fraction(bulk.rock_feo) if bulk.rock_feo else 0.0

    if bulk.envelope > 0.0:
        core_mass_earth = mass_earth * (1.0 - bulk.envelope)
        if not _inside(CLOSED_FORM_RANGES["envelope"], mass_earth=mass_earth, envelope=bulk.envelope, flux_earth=flux_earth,
                       age_gy=age_gy) or core_mass_earth > CLOSED_FORM_RANGES["dry"]["mass_earth"][1]:
            raise OutsideRange(f"{mass_earth!r} Earth masses with an envelope of {bulk.envelope!r} is outside Lopez and Fortney's range")
        if temperature_k is None:
            temperature_k = 278.3 * flux_earth ** 0.25
        gravity = G * mass_kg / (core_radius_earth(core_mass_earth) * R_EARTH) ** 2
        radius = envelope_radius_earth(mass_earth, bulk.envelope, flux_earth, age_gy, temperature_k, gravity)
        return _finish(radius * R_EARTH / 1000.0, mass_kg, "closed-envelope")

    if not prefer_solver:
        if _inside(SMALL_FIT_RANGE, mass_earth=mass_earth, cmf=cmf, water=water_like) and (
                mass_earth < CLOSED_FORM_RANGES["dry"]["mass_earth"][0] or (water_like > 0.0 and mass_earth < CLOSED_FORM_RANGES["water"]["mass_earth"][0])):
            # a body is never less dense than its grains with no pores: the fit, good to a few percent, is held under that radius
            grains = mixing_density_kg_m3(bulk.metal, bulk.rock, bulk.water, bulk.ice, 0.0)
            grain_radius_km = (3.0 * mass_kg / (4.0 * math.pi * grains)) ** (1.0 / 3.0) / 1000.0
            return _finish(min(small_body_radius_earth(mass_earth, cmf, water_like) * R_EARTH / 1000.0, grain_radius_km), mass_kg,
                           "closed-small")
        if water_like <= 0.0 and _inside(CLOSED_FORM_RANGES["dry"], mass_earth=mass_earth, cmf=cmf, x_fe=x_fe,
                                         core_scale=core_scale):
            return _finish(radius_fit_earth(mass_earth, cmf, x_fe, core_scale) * R_EARTH / 1000.0, mass_kg, "closed-dry")
        if water_like > 0.0 and _inside(CLOSED_FORM_RANGES["water"], mass_earth=mass_earth, water=water_like, cmf=cmf):
            return _finish(water_world_radius_earth(mass_earth, cmf, water_like) * R_EARTH / 1000.0, mass_kg, "closed-water")

    result = layered_structure(mass_earth, bulk.metal, bulk.rock, water_like, x_fe, core_scale)
    return _finish(result["radius"] * R_EARTH / 1000.0, mass_kg, "solver", result["moment_of_inertia"],
                   None if result["cmb_pressure"] is None else result["cmb_pressure"] / 1.0e9)


def x_fe_from_feo_mass_fraction(feo, si_per_cation=1.0):
    """The mantle's Fe/(Fe+Mg), molar, from its FeO mass fraction for a mantle with Si/(Fe+Mg) = `si_per_cation`
    (bisection; molar masses from `physics.elements`, part B 3.8)."""
    from planetgen.physics import elements

    mu_feo = elements.molar_mass("FeO") * 1000.0
    mu_mgo = elements.molar_mass("MgO") * 1000.0
    mu_sio2 = elements.molar_mass("SiO2") * 1000.0
    low, high = 0.0, 0.99
    for _ in range(60):
        x = 0.5 * (low + high)
        share = mu_feo * x / (mu_feo * x + mu_mgo * (1.0 - x) + si_per_cation * mu_sio2)
        if share < feo:
            low = x
        else:
            high = x
    return 0.5 * (low + high)
