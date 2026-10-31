# planetgen/physics/tidal_heating.py

"""
Tidal heating from the orbit's real eccentricity (GEN.737)
==========================================================

A body on an eccentric orbit is flexed as its distance from what it goes round changes, and the flexing is dissipated as heat. For a
body that keeps one face to its host (a moon, or a close planet) Peale, Cassen and Reynolds 1979 give the power as

    E = (21/2) (k2/Q) G Mh^2 n R^5 e^2 / a^6

with `Mh` the host's mass, `n` the mean motion, `R` the body's radius, `a` the semi-major axis and `e` the eccentricity stored with the
body (GEN.727); `k2/Q` is its Love number over its tidal quality factor, the rheology. The eccentricity is the orbit's own, so the
heating follows the orbit: a moon with a forced eccentricity (the Laplace resonance of Io, Europa and Ganymede) keeps it, a planet that
has been circularised (`orbit_elements.circularised`) has little left. The power is taken as constant over the body's age (a stand-in:
the eccentricity was larger for a planet that has since circularised, so its past heating was more).

Rheology of an icy body (a water mass fraction of `hydrosphere.ICY_MIN_WATER_FRACTION` or more): the ice shell's own (`icy_tide`, GEN.960), in equilibrium with the shell
it makes; the heat goes into the shell (`hydrosphere`), over its surface.

Rheology of a rocky body (GEN.955): its Love number and quality factor are `tidal_rheology.rock_response` of its own mass and radius, the
tide's frequency, its mantle's temperature and the viscosity it drew (`interior_heat.Interior`), so a cold mantle is stiff and barely
dissipates (Q up to `tidal_rheology.COLD_Q_CEILING`), a warm one near the solidus is flexed close to its relaxation time and dissipates
most, and a partly molten one is weak (`k2` rises toward the fluid 3/2). The Moon comes out at the lunar laser ranging's Q of 38 to 60 and the
solid Earth's tide at a few hundred, with no number fitted to either. The heat of the tide is the Peale formula with that `k2/Q` at the
mantle's temperature, so it is a function of the temperature and the box (`interior_heat.heat`) integrates it with the cooling: heat and
temperature feed back, and the box is bistable (a warmer mantle is tidally weaker once its melt passes the peak).

The **heat pipe**. A mantle whose melt reaches `tidal_rheology.CRITICAL_MELT_FRACTION` (the crystals no longer touch: the mush is a
suspension and the melt segregates) erupts: the melt carries the heat out at the surface as fast as the interior makes it (the heat pipe of
Io, O'Reilly and Davies 1981, about 2 W/m2 for centimetres a year of resurfacing; Moore 2001). The mantle then stays at
`heat_pipe_temperature_k`, the temperature of that melt fraction, and the heat it radiates is the radiogenic heat plus the tide it is flexed with
there. The pipe takes over where the radiogenic heat plus the tide at that temperature is more than the mantle's conduction can carry away
(the box would otherwise pass it), so no tidal power is chosen: the melt state and the fitted `k2/Q` of Io are no longer numbers of ours.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, hydrosphere, tidal_rheology

ROCK_RIGIDITY_PA = 7.0e10
"""float: The shear modulus of the solid mantle, Pa (PREM mantle average: 65 to 80 GPa). Source tier T2."""
GIANT_TIDAL_Q = 1.0e5
"""float: A gas giant's tidal quality factor (Jupiter's 3e4 to 1e6, Saturn's 1.8e4 from the motions of their moons). Source tier T3, tweak L0."""
GIANT_K2 = 0.3
"""float: And its Love number (Jupiter 0.49, Saturn 0.39, Neptune 0.13, Uranus 0.1). Source tier T3, tweak level L0."""
SOLIDUS_SURFACE_K = 1393.0
"""float: The peridotite solidus at the surface, K (Hirschmann 2000, recalled). Source tier T2, tweak level L0."""
SOLIDUS_PRESSURE_SCALE_GPA = 2.2
"""float: The pressure scale of the Simon-Glatzel curve `T_s = T_0 (1 + p / a)^b` that carries the solidus up from the surface to the core-mantle
boundary. Source tier T3, tweak level L2: fitted through two recalled points, 2,200 K at 10 GPa (Hirschmann 2000) and 4,180 K at 135 GPa (Fiquet et
al. 2010), where it gives Hirschmann's 1,750 K at 3 GPa without being fitted to it."""
SOLIDUS_EXPONENT = 0.267
"""float: And its exponent (see `SOLIDUS_PRESSURE_SCALE_GPA`). Source tier T3, tweak level L2."""
MELT_RANGE_K = 450.0
"""float: How far a mantle's liquidus lies above its solidus, K (a peridotite liquidus 400 to 500 K over its solidus, Hirschmann 2000, recalled).
Source tier T2, tweak level L0."""


def solidus_k(pressure_gpa):
    """The mantle's solidus at `pressure_gpa`, K: `SOLIDUS_SURFACE_K (1 + p / a)^b`, 1,750 K at 3 GPa, about 3,400 K at the middle of an Earth's mantle."""
    return SOLIDUS_SURFACE_K * (1.0 + pressure_gpa / SOLIDUS_PRESSURE_SCALE_GPA) ** SOLIDUS_EXPONENT


def liquidus_k(pressure_gpa):
    """The melt's liquidus there: all melt (`MELT_RANGE_K` over the solidus)."""
    return solidus_k(pressure_gpa) + MELT_RANGE_K


def heat_pipe_temperature_k(pressure_gpa):
    """The mantle temperature a melt-limited body holds, K: that of the melt fraction `tidal_rheology.CRITICAL_MELT_FRACTION` at which the crystals stop
    touching and the melt segregates and erupts (the heat pipe of Moore 2001). Derived from TR2, which replaces the chosen melt share of
    0.3 (TH7)."""
    return solidus_k(pressure_gpa) + tidal_rheology.CRITICAL_MELT_FRACTION * MELT_RANGE_K


@dataclass(frozen=True)
class Forcing:
    """
    What flexes a body, the tide-raising side of `E`.

    Attributes:
        host_mass_kg (float): The mass of what it orbits (its star for a planet, its planet for a moon).
        semi_major_axis_m (float): Its semi-major axis.
        eccentricity (float): Its eccentricity.
        moons (tuple): A planet's moons, `(mass_kg, distance_m, retrograde, orbital_rate_rad_s)` each (the rate signed, negative for a retrograde moon, from the planet's whole mass):
            their tides, at the planet's spin, flex it too (GEN.959). Empty for a moon.
        spin_rate_rad_s (float): The planet's spin, rad/s, negative for a retrograde one (what the moons' tides are raised against).
    """

    host_mass_kg: float
    semi_major_axis_m: float
    eccentricity: float
    moons: tuple = ()
    spin_rate_rad_s: float = 0.0


def moon_tides_of(body):
    """The `(mass_kg, distance_m, retrograde, orbital_rate_rad_s)` of each of a planet's moons and the planet's spin rate (rad/s, negative when retrograde) from the stored
    fields, or `((), 0.0)` for a moon, a planet with no moons or one whose day is not stored."""
    moons = getattr(body, "moons", None)
    day_hours = getattr(body, "rotation_period_hours", None)
    if not moons or getattr(body, "is_moon", False) or not day_hours:
        return (), 0.0
    tilt = getattr(body, "axial_tilt_deg", None) or 0.0
    spin = 2.0 * math.pi / (day_hours * 3600.0) * (-1.0 if tilt > 90.0 else 1.0)
    tides = []
    for moon in moons:
        distance_m = moon.distance * constants.AU_TO_M
        retrograde = (getattr(moon, "orbital_inclination_deg", 0.0) or 0.0) > 90.0
        rate = math.sqrt(constants.G * (body.mass + moon.mass) / distance_m ** 3) * (-1.0 if retrograde else 1.0)
        tides.append((moon.mass, distance_m, retrograde, rate))
    tides = tuple(tides)
    return tides, spin


def forcing_of(body):
    """The `Forcing` of a planet or moon from its stored orbit and moons, or `None` for a circle with no moon tide, a body with no host or no orbit."""
    host_mass = getattr(body, "primary_mass_kg", None)
    eccentricity = getattr(body, "eccentricity", None)
    distance_au = getattr(body, "distance", None)
    if not host_mass or not distance_au:
        return None
    moons, spin = moon_tides_of(body)
    if not moons and (not eccentricity or eccentricity <= 0.0):
        return None
    return Forcing(host_mass, distance_au * constants.AU_TO_M, eccentricity if eccentricity and eccentricity > 0.0 else 0.0, moons, spin)


def is_icy(body):
    """Whether `body` is an icy one (a water mass fraction of `hydrosphere.ICY_MIN_WATER_FRACTION` or more)."""
    return (getattr(body, "water_mass_fraction", None) or 0.0) >= hydrosphere.ICY_MIN_WATER_FRACTION


def mean_motion_rad_s(forcing):
    """The orbital rate of a body on `forcing`'s orbit, rad/s: the frequency of its eccentricity tide."""
    return math.sqrt(constants.G * forcing.host_mass_kg / forcing.semi_major_axis_m ** 3)


def peale_power_w(k2_over_q, forcing, radius_m):
    """The eccentricity tide's power, W, in a body of radius `radius_m` with the given `k2/Q` (the formula in the module's docstring)."""
    n = mean_motion_rad_s(forcing)
    return (10.5 * k2_over_q * constants.G * forcing.host_mass_kg ** 2 * n * radius_m ** 5
            * forcing.eccentricity ** 2 / forcing.semi_major_axis_m ** 6)


def rigid_k2(mass_kg, radius_m):
    """`k2` of a homogeneous, incompressible sphere of rock with rigidity `ROCK_RIGIDITY_PA`: `(3/2) / (1 + 19 mu / (2 rho g R))`, `rho g R`
    being `3 G M^2 / (4 pi R^4)`."""
    pressure = 3.0 * constants.G * mass_kg ** 2 / (4.0 * math.pi * radius_m ** 4)
    return 1.5 / (1.0 + 9.5 * ROCK_RIGIDITY_PA / pressure)
