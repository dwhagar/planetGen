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

Rheology of an icy body (a water mass fraction of `hydrosphere.ICY_MIN_WATER_FRACTION` or more): one `k2/Q` of 3e-4, the figure that
gives Europa's 1.4e11 W (`ICY_K2_OVER_Q`); the heat goes into its ice shell (`hydrosphere`), over its surface.

Rheology of a rocky body, in two states:

* Solid: `k2` from the rigid-body formula for a homogeneous sphere (Love 1911; `rigid_k2`), `k2 = (3/2) / (1 + 19 mu / (2 rho g R))` with the
  rigidity `mu` of rock, which gives the Moon's 0.021 (observed 0.024); `Q` of 100 (`constants.MOON_TIDAL_DISSIPATION_Q`, Goldreich and Soter
  1966, the same used to circularise it). The Moon comes out at 2e8 W (about 6e8 W at the Moon's own Q of 38).
* Melt-limited: when that solid-state heat is at least the body's own heat flow and would warm its mantle (the one-box model of
  `interior_heat`) to the temperature at which a share `HEAT_PIPE_MELT_FRACTION` of it is melt, the melt erupts and carries the heat out (the
  heat pipe of Io, O'Reilly and Davies 1981: a resurfacing of centimetres a year takes about 2 W/m2). The mantle stays at that temperature
  (`HEAT_PIPE_TEMPERATURE_K`) and the partly molten body flexes far more, `MOLTEN_K2_OVER_Q` (Io's: about 1e14 W at its 0.0041
  eccentricity, Peale et al. 1979, Moore 2003). No solid-state model gets there, because the heat leaves by volcanism, not conduction,
  which a one-box conduction model does not hold: this state and its two numbers are stand-ins (see `docs/design/research/gap-filling-plan.md`
  section 6.12).
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, hydrosphere

ROCK_RIGIDITY_PA = 7.0e10
"""float: The shear modulus of the solid mantle, Pa (PREM mantle average: 65 to 80 GPa). Source tier T2."""
ROCK_TIDAL_Q = constants.MOON_TIDAL_DISSIPATION_Q
"""float: The tidal quality factor of a solid rocky body (the same one that circularises it)."""
ICY_K2 = 0.03
"""float: An icy body's Love number; with `ROCK_TIDAL_Q` it gives `ICY_K2_OVER_Q`."""
ICY_K2_OVER_Q = ICY_K2 / ROCK_TIDAL_Q
"""float: An icy body's `k2/Q` of 3e-4: Europa's 1.4e11 W at its eccentricity of 0.009 (Peale, Cassen and Reynolds 1979 for Io, applied to the
ice shell by Ojakangas and Stevenson 1989). Source tier T3, tweak level L2 (tuned to the published figure)."""
GIANT_TIDAL_Q = 1.0e5
"""float: A gas giant's tidal quality factor (Jupiter's 3e4 to 1e6, Saturn's 1.8e4 from the motions of their moons). Source tier T3, tweak L0."""
GIANT_K2 = 0.3
"""float: And its Love number (Jupiter 0.49, Saturn 0.39, Neptune 0.13, Uranus 0.1). Source tier T3, tweak level L0."""
MOLTEN_K2_OVER_Q = 0.015
"""float: A partly molten rocky body's `k2/Q`: the value that gives Io the 1e14 W it radiates (Peale et al. 1979; Io's observed 0.9 to 1.5e14 W
over `k2 = 1.5`, `Q = 100`). Source tier T8 (fitted to the one body that shows it), tweak level L4."""
SOLIDUS_K = 1700.0
"""float: The mantle's temperature at which melt begins in the one-box model, K (the peridotite solidus, 1394 K at the surface and 1750 K at
3 GPa, Hirschmann 2000). Source tier T2."""
LIQUIDUS_K = 2200.0
"""float: And at which it is all melt, K (a peridotite liquidus 400 to 500 K over its solidus). Source tier T2."""
HEAT_PIPE_MELT_FRACTION = 0.3
"""float: The melt share at which the melt segregates and erupts fast enough to carry the heat out (the rheological transition lies at 0.3 to 0.4,
Arzi 1978; Moore 2003 puts Io's mantle near 0.2 to 0.3). Source tier T3, tweak level L3."""
HEAT_PIPE_TEMPERATURE_K = SOLIDUS_K + HEAT_PIPE_MELT_FRACTION * (LIQUIDUS_K - SOLIDUS_K)
"""float: The mantle temperature a melt-limited body holds, K."""


@dataclass(frozen=True)
class Forcing:
    """
    What a body orbits, the tide-raising side of `E`.

    Attributes:
        host_mass_kg (float): The mass of what it orbits (its star for a planet, its planet for a moon).
        semi_major_axis_m (float): Its semi-major axis.
        eccentricity (float): Its eccentricity.
    """

    host_mass_kg: float
    semi_major_axis_m: float
    eccentricity: float


def forcing_of(body):
    """The `Forcing` of a planet or moon from its stored orbit, or `None` for a circle, a body with no host or no orbit."""
    host_mass = getattr(body, "primary_mass_kg", None)
    eccentricity = getattr(body, "eccentricity", None)
    distance_au = getattr(body, "distance", None)
    if not host_mass or not eccentricity or not distance_au or eccentricity <= 0.0:
        return None
    return Forcing(host_mass, distance_au * constants.AU_TO_M, eccentricity)


def is_icy(body):
    """Whether `body` is an icy one (a water mass fraction of `hydrosphere.ICY_MIN_WATER_FRACTION` or more)."""
    return (getattr(body, "water_mass_fraction", None) or 0.0) >= hydrosphere.ICY_MIN_WATER_FRACTION


def response(body):
    """`(k2, Q)` of a body, its tidal Love number and quality factor: a gas giant's, an icy body's or, for rock, `rigid_k2` of its own mass and
    radius with `ROCK_TIDAL_Q` (the solid state; the melt-limited state of GEN.737 is a heat-flow state, not one of the spin's)."""
    if getattr(body, "body_type", None) == "g":
        return GIANT_K2, GIANT_TIDAL_Q
    if is_icy(body):
        return ICY_K2, ROCK_TIDAL_Q
    return rigid_k2(body.mass, body.radius * 1000.0), ROCK_TIDAL_Q


def peale_power_w(k2_over_q, forcing, radius_m):
    """The eccentricity tide's power, W, in a body of radius `radius_m` with the given `k2/Q` (the formula in the module's docstring)."""
    n = math.sqrt(constants.G * forcing.host_mass_kg / forcing.semi_major_axis_m ** 3)
    return (10.5 * k2_over_q * constants.G * forcing.host_mass_kg ** 2 * n * radius_m ** 5
            * forcing.eccentricity ** 2 / forcing.semi_major_axis_m ** 6)


def rigid_k2(mass_kg, radius_m):
    """`k2` of a homogeneous, incompressible sphere of rock with rigidity `ROCK_RIGIDITY_PA`: `(3/2) / (1 + 19 mu / (2 rho g R))`, `rho g R`
    being `3 G M^2 / (4 pi R^4)`."""
    pressure = 3.0 * constants.G * mass_kg ** 2 / (4.0 * math.pi * radius_m ** 4)
    return 1.5 / (1.0 + 9.5 * ROCK_RIGIDITY_PA / pressure)


def solid_power_w(forcing, mass_kg, radius_m):
    """A solid rocky body's tide, W."""
    return peale_power_w(rigid_k2(mass_kg, radius_m) / ROCK_TIDAL_Q, forcing, radius_m)


def molten_power_w(forcing, radius_m):
    """A melt-limited rocky body's tide, W."""
    return peale_power_w(MOLTEN_K2_OVER_Q, forcing, radius_m)


def icy_flux_w_m2(forcing, radius_km):
    """An icy body's tide over its surface, W/m2."""
    radius_m = radius_km * 1000.0
    return peale_power_w(ICY_K2_OVER_Q, forcing, radius_m) / (4.0 * math.pi * radius_m ** 2)
