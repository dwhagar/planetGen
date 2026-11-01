# planetgen/physics/king.py

"""
[R2359] models for globular clusters (GEN.159).

A King model is self-similar: for one central potential `W0` the density, the cumulative mass and the tidal radius
are one curve in units of the King radius `r0`, so ten tables (`W0` = 3 to 12, `tuning.KING_W0_RANGE`), each solved
once and resampled to `tuning.KING_TABLE_POINTS` log-spaced radii, serve every cluster: a cluster of mass `M` and
half-mass radius `r_h` scales its table by `r0 = r_h / (r_h in units of r0)`. A star's distance from the centre is
then one lookup of a random mass fraction (`radius_pc`, an inverse transform on a monotone array).

The solver integrates `W'' + 2 W'/x = -9 rho(W) / rho(W0)` outward from the centre (`W` falls from `W0` to 0 at the
tidal radius `x_t`), with `rho(W) = exp(W) erf(sqrt W) - sqrt(4 W / pi) (1 + 2 W / 3)`. It reproduces the textbook
concentrations `c = log10(x_t)` of 0.67, 1.03, 1.53, 2.12 and 2.74 for `W0` = 3, 5, 7, 9 and 12.
"""

import dataclasses
import functools
import math

import numpy as np
from scipy.integrate import solve_ivp

from planetgen import tuning

_X_FIRST = 1e-3
"""float: The first radius of a table, in King radii; the centre (mass 0) is added before it."""


def king_density(w):
    """The King density `rho(W)` for a dimensionless potential `W` (0 at and beyond the tidal edge)."""
    if w <= 0.0:
        return 0.0
    return math.exp(w) * math.erf(math.sqrt(w)) - math.sqrt(4.0 * w / math.pi) * (1.0 + 2.0 * w / 3.0)


@dataclasses.dataclass(frozen=True)
class KingTable:
    """One King model in King radii.

    Attributes:
        w0 (int): The central potential.
        concentration (float): `log10(r_t / r0)`.
        tidal_x (float): The tidal radius over `r0`.
        half_mass_x (float): The radius holding half the mass, over `r0`.
        x (numpy.ndarray): The radii over `r0`, from 0, log-spaced after the first.
        mass_fraction (numpy.ndarray): The mass inside each radius over the whole, rising from 0 to 1.
        density (numpy.ndarray): The density at each radius over the central density.
        mass_over_central (float): The whole mass over `4 pi r0^3 rho0`: `integral of x^2 rho / rho0`.
    """

    w0: int
    concentration: float
    tidal_x: float
    half_mass_x: float
    x: np.ndarray
    mass_fraction: np.ndarray
    density: np.ndarray
    mass_over_central: float

    def radius_x(self, fraction):
        """The radius, over `r0`, inside which `fraction` (0 to 1) of the mass lies."""
        return float(np.interp(fraction, self.mass_fraction, self.x))

    def mass_within_x(self, x):
        """The share of the mass inside `x` (over `r0`; 1 at and beyond the tidal radius)."""
        return float(np.interp(x, self.x, self.mass_fraction))

    def density_at_x(self, x):
        """The density at `x` over the central density (0 beyond the tidal radius)."""
        if x >= self.tidal_x:
            return 0.0
        return float(np.interp(x, self.x, self.density))


@functools.lru_cache(maxsize=None)
def king_table(w0):
    """The solved table of the King model with central potential `w0` (a whole number in `tuning.KING_W0_RANGE`)."""
    low, high = tuning.KING_W0_RANGE
    if w0 != int(w0) or not low <= w0 <= high:
        raise ValueError(f"a King table exists for W0 {low} to {high}, not {w0!r}")
    w0 = int(w0)
    rho0 = king_density(w0)
    first = 1e-6

    def rates(x, y):
        w, slope, _mass = y
        return [slope, -9.0 * king_density(w) / rho0 - 2.0 * slope / x, x * x * king_density(w) / rho0]

    def edge(_x, y):
        return y[0]

    edge.terminal = True
    edge.direction = -1
    # Near the centre W = W0 - 1.5 x^2.
    solution = solve_ivp(rates, [first, 1e5], [w0 - 1.5 * first ** 2, -3.0 * first, 0.0], events=edge,
                         dense_output=True, rtol=1e-10, atol=1e-12)
    tidal = float(solution.t_events[0][0])
    grid = np.concatenate(([0.0], np.geomspace(_X_FIRST, tidal, tuning.KING_TABLE_POINTS - 1)))
    inside = grid[1:]
    values = solution.sol(inside)
    whole = float(solution.sol(tidal)[2])
    mass = np.concatenate(([0.0], values[2] / whole))
    mass[-1] = 1.0
    density = np.concatenate(([1.0], [king_density(w) / rho0 for w in values[0]]))
    half = float(np.interp(0.5, mass, grid))
    return KingTable(w0=w0, concentration=math.log10(tidal), tidal_x=tidal, half_mass_x=half, x=grid,
                     mass_fraction=mass, density=density, mass_over_central=whole)


def scale_radius_pc(w0, half_mass_radius_pc):
    """The King radius `r0` (pc) of a cluster of half-mass radius `half_mass_radius_pc` and central potential `w0`."""
    return half_mass_radius_pc / king_table(w0).half_mass_x


def tidal_radius_pc(w0, half_mass_radius_pc):
    """The King model's own tidal radius (pc) of such a cluster."""
    return scale_radius_pc(w0, half_mass_radius_pc) * king_table(w0).tidal_x


def radius_pc(w0, half_mass_radius_pc, fraction):
    """The distance from the centre (pc) inside which `fraction` of a cluster's mass lies, by inverse transform:
    a random fraction gives a star's radius."""
    return scale_radius_pc(w0, half_mass_radius_pc) * king_table(w0).radius_x(fraction)


def mass_within_pc(w0, half_mass_radius_pc, r_pc):
    """The share of a cluster's mass within `r_pc` of its centre."""
    return king_table(w0).mass_within_x(r_pc / scale_radius_pc(w0, half_mass_radius_pc))


def central_density_per_pc3(w0, mass, half_mass_radius_pc):
    """The density at the centre, in the units of `mass` per cubic parsec."""
    r0 = scale_radius_pc(w0, half_mass_radius_pc)
    return mass / (4.0 * math.pi * r0 ** 3 * king_table(w0).mass_over_central)


def density_per_pc3(w0, mass, half_mass_radius_pc, r_pc):
    """The density at `r_pc` from the centre, in the units of `mass` per cubic parsec (0 beyond the tidal radius)."""
    table = king_table(w0)
    r0 = scale_radius_pc(w0, half_mass_radius_pc)
    return central_density_per_pc3(w0, mass, half_mass_radius_pc) * table.density_at_x(r_pc / r0)


def expected_systems(mass_solar, mean_system_mass_solar):
    """How many systems a cluster of this mass holds: its mass over the mean mass of its systems (the caller's, from the
    cluster's own age: `galaxy.components.Component.mean_system_mass_solar`)."""
    return mass_solar / mean_system_mass_solar


def jacobi_radius_pc(mass_solar, galactocentric_pc, galaxy_mass_solar):
    """The largest radius a cluster of `mass_solar` can keep at `galactocentric_pc` from the galactic centre (use the
    pericentre of an eccentric orbit): `R (M / (2 M_gal(<R)))^(1/3)`, where `galaxy_mass_solar` is the galaxy's own mass
    inside `R` (`galaxy.enclosed_mass.enclosed_mass_solar`, GEN.339: it was an isothermal halo of the Milky Way's 220 km/s)."""
    return galactocentric_pc * (mass_solar / (2.0 * galaxy_mass_solar)) ** (1.0 / 3.0)


def limit_w0(w0, half_mass_radius_pc, jacobi_pc):
    """The central potential to give a cluster that wants `w0`: lowered, a step at a time, until its King tidal radius
    fits inside `jacobi_pc`, but never below `tuning.KING_W0_RANGE`'s low end."""
    low, _high = tuning.KING_W0_RANGE
    while w0 > low and tidal_radius_pc(w0, half_mass_radius_pc) > jacobi_pc:
        w0 -= 1
    return w0


def nearest_w0(concentration):
    """The `W0` whose concentration `c = log10(r_t / r0)` is nearest `concentration` (a catalogue gives `c`)."""
    low, high = tuning.KING_W0_RANGE
    return min(range(low, high + 1), key=lambda w0: abs(king_table(w0).concentration - concentration))
