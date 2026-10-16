# planetgen/physics/component_profiles.py

"""
The density profiles of a galaxy's components (GEN.237), one table per profile, in the units the sector gate reads.

`galaxy_components` rows name one of three profiles (`galaxy.components.PROFILES`): `king` for a globular cluster
(`physics/king.py`, a table per central potential), and the two a satellite galaxy uses:

* `exponential`: density `exp(-r / a)`, for the Magellanic Clouds and the Sagittarius remnant. Half its mass lies inside
  `2.674 a` before the cut.
* `plummer`: density `(1 + (r / a)^2)^(-5/2)`, for the dwarf spheroidals. Half its mass lies inside `1.305 a`.

Both are cut off at a tidal radius, `tuning.SATELLITE_TIDAL_OVER_HALF_MASS` times the half-mass radius of the uncut
profile, and renormalised so the cut profile holds the whole mass (the cut loses 1e-5 of an exponential and 2% of a
Plummer sphere). Each table has the fields of `king.KingTable` that the gate uses (`x`, `density`, `mass_fraction`,
`tidal_x`, `half_mass_x`, `mass_over_central`), so `ComponentField` treats the three alike; this module is the one place
that turns a `Component` into its table, its scale radius and its central density.
"""

import dataclasses
import functools
import math

import numpy as np

from planetgen import tuning
from planetgen.physics import king

_X_FIRST = 1e-3

_UNCUT_HALF_MASS_X = {"exponential": 2.674, "plummer": 1.305}
"""dict: The radius, in scale radii, holding half the mass of the uncut profile (used to place the cut)."""


@dataclasses.dataclass(frozen=True)
class ProfileTable:
    """One non-King profile in scale radii (the fields `king.KingTable` has, less `w0`).

    Attributes:
        profile (str): `exponential` or `plummer`.
        concentration (float): `log10(tidal_x)`.
        tidal_x (float): The cut radius over the scale radius.
        half_mass_x (float): The radius holding half the (cut) mass, over the scale radius.
        x (numpy.ndarray): The radii over the scale radius, from 0, log-spaced after the first.
        mass_fraction (numpy.ndarray): The mass inside each radius over the whole, rising from 0 to 1.
        density (numpy.ndarray): The density at each radius over the central density.
        mass_over_central (float): The whole mass over `4 pi a^3 rho0`.
    """

    profile: str
    concentration: float
    tidal_x: float
    half_mass_x: float
    x: np.ndarray
    mass_fraction: np.ndarray
    density: np.ndarray
    mass_over_central: float

    def radius_x(self, fraction):
        """The radius, over the scale radius, inside which `fraction` (0 to 1) of the mass lies."""
        return float(np.interp(fraction, self.mass_fraction, self.x))

    def mass_within_x(self, x):
        """The share of the mass inside `x` (1 at and beyond the cut)."""
        return float(np.interp(x, self.x, self.mass_fraction))

    def density_at_x(self, x):
        """The density at `x` over the central density (0 beyond the cut)."""
        if x >= self.tidal_x:
            return 0.0
        return float(np.interp(x, self.x, self.density))


def _density_and_mass(profile, x):
    """`(density over the central density, the integral of x^2 density from 0 to x)`: the mass inside `x` over
    `4 pi a^3 rho0`, in closed form."""
    if profile == "exponential":
        return np.exp(-x), 2.0 - np.exp(-x) * (x * x + 2.0 * x + 2.0)
    return (1.0 + x * x) ** -2.5, x ** 3 / (3.0 * (1.0 + x * x) ** 1.5)


@functools.lru_cache(maxsize=None)
def profile_table(profile):
    """The table of `profile` (`exponential` or `plummer`)."""
    if profile not in _UNCUT_HALF_MASS_X:
        raise ValueError(f"no table for profile {profile!r}; expected one of {sorted(_UNCUT_HALF_MASS_X)}")
    tidal = tuning.SATELLITE_TIDAL_OVER_HALF_MASS * _UNCUT_HALF_MASS_X[profile]
    grid = np.concatenate(([0.0], np.geomspace(_X_FIRST, tidal, tuning.KING_TABLE_POINTS - 1)))
    density, mass = _density_and_mass(profile, grid)
    whole = float(mass[-1])
    mass_fraction = mass / whole
    mass_fraction[-1] = 1.0
    half = float(np.interp(0.5, mass_fraction, grid))
    return ProfileTable(profile=profile, concentration=math.log10(tidal), tidal_x=tidal, half_mass_x=half, x=grid,
                        mass_fraction=mass_fraction, density=density, mass_over_central=whole)


def table_for(component):
    """The table of `component`'s profile: the King table of its central potential, or the exponential or Plummer
    table."""
    if component.profile == "king":
        return king.king_table(int(component.king_w0))
    return profile_table(component.profile)


def scale_radius_pc(component):
    """The component's scale radius (pc): its half-mass radius over the table's half-mass radius in scale radii."""
    return component.half_mass_radius_pc / table_for(component).half_mass_x


def cut_radius_pc(profile, half_mass_radius_pc, king_w0=None):
    """The tidal radius (pc) a component of this profile and half-mass radius has: the King model's own, or the cut
    of an exponential or Plummer profile."""
    if profile == "king":
        return king.tidal_radius_pc(int(king_w0), half_mass_radius_pc)
    table = profile_table(profile)
    return half_mass_radius_pc / table.half_mass_x * table.tidal_x


def central_density_per_pc3(component):
    """The density at the centre, in the units of the component's mass per cubic parsec."""
    r0 = scale_radius_pc(component)
    return component.mass_solar / (4.0 * math.pi * r0 ** 3 * table_for(component).mass_over_central)
