# planetgen/physics/neutron_star_structure.py

"""
The structure of a neutron star from an equation of state (GEN.968, parent GEN.743).

The Tolman-Oppenheimer-Volkoff equations of a static, spherical star in general relativity, integrated over the
twenty-one tabulated cold-matter equations of state of `datasets.neutron_star_matter` [R2834]:

    dP/dr = -(e + P) (m + 4 pi r^3 P) / (r (r - 2 m)),        dm/dr = 4 pi r^2 e,

in geometrised units (G = c = 1) with lengths in kilometres, so that a pressure and an energy density are in km^-2 and a
mass is the length G M / c^2. The integration runs in the logarithm of the pressure, from the centre (a series start a
fraction of a part per million off the central pressure) to the lowest pressure of the table, which is the surface, and
the equation of state between two rows is the PCHIP interpolant of log e against log P and of log n against log P (no
overshoot, so the sound speed stays what the table says). A star is fixed by its central pressure; the stable branch of
an equation of state is the central pressures up to the one that holds the most mass (`maximum_mass`), and above it the
star collapses to a black hole.

What the module answers by integration, which needs the downloaded tables and scipy: a star's radius, mass and layers (`structure`),
the mass-radius curve of one equation of state (`integrated_curve`) and the heaviest star it holds (`integrated_maximum_mass_solar`).
What generation reads is the reduced table those make (`scripts/build_neutron_star_table.py` writes `data/neutron_star_curves.json`: 25
stars of each equation of state): the largest mass an equation of state supports (`maximum_mass_solar`), the radius and the crust of a
star of a given mass (`radius_km`, `crust_thickness_km`), the equations of state that hold a star (`supporting`) and the band the
twenty-one span at a mass (`radius_band_km`). `physics.neutron_star` (GEN.969) draws a star's radius from them, weighted by the [R3075] constraints; the crust thickness (`crust_thickness_km`) is not read by generation yet.

Layers. The crust ends at the density where the matter is uniform (crust-core transition, `CRUST_CORE_DENSITY_FM3`,
a stand-in: the tables do not mark it, the literature range is 0.07 to 0.10 fm^-3) and its outer part ends at neutron
drip (`NEUTRON_DRIP_DENSITY_FM3`, 2.6e-4 fm^-3, the baryon density the BSk22 outer crust of every table ends at).
"""

import functools
import json
import math
import os
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar

from planetgen.datasets import neutron_star_matter
from planetgen.physics import constants, unit_registry as units

_msun_to_kg = units.converter("Msun", "kg")

MEV_FM3_PER_KM2 = 1.0 / (1.602176634e-13 * 1e45 * constants.G / constants.SPEED_OF_LIGHT_M_S ** 4 * 1e6)
"""float: One km^-2 (the geometrised pressure) in MeV fm^-3. G/c^4 turns a J/m^3 into m^-2."""

SOLAR_MASS_KM = _msun_to_kg(constants.G) / constants.SPEED_OF_LIGHT_M_S ** 2 / 1e3
"""float: The length G M / c^2 of a solar mass, km (1.4766)."""

CRUST_CORE_DENSITY_FM3 = 0.08
"""float: The baryon density, fm^-3, at the bottom of the crust. [Stand-in, GEN.725 NS1: the tables do not mark the
transition; 0.07 to 0.10 is the range of the compressible-liquid-drop and RMF calculations, [R2834: Table 1] and
their references.]"""

NEUTRON_DRIP_DENSITY_FM3 = 2.6e-4
"""float: The baryon density, fm^-3, at neutron drip, the bottom of the outer crust of the BSk22 crust all tables share
([R0722], Living Rev. Relativ. 11, 10)."""

_START_FRACTION = 1e-7
_RESAMPLE_POINTS = 6000
_RTOL = 1e-6
"""float: The integrator's relative tolerance; the radius and mass of a star agree with a tolerance of 1e-10 to about 1e-5."""
_CENTRAL_PRESSURE_STEPS = 24
_CURVE_MASSES = 25
_CURVE_START_DENSITY_FM3 = 0.2
"""float: The central baryon density, fm^-3, where a tabulated curve starts (0.15 to 0.5 solar masses). Below it the stable
branch continues down to stars of a tenth of a solar mass and tens of kilometres, which no collapse forms."""


@dataclass(frozen=True)
class Structure:
    """One star: its fixed central pressure and what the integration gave."""

    eos_index: int
    central_pressure_mev_fm3: float
    central_baryon_density_fm3: float
    mass_solar: float
    radius_km: float
    crust_core_radius_km: float
    drip_radius_km: float
    crust_mass_solar: float
    """float: The mass of the layer above the crust-core transition; zero when the star's centre is below it."""

    @property
    def crust_thickness_km(self):
        """The thickness of the crust, km."""
        return self.radius_km - self.crust_core_radius_km

    @property
    def outer_crust_thickness_km(self):
        """The thickness of the layer above neutron drip, km."""
        return self.radius_km - self.drip_radius_km

    @property
    def compactness(self):
        """G M / (R c^2), dimensionless (0.5 is a black hole, 4/9 the Buchdahl limit of a star of constant density)."""
        return self.mass_solar * SOLAR_MASS_KM / self.radius_km


class _Eos:
    """One equation of state as interpolants of the logarithms, in geometrised units."""

    def __init__(self, table):
        self.index = table.index
        pressure = np.log(table.pressure_mev_fm3 / MEV_FM3_PER_KM2)
        if np.any(np.diff(pressure) <= 0):
            raise ValueError(f"EOS{table.index}: the pressure must rise with the density")
        self.log_p_min, self.log_p_max = pressure[0], pressure[-1]
        energy = PchipInterpolator(pressure, np.log(table.energy_density_mev_fm3 / MEV_FM3_PER_KM2))
        density = PchipInterpolator(pressure, np.log(table.baryon_density_fm3))
        self._log_p_of_density = PchipInterpolator(np.log(table.baryon_density_fm3), pressure)
        self.density_max_fm3 = float(table.baryon_density_fm3[-1])
        # the interpolants resampled finely: the integrator calls them tens of thousands of times a star
        self._grid = np.linspace(self.log_p_min, self.log_p_max, _RESAMPLE_POINTS)
        self._log_energy = energy(self._grid)
        self._log_density = density(self._grid)

    def energy(self, log_p):
        """The energy density, km^-2, at ln(pressure in km^-2)."""
        return math.exp(np.interp(log_p, self._grid, self._log_energy))

    def density_fm3(self, log_p):
        """The baryon density, fm^-3."""
        return math.exp(np.interp(log_p, self._grid, self._log_density))

    def log_pressure(self, density_fm3):
        """ln(pressure, km^-2) at a baryon density in fm^-3."""
        return float(self._log_p_of_density(math.log(density_fm3)))


@functools.lru_cache(maxsize=1)
def _eos_set():
    return {table.index: _Eos(table) for table in neutron_star_matter.equations_of_state()}


def _eos(eos):
    """The `_Eos` of an equation-of-state number (1 to 21); an `_Eos` stands for itself."""
    if isinstance(eos, _Eos):
        return eos
    try:
        return _eos_set()[eos]
    except KeyError:
        raise ValueError(f"equation of state {eos!r}: the compilation holds 1 to {neutron_star_matter.EOS_COUNT}") from None


def _integrate(eos, log_p_centre, dense=False):
    """Integrates the TOV equations in ln P from the centre; returns the solver's solution over (r, m), with its
    interpolant `sol` when `dense`."""
    p_c = math.exp(log_p_centre)
    e_c = eos.energy(log_p_centre)
    r0 = math.sqrt(3.0 * _START_FRACTION / (2.0 * math.pi * (e_c + p_c) * (e_c + 3.0 * p_c)) * p_c)
    s0 = log_p_centre + math.log1p(-_START_FRACTION)
    m0 = 4.0 / 3.0 * math.pi * r0 ** 3 * e_c

    def rhs(s, y):
        r, m = y
        p = math.exp(s)
        e = eos.energy(s)
        slope = -p * r * (r - 2.0 * m) / ((e + p) * (m + 4.0 * math.pi * r ** 3 * p))
        return [slope, 4.0 * math.pi * r * r * e * slope]

    return solve_ivp(rhs, (s0, eos.log_p_min), [r0, m0], method="DOP853", rtol=_RTOL, atol=[1e-8, 1e-8],
                     dense_output=dense)


def structure(eos, central_pressure_mev_fm3):
    """The `Structure` of the star of equation of state `eos` (1 to 21) whose central pressure is given, MeV fm^-3."""
    eos = _eos(eos)
    log_p = math.log(central_pressure_mev_fm3 / MEV_FM3_PER_KM2)
    if not eos.log_p_min < log_p < eos.log_p_max:
        raise ValueError(f"EOS{eos.index}: a central pressure of {central_pressure_mev_fm3} MeV/fm^3 is outside the table")
    solution = _integrate(eos, log_p, dense=True)
    radius, mass = solution.y[0, -1], solution.y[1, -1]

    def reaches(density_fm3):
        """The radius and enclosed mass where the density falls to `density_fm3`, or None if the star's centre is below."""
        if eos.density_fm3(log_p) <= density_fm3:
            return None
        s = eos.log_pressure(density_fm3)
        r, m = solution.sol(s)
        return float(r), float(m)

    crust = reaches(CRUST_CORE_DENSITY_FM3)
    drip = reaches(NEUTRON_DRIP_DENSITY_FM3)
    return Structure(
        eos_index=eos.index,
        central_pressure_mev_fm3=central_pressure_mev_fm3,
        central_baryon_density_fm3=eos.density_fm3(log_p),
        mass_solar=mass / SOLAR_MASS_KM,
        radius_km=radius,
        crust_core_radius_km=crust[0] if crust else 0.0,
        drip_radius_km=drip[0] if drip else 0.0,
        crust_mass_solar=(mass - crust[1]) / SOLAR_MASS_KM if crust else 0.0)


def _mass_solar(eos, log_p):
    return _integrate(eos, log_p).y[1, -1] / SOLAR_MASS_KM


@functools.lru_cache(maxsize=None)
def _maximum_mass(index):
    """(log pressure, mass in solar masses) of the heaviest star of the equation of state, or the table's end if the
    mass is still rising there."""
    eos = _eos(index)
    low = eos.log_pressure(_CURVE_START_DENSITY_FM3)
    grid = np.linspace(low, eos.log_p_max - 1e-6, _CENTRAL_PRESSURE_STEPS)
    masses = np.array([_mass_solar(eos, s) for s in grid])
    top = int(np.argmax(masses))
    if top == len(grid) - 1:
        return float(grid[-1]), float(masses[-1])
    found = minimize_scalar(lambda s: -_mass_solar(eos, s), bounds=(grid[max(top - 1, 0)], grid[top + 1]),
                            method="bounded", options={"xatol": 1e-6})
    return float(found.x), float(-found.fun)


def integrated_maximum_mass_solar(eos):
    """The largest mass, solar masses, the equation of state (1 to 21) supports, by integration; a heavier star collapses."""
    return _maximum_mass(_eos(eos).index)[1]


@functools.lru_cache(maxsize=None)
def _curve(index):
    eos = _eos(index)
    log_p_max = _maximum_mass(index)[0]
    low = eos.log_pressure(_CURVE_START_DENSITY_FM3)
    masses, radii = [], []
    for s in np.linspace(low, log_p_max, 2 * _CURVE_MASSES):
        solution = _integrate(eos, float(s))
        masses.append(solution.y[1, -1] / SOLAR_MASS_KM)
        radii.append(solution.y[0, -1])
    masses, radii = np.asarray(masses), np.asarray(radii)
    stable = np.maximum.accumulate(masses) == masses
    return masses[stable], radii[stable]


def integrated_curve(eos):
    """`(masses in solar masses, radii in km)`: the stable branch of the equation of state (1 to 21), rising in mass, from
    the star of density 0.1 fm^-3 at its centre to the heaviest."""
    masses, radii = _curve(_eos(eos).index)
    return masses.copy(), radii.copy()


def curve_with_crust(eos, count):
    """`(masses, radii, crust thicknesses, crust masses)` of `count` stars of the equation of state (1 to 21), rising in mass from the star
    of density `_CURVE_START_DENSITY_FM3` at its centre to the heaviest (the arrays `scripts/build_neutron_star_table.py` writes)."""
    eos = _eos(eos)
    log_p_max = _maximum_mass(eos.index)[0]
    low = eos.log_pressure(_CURVE_START_DENSITY_FM3)
    stars = [structure(eos, math.exp(float(s)) * MEV_FM3_PER_KM2) for s in np.linspace(low, log_p_max, count)]
    return (np.array([star.mass_solar for star in stars]), np.array([star.radius_km for star in stars]),
            np.array([star.crust_thickness_km for star in stars]), np.array([star.crust_mass_solar for star in stars]))


# --- the reduced table generation reads -----------------------------------------------------------------------------------------

TABLE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "neutron_star_curves.json")
"""str: The mass-radius curves and crusts of the 21 equations of state, integrated by `scripts/build_neutron_star_table.py`."""


@functools.lru_cache(maxsize=1)
def _table():
    with open(TABLE, encoding="utf-8") as handle:
        document = json.load(handle)
    return {int(index): {"maximum_mass_solar": model["maximum_mass_solar"],
                         "mass_solar": np.array(model["mass_solar"]), "radius_km": np.array(model["radius_km"]),
                         "crust_thickness_km": np.array(model["crust_thickness_km"]),
                         "crust_mass_solar": np.array(model["crust_mass_solar"])}
            for index, model in document["models"].items()}


def _model(eos):
    try:
        return _table()[eos]
    except KeyError:
        raise ValueError(f"equation of state {eos!r}: the table holds 1 to {neutron_star_matter.EOS_COUNT}") from None


def maximum_mass_solar(eos):
    """The largest mass, solar masses, the equation of state (1 to 21) supports; a heavier star collapses."""
    return _model(eos)["maximum_mass_solar"]


def radius_km(eos, mass_solar):
    """The radius, km, of a star of this mass of equation of state `eos` (1 to 21), or None if the equation of state does
    not support so heavy a star or the mass is below the curve (0.15 to 0.5 solar masses, by equation of state)."""
    model = _model(eos)
    masses = model["mass_solar"]
    if not masses[0] <= mass_solar <= masses[-1]:
        return None
    return float(np.interp(mass_solar, masses, model["radius_km"]))


def crust_thickness_km(eos, mass_solar):
    """The thickness, km, of the crust (the layer above `CRUST_CORE_DENSITY_FM3`) of a star of this mass, or None as for `radius_km`."""
    model = _model(eos)
    masses = model["mass_solar"]
    if not masses[0] <= mass_solar <= masses[-1]:
        return None
    return float(np.interp(mass_solar, masses, model["crust_thickness_km"]))


def supporting(mass_solar):
    """The equations of state (1 to 21) that support a star of this mass, in order."""
    return tuple(index for index in range(1, neutron_star_matter.EOS_COUNT + 1) if radius_km(index, mass_solar) is not None)


def radius_band_km(mass_solar):
    """`(lowest, median, highest)` radius, km, of a star of this mass over the equations of state that support it, or
    None if none does."""
    radii = [radius_km(index, mass_solar) for index in supporting(mass_solar)]
    if not radii:
        return None
    return min(radii), float(np.median(radii)), max(radii)
