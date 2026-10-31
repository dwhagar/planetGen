# planetgen/physics/radiative_convective.py

"""
A one-dimensional radiative-convective model of an air: the surface temperature that balances the absorbed sunlight (GEN.715, GEN.717,
GEN.718; docs/design/composition/infrared-opacity.md)
=======================================================================================================================================

The air of a body is cut into `LAYERS` layers equally spaced in ln P between a model top (a ten-thousandth of the surface pressure)
and the surface. A temperature profile is found that

* follows the adiabat from the surface up to the tropopause (the dry adiabat of the mixture with its heat capacity at each
  temperature, the water pseudo-adiabat once the vapour saturates, `Column.adiabat`), and
* is in radiative equilibrium above it (every layer takes up as much infrared as it gives, plus any sunlight deposited in it), and
* sends out to space exactly the absorbed flux.

The infrared flux is solved band by band and g point by g point (`line_opacity.layer_depth`) with the two-stream form
for a source linear in optical depth between level temperatures (diffusivity factor 1.66). The tropopause is where the radiative
equilibrium temperature meets the adiabat. The water vapour has Manabe and Wetherald's (1967) relative humidity
profile, `RH = RH_s (p/p_s - 0.02) / 0.98`, with the surface humidity the mixture's own vapour pressure over the saturation pressure
at the surface (a surface that is saturated holds the saturation; a drier one holds what the mixture gives), and no vapour
above the coldest level below it (the cold trap).

Why this and not the grey closure of GEN.703: the closure needs one number for the whole column and cannot say what a 15 micron
band does when water vapour thins aloft (GEN.717) or what a haze that absorbs sunlight high up does to the ground (GEN.718).
Here the vapour profile is a profile, and `Column.sunlight_aloft` puts a share of the absorbed flux where a haze takes it.

Limits (stated, registry rows IR1 to IR8 in docs/design/research/gap-filling-plan.md section 6): the only infrared cloud is the sulfuric
acid one the caller passes (`cloud`), scaled for scattering, and no condensation of anything but water; the humidity profile is Manabe and Wetherald's, not a computed one; sunlight is deposited by the
caller, never absorbed by the gases here; the stratosphere is clear-sky radiative equilibrium without ozone heating
(a gas without a table has stand-in bands); the tables hold 100 K as the lowest temperature and 100 bar as the highest pressure, beyond
which a lookup takes the edge; the model is a global mean (one column, the absorbed flux already averaged over the sphere).

Cost: a call takes one to ten seconds (about 15 surface temperatures tried, each with a handful of Newton steps for the stratosphere,
each step a few dozen flux sweeps). The generator needs it for few bodies; a table of results is the way to many (a follow-up).

Use:

    solution = radiative_convective.solve(partial_pressures_pa, gravity_m_s2, absorbed_flux_w_m2)
    solution.surface_temperature_k, solution.tropopause_pressure_pa
"""

import functools
import math
import re
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from planetgen.physics import (
    constants,
    gases,
    greenhouse,
    cloud_infrared,
    line_opacity,
    mixture,
    radiative_kernel,
    saturation,
    species_thermo,
)
from planetgen.physics.air_column import AVOGADRO

LAYERS = 40
"""int: The number of layers (the surface temperature moves by about a kelvin between 30 and 60 layers)."""
TOP_FRACTION = 1.0e-4
"""float: The model top as a share of the surface pressure."""
DIFFUSIVITY = 1.66
"""float: The diffusivity factor of the two-stream flux (Elsasser 1942): the slant path that stands in for the angular integral of the radiance."""
HUMIDITY_TOP = 0.02
"""float: STAND-IN (IR7, tier T1 on Earth, T6 on any other world, L0): Manabe and Wetherald's (1967) `p/p_s` below which the relative humidity is zero; their profile
(relative humidity linear in pressure from the surface value to zero here) is Earth's, applied to every world with water vapour in its air."""
MAX_WATER_FRACTION = 0.5
"""float: STAND-IN (IR7, tier T7, L0): the most water vapour a layer holds, as a mole fraction (a steam air is not modelled by this profile)."""
WATER_TABLE_TOP_K = min(640.0, saturation.species("H2O").liquid.t_max_k)
"""float: Where the water table ends and above which the vapour is taken as unable to condense: 640 K, just under water's critical temperature (the
`saturation` module's, 647.1 K), so the latent heat, which diverges at the critical point, stays finite on the grid. A numerical limit of this
table, not a saturation constant: the curve itself is `saturation`'s."""
METHANE_TABLE_TOP_K = 188.0
"""float: Where the methane curve ends and above which methane is taken as unable to condense: 188 K, just under its critical temperature (190.56 K), so the latent
heat, which diverges at the critical point, stays finite (GEN.1030)."""

SUNLIGHT_WIDTH_LN_P = 1.0
"""float: STAND-IN (IR8, tier T7, L0): the width of the layer of sunlight deposited aloft, in ln P (a Gaussian: a haze or an upper methane band is a few
scale heights deep). The share and the pressure are inputs the caller supplies; nothing here derives them from the aerosols."""
CLOUD_LONGWAVE_SHARE = 0.13
"""float: STAND-IN (IR12, tier T9, L3): the share of the clear-sky outgoing infrared flux the cloud of an air whose surface humidity is `EARTH_HUMIDITY` keeps in.
The model is a clear-sky one and the albedo of a body is its all-sky one, so without this a world that keeps the Earth's cloud (and its cloud shortwave
cooling) would lose the cloud's longwave warming: 26.5 W/m2 of the 266 W/m2 the clear-sky Earth sends out (CERES EBAF, 0.10; that leaves the Earth twin of the
fixed point at 284.6 K, the clear-sky model's own bias of 3 K, and 0.13 puts it at 287.9 K, the measured 288 K, which is the value held). The sky sends out
`(1 - share) * clear` and the cloud scales with the surface humidity (to 1.3 times at most, a saturated surface) and is nil without water vapour. A cloud
that is a share of the flux, not a layer in the column, masks the carbon dioxide forcing and the vapour feedback by the same tenth, as the cloudy
sky of the Earth does, and keeps the feedbacks of the clear column whole (a grey layer over the sky hid all of the vapour's)."""
EARTH_HUMIDITY = 0.77
"""float: The mean relative humidity of the Earth's surface air (Manabe and Wetherald 1967), the humidity at which `CLOUD_LONGWAVE_SHARE` holds (tier T1)."""
FINITE_DIFFERENCE_K = 0.5
"""float: The step of the Jacobian of the stratosphere, K."""
STRATOSPHERE_STEPS = 25
STRATOSPHERE_TOLERANCE = 1.0e-3
"""float: The stratosphere is settled when no layer's heating is off by more than this share of the absorbed flux."""
MAX_STEP_K = 40.0
"""float: The largest change of a level temperature in one Newton step."""
SURFACE_TOLERANCE_K = 0.01
"""float: The surface temperature is found to this."""
POLISH_STEPS = 3
NULL_MODE_CUT = 1.0e-4
"""float: Singular values of the stratosphere's Jacobian under this share of the largest are dropped in a Newton step. In a thin layer the heating depends on the sum of its two
level temperatures only, so a zig-zag from level to level changes nothing and the minimum-norm step leaves it where it was."""
SMOOTHING_PASSES = 3
"""int: The guess a Newton solve starts from is the last profile smoothed this many times by [1/4, 1/2, 1/4], so that no zig-zag is carried over."""
CONSISTENT_PASSES = 4
"""int: At most this many passes between the optics and the profile at one surface temperature."""
CONSISTENT_K = 0.15
"""float: The optics and the profile are consistent when the levels move by less than this (K) between passes."""
LOWEST_TEMPERATURE_FRACTION = 0.3
"""float: No level is colder than this share of the skin temperature."""
HIGHEST_SURFACE_K = 4000.0
"""float: The search for the surface temperature gives up here (the thermodynamic fits stop at 4000 K)."""
BALANCE_TOLERANCE = 0.02
"""float: The polishing stops when the outgoing flux is within this (W/m2) of the absorbed one."""
CONVERGED_W_M2 = 0.25
CONVERGED_SHARE = 0.005
"""float: A solution counts as converged when the imbalance is under the larger of `CONVERGED_W_M2` and this share of the absorbed flux."""

SKIN_NOTE = "skin temperature: the radiative equilibrium of a thin grey layer, (F / 2 sigma)^(1/4)"

SIGMA = constants.STEFAN_BOLTZMANN_CONSTANT
WATER_MOLAR_MASS = gases.gas("H2O").molar_mass
R = constants.R


@dataclass(frozen=True)
class Solution:
    """
    The radiative-convective state of one air.

    Attributes:
        surface_temperature_k (float): The temperature of the surface (and of the air at it).
        level_pressure_pa (ndarray): The pressure of each level, top to surface.
        level_temperature_k (ndarray): The temperature at each level.
        tropopause_pressure_pa (float): The pressure where the profile leaves the adiabat for radiative equilibrium.
        tropopause_temperature_k (float): The temperature there.
        outgoing_w_m2 (float): The infrared flux to space of the solution (equal to the absorbed flux when `converged`).
        water_fraction (ndarray): The water vapour mole fraction of each layer.
        converged (bool): The balance was met to within `CONVERGED_W_M2` or `CONVERGED_SHARE` of the absorbed flux.
    """

    surface_temperature_k: float
    level_pressure_pa: np.ndarray
    level_temperature_k: np.ndarray
    tropopause_pressure_pa: float
    tropopause_temperature_k: float
    outgoing_w_m2: float
    water_fraction: np.ndarray
    converged: bool


# ---- water and heat capacity ----

@functools.lru_cache(maxsize=1)
def _water_curve():
    """`(temperature grid, ln of the saturation pressure, latent heat J/mol)` of water from 30 K to `WATER_TABLE_TOP_K`."""
    grid = np.arange(30.0, WATER_TABLE_TOP_K + 0.25, 0.5)
    log_pressure = np.log([saturation.pressure_pa("H2O", float(t)) for t in grid])
    latent = R * grid ** 2 * np.gradient(log_pressure, grid)
    return grid, log_pressure, latent


def saturation_pressure_pa(temperature_k):
    """Water's saturation pressure at `temperature_k` (scalar or array), Pa; infinite above `WATER_TABLE_TOP_K`."""
    grid, log_pressure, _latent = _water_curve()
    temperature = np.asarray(temperature_k, dtype=float)
    value = np.exp(np.interp(np.maximum(temperature, grid[0]), grid, log_pressure))
    return np.where(temperature > grid[-1], np.inf, value)


def latent_heat_j_mol(temperature_k):
    """Water's latent heat at `temperature_k`, J/mol, from the slope of its saturation curve (Clausius and Clapeyron)."""
    grid, _log, latent = _water_curve()
    return np.interp(temperature_k, grid, latent)


@functools.lru_cache(maxsize=1)
def _methane_curve():
    """`(temperature grid, ln of the saturation pressure, latent heat J/mol)` of methane from 30 K to `METHANE_TABLE_TOP_K` (GEN.1030)."""
    grid = np.arange(30.0, METHANE_TABLE_TOP_K + 0.25, 0.5)
    log_pressure = np.log([saturation.pressure_pa("CH4", float(t)) for t in grid])
    latent = R * grid ** 2 * np.gradient(log_pressure, grid)
    return grid, log_pressure, latent


def methane_saturation_pressure_pa(temperature_k):
    """Methane's saturation pressure at `temperature_k` (scalar or array), Pa; infinite above `METHANE_TABLE_TOP_K`."""
    grid, log_pressure, _latent = _methane_curve()
    temperature = np.asarray(temperature_k, dtype=float)
    value = np.exp(np.interp(np.maximum(temperature, grid[0]), grid, log_pressure))
    return np.where(temperature > grid[-1], np.inf, value)


def methane_latent_heat_j_mol(temperature_k):
    """Methane's latent heat at `temperature_k`, J/mol, from the slope of its saturation curve (Clausius and Clapeyron)."""
    grid, _log, latent = _methane_curve()
    return np.interp(temperature_k, grid, latent)


_CP_GRID = np.arange(40.0, 4001.0, 20.0)

H2_ROTATION_TEMPERATURE_K = 87.55
"""float: The rotational temperature of molecular hydrogen, hc B / k with B = 60.853 cm-1 (Herzberg 1950; NIST Chemistry WebBook) and hc/k = 1.43878 cm K."""
EQUILIBRIUM_HYDROGEN_BELOW_K = 200.0
"""float: Below this temperature the heat capacity of hydrogen is that of equilibrium ortho and para forms (GEN.1030); the NASA fit of `species_thermo` starts here, where the two agree to 1 percent."""


def equilibrium_hydrogen_cp_over_r(temperature_k):
    """
    The heat capacity at constant pressure over R of molecular hydrogen whose ortho and para forms stay in equilibrium (translation 5/2 plus the rotational
    fluctuation, levels of weight (2J + 1) with a nuclear spin factor of 3 for odd J and 1 for even J). Cold hydrogen is far from the 7/2 of a hot one: it
    reaches 4.4 at 60 K, which makes the lapse rate of a giant's troposphere there about a fifth shallower (Gierasch 1985).
    """
    temperature = np.asarray(temperature_k, dtype=float)[..., None]
    level = np.arange(40)
    energy = level * (level + 1) * H2_ROTATION_TEMPERATURE_K
    weight = (2 * level + 1) * np.where(level % 2 == 1, 3.0, 1.0)
    boltzmann = weight * np.exp(-energy / temperature)
    total = boltzmann.sum(axis=-1)
    mean = (boltzmann * energy).sum(axis=-1) / total
    spread = (boltzmann * energy ** 2).sum(axis=-1) / total - mean ** 2
    return 2.5 + spread / temperature[..., 0] ** 2


def _atoms(formula):
    """The number of atoms in a molecule written without brackets (`"HC3N"` is 5)."""
    return sum(int(count or 1) for _symbol, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula))


@functools.cache
def _cp_over_r(formula):
    """
    The heat capacity at constant pressure over R of a gas on `_CP_GRID`, from the first source the gas has (GEN.852, row IR9):

    1. its NASA Glenn fit (`species_thermo.gas`; the cold end of the fit stands in below it, the hot end above it; for hydrogen, below
       `EQUILIBRIUM_HYDROGEN_BELOW_K` the heat capacity of ortho and para forms in equilibrium instead, GEN.1030, row IR15);
    2. the heat capacities the gas table gives at 100, 300 and 1000 K (`gases.Gas.cp_j_mol_k`), interpolated linearly in temperature and
       held at the end values beyond them (a polyatomic gas gains heat capacity above 1000 K, so this is low there);
    3. a rigid-rotor count of the atoms: 5/2 for one, 7/2 for two, and 4 + (3N - 6) / 2 for a molecule of N atoms, half its vibrations awake.
    """
    try:
        thermo = species_thermo.gas(formula)
    except KeyError:
        pass
    else:
        low, high = thermo.range_k
        values = np.array([thermo.cp_over_r(min(max(t, low), high)) for t in _CP_GRID])
        if formula == "H2":
            join = float(equilibrium_hydrogen_cp_over_r(EQUILIBRIUM_HYDROGEN_BELOW_K))
            cold = _CP_GRID < EQUILIBRIUM_HYDROGEN_BELOW_K
            values[cold] = equilibrium_hydrogen_cp_over_r(_CP_GRID[cold]) + (thermo.cp_over_r(max(low, EQUILIBRIUM_HYDROGEN_BELOW_K)) - join)
        return values
    known = sorted((float(t), v / R) for t, v in gases.gas(formula).cp_j_mol_k.items() if v is not None)
    if known:
        return np.interp(_CP_GRID, [t for t, _ in known], [v for _, v in known])
    atoms = _atoms(formula)
    return np.full_like(_CP_GRID, {1: 2.5, 2: 3.5}.get(atoms, 4.0 + 0.5 * (3 * atoms - 6)))


# ---- the column ----

class Column:
    """An air cut into layers: its geometry, its composition and the pieces of the thermodynamics the solver needs."""

    def __init__(self, partial_pressures_pa, gravity_m_s2, layers=LAYERS, water_fraction=None, cloud=None):
        self.cloud = cloud
        self.fixed_water = None if water_fraction is None else np.asarray(water_fraction, dtype=float)
        fractions = mixture.mole_fractions(partial_pressures_pa)
        self.gravity = gravity_m_s2
        self.surface_pressure_pa = mixture.total_pressure_pa(partial_pressures_pa)
        self.layers = layers
        top = self.surface_pressure_pa * TOP_FRACTION
        self.level_pressure = np.exp(np.linspace(math.log(top), math.log(self.surface_pressure_pa), layers + 1))
        self.layer_pressure = np.sqrt(self.level_pressure[1:] * self.level_pressure[:-1])
        self.layer_mass = np.diff(self.level_pressure) / gravity_m_s2
        self.water_surface = fractions.get("H2O", 0.0)
        dry_total = 1.0 - self.water_surface
        self.dry = {formula: fraction / dry_total for formula, fraction in fractions.items() if formula != "H2O"} if dry_total > 0.0 else {}
        self.dry_molar_mass = sum(fraction * gases.gas(formula).molar_mass for formula, fraction in self.dry.items())
        self.cp_dry = R * sum(fraction * _cp_over_r(formula) for formula, fraction in self.dry.items()) if self.dry else np.zeros_like(_CP_GRID)
        self.cp_water = R * _cp_over_r("H2O")
        # the condensable of the adiabat (GEN.1030): water when the air holds any; else methane, whose condensation makes the lapse rate of the cold
        # hydrogen giants about half the dry one (Uranus and Neptune near 1.3 bar). Every other gas is the background it rides in.
        methane = fractions.get("CH4", 0.0) if self.water_surface <= 0.0 else 0.0
        if self.water_surface > 0.0:
            self.condensable_surface, self.background_molar_mass, self.cp_background = self.water_surface, self.dry_molar_mass, self.cp_dry
            self.cp_condensable = self.cp_water
            self._saturation, self._latent = saturation_pressure_pa, latent_heat_j_mol
        elif 0.0 < methane < 1.0:
            others = {formula: fraction / (1.0 - methane) for formula, fraction in fractions.items() if formula != "CH4"}
            self.condensable_surface = methane
            self.background_molar_mass = sum(fraction * gases.gas(formula).molar_mass for formula, fraction in others.items())
            self.cp_background = R * sum(fraction * _cp_over_r(formula) for formula, fraction in others.items())
            self.cp_condensable = R * _cp_over_r("CH4")
            self._saturation, self._latent = methane_saturation_pressure_pa, methane_latent_heat_j_mol
        else:
            self.condensable_surface, self.background_molar_mass, self.cp_background = 0.0, 0.0, self.cp_dry
            self.cp_condensable = self.cp_water
            self._saturation, self._latent = saturation_pressure_pa, latent_heat_j_mol
        self.prepared = line_opacity.prepare(partial_pressures_pa)

    # -- the adiabat --

    def _dlnT_dlnP(self, temperature, pressure, condensing=None):
        """d ln T / d ln P along the adiabat: dry or unsaturated, or the pseudo-adiabat of the condensable (water, else methane) once its vapour saturates
        (`condensing`: take that branch whatever the state, so a step that starts on one branch uses it at its midpoint too; None decides from the state)."""
        x = self.condensable_surface
        if x <= 0.0:
            return R / np.interp(temperature, _CP_GRID, self.cp_dry)
        cp_dry = np.interp(temperature, _CP_GRID, self.cp_background)
        cp_water = np.interp(temperature, _CP_GRID, self.cp_condensable)
        saturated = float(self._saturation(temperature))
        if not (x * pressure >= saturated if condensing is None else condensing):
            return R / ((1.0 - x) * cp_dry + x * cp_water)
        latent = float(self._latent(temperature))
        if self.background_molar_mass <= 0.0 or pressure <= saturated * (1.0 + 1.0e-9):
            return min(R * temperature / latent, R / cp_water)
        r = saturated / (pressure - saturated)
        numerator = (1.0 + r) * (R + r * latent / temperature)
        denominator = (cp_dry + r * cp_water - r * latent / temperature + r * latent / temperature
                       + r * (1.0 + r) * latent ** 2 / (R * temperature ** 2))
        return numerator / denominator

    def _saturated(self, temperature, pressure):
        """Whether the surface's vapour mixing ratio of the condensable saturates at this temperature and pressure."""
        return self.condensable_surface > 0.0 and self.condensable_surface * pressure >= float(self._saturation(temperature))

    def _advance(self, temperature, log_pressure, step):
        """One second-order (midpoint) step of ln T along the adiabat, with the branch (dry or saturated) taken at the start."""
        branch = self._saturated(temperature, math.exp(log_pressure))
        slope = self._dlnT_dlnP(temperature, math.exp(log_pressure), branch)
        midpoint = temperature * math.exp(0.5 * step * slope)
        slope = self._dlnT_dlnP(midpoint, math.exp(log_pressure + 0.5 * step), branch)
        return temperature * math.exp(step * slope)

    def adiabat(self, surface_temperature_k):
        """
        The adiabat's temperature at every level but the surface (top first), integrating up from the surface with second-order steps in
        ln P. The slope has a corner where the vapour saturates; a step that crosses it is cut at the crossing (bisection) so the
        temperature is a smooth function of the surface temperature.
        """
        temperature = float(surface_temperature_k)
        log_pressure = math.log(self.surface_pressure_pa)
        result = np.empty(self.layers)
        for index in range(self.layers - 1, -1, -1):
            target = math.log(self.level_pressure[index])
            steps = max(1, int(abs(log_pressure - target) / 0.05) + 1)
            step = (target - log_pressure) / steps
            for _ in range(steps):
                remaining = step
                for _crossing in range(4):
                    before = self._saturated(temperature, math.exp(log_pressure))
                    trial = self._advance(temperature, log_pressure, remaining)
                    if self._saturated(trial, math.exp(log_pressure + remaining)) == before:
                        temperature, log_pressure = trial, log_pressure + remaining
                        break
                    low, high = 0.0, 1.0
                    for _ in range(30):
                        middle = 0.5 * (low + high)
                        moved = self._advance(temperature, log_pressure, remaining * middle)
                        if self._saturated(moved, math.exp(log_pressure + remaining * middle)) == before:
                            low = middle
                        else:
                            high = middle
                    moved = self._advance(temperature, log_pressure, remaining * high)
                    temperature, log_pressure = moved, log_pressure + remaining * high
                    remaining *= 1.0 - high
                    if remaining == 0.0:
                        break
            result[index] = temperature
        return result

    # -- the water --

    def water_profile(self, layer_temperature, surface_temperature):
        """The water vapour mole fraction of each layer: Manabe and Wetherald's humidity profile, no more than the cold trap leaves."""
        if self.fixed_water is not None:
            return self.fixed_water
        if not self.dry:
            return np.ones(self.layers)
        if self.water_surface <= 0.0:
            return np.zeros(self.layers)
        saturated_surface = float(saturation_pressure_pa(surface_temperature))
        shape = np.clip((self.layer_pressure / self.surface_pressure_pa - HUMIDITY_TOP) / (1.0 - HUMIDITY_TOP), 0.0, None)
        if math.isfinite(saturated_surface):
            humidity = min(1.0, self.water_surface * self.surface_pressure_pa / saturated_surface)
            saturated = saturation_pressure_pa(layer_temperature)
            fraction = np.where(np.isfinite(saturated), humidity * shape * saturated / self.layer_pressure, self.water_surface * shape)
        else:
            fraction = self.water_surface * shape
        fraction = np.minimum(fraction, MAX_WATER_FRACTION)
        return np.minimum.accumulate(fraction[::-1])[::-1]

    # -- the optics --

    def optical_depth(self, level_temperature):
        """The optical depth of every layer (layer, band, g) and the water profile, for the temperatures at the levels (top first, surface last)."""
        layer_temperature = 0.5 * (level_temperature[:-1] + level_temperature[1:])
        water = self.water_profile(layer_temperature, level_temperature[-1])
        mean_molar_mass = (1.0 - water) * self.dry_molar_mass + water * WATER_MOLAR_MASS
        molecules_cm2 = self.layer_mass / mean_molar_mass * AVOGADRO * 1.0e-4
        depth = line_opacity.layer_depth(self.prepared, layer_temperature, self.layer_pressure, water, molecules_cm2)
        if self.cloud is not None:
            tau_550, radius_um = self.cloud
            depth = depth + cloud_infrared.layer_depth(self.layer_pressure, self.layer_mass, layer_temperature, tau_550, radius_um)[:, :, None]
        return depth, water

    def cloud_share(self, surface_temperature):
        """The share of the clear-sky outgoing flux the cloud keeps in (stand-in IR12): `CLOUD_LONGWAVE_SHARE` at the humidity of the Earth, in proportion to the
        surface humidity, nil where the air holds no water vapour."""
        if self.water_surface <= 0.0:
            return 0.0
        if not self.dry:
            return CLOUD_LONGWAVE_SHARE / EARTH_HUMIDITY
        saturated = float(saturation_pressure_pa(surface_temperature))
        humidity = min(1.0, self.water_surface * self.surface_pressure_pa / saturated) if math.isfinite(saturated) else 1.0
        return CLOUD_LONGWAVE_SHARE * humidity / EARTH_HUMIDITY


# ---- the flux ----

def sweep(level_temperature, depth, surface_temperature):
    """
    The upward and downward infrared flux at every level, W/m2: `level_temperature` is (profile, level), `depth` the layers' optical
    depths (layer, band, g), `surface_temperature` (profile). The source is linear in optical depth between level temperatures. The
    layer-by-layer sweep is `radiative_kernel.sweep`.
    """
    emission = line_opacity.band_emission(level_temperature)                  # (profile, level, band)
    ground = line_opacity.band_emission(surface_temperature)                  # (profile, band)
    return radiative_kernel.sweep(np.ascontiguousarray(emission), np.ascontiguousarray(ground), np.ascontiguousarray(depth),
                                  line_opacity.G_WEIGHTS, DIFFUSIVITY)


def _upper_levels(layer_temperature):
    """The temperatures at the levels above the surface from layer temperatures: the mean of neighbours, the top one extrapolated."""
    count, layers = layer_temperature.shape
    levels = np.empty((count, layers + 1))
    levels[:, 1:-1] = 0.5 * (layer_temperature[:, 1:] + layer_temperature[:, :-1])
    levels[:, 0] = layer_temperature[:, 0] - 0.5 * (layer_temperature[:, 1] - layer_temperature[:, 0])
    return levels


# ---- the solver ----

class _Solver:
    """Finds the surface temperature whose radiative-convective profile sends the absorbed flux to space."""

    def __init__(self, column, absorbed_w_m2, heating_w_m2=None):
        self.column = column
        self.absorbed = absorbed_w_m2
        self.heating = np.zeros(column.layers) if heating_w_m2 is None else np.asarray(heating_w_m2, dtype=float)
        self.effective = (absorbed_w_m2 / SIGMA) ** 0.25
        self.skin = (absorbed_w_m2 / (2.0 * SIGMA)) ** 0.25
        self.floor = LOWEST_TEMPERATURE_FRACTION * self.skin
        self.adiabats = {}
        self.jacobians = {}
        self.stratosphere = np.full(column.layers, self.skin)
        self.tropopause = None
        self.last = None

    def adiabat(self, surface_temperature):
        key = round(surface_temperature, 7)
        if key not in self.adiabats:
            if len(self.adiabats) > 300:
                self.adiabats.clear()
            self.adiabats[key] = self.column.adiabat(surface_temperature)
        return self.adiabats[key]

    def _flux(self, levels, depth):
        up, down = sweep(levels, depth, levels[:, -1])
        net = up - down
        return net, net[:, 1:] - net[:, :-1] + self.heating[None, :]

    def _equilibrium(self, adiabat, surface, m, depth, start):
        """Solves the radiative equilibrium of the `m` top levels with the others on the adiabat: Newton's method with a finite-difference
        Jacobian that is made again only when a step has not cut the residual in half."""
        layers = self.column.layers
        base = np.empty(layers + 1)
        base[:layers] = adiabat
        base[layers] = surface
        x = start[:m].copy()
        jacobian = self.jacobians.get(m)
        previous = None
        scale = max(self.absorbed, 1.0e-9)
        for _ in range(STRATOSPHERE_STEPS):
            levels = base.copy()
            levels[:m] = x
            net, heating = self._flux(levels[None, :], depth)
            residual = heating[0, :m]
            size = float(np.abs(residual).max()) if m else 0.0
            if m == 0 or size < STRATOSPHERE_TOLERANCE * scale:
                break
            if jacobian is None or (previous is not None and size > 0.5 * previous) or (previous is None and size > 5.0e-3 * scale):
                batch = np.tile(levels, (m + 1, 1))
                for j in range(m):
                    batch[j + 1, j] += FINITE_DIFFERENCE_K
                _net, heated = self._flux(batch, depth)
                jacobian = ((heated[1:, :m] - residual[None, :]) / FINITE_DIFFERENCE_K).T
                self.jacobians[m] = jacobian
            previous = size
            step = -np.linalg.lstsq(jacobian, residual, rcond=NULL_MODE_CUT)[0]
            if not np.all(np.isfinite(step)):
                step = np.zeros_like(step)
            biggest = float(np.abs(step).max())
            if biggest > MAX_STEP_K:
                step *= MAX_STEP_K / biggest
            x = np.maximum(x + step, self.floor)
            if biggest < 1.0e-3:
                break
        levels = base.copy()
        levels[:m] = x
        net, heating = self._flux(levels[None, :], depth)
        return x, levels, heating[0], net[0]

    def _profile(self, surface_temperature, adiabat, m, depth, start):
        """The radiative-convective profile for fixed optics: the tropopause level found by moving `m` until the layer below it is convective
        and the one above it is not, then the equilibrium above it. Returns `(m, x, levels)`."""
        layers = self.column.layers
        tried = {}
        accepted = False
        for _ in range(10):
            x, levels, heating, _net = self._equilibrium(adiabat, surface_temperature, m, depth, start)
            below = heating[m] if m < layers else -1.0
            too_stable = m < layers - 1 and below > 1.0e-3 * self.absorbed
            too_unstable = m > 1 and x[m - 1] < adiabat[m - 1] - 0.05
            violation = max(below / self.absorbed if too_stable else 0.0,
                            (adiabat[m - 1] - x[m - 1]) / 50.0 if too_unstable else 0.0)
            tried[m] = violation
            if not too_stable and not too_unstable:
                accepted = True
                break
            if too_stable and not too_unstable:
                new = m + 1
            elif too_unstable and not too_stable:
                new = m - 1
            else:
                new = m + 1 if below / self.absorbed > (adiabat[m - 1] - x[m - 1]) / 50.0 else m - 1
            if new in tried:
                break
            if new > m:
                start = np.concatenate([x, [max(adiabat[m], self.skin)], np.full(layers - m - 1, self.skin)])[:layers]
            else:
                start = np.concatenate([x[:new], np.full(layers - new, self.skin)])
            m = new
        if not accepted:
            m = min(tried, key=tried.get)
            x, levels, _heating, _net = self._equilibrium(adiabat, surface_temperature, m, depth, start)
        return m, x, levels

    def outgoing(self, surface_temperature):
        """The infrared flux to space of the best profile at one surface temperature. The optics come from the profile and the profile from
        the optics, so the two are passed back and forth until the levels move by less than `CONSISTENT_K`."""
        column = self.column
        layers = column.layers
        adiabat = self.adiabat(surface_temperature)
        if self.tropopause is None:
            warmer = adiabat > self.skin
            self.tropopause = max(1, int(np.argmax(warmer)) if warmer.any() else layers)
        m = min(self.tropopause, layers)
        levels = np.empty(layers + 1)
        levels[:layers] = adiabat
        levels[layers] = surface_temperature
        smooth = self.stratosphere.copy()
        for _ in range(SMOOTHING_PASSES):
            smooth[1:-1] = 0.25 * smooth[:-2] + 0.5 * smooth[1:-1] + 0.25 * smooth[2:]
        levels[:m] = np.maximum(smooth[:m], adiabat[:m])
        start = np.maximum(smooth, self.skin * 0.5)
        for _ in range(CONSISTENT_PASSES):
            depth, water = column.optical_depth(levels)
            m, x, new = self._profile(surface_temperature, adiabat, m, depth, start)
            new = np.maximum(np.nan_to_num(new, nan=self.skin), self.floor)
            change = float(np.abs(new - levels).max())
            levels = new
            start = np.concatenate([x, np.full(layers - m, self.skin)])
            if change < CONSISTENT_K:
                break
        self.tropopause = m
        self.stratosphere[:m] = levels[:m]
        depth, water = column.optical_depth(levels)
        up, down = sweep(levels[None, :], depth, levels[None, -1])
        self.last = (surface_temperature, levels, m, water)
        return float((up - down)[0, 0]) * (1.0 - column.cloud_share(surface_temperature))

    def residual(self, surface_temperature):
        return self.outgoing(surface_temperature) - self.absorbed

    def solve(self):
        low = self.effective
        high = self.effective * 1.3
        f_low = self.residual(low)
        f_high = self.residual(high)
        count = 0
        while f_high < 0.0 and count < 25 and high < HIGHEST_SURFACE_K:
            low, f_low = high, f_high
            high = min(self.effective + (high - self.effective) * 1.5 + 3.0, HIGHEST_SURFACE_K)
            f_high = self.residual(high)
            count += 1
        count = 0
        while f_low > 0.0 and count < 25:
            high, f_high = low, f_low
            low = max(self.effective * 0.5, low - (high - low) * 1.5 - 3.0)
            f_low = self.residual(low)
            count += 1
            if low <= self.effective * 0.5:
                break
        if f_low <= 0.0 <= f_high:
            try:
                surface = brentq(self.residual, low, high, xtol=SURFACE_TOLERANCE_K, rtol=1.0e-6, maxiter=40)
            except ValueError:
                surface = low if abs(f_low) < abs(f_high) else high
        else:
            surface = low if abs(f_low) < abs(f_high) else high
        balance = self.residual(surface)
        for _ in range(POLISH_STEPS):
            if abs(balance) < BALANCE_TOLERANCE:
                break
            other = self.residual(surface + 0.2)
            slope = (other - balance) / 0.2
            if not slope > 0.0:
                break
            surface = max(self.effective * 0.5, surface - balance / slope)
            balance = self.residual(surface)
        return surface, balance


def sunlight_aloft(column, absorbed_w_m2, share, pressure_pa):
    """
    The heating of each layer, W/m2, of a share of the absorbed flux deposited aloft: a Gaussian in ln P of width `SUNLIGHT_WIDTH_LN_P`
    centred on `pressure_pa` (a haze or an upper methane band takes sunlight there instead of the ground).
    """
    if not 0.0 <= share < 1.0:
        raise ValueError("the share of sunlight deposited aloft is at least 0 and under 1")
    if share == 0.0:
        return np.zeros(column.layers)
    weight = np.exp(-0.5 * ((np.log(column.layer_pressure) - math.log(pressure_pa)) / SUNLIGHT_WIDTH_LN_P) ** 2)
    return share * absorbed_w_m2 * weight / weight.sum()


def solve(partial_pressures_pa, gravity_m_s2, absorbed_flux_w_m2, aloft_share=0.0, aloft_pressure_pa=None, layers=LAYERS,
          water_fraction=None, cloud=None):
    """
    The radiative-convective state of an air under an absorbed flux.

    Args:
        partial_pressures_pa (dict): Gas-table formula to partial pressure at the surface, Pa. The water vapour is the surface's own.
        gravity_m_s2 (float): Surface gravity.
        absorbed_flux_w_m2 (float): The globally averaged absorbed flux (`greenhouse.absorbed_flux_w_m2`), all of it reaching the ground
            except the `aloft_share`.
        aloft_share (float): The share of it taken high up, by a haze or upper-air absorbers (0 to under 1); the rest heats the surface.
        aloft_pressure_pa (float or None): The pressure where that share is taken (required when `aloft_share` is above 0).
        layers (int): The number of layers.
        water_fraction (ndarray or None): Hold the water vapour mole fraction of each layer (top first, `layers` long) at these values instead of
            the humidity profile (a previous solution's `water_fraction`: the experiment of a change with the vapour fixed).
        cloud (tuple or None): `(tau_550, radius_um)` of a sulfuric acid cloud (`sulfuric_cloud.of_draw`): its depth at 550 nm and droplet radius.
            Its infrared opacity (`cloud_infrared`) sits in the layers between the level cooler than the cloud base temperature and its top.

    Returns:
        Solution: Its `converged` is False when the search for the surface temperature did not reach a balance (for example when the sunlight
        aloft is all the flux there is).

    Raises:
        ValueError: The mixture is empty, or a number is not positive and finite.
    """
    for name, value in (("gravity_m_s2", gravity_m_s2), ("absorbed_flux_w_m2", absorbed_flux_w_m2)):
        if not (value > 0.0 and math.isfinite(value)):
            raise ValueError(f"{name} must be a positive number, not {value!r}")
    column = Column(partial_pressures_pa, gravity_m_s2, layers, water_fraction, cloud)
    heating = None
    if aloft_share > 0.0:
        if aloft_pressure_pa is None or not aloft_pressure_pa > 0.0:
            raise ValueError("aloft_pressure_pa is needed when aloft_share is above 0")
        heating = sunlight_aloft(column, absorbed_flux_w_m2, aloft_share, aloft_pressure_pa)
    effective = (absorbed_flux_w_m2 / SIGMA) ** 0.25
    if not column.prepared.absorbers:
        levels = np.full(layers + 1, effective)
        return Solution(effective, column.level_pressure, levels, float(column.level_pressure[0]), effective, absorbed_flux_w_m2,
                        np.zeros(layers), True)
    solver = _Solver(column, absorbed_flux_w_m2, heating)
    surface, balance = solver.solve()
    _temperature, levels, m, water = solver.last
    return Solution(surface_temperature_k=float(surface), level_pressure_pa=column.level_pressure, level_temperature_k=levels,
                    tropopause_pressure_pa=float(column.level_pressure[m]), tropopause_temperature_k=float(levels[m]),
                    outgoing_w_m2=float(absorbed_flux_w_m2 + balance), water_fraction=water, converged=abs(balance) <= max(CONVERGED_W_M2, CONVERGED_SHARE * absorbed_flux_w_m2))


@dataclass(frozen=True)
class Warming:
    """
    The greenhouse warming of one air in the repo's variable.

    Attributes:
        effective_temperature_k (float): `T_e`, the airless grey temperature of the absorbed flux.
        surface_temperature_k (float): `T_s` of the radiative-convective model.
        factor (float): `G = (T_s / T_e)^4 - 1`, the repo's greenhouse variable (never below 0).
        tropopause_pressure_pa (float): Where the profile leaves the adiabat for radiative equilibrium.
        tropopause_temperature_k (float): The temperature there.
        converged (bool): The model balanced the absorbed flux.
        solution (Solution): The whole profile.
    """

    effective_temperature_k: float
    surface_temperature_k: float
    factor: float
    tropopause_pressure_pa: float
    tropopause_temperature_k: float
    converged: bool
    solution: Solution


def warming(partial_pressures_pa, gravity_m_s2, absorbed_flux_w_m2, aloft_share=0.0, aloft_pressure_pa=None):
    """
    The warming of an air under an absorbed flux (the globally averaged W/m2 of `greenhouse.absorbed_flux_w_m2`) from `solve`, in the
    variable `atmosphere_state.greenhouse_warming` gives. A call takes one to several seconds (up to half a minute for a dense steam
    or carbon dioxide air), so the generator does not call it yet: it stays on the fast closure until a lookup table of the model's
    answers stands in front of it (filed with TODO).

    Raises:
        ValueError: As `solve`.
    """
    solution = solve(partial_pressures_pa, gravity_m_s2, absorbed_flux_w_m2, aloft_share, aloft_pressure_pa)
    effective = greenhouse.effective_temperature_k(absorbed_flux_w_m2)
    surface = solution.surface_temperature_k
    return Warming(effective_temperature_k=effective, surface_temperature_k=surface,
                   factor=max(greenhouse.greenhouse_factor(surface, effective), 0.0), tropopause_pressure_pa=solution.tropopause_pressure_pa,
                   tropopause_temperature_k=solution.tropopause_temperature_k, converged=solution.converged, solution=solution)
