# planetgen/physics/infrared_opacity.py

"""
The grey infrared optical depth tau0 of an air, built from its mixture on a wavenumber grid (GEN.703,
docs/design/composition/atmosphere-composition-parts/03-infrared-greenhouse-thermal.md sections 2, 3 and 6.2)
==========================================================================================================

`atmosphere_state` used one grey mass absorption coefficient per gas. A trace gas on a cold world then did almost nothing
(forty parts per million of carbon dioxide gave a thousandth of a unit) and the temperature never entered. Here the optical depth
is built the way route 1 of part 03 section 6.2 defines it, from a spectrum:

    tau0 = -ln( sum_k w_k(T) exp(-tau_k) )

the equivalent grey thickness of Koll (2019, arXiv:1907.13145, Eq. 13), `w_k(T)` the share of a blackbody at the surface
temperature `T` that lies in the bin `k` of a 10 cm-1 grid (0 to 6000 cm-1, with a last bin to infinity), `tau_k` the column optical
depth in it. A bin that lets light through (a window) holds `tau0` down however thick the rest is, and a band that sits under the
Planck peak of the surface (carbon dioxide at 15 um on a temperate world) matters far more than one far from it. The pieces of `tau_k`:

* **Bands** (gas table, `gas(f).infrared.bands`): each band of an absorbing gas is a plateau (its range, or its centre) of optical
  depth `kappa0 u` with exponential wings, `tau(nu) = kappa0 u exp(-d / WING_E_FOLD_CM1)`, `u` the gas's column, kg/m2, and `d` the
  distance from the plateau. The width of a saturated band grows as `2 WING_E_FOLD_CM1 ln(kappa0 u)`, the logarithm in the measured
  forcing of carbon dioxide (`F = 5.32 ln(C/C0)`, part 03 section 2.2, Byrne and Goldblatt 2014a, READ): the e-fold is set from that
  law (about 10 cm-1), and `kappa0` enters only through a logarithm, so a factor of a few in it moves nothing much.
* **Far wings** (`FAR_WING_M2_KG`): the pressure-broadened opacity between the bands, `kappa_w u (P / P_REF)` within a reach of
  600 cm-1 times `P / P_REF` of a band. Koll 2019's statements are its only source (a pure H2O air is thick above about 1e-2 bar, a pure CO2 air is of order
  unity at 1 bar, CO is thin even at 1 bar). It is what makes a dense CO2 or steam air a runaway-like greenhouse, and it grows as
  the square of the pressure (the exponent `PRESSURE_EXPONENT` of the closure is 2).
* **Water vapour** (GEN.812): the partial pressure of H2O is the surface's, and the vapour thins with height (`vapour_profile`: the
  humidity falls linearly with pressure, Manabe and Wetherald 1967; the temperature follows the troposphere; the mixing ratio never
  exceeds the surface's). Well mixed, an Earth twin at 0.77 humidity held three times the water the Earth does and came out at 319 K, not
  288 K. The H2O bands and far wing take a column that gives the same pressure-weighted opacity as the profile; the continuum is summed
  over the layers.
* **Water vapour continuum** (part 03 section 3.2, MT_CKD 4.3, `WATER_CONTINUUM`): the self and foreign coefficients at 200 to
  3600 cm-1, `tau(nu) = sum over layers of (sigma_self(nu, T) x_w + sigma_for(nu) (1 - x_w)) x_w n^2 dz / n_ref`; it closes the 8 to 12 um
  window on a humid world and is why a warm ocean world runs away.
* **Clouds** (`CLOUD_LONGWAVE_TAU0`, a stand-in): the module is a clear-sky greenhouse (outgoing 266 W/m2 at 288 K for the Earth, the
  measured clear-sky value), but the albedo of a body is its all-sky one, so one grey depth for the longwave effect of the cloud is added,
  scaled with the surface humidity and nil without liquid water.
* **Collision-induced pairs** (part 03 section 3.1, HITRAN, `CIA_PAIRS`): `kP(T) x1 x2 n0^2 H / 2` is the optical depth if the
  pair absorbed like its peak everywhere in the file's range; the saturating `ln(1 + x)` turns it into the equivalent thickness (it
  reproduces the four cases near the Planck peak in the table of part 03 section 3.1 to 25 percent and is within a factor of 3 for
  the cold ones; it is a lower bound where the file stops short of the temperature).

Limits (stated; the opacities are STAND-INS, labelled in the gap plan registry, until real line data are in): the repo holds band
centres only, no line strengths, so `kappa0` and the wing coefficients are class values; the bands of a gas without a listed band get
a plateau between 600 and 1500 cm-1 (the fingerprint region); condensation aloft is not modelled beyond the vapour profile (no cloud physics, only the grey cloud depth); the
hot bands of carbon dioxide count a hundredth as much as its fundamentals; no temperature dependence of the band strengths; the
pressure exponent of the closure is 2 for every air and the lapse scaling alpha runs 0.6 (condensable vapour supplies all the
opacity) to 0.8 (none does) and 0.85 for an air of hydrogen and helium (part 03 section 5.3: Earth 0.6, Venus 0.8, Titan 0.77,
Jupiter 0.85). Real k-tables of the gases (correlated-k, offline with ExoMol or HITRAN data) would replace the band and wing
numbers one for one; nothing else here changes.

Use:

    opacity = infrared_opacity.optical_depth(mixture, temperature_k, gravity_m_s2)
    opacity.tau0, opacity.by_part, opacity.by_gas
    # atmosphere_state.greenhouse_warming solves the surface temperature with tau0(T_s) (the fixed point)
"""

import functools
import math
from dataclasses import dataclass

import numpy as np

from planetgen.physics import constants, gases, mixture, saturation

SECOND_RADIATION_CM_K = (2.0 * math.pi * constants.REDUCED_PLANCK * constants.SPEED_OF_LIGHT_M_S / constants.BOLTZMANN) * 100.0
"""float: `h c / k` in cm K, the Planck exponent `x = c2 nu / T` for a wavenumber in cm-1."""

BIN_CM1 = 10.0
"""float: The width of a grid bin, cm-1."""
TOP_CM1 = 6000.0
"""float: The top of the grid; the last bin runs from here to infinity."""
EDGES_CM1 = np.arange(0.0, TOP_CM1 + BIN_CM1, BIN_CM1)
"""ndarray: Bin edges, cm-1."""
CENTRES_CM1 = np.concatenate([(EDGES_CM1[:-1] + EDGES_CM1[1:]) / 2.0, [TOP_CM1 + 200.0]])
"""ndarray: Bin centres, cm-1, the tail bin placed above the top."""

PRESSURE_REFERENCE_PA = 1.0e5
"""float: `P_REF`, one bar: the pressure the far-wing coefficients are stated at."""
WING_E_FOLD_CM1 = 10.0
"""float: The e-folding distance of a band's wing. Set from the logarithm of the carbon dioxide forcing law (Byrne and Goldblatt 2014a,
`F = 5.32 ln(C/C0)`): a Planck flux difference of about 0.28 W/m2 per cm-1 at 15 um on a temperate world gives about 19 cm-1 of
saturated width per e-fold of the column, 10 cm-1 on each side."""
BAND_KAPPA0_M2_KG = 1.0e3
"""float: STAND-IN: the plateau opacity of a fundamental band of a gas the table flags `"ir"`, m2/kg (a strong band, line-averaged; enters
the width only as `ln(kappa0 u)`)."""
WEAK_KAPPA0_M2_KG = 1.0e1
"""float: STAND-IN: the same for a gas the table flags `"weak"`."""
HOT_BAND_FACTOR = 1.0e-6
"""float: STAND-IN: a hot band (label contains `hot`) counts this much of a fundamental's plateau. The 15 um fundamental alone gives
the 3.7 W/m2 of a doubling of CO2 that the forcing law states; a hot band or Fermi dyad that saturates at Earth's column adds a
further 1 to 3 W/m2 each, the line-averaged opacity of such a sparse band is orders of magnitude under a fundamental's."""
DYAD_BAND_FACTOR = 1.0e-5
"""float: STAND-IN: a Fermi dyad (label contains `dyad`) counts this much of a fundamental's plateau (weak in the infrared)."""
GENERIC_BAND_CM1 = (600.0, 1500.0)
"""tuple: STAND-IN: the plateau of a gas flagged `"ir"` whose table lists no band."""

FAR_WING_M2_KG = {
    "H2O": 1.0e-2,
    "CO2": 3.0e-4,
    "CO": 1.0e-5,
}
"""
dict: The far-wing mass opacity of a gas at one bar, m2 per kg of the gas. STAND-IN, read from Koll 2019's statement of how thick a pure
gas is (arXiv:1907.13145, text of his Fig. 4, part 03 section 6.2): pure H2O air thick above about 1e-2 bar, pure CO2 of order unity
at 1 bar (taken as 3), CO thin even at 1 bar. Not tuned to Earth or Venus.
"""
GENERIC_FAR_WING_M2_KG = 1.0e-4
"""float: STAND-IN: a gas flagged `"ir"` with no entry above: weaker than CO2 (a rank of 16 in the saturated metric of part 03; CH4 is 27)."""
WING_REACH_CM1 = 600.0
"""float: STAND-IN: how far from a band the far wing reaches at one bar; the reach grows with the pressure (line widths do), so a dense
air has wings across the whole spectrum and a thin one only beside its bands."""
WEAK_FAR_WING_M2_KG = 1.0e-5
"""float: STAND-IN: a gas flagged `"weak"`: as thin as CO."""

LAPSE_SCALING_CONDENSABLE = 0.6
LAPSE_SCALING_DRY = 0.8
LAPSE_SCALING_HYDROGEN = 0.85
"""float: Robinson and Catling's alpha (0.6 to 0.9) for an air whose opacity is all condensable vapour, none of it, and an H2/He one."""
PRESSURE_EXPONENT = 2.0
"""float: The exponent n of the closure: pressure-broadened lines and pressure-induced absorption."""

WATER_CONTINUUM = (
    # nu cm-1, sigma_self (cm2 per molecule at n_h2o = n_ref, 296 K), sigma_foreign (same, n_for = n_ref), self temperature exponent
    (200.0, 2.7e-20, 2.0e-21, 2.76),
    (400.0, 6.7e-21, 7.0e-23, 1.63),
    (600.0, 1.4e-21, 6.8e-24, 3.06),
    (800.0, 3.5e-22, 1.0e-24, 5.06),
    (1000.0, 1.3e-22, 2.4e-25, 5.64),
    (1200.0, 8.4e-23, 5.0e-25, 6.11),
    (2500.0, 8.5e-24, 1.8e-27, 6.21),
    (3600.0, 1.4e-21, 1.2e-22, 2.42),
)
"""tuple: The MT_CKD 4.3 water vapour continuum of part 03 section 3.2 (DATA, computed there from the AER file at 296 K)."""
WATER_CONTINUUM_REFERENCE_K = 296.0
WATER_CONTINUUM_REFERENCE_DENSITY_CM3 = 1013.0e2 / (constants.BOLTZMANN * WATER_CONTINUUM_REFERENCE_K) * 1.0e-6
"""float: `n_ref`, the number density at 1013 mbar and 296 K, per cm3 (the file's reference)."""

CIA_PAIRS = {
    # pair (sorted): ((T K, kP cm5 per molecule2) ..., lowest and highest wavenumber of the file, cm-1)
    ("N2", "N2"): (((100.0, 1.7e-45), (150.0, 6.8e-46), (200.0, 3.9e-46), (250.0, 2.6e-46), (300.0, 1.9e-46), (400.0, 1.5e-46)), (0.0, 5000.0)),
    ("H2", "N2"): (((100.0, 7.6e-45), (150.0, 9.7e-45), (200.0, 9.7e-45), (250.0, 8.8e-45), (300.0, 7.6e-45), (400.0, 5.9e-45)), (0.0, 5100.0)),
    ("H2", "H2"): (((200.0, 3.0e-45), (250.0, 3.1e-45), (300.0, 3.0e-45), (400.0, 2.7e-45), (600.0, 2.4e-45), (800.0, 2.3e-45),
                    (1000.0, 2.8e-45), (1500.0, 3.5e-45), (2000.0, 3.9e-45), (3000.0, 2.6e-45)), (20.0, 10000.0)),
    ("H2", "He"): (((200.0, 1.1e-45), (250.0, 1.1e-45), (300.0, 1.2e-45), (400.0, 1.2e-45), (600.0, 1.1e-45), (800.0, 1.3e-45),
                    (1000.0, 1.9e-45), (1500.0, 3.0e-45), (2000.0, 3.6e-45), (3000.0, 3.2e-45)), (20.0, 20000.0)),
    ("CO2", "CO2"): (((200.0, 4.1e-45), (250.0, 3.0e-45), (300.0, 3.0e-45), (400.0, 3.5e-45), (600.0, 3.4e-45), (800.0, 8.0e-47)), (1.0, 4500.0)),
    ("CO2", "H2"): (((100.0, 6.7e-44), (200.0, 3.8e-44), (300.0, 2.7e-44), (400.0, 2.0e-44), (600.0, 1.2e-44)), (5.0, 1500.0)),
    ("CH4", "CO2"): (((100.0, 1.2e-43), (200.0, 1.8e-44), (300.0, 7.2e-45), (400.0, 3.5e-45), (600.0, 1.5e-45)), (0.0, 720.0)),
    ("CH4", "N2"): (((100.0, 1.4e-44), (150.0, 6.8e-45), (200.0, 4.2e-45), (250.0, 3.0e-45), (300.0, 2.3e-45), (400.0, 1.6e-45)), (0.0, 800.0)),
    ("CH4", "H2"): (((111.0, 1.5e-44), (144.0, 1.9e-44), (186.0, 1.9e-44), (240.0, 1.7e-44), (310.0, 1.5e-44), (400.0, 1.1e-44)), (0.02, 1946.0)),
}
"""dict: The HITRAN collision-induced pairs of part 03 section 3.1 (DATA, computed there): `kP(T)` at the file's temperatures."""

VAPOUR_LAYERS = 96
"""int: Layers of equal depth in `ln p` that the water vapour profile is summed over (a midpoint sum; 96 and 192 agree to within 0.5 percent of
the continuum depth)."""
HUMIDITY_TOP = 0.02
"""float: The fraction of the surface pressure at which the relative humidity reaches zero: Manabe and Wetherald (1967) give
`h(p) = h_s (p / p_s - 0.02) / 0.98`, a humidity that falls linearly with pressure (tier: verified for the Earth's mean state, tweak level L0;
carried to other worlds it is a STAND-IN, GEN.725 row GH1)."""

CLOUD_LONGWAVE_TAU0 = 0.35
"""float: STAND-IN (GEN.725 row GH2): the grey longwave optical depth of the cloud an air holds at the Earth's mean surface humidity
(`EARTH_HUMIDITY`); the cloud scales with the surface humidity and is nil without water vapour. Set so that the Earth's clear-sky greenhouse
(what the bands, wings and continuum give, and the module's 266 W/m2 at 288 K matches the measured clear-sky outgoing flux) loses the
cloud longwave effect of the measured all-sky Earth, 26.5 W/m2 (CERES EBAF; tier: educated guess from the needed output, tweak level L2)."""
EARTH_HUMIDITY = 0.77
"""float: The mean relative humidity of the Earth's surface air (Manabe and Wetherald 1967), the humidity at which `CLOUD_LONGWAVE_TAU0` holds."""

ITERATIONS = 40
"""int: At most this many steps in the surface temperature fixed point."""
UNDAMPED_STEPS = 10
"""int: The steps taken straight (the fixed point contracts by about a tenth a step); later ones average with the last."""
TOLERANCE_K = 1.0e-3
"""float: The fixed point stops when the surface temperature changes by less than this."""


@dataclass(frozen=True)
class Opacity:
    """
    The infrared optical depth of one air at one surface temperature.

    Attributes:
        tau0 (float): The equivalent grey thickness, `-ln(sum w_k exp(-tau_k))`.
        by_part (dict): `tau0` of the air with only that part kept: `"bands"`, `"far_wing"`, `"water_continuum"`, `"collision_induced"`
            (`"clouds"` is the cloud's own `tau0`, which is added to the rest).
        by_gas (dict): `tau0` of the air with only that gas's bands and far wing kept (continuum and pairs are listed under their own part).
        condensable_share (float): The share of the sum of `by_gas` that is H2O, 0 to 1.
        hydrogen_dominated (bool): H2 and He make up more than half of the molecules.
    """

    tau0: float
    by_part: dict
    by_gas: dict
    condensable_share: float
    hydrogen_dominated: bool


@dataclass(frozen=True)
class Vapour:
    """
    The water vapour of an air from the surface up to where the humidity runs out, in `VAPOUR_LAYERS` layers of equal `ln p`.

    Attributes:
        pressure_pa (ndarray): The pressure at the middle of each layer, Pa, falling from the surface.
        log_step (float): The depth of a layer in `ln p`.
        temperature_k (ndarray): The temperature of each layer, K (the troposphere's `T_s (p / p_s)^beta`).
        fraction (ndarray): The mole fraction of water in each layer, at most the surface's.
        surface_fraction (float): The mole fraction of water at the surface.
        cloud_cover (float): The cloud the air holds, as a multiple of the Earth's: the surface humidity over `EARTH_HUMIDITY`, nil where the
            surface is above water's critical temperature (no liquid, no cloud of it).
    """

    pressure_pa: np.ndarray
    log_step: float
    temperature_k: np.ndarray
    fraction: np.ndarray
    surface_fraction: float
    cloud_cover: float


_PLANCK_TERMS = 40
"""int: The terms of the Planck fraction's series for a large `x`."""
_PLANCK_N = np.arange(1, _PLANCK_TERMS + 1).astype(float)[:, None]
_PLANCK_N2 = (np.arange(1, _PLANCK_TERMS + 1) ** 2).astype(float)[:, None]
_PLANCK_N3 = (np.arange(1, _PLANCK_TERMS + 1) ** 3).astype(float)[:, None]
_PLANCK_NEGLIGIBLE = 45.0
"""float: Term `n` of the series is `exp(-(n - 1) x)` times the first at most, so from `(n - 1) x` of 45 on (3e-20) it is under half the
last bit of the sum and adding it changes nothing: a sum that stops there has the bits of one that goes on to the 40th term."""


_PLANCK_STEPS = (2.0, 6.0, 20.0, 50.0)
"""tuple: The `x` from which a series of fewer terms is as good as the 40 (`_PLANCK_NEGLIGIBLE`), the smallest `x` of a group fixing its terms."""


def _planck_sum(large, count):
    """The series of the Planck fraction at `large` (a 1-d array of `x` of 2 and over), its first `count` terms added in order."""
    n, n2, n3 = _PLANCK_N[:count], _PLANCK_N2[:count], _PLANCK_N3[:count]
    terms = np.exp(-n * large) / n * (large ** 3 + 3.0 * large ** 2 / n + 6.0 * large / n2 + 6.0 / n3)
    return np.add.accumulate(terms, axis=0)[-1]


def planck_fraction_below(x):
    """The share of a blackbody's power below the dimensionless frequency `x = h c nu / (k T)`, array in, array out."""
    x = np.asarray(x, dtype=float)
    shape = x.shape
    x = x.reshape(-1)
    result = np.empty_like(x)
    below = x < 2.0
    small = x[below]
    result[below] = 15.0 / math.pi ** 4 * (small ** 3 / 3.0 - small ** 4 / 8.0 + small ** 5 / 60.0 - small ** 7 / 5040.0
                                           + small ** 9 / 272160.0 - small ** 11 / 13305600.0)
    large = x[~below]
    if large.size:
        # The terms are formed at once and added in order (`add.accumulate` adds one after the other where `add.reduce` adds pairwise), the
        # sum having the bits of a loop that adds `n = 1, 2, ...` to a running total (PERF.201); each `x` has only as many terms as it needs.
        total = np.empty_like(large)
        done = np.zeros(large.shape, dtype=bool)
        for least in reversed(_PLANCK_STEPS):
            group = ~done & (large >= least)
            if group.any():
                total[group] = _planck_sum(large[group], 1 + math.ceil(_PLANCK_NEGLIGIBLE / least))
                done |= group
        rest = ~done                                                 # a NaN: no comparison holds, the whole sum is made
        if rest.any():
            total[rest] = _planck_sum(large[rest], _PLANCK_TERMS)
        result[~below] = 1.0 - 15.0 / math.pi ** 4 * total
    return result.reshape(shape)


def planck_weights(temperature_k):
    """The share of a blackbody at `temperature_k` in each grid bin (they add up to 1; the last bin runs to infinity)."""
    below = planck_fraction_below(EDGES_CM1 * SECOND_RADIATION_CM_K / temperature_k)
    return np.concatenate([np.diff(below), [1.0 - below[-1]]])


def _plateaus(formula):
    """`[(low, high, factor)]`, cm-1: where a gas's bands sit, from the gas table (a centre is a zero-width plateau)."""
    row = gases.gas(formula)
    if row.infrared is None or not row.infrared.bands:
        if row.ir == "ir":
            return [(*GENERIC_BAND_CM1, 1.0)]
        return []
    plateaus = []
    for band in row.infrared.bands:
        factor = HOT_BAND_FACTOR if "hot" in band.label else DYAD_BAND_FACTOR if "dyad" in band.label else 1.0
        low, high = (band.wavenumber_cm1, band.wavenumber_cm1) if band.wavenumber_cm1 is not None else band.range_cm1
        plateaus.append((float(low), float(high), factor))
    return plateaus


def _band_depth(formula, column_kg_m2, pressure_scaling):
    """The optical depth of a gas's bands in every bin: the plateau opacity grows with the pressure (line widths do)."""
    kappa0 = (BAND_KAPPA0_M2_KG if gases.gas(formula).ir == "ir" else WEAK_KAPPA0_M2_KG) * pressure_scaling
    depth = np.zeros_like(CENTRES_CM1)
    for factor, falloff in _band_falloffs(formula):
        depth += factor * kappa0 * column_kg_m2 * falloff
    return depth


@functools.lru_cache(maxsize=None)
def _band_falloffs(formula):
    """`[(factor, falloff)]` for a gas's plateaus: the shape each one has over the bins, which depends on neither the column nor the
    pressure (PERF.201: worked out once per gas; the arrays are shared and cannot be written to)."""
    shapes = []
    for low, high, factor in _plateaus(formula):
        distance = np.maximum(np.maximum(low - CENTRES_CM1, CENTRES_CM1 - high), 0.0)
        falloff = np.exp(-distance / WING_E_FOLD_CM1)
        falloff.flags.writeable = False
        shapes.append((factor, falloff))
    return tuple(shapes)


def _wing_depth(formula, column_kg_m2, pressure_scaling):
    """The far-wing optical depth of a gas in every bin: `kappa_w u (P / P_REF)` within `WING_REACH_CM1 (P / P_REF)` of a band."""
    kappa = far_wing_m2_kg(formula)
    depth = np.zeros_like(CENTRES_CM1)
    if kappa <= 0.0:
        return depth
    reach = WING_REACH_CM1 * pressure_scaling
    for low, high, _factor in _plateaus(formula):
        near = (CENTRES_CM1 >= low - reach) & (CENTRES_CM1 <= high + reach)
        depth[near] = kappa * column_kg_m2 * pressure_scaling
    return depth


def far_wing_m2_kg(formula):
    """The far-wing mass opacity of a gas at one bar, m2/kg: its entry in `FAR_WING_M2_KG`, else the table's flag (`"ir"` generic,
    `"weak"` thin, anything else none)."""
    if formula in FAR_WING_M2_KG:
        return FAR_WING_M2_KG[formula]
    flag = gases.gas(formula).ir
    if flag == "ir":
        return GENERIC_FAR_WING_M2_KG
    if flag == "weak":
        return WEAK_FAR_WING_M2_KG
    return 0.0


def vapour_profile(partial_pressures_pa, temperature_k):
    """
    The water vapour of an air by height, or `None` for an air with none.

    The relative humidity falls from the surface's with pressure (Manabe and Wetherald 1967, `HUMIDITY_TOP`), the temperature falls as
    the troposphere's `T_s (p / p_s)^beta` (`beta = alpha (gamma - 1) / gamma`, `alpha` that of an air whose opacity is condensable
    vapour, `LAPSE_SCALING_CONDENSABLE`), and the vapour stands at that humidity of its saturation pressure at the layer's
    temperature, never above the surface's mole fraction. The surface humidity is the surface vapour pressure over its saturation
    pressure at `temperature_k` (1 above the critical temperature, where there is none, where the saturation pressure has fallen to nothing, and for a surface that holds more than it can).
    An air whose vapour is well under saturation keeps its mixing ratio until the humidity or the cold takes it; one at the
    saturation limit loses it as fast as the saturation pressure falls.
    """
    if partial_pressures_pa.get("H2O", 0.0) <= 0.0:
        return None
    gamma = mixture.heat_capacity_ratio(partial_pressures_pa, temperature_k)
    total = mixture.total_pressure_pa(partial_pressures_pa)
    surface_fraction = partial_pressures_pa["H2O"] / total
    surface_saturation = saturation.pressure_pa("H2O", temperature_k)
    humidity = min(partial_pressures_pa["H2O"] / surface_saturation, 1.0) if 0.0 < surface_saturation < math.inf else 1.0
    step = -math.log(HUMIDITY_TOP) / VAPOUR_LAYERS
    ratio = np.exp(-(np.arange(VAPOUR_LAYERS) + 0.5) * step)
    exponent = LAPSE_SCALING_CONDENSABLE * (gamma - 1.0) / gamma
    layer_k = temperature_k * ratio ** exponent
    held = np.array([humidity * (r - HUMIDITY_TOP) / (1.0 - HUMIDITY_TOP) * saturation.pressure_pa("H2O", t) / (total * r)
                     for r, t in zip(ratio, layer_k)])
    return Vapour(pressure_pa=total * ratio, log_step=step, temperature_k=layer_k,
                  fraction=np.minimum(held, surface_fraction), surface_fraction=surface_fraction,
                  cloud_cover=humidity / EARTH_HUMIDITY if math.isfinite(surface_saturation) else 0.0)


def _vapour_pressure_weight(vapour):
    """How much of its column a gas that stood at the surface's mole fraction all the way up would have to pass for the vapour's
    pressure-broadened opacity to be the same: `sum q p^2 / (q_s sum p^2)` over the layers (1 for a well-mixed gas)."""
    return float(np.sum(vapour.fraction * vapour.pressure_pa ** 2) / (vapour.surface_fraction * np.sum(vapour.pressure_pa ** 2)))


def _interpolate(points, x):
    """Linear in x of `[(x, y)]`, flat outside the points."""
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return float(np.interp(x, xs, ys))


def pair_coefficient(pair, temperature_k):
    """`kP` of a collision-induced pair, cm5 per molecule2, at `temperature_k` (linear between the file's temperatures, the nearest
    one outside them), or `None` for a pair the file does not hold."""
    entry = CIA_PAIRS.get(tuple(sorted(pair)))
    if entry is None:
        return None
    return _interpolate(entry[0], temperature_k)


@functools.lru_cache(maxsize=1)
def _continuum_on_grid():
    """The water continuum's table on the grid bins it covers, which do not depend on the temperature: the mask of those bins, the
    self coefficient at the reference temperature, its temperature exponent and the foreign coefficient (PERF.201: worked out once)."""
    nus = [row[0] for row in WATER_CONTINUUM]
    log_self = [math.log(row[1]) for row in WATER_CONTINUUM]
    log_for = [math.log(row[2]) for row in WATER_CONTINUUM]
    exponent = [row[3] for row in WATER_CONTINUUM]
    inside = (CENTRES_CM1 >= nus[0]) & (CENTRES_CM1 <= nus[-1])
    nu = CENTRES_CM1[inside]
    return inside, np.exp(np.interp(nu, nus, log_self)), np.interp(nu, nus, exponent), np.exp(np.interp(nu, nus, log_for))


def _water_continuum_depth(vapour, molar_mass_kg_mol, gravity_m_s2):
    """The continuum's optical depth in every bin, summed over the layers of the vapour profile (zero without water vapour)."""
    depth = np.zeros_like(CENTRES_CM1)
    if vapour is None:
        return depth
    q = vapour.fraction[:, None]
    density_cm3 = (vapour.pressure_pa / (constants.BOLTZMANN * vapour.temperature_k) * 1.0e-6)[:, None]
    height_cm = (constants.R * vapour.temperature_k / (molar_mass_kg_mol * gravity_m_s2) * vapour.log_step * 100.0)[:, None]
    inside, sigma_self_reference, exponent, sigma_for = _continuum_on_grid()
    cold = (WATER_CONTINUUM_REFERENCE_K / vapour.temperature_k)[:, None]
    sigma_self = sigma_self_reference * cold ** exponent
    # tau = (sigma_self (n_w / n_ref) + sigma_for (n_for / n_ref)) n_w dz per layer: n_w = q n, n_for = (1 - q) n, dz = H d(ln p)
    reference = WATER_CONTINUUM_REFERENCE_DENSITY_CM3
    depth[inside] = (((sigma_self * q + sigma_for * (1.0 - q)) * q * density_cm3 ** 2 / reference) * height_cm).sum(axis=0)
    return depth


def _pair_depth(fractions, temperature_k, pressure_pa, scale_height_m):
    """The collision-induced pairs' optical depth in every bin: the pairs' thin depths added over the range of each file, then `ln(1 + x)`."""
    depth = np.zeros_like(CENTRES_CM1)
    density_cm3 = pressure_pa / (constants.BOLTZMANN * temperature_k) * 1.0e-6
    height_cm = scale_height_m * 100.0
    for pair, (_points, (low, high)) in CIA_PAIRS.items():
        x1 = fractions.get(pair[0], 0.0)
        x2 = fractions.get(pair[1], 0.0)
        if x1 <= 0.0 or x2 <= 0.0:
            continue
        # the HITRAN convention: kP n1 n2 per unit length, and a pair of unlike gases counts both orders once
        thin = pair_coefficient(pair, temperature_k) * x1 * x2 * density_cm3 ** 2 * height_cm / 2.0
        depth[(CENTRES_CM1 >= low) & (CENTRES_CM1 <= high)] += thin
    return np.log1p(depth)


def _equivalent(weights, depth):
    """`-ln(sum w exp(-tau))`, taken about the thinnest bin that has weight so a very thick air cannot underflow."""
    present = weights > 0.0
    floor = float(np.min(depth[present]))
    return max(floor - math.log(float(np.sum(weights[present] * np.exp(-(depth[present] - floor))))), 0.0)


def optical_depth(partial_pressures_pa, temperature_k, gravity_m_s2):
    """
    The grey infrared optical depth of an air at a surface temperature.

    Raises:
        ValueError: The mixture is empty, or the temperature or gravity is not positive.
    """
    if not (temperature_k > 0.0 and math.isfinite(temperature_k)) or not (gravity_m_s2 > 0.0 and math.isfinite(gravity_m_s2)):
        raise ValueError("temperature_k and gravity_m_s2 must be positive numbers")
    total = mixture.total_pressure_pa(partial_pressures_pa)
    fractions = mixture.mole_fractions(partial_pressures_pa)
    molar_mass = mixture.mean_molar_mass_kg_mol(partial_pressures_pa)
    scale_height = constants.R * temperature_k / (molar_mass * gravity_m_s2)
    weights = planck_weights(temperature_k)
    columns = {formula: fraction * gases.gas(formula).molar_mass * total / (molar_mass * gravity_m_s2)
               for formula, fraction in fractions.items()}
    vapour = vapour_profile(partial_pressures_pa, temperature_k)
    if vapour is not None:
        # the vapour's bands and wings see the pressure of the layer it stands in, so its column counts as much as the well-mixed one that
        # gives the same pressure-weighted opacity
        columns["H2O"] *= _vapour_pressure_weight(vapour)
    scaling = total / PRESSURE_REFERENCE_PA

    bands = np.zeros_like(CENTRES_CM1)
    wing = np.zeros_like(CENTRES_CM1)
    by_gas = {}
    for formula, column in columns.items():
        gas_bands = _band_depth(formula, column, scaling)
        gas_wing = _wing_depth(formula, column, scaling)
        if gas_wing.any() or gas_bands.any():
            by_gas[formula] = _equivalent(weights, gas_bands + gas_wing)
        bands += gas_bands
        wing += gas_wing
    continuum = _water_continuum_depth(vapour, molar_mass, gravity_m_s2)
    pairs = _pair_depth(fractions, temperature_k, total, scale_height)
    clouds = CLOUD_LONGWAVE_TAU0 * vapour.cloud_cover if vapour is not None else 0.0
    by_part = {"bands": _equivalent(weights, bands), "far_wing": _equivalent(weights, wing), "water_continuum": _equivalent(weights, continuum),
               "collision_induced": _equivalent(weights, pairs), "clouds": clouds}
    tau0 = _equivalent(weights, bands + wing + continuum + pairs) + clouds
    absorber = sum(by_gas.values())
    share = by_gas.get("H2O", 0.0) / absorber if absorber > 0.0 else 0.0
    hydrogen = fractions.get("H2", 0.0) + fractions.get("He", 0.0) > 0.5
    return Opacity(tau0=tau0, by_part=by_part, by_gas=by_gas, condensable_share=share, hydrogen_dominated=hydrogen)


def lapse_scaling(opacity):
    """Robinson and Catling's alpha for an air: 0.6 when condensable vapour supplies all the opacity, 0.8 when none does, 0.85 for H2/He."""
    if opacity.hydrogen_dominated:
        return LAPSE_SCALING_HYDROGEN
    return LAPSE_SCALING_DRY + (LAPSE_SCALING_CONDENSABLE - LAPSE_SCALING_DRY) * opacity.condensable_share
