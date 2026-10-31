# planetgen/physics/line_opacity.py

"""
The infrared optical depth of a layer of air, band by band, from line-list tables (GEN.715, GEN.717, GEN.718;
docs/design/composition/line-opacity.md; docs/design/composition/atmosphere-composition-parts/03-infrared-greenhouse-thermal.md
sections 2, 3 and 6)
================================================================================================================================

GEN.703 built one grey optical depth per air from band centres, with the strength of each band a stand-in. This module replaces
that: the optical depth of one layer (a pressure, a temperature and the molecules above one square centimetre) on 130 wavenumber
bands from 0 to 12000 cm-1, each band divided into 10 g points of a correlated-k distribution, so a band that is black in its
line centres and open between them keeps its windows. `radiative_convective` uses it layer by layer.

* **Line absorption** of water vapour, carbon dioxide and methane: the ExoMolOP correlated-k tables (`datasets/ir_ktables.py`,
  POKAZATEL, UCL-4000 and YT34to10 line lists), read at the layer's temperature (linear in 1/T) and pressure (linear in ln P) and
  multiplied by the molecules of the gas. Gases overlap by the random-overlap method (`overlap`): the k-distributions are
  combined as if the lines of different gases were uncorrelated, then sorted again by strength.
* **The water vapour continuum** (MT_CKD 4.3, part 03 section 3.2, `WATER_CONTINUUM`): self and foreign coefficients that close
  the 8 to 12 micron window of a humid air.
* **Collision-induced absorption** (`datasets/hitran_cia.py`): `k(T) n1 n2` per unit length for the pairs the HITRAN files hold
  (N2-N2, O2-O2, N2-O2, CO2-CO2, N2-CH4, N2-H2, CO2-CH4, CO2-H2, H2-H2, CH4-H2, H2-He), the mean over each band.
* **Gases without a table** (stand-ins, `STANDIN_*` constants, registry rows IR1 and IR2 of the gap plan): the gas table lists
  the centres of the bands of every infrared-active gas and nothing more, so such a gas gets a plateau of opacity at each centre
  with exponential wings (the form of GEN.703), turned into a k-distribution by sorting the opacity inside a band. The plateau
  strength and the wing reach are class values, not measurements.
* **The far wings of the gases** (`physics/co2_far_wing.py`, GEN.1019, GEN.1020): the tables cut every line at 25 cm-1 from its centre, so the
  opacity of a dense carbon dioxide air between its bands is missing (Venus' surface is 92 bar of carbon dioxide). The wing of every
  carbon dioxide band (Perrin and Hartmann's sub-Lorentzian factor) and the wings of the minor gases the tables hold (water, sulfur dioxide,
  carbon monoxide) in a carbon dioxide air are added, with the band strengths read off each table at its lowest pressure. It is skipped
  where the carbon dioxide partial pressure is under `WING_MIN_CO2_PA` (it goes as that pressure times the total).
* **The sulfuric acid cloud** (`physics/cloud_infrared.py`, GEN.1020): `radiative_convective` adds the cloud's infrared depth to the
  layers of the column; this module is not involved.

If a table dataset is missing and the run disregards missing data, the gases it would hold take the stand-in path and the pairs
are left out (`fetch.stand_ins_used` lists it); any other run raises `fetch.DataMissing`.

Use:

    prepared = line_opacity.prepare(partial_pressures_pa)
    tau = line_opacity.layer_depth(prepared, temperature_k, pressure_pa, water_fraction, molecules_cm2)    # (layer, band, g)
"""

import functools
import math
from dataclasses import dataclass, field

import numpy as np

from planetgen.datasets import fetch, hitran_cia, ir_ktables
from planetgen.physics import co2_far_wing, constants, gases, mixture
from planetgen.physics.air_column import AVOGADRO

SECOND_RADIATION_CM_K = (2.0 * math.pi * constants.REDUCED_PLANCK * constants.SPEED_OF_LIGHT_M_S / constants.BOLTZMANN) * 100.0
"""float: `h c / k` in cm K, the Planck exponent `x = c2 nu / T` for a wavenumber in cm-1."""

EDGES_CM1 = hitran_cia.EDGES_CM1
"""ndarray: The 131 edges, cm-1, of the 130 bands (0 to 200 cm-1, then the k-tables' 129)."""
CENTRES_CM1 = 0.5 * (EDGES_CM1[1:] + EDGES_CM1[:-1])
"""ndarray: The band centres, cm-1."""
BANDS = len(CENTRES_CM1)
G_POINTS, G_WEIGHTS, G_EDGES = ir_ktables.g_grid()
"""ndarray: The g points of the k-distribution, their weights (adding to 1) and the 11 edges of the intervals."""
GRID_POINTS = len(G_POINTS)

PRESSURE_REFERENCE_PA = 1.0e5
"""float: `P_REF`, one bar: the pressure the stand-in opacities are stated at."""

TABLE_GASES = tuple(ir_ktables.GASES)
"""tuple: The gases with a correlated-k table."""

OVERLAP_SMALL = 0.05
"""float: A gas whose largest optical depth in a band is under this adds only its mean to the other's distribution (no re-sort)."""

STANDIN_SUBSAMPLES = 25
"""int: How many wavenumbers a band is sampled at to turn a stand-in opacity into a distribution."""
STANDIN_WING_E_FOLD_CM1 = 10.0
"""float: STAND-IN (IR1, tier T6, L1): the e-folding distance of a band's wing. Set in GEN.703 from the logarithm of the carbon dioxide
forcing law (Byrne and Goldblatt 2014a, `F = 5.32 ln(C/C0)`): about 19 cm-1 of saturated width per e-fold of the column, 10 cm-1 on each side."""
STANDIN_KAPPA0_M2_KG = 1.0e3
"""float: STAND-IN (IR1, tier T7, L0): the plateau opacity of a fundamental band of a gas the table flags `"ir"`, m2/kg at one bar (a strong
band; it enters the width only as a logarithm)."""
STANDIN_WEAK_KAPPA0_M2_KG = 1.0e1
"""float: STAND-IN (IR1, tier T7, L0): the same for a gas the table flags `"weak"`."""
STANDIN_HOT_BAND_FACTOR = 1.0e-6
"""float: STAND-IN (IR1, tier T7, L1): a hot band (label contains `hot`) counts this much of a fundamental's plateau."""
STANDIN_DYAD_BAND_FACTOR = 1.0e-5
"""float: STAND-IN (IR1, tier T7, L1): a Fermi dyad (label contains `dyad`) counts this much of a fundamental's plateau."""
STANDIN_GENERIC_BAND_CM1 = (600.0, 1500.0)
"""tuple: STAND-IN (IR1, tier T7, L0): the plateau of a gas flagged `"ir"` whose table lists no band."""
STANDIN_FAR_WING_M2_KG = {"CO": 1.0e-5}
"""dict: STAND-IN (IR2, tier T5, L0): the far-wing mass opacity at one bar of a gas without a table, m2/kg, read from Koll 2019's statement (CO
thin even at one bar). A gas flagged `"ir"` without an entry takes `STANDIN_GENERIC_WING_M2_KG`, one flagged `"weak"`
`STANDIN_WEAK_WING_M2_KG`."""
STANDIN_GENERIC_WING_M2_KG = 1.0e-4
STANDIN_WEAK_WING_M2_KG = 1.0e-5
STANDIN_WING_REACH_CM1 = 600.0
"""float: STAND-IN (IR2, tier T7, L0): how far from a band the far wing reaches at one bar; it grows with the pressure (line widths do)."""

WING_MIN_CO2_PA = 1.0e3
"""float: A layer whose carbon dioxide partial pressure is under this adds no far wing (a numerical floor: the wing goes as that pressure times
the total, 1e-9 of Venus's at Earth's air; it saves the cost of the wing for every layer of a thin air)."""

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
"""tuple: The MT_CKD 4.3 water vapour continuum of part 03 section 3.2 (DATA, computed there from the AER file at 296 K), eight nodes (IR4)."""
WATER_CONTINUUM_REFERENCE_K = 296.0
WATER_CONTINUUM_REFERENCE_DENSITY_CM3 = 1013.0e2 / (constants.BOLTZMANN * WATER_CONTINUUM_REFERENCE_K) * 1.0e-6
"""float: `n_ref`, the number density at 1013 mbar and 296 K, per cm3 (the file's reference)."""

# STAND-IN (IR5, tier T6, L1): the k-tables start at 200 cm-1 (50 micron); the band from 0 to 200 cm-1 takes the first band's distribution
# (the water rotation band is as strong below 200 cm-1 as it is just above).


@dataclass(frozen=True)
class KTable:
    """One gas's reduced table: `logk` is log10 of k (cm2 per molecule), (band, temperature, pressure, g)."""

    logk: np.ndarray
    inverse_temperature: np.ndarray
    log_pressure_first: float
    log_pressure_step: float
    n_pressures: int


@dataclass(frozen=True)
class Prepared:
    """
    An air prepared for `layer_depth`.

    Attributes:
        dry (dict): Each gas's mole fraction of the dry air (water vapour left out; they add up to 1).
        has_water (bool): The mixture holds water vapour.
        tabled (tuple): The gases (water included when present) that have a k-table.
        standin (tuple): The infrared-active gases without a table (the stand-in path).
        pairs (tuple): `(key, first, second)` of the collision-induced pairs the mixture holds.
        absorbers (bool): Anything in the air absorbs in the infrared at all.
        tables (dict): The tables of `tabled`, formula to `KTable`.
        table_temperatures (ndarray): Their temperature axis, K.
        pair_tables (dict): The pairs' `{key: (temperatures, k)}`.
    """

    dry: dict
    has_water: bool
    tabled: tuple
    standin: tuple
    pairs: tuple
    absorbers: bool
    tables: dict = field(default_factory=dict, compare=False)
    table_temperatures: np.ndarray = field(default=None, compare=False)
    pair_tables: dict = field(default_factory=dict, compare=False)


def planck_fraction_below(x):
    """The share of a blackbody's power below the dimensionless frequency `x = h c nu / (k T)`, array in, array out."""
    x = np.asarray(x, dtype=float)
    small = np.minimum(x, 2.0)
    series = 15.0 / math.pi ** 4 * (small ** 3 / 3.0 - small ** 4 / 8.0 + small ** 5 / 60.0 - small ** 7 / 5040.0
                                    + small ** 9 / 272160.0 - small ** 11 / 13305600.0)
    large = np.maximum(x, 2.0)
    total = np.zeros_like(large)
    for n in range(1, 41):
        total += np.exp(-n * large) / n * (large ** 3 + 3.0 * large ** 2 / n + 6.0 * large / n ** 2 + 6.0 / n ** 3)
    return np.where(x < 2.0, series, 1.0 - 15.0 / math.pi ** 4 * total)


_PLANCK_X = np.linspace(0.0, 60.0, 60001)
_PLANCK_FRACTION = planck_fraction_below(_PLANCK_X)


_EMISSION_STEP_K = 0.25
_EMISSION_FIRST_K = 5.0
_EMISSION_TOP_K = 4100.0


@functools.lru_cache(maxsize=1)
def _emission_table():
    """The blackbody power in each band on a fine temperature grid, W/m2, (temperature, band), so a lookup is a gather and one linear interpolation."""
    grid = np.arange(_EMISSION_FIRST_K, _EMISSION_TOP_K + _EMISSION_STEP_K, _EMISSION_STEP_K)
    x = SECOND_RADIATION_CM_K * EDGES_CM1[None, :] / grid[:, None]
    return constants.STEFAN_BOLTZMANN_CONSTANT * grid[:, None] ** 4 * np.diff(np.interp(x, _PLANCK_X, _PLANCK_FRACTION), axis=1)


def band_emission(temperature_k):
    """
    The power a blackbody radiates in each of the 130 bands, W/m2: `temperature_k` is an array of any shape, the result has a last axis of
    bands. Linear between the table's temperatures (every quarter kelvin from 5 to 4100 K; the ends hold outside).
    """
    table = _emission_table()
    position = np.clip((np.asarray(temperature_k, dtype=float) - _EMISSION_FIRST_K) / _EMISSION_STEP_K, 0.0, len(table) - 1 - 1.0e-9)
    index = position.astype(int)
    weight = (position - index)[..., None]
    return table[index] * (1.0 - weight) + table[index + 1] * weight


def reduce_distribution(k, weights, edges):
    """A k-distribution reduced to the g `edges`: `k` (..., n) with (n,) `weights` adding to 1, sorted by k, cut at `edges`, each interval's
    mean (..., len(edges) - 1)."""
    return ir_ktables.reduce_distribution(k, weights, edges)


# ---- the k-tables ----

def _tabled():
    """`(tables, temperatures)` from the downloaded k-tables, or `({}, None)` when the dataset is missing and the run disregards that."""
    located = fetch.locate(ir_ktables.DATASET)
    if located.use_stand_in("infrared line opacity from the plateau-and-wing band stand-ins"):
        return {}, None
    return _build_tables(located.path(ir_ktables.DERIVED))


@functools.lru_cache(maxsize=2)
def _build_tables(path):
    stored = ir_ktables.load(path)
    if not np.array_equal(stored["edges_cm1"], ir_ktables.BAND_EDGES_CM1):
        raise ValueError("the downloaded k-tables were made for other bands than this module's; run the setup step again")
    pressure = stored["pressure_bar"]
    steps = np.diff(np.log(pressure))
    if not np.allclose(steps, steps[0]):
        raise ValueError("the k-table pressures must be equally spaced in ln P")
    temperature = stored["temperature_k"]
    tables = {}
    for gas in TABLE_GASES:
        logk = np.log10(np.maximum(stored[f"k_{gas}"], 1.0e-45)).astype(np.float32)
        tables[gas] = KTable(logk, 1.0 / temperature, float(np.log(pressure[0])), float(steps[0]), len(pressure))
    return tables, temperature


def _lookup(table, temperatures, temperature, pressure_bar):
    """k (cm2 per molecule) of one gas at layer temperatures and pressures: shape (layer, band, g), the 0 to 200 cm-1 band first (the first band's)."""
    temperature = np.clip(np.asarray(temperature, dtype=float), temperatures[0], temperatures[-1])
    index = np.clip(np.searchsorted(temperatures, temperature) - 1, 0, len(temperatures) - 2)
    inverse = table.inverse_temperature
    weight_t = np.clip((1.0 / temperature - inverse[index + 1]) / (inverse[index] - inverse[index + 1]), 0.0, 1.0)
    position = np.clip((np.log(np.maximum(pressure_bar, 1.0e-30)) - table.log_pressure_first) / table.log_pressure_step,
                       0.0, table.n_pressures - 1 - 1.0e-9)
    ip = position.astype(int)
    weight_p = position - ip
    logk = table.logk
    mix = (logk[:, index, ip, :] * (weight_t * (1 - weight_p))[None, :, None]
           + logk[:, index + 1, ip, :] * ((1 - weight_t) * (1 - weight_p))[None, :, None]
           + logk[:, index, ip + 1, :] * (weight_t * weight_p)[None, :, None]
           + logk[:, index + 1, ip + 1, :] * ((1 - weight_t) * weight_p)[None, :, None])
    k = (10.0 ** mix.astype(float)).transpose(1, 0, 2)
    return np.concatenate([k[:, :1], k], axis=1)


# ---- collision-induced pairs ----

def _pair_tables():
    """`{pair key: (temperatures, k)}` from the downloaded files, or `{}` when missing and disregarded."""
    located = fetch.locate(hitran_cia.DATASET)
    if located.use_stand_in("collision-induced absorption left out"):
        return {}
    return _build_pairs(located.path(hitran_cia.DERIVED))


@functools.lru_cache(maxsize=2)
def _build_pairs(path):
    stored = hitran_cia.load(path)
    return {pair: (stored[f"temperature_{pair}"], stored[f"k_{pair}"]) for pair in hitran_cia.FILES}


def _fraction(prepared, formula, water_fraction):
    if formula == "H2O":
        return water_fraction
    return prepared.dry.get(formula, 0.0) * (1.0 - water_fraction)


def _pair_depth(prepared, temperature, water_fraction, molecules_cm2, density_cm3):
    """Collision-induced optical depth (layer, band): `k(T) x1 x2 n N` summed over the pairs of the air (`n` the number density, `N` the molecules
    above a square centimetre in the layer)."""
    depth = np.zeros((len(temperature), BANDS))
    for key, first, second in prepared.pairs:
        weight = _fraction(prepared, first, water_fraction) * _fraction(prepared, second, water_fraction) * density_cm3 * molecules_cm2
        if not np.any(weight > 0.0):
            continue
        temperatures, k = prepared.pair_tables[key]
        coefficient = np.empty((len(temperature), BANDS))
        for band in range(BANDS):
            coefficient[:, band] = np.interp(temperature, temperatures, k[:, band])
        depth += coefficient * weight[:, None]
    return depth


# ---- water continuum ----

def _continuum_depth(temperature, pressure_pa, water_fraction, molecules_cm2):
    """The MT_CKD continuum's optical depth (layer, band); zero outside 200 to 3600 cm-1."""
    nus = np.array([row[0] for row in WATER_CONTINUUM])
    log_self = np.log([row[1] for row in WATER_CONTINUUM])
    log_foreign = np.log([row[2] for row in WATER_CONTINUUM])
    exponent = np.array([row[3] for row in WATER_CONTINUUM])
    inside = (CENTRES_CM1 >= nus[0]) & (CENTRES_CM1 <= nus[-1])
    nu = np.where(inside, CENTRES_CM1, nus[0])
    sigma_self = (np.exp(np.interp(nu, nus, log_self))[None, :]
                  * (WATER_CONTINUUM_REFERENCE_K / temperature[:, None]) ** np.interp(nu, nus, exponent)[None, :])
    sigma_foreign = np.exp(np.interp(nu, nus, log_foreign))[None, :]
    density = (pressure_pa / (constants.BOLTZMANN * temperature) * 1.0e-6)[:, None]
    water = (water_fraction * molecules_cm2)[:, None]
    depth = (water * (sigma_self * water_fraction[:, None] + sigma_foreign * (1.0 - water_fraction[:, None]))
             * density / WATER_CONTINUUM_REFERENCE_DENSITY_CM3)
    return np.where(inside[None, :], depth, 0.0)


# ---- stand-in bands ----

def _plateaus(formula):
    """`[(low, high, factor)]`, cm-1: where a gas's bands sit, from the gas table (a centre is a zero-width plateau)."""
    row = gases.gas(formula)
    if row.infrared is None or not row.infrared.bands:
        if row.ir == "ir":
            return [(*STANDIN_GENERIC_BAND_CM1, 1.0)]
        return []
    plateaus = []
    for band in row.infrared.bands:
        factor = STANDIN_HOT_BAND_FACTOR if "hot" in band.label else STANDIN_DYAD_BAND_FACTOR if "dyad" in band.label else 1.0
        low, high = (band.wavenumber_cm1, band.wavenumber_cm1) if band.wavenumber_cm1 is not None else band.range_cm1
        plateaus.append((float(low), float(high), factor))
    return plateaus


def _wing_mass_opacity(formula):
    """The far-wing mass opacity of a gas without a table at one bar, m2/kg."""
    if formula in STANDIN_FAR_WING_M2_KG:
        return STANDIN_FAR_WING_M2_KG[formula]
    flag = gases.gas(formula).ir
    if flag == "ir":
        return STANDIN_GENERIC_WING_M2_KG
    if flag == "weak":
        return STANDIN_WEAK_WING_M2_KG
    return 0.0


@functools.cache
def _subsamples():
    """The wavenumbers a stand-in band is sampled at, (band, sub-sample), cm-1."""
    fraction = (np.arange(STANDIN_SUBSAMPLES) + 0.5) / STANDIN_SUBSAMPLES
    return EDGES_CM1[:-1, None] + np.diff(EDGES_CM1)[:, None] * fraction[None, :]


def _standin_distribution(formula, pressure_pa):
    """A stand-in gas's mass opacity (layer, band, g), m2/kg: the plateaus (growing with the pressure) and far wings, sorted into a distribution."""
    scaling = np.asarray(pressure_pa, dtype=float) / PRESSURE_REFERENCE_PA
    nu = _subsamples()
    kappa0 = STANDIN_KAPPA0_M2_KG if gases.gas(formula).ir == "ir" else STANDIN_WEAK_KAPPA0_M2_KG
    wing = _wing_mass_opacity(formula)
    opacity = np.zeros((len(scaling), BANDS, STANDIN_SUBSAMPLES))
    for low, high, factor in _plateaus(formula):
        distance = np.maximum(np.maximum(low - nu, nu - high), 0.0)
        opacity += (factor * kappa0 * np.exp(-distance / STANDIN_WING_E_FOLD_CM1))[None] * scaling[:, None, None]
        if wing > 0.0:
            reach = STANDIN_WING_REACH_CM1 * scaling
            near = (nu[None] >= low - reach[:, None, None]) & (nu[None] <= high + reach[:, None, None])
            opacity += np.where(near, wing * scaling[:, None, None], 0.0)
    return reduce_distribution(opacity, np.full(STANDIN_SUBSAMPLES, 1.0 / STANDIN_SUBSAMPLES), G_EDGES)


def _wing_depth(prepared, temperature, pressure, water, molecules):
    """The far wings of the carbon dioxide and the minor gases in a carbon dioxide air (layer, band), grey across the g points."""
    depth = np.zeros((len(temperature), BANDS))
    co2 = prepared.dry.get("CO2", 0.0) * (1.0 - water)
    active = (co2 * pressure >= WING_MIN_CO2_PA) & (co2 > 0.0)
    if "CO2" not in prepared.tabled or not active.any():
        return depth
    for gas in ("CO2", *(g for g in prepared.tabled if g in co2_far_wing.number("minor_gamma_cm1_atm"))):
        x = water if gas == "H2O" else prepared.dry.get(gas, 0.0) * (1.0 - water)
        mask = active & (x > 0.0)
        if not mask.any():
            continue
        k = _lookup(prepared.tables[gas], prepared.table_temperatures, temperature[mask], np.full(int(mask.sum()), co2_far_wing.REFERENCE_PRESSURE_BAR))
        strengths = co2_far_wing.band_strengths((k * G_WEIGHTS).sum(axis=-1), EDGES_CM1)
        if gas == "CO2":
            cross = co2_far_wing.cross_section_cm2(strengths, EDGES_CM1, temperature[mask], pressure[mask], co2[mask])
        else:
            cross = co2_far_wing.minor_gas_cross_section_cm2(gas, strengths, EDGES_CM1, temperature[mask], pressure[mask], co2[mask])
        depth[mask] += cross * (x[mask] * molecules[mask])[:, None]
    return depth


# ---- overlap ----

def resort(first, second):
    """The random overlap of two optical-depth distributions (..., g): all `g x g` sums with their joint weights, sorted, cut at the g edges."""
    count = first.shape[-1]
    combined = (first[..., :, None] + second[..., None, :]).reshape(first.shape[:-1] + (count * count,))
    weights = (G_WEIGHTS[:, None] * G_WEIGHTS[None, :]).reshape(-1)
    return reduce_distribution(combined, weights, G_EDGES)


def overlap(first, second):
    """The optical depth (layer, band, g) of two gases that absorb together, by random overlap: where one is thin it adds its mean."""
    mean_first = (first * G_WEIGHTS).sum(-1)
    mean_second = (second * G_WEIGHTS).sum(-1)
    second_small = second[..., -1] < OVERLAP_SMALL
    first_small = (first[..., -1] < OVERLAP_SMALL) & ~second_small
    result = np.where(second_small[..., None], first + mean_second[..., None],
                      np.where(first_small[..., None], second + mean_first[..., None], 0.0))
    both = ~(second_small | first_small)
    if both.any():
        result[both] = resort(first[both], second[both])
    return result


# ---- the air and its layers ----

def prepare(partial_pressures_pa):
    """
    Prepares an air for `layer_depth`.

    Raises:
        ValueError: The mixture is empty.
        fetch.DataMissing: A table dataset is missing and the run does not disregard that.
    """
    fractions = mixture.mole_fractions(partial_pressures_pa)
    water = fractions.get("H2O", 0.0)
    dry_total = 1.0 - water
    dry = {formula: fraction / dry_total for formula, fraction in fractions.items() if formula != "H2O"} if dry_total > 0.0 else {}
    tables, temperatures = _tabled()
    tabled = tuple(formula for formula in TABLE_GASES if tables and fractions.get(formula, 0.0) > 0.0)
    standin = tuple(formula for formula, fraction in fractions.items()
                    if formula not in tabled and fraction > 0.0 and gases.gas(formula).is_infrared_active)
    pair_tables = _pair_tables()
    pairs = tuple((key, *key.split("-")) for key in pair_tables if all(fractions.get(formula, 0.0) > 0.0 for formula in key.split("-")))
    return Prepared(dry=dry, has_water=water > 0.0, tabled=tabled, standin=standin, pairs=pairs,
                    absorbers=bool(tabled or standin or pairs), tables=tables, table_temperatures=temperatures,
                    pair_tables=pair_tables)


def layer_depth(prepared, temperature_k, pressure_pa, water_fraction, molecules_cm2):
    """
    The infrared optical depth of a stack of layers on the bands and g points.

    Args:
        prepared (Prepared): From `prepare`.
        temperature_k (ndarray): Each layer's temperature (layers top to bottom), K.
        pressure_pa (ndarray): Each layer's pressure, Pa.
        water_fraction (ndarray): Each layer's mole fraction of water vapour (the dry gases share the rest in `prepared.dry`'s proportions).
        molecules_cm2 (ndarray): The molecules above one square centimetre in each layer.

    Returns:
        ndarray: Optical depth (layer, 130 bands, 10 g points), the g points sorted by strength.
    """
    temperature = np.asarray(temperature_k, dtype=float)
    pressure = np.asarray(pressure_pa, dtype=float)
    water = np.asarray(water_fraction, dtype=float)
    molecules = np.asarray(molecules_cm2, dtype=float)
    layers = len(temperature)
    if not prepared.absorbers:
        return np.zeros((layers, BANDS, GRID_POINTS))
    parts = []
    for gas in prepared.tabled:
        x = water if gas == "H2O" else prepared.dry.get(gas, 0.0) * (1.0 - water)
        if np.any(x > 0.0):
            k = _lookup(prepared.tables[gas], prepared.table_temperatures, temperature, pressure / 1.0e5)
            parts.append(k * (x * molecules)[:, None, None])
    for gas in prepared.standin:
        x = water if gas == "H2O" else prepared.dry.get(gas, 0.0) * (1.0 - water)
        mass_kg_m2 = x * molecules * 1.0e4 * gases.gas(gas).molar_mass / AVOGADRO
        parts.append(_standin_distribution(gas, pressure) * mass_kg_m2[:, None, None])
    depth = np.zeros((layers, BANDS, GRID_POINTS))
    if parts:
        parts.sort(key=lambda part: -float(part.sum()))
        depth = parts[0]
        for part in parts[1:]:
            depth = overlap(depth, part)
    grey = np.zeros((layers, BANDS))
    if prepared.has_water:
        grey += _continuum_depth(temperature, pressure, water, molecules)
    if prepared.pairs:
        density = pressure / (constants.BOLTZMANN * temperature) * 1.0e-6
        grey += _pair_depth(prepared, temperature, water, molecules, density)
    grey += _wing_depth(prepared, temperature, pressure, water, molecules)
    return depth + grey[:, :, None]
