# planetgen/physics/greenhouse_table.py

"""
The greenhouse of an air, read from a table fitted to the radiative-convective model (GEN.850;
docs/design/composition/greenhouse-table.md)
===================================================================================================================

`radiative_convective.warming` finds the surface temperature of an air from the ExoMolOP correlated-k tables and the HITRAN
collision-induced pairs, and takes one to a few seconds (a minute for a dense air). The surface-air fixed point asks for the greenhouse of
every body five to ten times, so the generator reads this table instead: a small network, fitted by `scripts/fit_greenhouse_table.py` to
the model's own answers on tens of thousands of airs (the airs the generator makes and perturbed copies of them), whose weights are the
file `physics/data/greenhouse_table.json`. It is a stand-in for the model and carries a registry row (GEN.725, fill id IR11) with the error it
was measured to have; the model remains the truth (`tests/test_greenhouse_table.py` holds the table to it where the data are on the machine).

    ln(T_s / T_e) = net(features)         T_e = (F / sigma)^(1/4), F the absorbed flux

The features (`FEATURE_NAMES`) are the logarithms of the effective temperature, the gravity and the surface pressure; of the mole
fractions of the three tabled gases (H2O, CO2, CH4) and of the gases that pair in collisions (N2, O2, H2, He); of the shares of the other infrared-active gases
(those the gas table flags `"ir"`, and those flagged `"weak"`, which have only stand-in bands); the mean molar mass and the heat-capacity ratio at
`T_e`; and the share of the absorbed sunlight that a haze takes aloft with the pressure it is taken at (the inputs of
`radiative_convective.solve`, GEN.851). A feature outside the range the table was fitted on is held at its edge, the way the k-tables themselves hold
their edge value beyond 100 bar and 100 K; the table was fitted up to a surface pressure of `PRESSURE_MAX_PA`, above which the model itself does not
converge.

Use:

    warming = greenhouse_table.warming(partial_pressures_pa, gravity_m_s2, absorbed_flux_w_m2)
    warming.surface_temperature_k, warming.factor      # the factor is G = (T_s / T_e)^4 - 1, the repo's greenhouse variable
"""

import functools
import hashlib
import json
import math
import os
from dataclasses import dataclass

import numpy as np

from planetgen.datasets import hitran_cia, ir_ktables
from planetgen.physics import constants, elements, gases, mixture
from planetgen.util import detmath

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "greenhouse_table.json")
"""str: The fitted table (`scripts/fit_greenhouse_table.py`)."""


TABLED_GASES = ("H2O", "CO2", "CH4")
"""tuple: The gases with a correlated-k table (`datasets/ir_ktables.py`)."""
PAIRING_GASES = ("N2", "O2", "H2", "He")
"""tuple: The gases whose collision-induced pairs with each other and with the tabled gases the model reads (`datasets/hitran_cia.py`)."""
FEATURE_NAMES = (
    "ln_effective_temperature", "ln_gravity", "ln_pressure",
    *("ln_fraction_" + formula for formula in TABLED_GASES + PAIRING_GASES),
    "ln_fraction_ir", "ln_fraction_weak", "ln_mean_molar_mass", "heat_capacity_ratio", "aloft_share", "ln_aloft_pressure_share",
)
"""tuple: The columns of `features`."""

FRACTION_FLOOR = 1.0e-9
"""float: Added to a mole fraction under its logarithm, so an absent gas is a finite feature."""
THIN_AIR_PA = 10.0
"""float: STAND-IN (IR11, tier T5, L0): below this surface pressure the log of the warming is scaled by `P / THIN_AIR_PA`: an air this thin is optically thin in every band, so
its greenhouse goes with its column, which the network, fitted on airs of 10 Pa and more, would otherwise hold at its edge."""
ABSORBER_SHARE_MIN = 0.01
"""float: STAND-IN (IR11, tier T5, L0): below this mole share of gases that absorb or collide in the infrared (every gas the gas table does not flag `none`: the noble gases) the log of
the warming is scaled by the share over this: an air of noble gases has no greenhouse, and the network, fitted on airs that always hold some absorber, does not know it."""
PRESSURE_MAX_PA = 1.0e8
"""float: The surface pressure the table was fitted up to (1000 bar; the k-tables stop at 100 bar, and the model does not converge beyond about 1e8 Pa)."""
HIGHEST_SURFACE_K = 4000.0
"""float: The hottest surface the table gives (the model's own ceiling, `radiative_convective.HIGHEST_SURFACE_K`)."""
DEFAULT_ALOFT_SHARE_OF_PRESSURE = 0.1
"""float: The pressure aloft, as a share of the surface pressure, of an air with no sunlight taken aloft (any value: the share is zero)."""


@dataclass(frozen=True)
class Warming:
    """
    The greenhouse warming of one air.

    Attributes:
        effective_temperature_k (float): `T_e`, the airless grey temperature of the absorbed flux.
        surface_temperature_k (float): `T_s` of the table.
        factor (float): `G = (T_s / T_e)^4 - 1`, the repo's greenhouse variable (never below 0).
    """

    effective_temperature_k: float
    surface_temperature_k: float
    factor: float


@dataclass(frozen=True)
class Table:
    """The fitted network: standardisation, the feature range, the layers and what the fit measured."""

    mean: np.ndarray
    scale: np.ndarray
    low: np.ndarray
    high: np.ndarray
    weights: tuple
    biases: tuple
    summary: dict


def data_digest():
    """A digest of the downloads the model reads: the sha256 of every file of the k-table and collision-induced datasets and the version of their
    reduction. The table records the digest it was fitted under, and `table` refuses one fitted to other data."""
    parts = [f"{dataset.key}:{dataset.derive_version}" for dataset in (ir_ktables.DATASET, hitran_cia.DATASET)]
    parts.extend(piece.sha256 for dataset in (ir_ktables.DATASET, hitran_cia.DATASET) for piece in dataset.pieces)
    return hashlib.sha256("\n".join(parts).encode("utf-8")).hexdigest()


@functools.cache
def table():
    """The fitted table, read once."""
    with open(DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    if tuple(data["features"]) != FEATURE_NAMES:
        raise elements.DataError(f"{DATA_FILE}: the features are not the module's; fit the table again (scripts/fit_greenhouse_table.py)")
    if data["data_digest"] != data_digest():
        raise elements.DataError(f"{DATA_FILE}: the table was fitted to other k-tables or collision-induced files than the catalogue's; "
                                 "fit it again (scripts/fit_greenhouse_table.py)")
    return Table(mean=np.array(data["mean"]), scale=np.array(data["scale"]), low=np.array(data["low"]), high=np.array(data["high"]),
                 weights=tuple(np.array(layer) for layer in data["weights"]), biases=tuple(np.array(layer) for layer in data["biases"]),
                 summary=data["summary"])


def features(partial_pressures_pa, gravity_m_s2, absorbed_flux_w_m2, aloft_share=0.0, aloft_pressure_pa=None):
    """The inputs of the network for an air (see the module docstring), a vector in the order of `FEATURE_NAMES`."""
    total = mixture.total_pressure_pa(partial_pressures_pa)
    fractions = {formula: pressure / total for formula, pressure in partial_pressures_pa.items()}
    effective = (absorbed_flux_w_m2 / constants.STEFAN_BOLTZMANN_CONSTANT) ** 0.25
    others = [(formula, fraction) for formula, fraction in fractions.items() if formula not in TABLED_GASES and formula not in PAIRING_GASES]
    infrared = sum(fraction for formula, fraction in others if gases.gas(formula).ir == "ir")
    weak = sum(fraction for formula, fraction in others if gases.gas(formula).ir == "weak")
    molar_mass = sum(fraction * gases.gas(formula).molar_mass for formula, fraction in fractions.items())
    if aloft_pressure_pa is None:
        aloft_pressure_pa = total * DEFAULT_ALOFT_SHARE_OF_PRESSURE
    row = [math.log(effective), math.log(gravity_m_s2), math.log(total)]
    row.extend(math.log(fractions.get(formula, 0.0) + FRACTION_FLOOR) for formula in TABLED_GASES + PAIRING_GASES)
    row.extend((math.log(infrared + FRACTION_FLOOR), math.log(weak + FRACTION_FLOOR), math.log(molar_mass),
                mixture.heat_capacity_ratio(partial_pressures_pa, effective), aloft_share, math.log(aloft_pressure_pa / total)))
    return np.array(row)


def log_ratio(vector, fitted=None):
    """`ln(T_s / T_e)` the network gives for a feature vector (held inside the fitted range), never below 0."""
    fitted = fitted or table()
    hidden = (np.clip(vector, fitted.low, fitted.high) - fitted.mean) / fitted.scale
    for weights, biases in zip(fitted.weights[:-1], fitted.biases[:-1]):
        hidden = np.tanh(detmath.matmul(hidden, weights) + biases)
    return max(float((detmath.matmul(hidden, fitted.weights[-1]) + fitted.biases[-1])[0]), 0.0)


def warming(partial_pressures_pa, gravity_m_s2, absorbed_flux_w_m2, aloft_share=0.0, aloft_pressure_pa=None):
    """
    The warming of an air under an absorbed flux (the globally averaged W/m2 of `greenhouse.absorbed_flux_w_m2`): the table's reading of what
    `radiative_convective.warming` finds. `aloft_share` of the flux is deposited at `aloft_pressure_pa` instead of the ground (a haze's sunlight).

    Raises:
        ValueError: The mixture is empty, or the flux or gravity is not positive, or the share is outside [0, 1).
    """
    for name, value in (("absorbed_flux_w_m2", absorbed_flux_w_m2), ("gravity_m_s2", gravity_m_s2)):
        if not (value > 0.0 and math.isfinite(value)):
            raise ValueError(f"{name} must be a positive number, not {value!r}")
    if not 0.0 <= aloft_share < 1.0:
        raise ValueError(f"aloft_share must be in [0, 1), not {aloft_share!r}")
    effective = (absorbed_flux_w_m2 / constants.STEFAN_BOLTZMANN_CONSTANT) ** 0.25
    ratio = log_ratio(features(partial_pressures_pa, gravity_m_s2, absorbed_flux_w_m2, aloft_share, aloft_pressure_pa))
    total = mixture.total_pressure_pa(partial_pressures_pa)
    absorbing = sum(pressure for formula, pressure in partial_pressures_pa.items() if gases.gas(formula).ir != "none") / total
    ratio *= min(1.0, total / THIN_AIR_PA, absorbing / ABSORBER_SHARE_MIN)
    surface = min(effective * math.exp(ratio), HIGHEST_SURFACE_K) if effective < HIGHEST_SURFACE_K else effective
    return Warming(effective_temperature_k=effective, surface_temperature_k=surface,
                   factor=max((surface / effective) ** 4 - 1.0, 0.0))
