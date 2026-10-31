# planetgen/physics/atmosphere_retention.py

"""
When a body can hold an atmosphere (GEN.355, docs/design/habitability/atmospheres-retention-and-classes.md section 3, registry 6.74)
===================================================================================================================================

The answer the design note gives, as functions. Three gates and one loss:

1. **The cosmic shoreline** (`shoreline_ratio`, `verdict_of`, `pressure_cap_pa`): bodies with and without air separate along the
   cumulative insolation over the fourth power of the escape velocity (Zahnle and Catling 2017). The ratio `r = I / (K v^4)` is under 0.7 for
   a body that holds a dense atmosphere and over 30 for one that holds none; `K` is the middle of the window the 13 bodies of the design
   table allow (3.0e-4 to 6.6e-4). `I` is the insolation in Earth's, or for a star that stays active for long the insolation times the star's
   cumulative XUV over the Sun's to 4.5 Gyr. The shoreline caps the pressure and does not set it.
2. **Which gases are held** (`held_species`, `minimum_molar_mass_kg_mol`): a gas stays when its Jeans parameter at the exobase is 30 or more.
   The exobase temperature is the thermosphere model's (`thermosphere.exobase_temperature_k`), not a second one-line fit.
3. **The cold side** (`cold_side_k`): the temperature under which a gas at a partial pressure is ice, from the vapour curves of
   `saturation` (so Titan holds its nitrogen as vapour at 94 K and Pluto only a pascal of it).
4. **Energy-limited loss** (`energy_limited_loss_fraction`): the share of a planet's mass the absorbed XUV could lift out of its well over the
   star's life, an upper bound on the loss of an envelope or an ocean's hydrogen.

The pressure stage of the generator takes its shoreline from the logistic of rows d of the atmosphere recipes (GEN.511); this module is the
library that answers the design question and is tested on the 13 bodies the question was written from. It changes no generated output.

Use:

    r = shoreline_ratio(escape_velocity_km_s(mass_earth, radius_earth), insolation_earth)
    verdict_of(r)                       # "dense", "moderate", "tenuous", "trace" or "none"
    held = held_species(mass_earth, radius_earth, teq_k, ("N2", "CO2", "H2"))
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, gases, saturation, thermosphere, unit_registry as units

_mearth_to_kg = units.converter("Mearth", "kg")

ESCAPE_VELOCITY_EARTH_KM_S = 11.186
"""float: Earth's escape velocity, km/s (the design table's v_esc is this times the square root of mass over radius, in Earth's)."""

SHORELINE_K = 3.85e-4
"""float: The constant of `S = K v^4` (Earth's insolation and km/s): the geometric middle of the window 3.0e-4 to 6.6e-4 that Titan, Pluto and
Triton (lower bound) and Ganymede (upper bound) leave, fitted to 13 Solar System bodies and contradicted by 55 Cnc e and TOI-561 b (row SH1)."""

SHORELINE_K_WINDOW = (3.0e-4, 6.6e-4)
"""tuple: The window `SHORELINE_K` is the middle of (row SH1)."""

VERDICT_EDGES = ((0.7, "dense"), (1.0, "moderate"), (2.5, "tenuous"), (30.0, "trace"))
"""tuple: `(upper edge of r, verdict)`; over the last edge the verdict is "none" (row SH2)."""

PRESSURE_CAPS_PA = (None, 1.0e4, 1.0e3, 1.0, 0.0)
"""tuple: The cap on the surface pressure, Pa, in each verdict of `VERDICT_EDGES` and then "none": no cap (dense), 0.1 bar, 10 mbar, 10 micro-bar
(1 Pa), none (row SH3)."""

JEANS_KEEP = 30.0
"""float: A gas of this Jeans parameter or more at the exobase stays for the age of a body (Catling and Kasting 2017); one under 15 is gone within
about 1 Gyr (row SH4)."""

EARTH_TEQ_K = 278.6
"""float: The equilibrium temperature of a body at 1 Earth insolation with no albedo."""

_AVOGADRO = 6.02214076e23
_GY_S = 1.0e9 * constants.SECONDS_PER_YEAR


@dataclass(frozen=True)
class Retention:
    """
    Attributes:
        escape_velocity_km_s (float): The body's.
        ratio (float): The shoreline ratio `r`.
        verdict (str): `verdict_of(ratio)`.
        pressure_cap_pa (float or None): The cap the verdict puts on the surface pressure; `None` for none.
        equilibrium_k (float): The equilibrium temperature at the albedo.
        exobase_k (float): The exobase temperature of the thermosphere model.
        held (tuple): The formulas of `species` that pass the Jeans gate.
    """
    escape_velocity_km_s: float
    ratio: float
    verdict: str
    pressure_cap_pa: object
    equilibrium_k: float
    exobase_k: float
    held: tuple


def escape_velocity_km_s(mass_earth, radius_earth):
    """The escape velocity of a body, km/s."""
    return ESCAPE_VELOCITY_EARTH_KM_S * math.sqrt(mass_earth / radius_earth)


def shoreline_ratio(escape_velocity, insolation_earth, k=SHORELINE_K):
    """The shoreline ratio `r = I / (K v^4)`: under 1 is on the side of the shoreline where air stays. `insolation_earth` is the cumulative
    insolation in Earth's, or for an active star its XUV-weighted figure (the star's cumulative XUV over the Sun's to 4.5 Gyr, times S)."""
    return insolation_earth / (k * escape_velocity ** 4)


def verdict_of(ratio):
    """`"dense"`, `"moderate"`, `"tenuous"`, `"trace"` or `"none"` for a shoreline ratio."""
    for edge, verdict in VERDICT_EDGES:
        if ratio < edge:
            return verdict
    return "none"


def pressure_cap_pa(ratio):
    """The most surface pressure the shoreline lets a body hold, Pa: `None` for no cap (a dense atmosphere is capped by its supply)."""
    for (edge, _), cap in zip(VERDICT_EDGES, PRESSURE_CAPS_PA):
        if ratio < edge:
            return cap
    return PRESSURE_CAPS_PA[-1]


def equilibrium_temperature_k(insolation_earth, albedo=0.3):
    """The equilibrium temperature of a body, K, in the insolation (Earth's) with the albedo."""
    return EARTH_TEQ_K * (insolation_earth * (1.0 - albedo)) ** 0.25


def jeans_parameter(molar_mass_kg_mol, escape_velocity, exobase_k):
    """The Jeans parameter `m v^2 / (2 k T)` of a gas at the exobase, for the escape velocity in km/s."""
    return molar_mass_kg_mol / _AVOGADRO * (escape_velocity * 1.0e3) ** 2 / (2.0 * constants.BOLTZMANN * exobase_k)


def minimum_molar_mass_kg_mol(escape_velocity, exobase_k):
    """The lightest gas a body keeps, kg/mol: the one whose Jeans parameter is `JEANS_KEEP`."""
    return 2.0 * constants.BOLTZMANN * exobase_k * JEANS_KEEP * _AVOGADRO / (escape_velocity * 1.0e3) ** 2


def exobase_k(insolation_earth, xuv_flux_earth, albedo=0.3, co2_fraction=0.0):
    """The exobase temperature, K: the thermosphere model's, from the mesopause at the equilibrium temperature and the XUV flux (Earth's)."""
    mesopause = thermosphere.mesopause_temperature_k(equilibrium_temperature_k(insolation_earth, albedo))
    return thermosphere.exobase_temperature_k(xuv_flux_earth, mesopause, co2_fraction)


def held_species(escape_velocity, exobase, species):
    """The formulas of `species` whose Jeans parameter at `exobase` (K) is `JEANS_KEEP` or more, for the escape velocity in km/s."""
    return tuple(formula for formula in species if jeans_parameter(gases.gas(formula).molar_mass, escape_velocity, exobase) >= JEANS_KEEP)


def cold_side_k(formula, partial_pressure_pa):
    """The temperature under which `formula` at `partial_pressure_pa` is ice and not vapour, K (`None` over its critical pressure)."""
    return saturation.condensation_temperature_k(formula, partial_pressure_pa)


def energy_limited_loss_fraction(mass_earth, radius_earth, cumulative_xuv_j_m2, efficiency=thermosphere.ENERGY_LIMITED_EFFICIENCY):
    """The share of the planet's mass the absorbed XUV could lift out of its well by its age, `eps pi R^3 E / (G M^2)`: an upper bound
    (every absorbed photon lifts gas) on the loss of an envelope, and the clock of an ocean's hydrogen."""
    radius = radius_earth * constants.EARTH_RADIUS_KM * 1.0e3
    mass = _mearth_to_kg(mass_earth)
    return efficiency * math.pi * radius ** 3 * cumulative_xuv_j_m2 / (constants.G * mass ** 2)


def retention(mass_earth, radius_earth, insolation_earth, xuv_index=None, xuv_flux_earth=1.0, albedo=0.3, co2_fraction=0.0,
              species=("H2", "He", "CH4", "H2O", "N2", "CO", "O2", "Ar", "CO2", "SO2"), k=SHORELINE_K):
    """
    The `Retention` of a body: its shoreline ratio and verdict, the pressure the verdict allows and the gases it holds.

    Args:
        mass_earth, radius_earth (float): The body's.
        insolation_earth (float): The cumulative insolation in Earth's (the present one for a Sun-like star).
        xuv_index (float, optional): The star's cumulative XUV over the Sun's to 4.5 Gyr; `None` takes the insolation alone.
        xuv_flux_earth (float): The XUV flux at the body now, in Earth's present.
        albedo, co2_fraction (float): For the equilibrium and exobase temperatures.
        species (tuple): The formulas to test.
        k (float): The shoreline constant.
    """
    velocity = escape_velocity_km_s(mass_earth, radius_earth)
    ratio = shoreline_ratio(velocity, insolation_earth * (xuv_index if xuv_index is not None else 1.0), k)
    exobase = exobase_k(insolation_earth, xuv_flux_earth, albedo, co2_fraction)
    return Retention(velocity, ratio, verdict_of(ratio), pressure_cap_pa(ratio), equilibrium_temperature_k(insolation_earth, albedo), exobase,
                     held_species(velocity, exobase, species))
