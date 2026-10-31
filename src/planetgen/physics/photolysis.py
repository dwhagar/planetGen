# planetgen/physics/photolysis.py

"""
The photolysis lifetime table (GEN.496): how long a molecule lasts against the far ultraviolet of its star.

A photolysis lifetime is a constant of the species, not a property of a body. `data/photolysis_lifetimes.json` holds,
for each species of the table, its lifetime at 1 AU from the quiet Sun (photodissociation plus photoionisation,
optically thin and unshielded), the range the sources give, the factor by which the active Sun shortens it where one
was read, and the species-specific ratio of its rate around a 4000 K and a 10000 K star to the Sun's at equal
91.2 to 200 nm flux. The rate in a given place is

    beta = (1 / lifetime) * S_i(T_eff) * A * F

with `F` the star's 91.2 to 200 nm flux at the place relative to the quiet Sun at 1 AU (`1 / r^2` for the Sun at `r`
AU: halving the distance quadruples the rate), `S_i` the ratio interpolated linearly in the logarithm of the rate and
of the temperature between the three templates (4000 K, the Sun at 5772 K, 10000 K) and `A` the active-Sun factor.

* **One scalar for "stellar UV" is wrong by orders of magnitude.** Around a 4000 K star the rate of CO falls to 3.5e-4
  of the Sun's, of N2 to 2.9e-5, CH4 to 0.013 and HCN to 0.027, while NH3 (0.88) and CH2O (0.78) barely change; around
  a 10000 K star CO2 rises 6.5 times and NH3 falls to 0.41.
* **Outside 4000 to 10000 K the ratio is the end template's** (`outside_templates` says so): the sources have no
  spectrum there, so a star cooler than 4000 K is under-corrected for CO, N2, CH4 and HCN. It is an open item.
* **The active-Sun factor is only for the species whose value was read** (H2O 1.68, CO2 2.35, HCN 2.1 from a search
  snippet) or computed (NH3 1.07 from measured solar spectra, `scripts/nh3_activity_factor.py`, GEN.794; a measured cycle
  maximum, milder than the Huebner and Mukherjee active Sun of H2O and CO2). For the others `active_sun=True` raises `NoData`
  rather than assume one.
* **Shielding is not in the lifetimes**: CO and N2 self-shield in a dense inner coma and in a thick air column, and the lifetimes
  in an atmosphere are longer (docs/design/composition/atmosphere-composition-parts/04-photochemistry-lifetimes-escape.md). For the
  species with a read mean UV cross section (H2O, NH3, CH4, CO2, H2S; Hu, Seager and Bains 2012 table 6) the column that absorbs the light
  is one optical depth, `shielded_column_mol_m2`, and `photon_limited_flux_mol_m2_s` is the most the light can break in a thick air
  (GEN.982). It is the gross loss: the products that recombine are not taken off.
* **Source values** are labelled in the file: the quiet-Sun lifetime of H2O, CO, CO2, HCN, CH3OH and CH2O is the
  comet-community value, of CH4, C2H6, H2S and N2 the Heays 2017 total (a factor 2 for C2H6, H2S and N2), and NH3
  6.7e3 s (1.5e-4 per s), not the 5.6e3 s that is a mean of quiet and active Sun.

The functions are pure.
"""

import json
import math
import os
from dataclasses import dataclass
from typing import Optional, Tuple

from planetgen.physics import elements

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "photolysis_lifetimes.json")

TEMPLATE_TEMPERATURES_K = (4000.0, 5772.0, 10000.0)
"""tuple: The effective temperatures of the three spectral templates; the middle one is the Sun's (ratio 1)."""


class UnknownSpecies(KeyError):
    """A species with no row in the photolysis table."""


class NoData(ValueError):
    """The table has no value for what was asked (an active-Sun factor that was never read)."""


@dataclass(frozen=True)
class Species:
    """
    One row of the table.

    Attributes:
        formula (str): The gas-table formula (formaldehyde is `CH2O`).
        lifetime_s (float): Lifetime at 1 AU, quiet Sun, seconds (dissociation plus ionisation).
        lifetime_range_s (tuple): The lowest and highest quiet-Sun value the sources give.
        beta_diss_per_s, beta_ion_per_s (float): Heays 2017 rates, quiet Sun, 1 AU.
        activity_factor (float or None): How many times shorter the lifetime is for the active Sun, if read.
        ratio_4000k, ratio_10000k (float): The dissociation rate relative to the Sun's at equal 91.2 to 200 nm flux.
        lifetime_source, status (str): Where the lifetime comes from, and how it was got.
        uv_cross_section_cm2 (float or None): The mean UV cross section (Hu, Seager and Bains 2012 table 6), `None` when none was read.
        uv_cross_section_source (str): Where it comes from.
    """

    formula: str
    lifetime_s: float
    lifetime_range_s: Tuple[float, float]
    beta_diss_per_s: float
    beta_ion_per_s: float
    activity_factor: Optional[float]
    ratio_4000k: float
    ratio_10000k: float
    lifetime_source: str
    status: str
    uv_cross_section_cm2: Optional[float] = None
    uv_cross_section_source: str = ""


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)["species"]
    table = {}
    for formula, row in data.items():
        low, high = row["lifetime_range_s"]
        if not (0 < low <= row["lifetime_s"] <= high):
            raise ValueError(f"photolysis_lifetimes.json: {formula} lifetime lies outside its own range")
        if row["ratio_4000k"] <= 0 or row["ratio_10000k"] <= 0:
            raise ValueError(f"photolysis_lifetimes.json: {formula} has a spectral ratio that is not above zero")
        if row["activity_factor"] is not None and row["activity_factor"] < 1:
            raise ValueError(f"photolysis_lifetimes.json: {formula} activity factor is below 1")
        table[formula] = Species(
            formula, row["lifetime_s"], (low, high), row["beta_diss_per_s"], row["beta_ion_per_s"],
            row["activity_factor"], row["ratio_4000k"], row["ratio_10000k"], row["lifetime_source"], row["status"],
            row.get("uv_cross_section_cm2"), row.get("uv_cross_section_source", ""),
        )
    return table


_TABLE = _load()


def species(formula):
    """The row of `formula`; raises `UnknownSpecies`."""
    try:
        return _TABLE[formula]
    except KeyError:
        raise UnknownSpecies(formula) from None


def formulas():
    """The formulas of the table, in file order."""
    return tuple(_TABLE)


def outside_templates(teff_k):
    """True when `teff_k` lies outside the 4000 to 10000 K the spectral templates cover."""
    return not TEMPLATE_TEMPERATURES_K[0] <= teff_k <= TEMPLATE_TEMPERATURES_K[-1]


def spectral_ratio(formula, teff_k):
    """
    The species' dissociation rate around a star of `teff_k` relative to the Sun's, at equal 91.2 to 200 nm flux:
    log-linear between the 4000 K, solar and 10000 K templates, the end template's value beyond them.
    """
    row = species(formula)
    if not teff_k > 0:
        raise ValueError("a temperature must be positive")
    cold, sun, hot = TEMPLATE_TEMPERATURES_K
    t = min(max(teff_k, cold), hot)
    if t <= sun:
        lo_t, hi_t, lo_r, hi_r = cold, sun, row.ratio_4000k, 1.0
    else:
        lo_t, hi_t, lo_r, hi_r = sun, hot, 1.0, row.ratio_10000k
    f = (math.log(t) - math.log(lo_t)) / (math.log(hi_t) - math.log(lo_t))
    return math.exp(math.log(lo_r) + f * (math.log(hi_r) - math.log(lo_r)))


def rate_per_s(formula, flux_ratio=1.0, teff_k=TEMPLATE_TEMPERATURES_K[1], active_sun=False):
    """
    The photolysis rate, per second, of `formula`.

    Args:
        formula (str): A species of the table.
        flux_ratio (float): The star's 91.2 to 200 nm flux at the place relative to the quiet Sun at 1 AU (`1 / r**2`
            for the Sun at `r` AU).
        teff_k (float): The star's effective temperature, kelvin.
        active_sun (bool): Apply the species' active-Sun factor (the star is in an active phase).

    Raises:
        NoData: If `active_sun` and no factor was read for the species.
    """
    row = species(formula)
    if not flux_ratio >= 0:
        raise ValueError("a flux ratio cannot be negative")
    activity = 1.0
    if active_sun:
        if row.activity_factor is None:
            raise NoData(f"no active-Sun factor was read for {formula}")
        activity = row.activity_factor
    return flux_ratio * activity * spectral_ratio(formula, teff_k) / row.lifetime_s


def lifetime_s(formula, flux_ratio=1.0, teff_k=TEMPLATE_TEMPERATURES_K[1], active_sun=False):
    """The photolysis lifetime, seconds: `1 / rate_per_s`; infinite with no flux."""
    rate = rate_per_s(formula, flux_ratio, teff_k, active_sun)
    return math.inf if rate == 0 else 1.0 / rate


def shielded_column_mol_m2(formula):
    """
    The column of `formula` that absorbs its light, mol per m2 (GEN.982): the one that has an optical depth of 1 at the species' mean UV cross
    section `sigma`, `1 / sigma` molecules per cm2. Gas below it sees none of the far ultraviolet, so a thick air loses the species only from this column
    however much of it there is (the critical-pressure rule of Hu, Seager and Bains 2012: a gas is photolysed above the pressure where its
    column is one optical depth).

    Raises:
        NoData: If no mean cross section was read for the species.
    """
    sigma = species(formula).uv_cross_section_cm2
    if sigma is None:
        raise NoData(f"no mean UV cross section was read for {formula}")
    return 1.0e4 / (sigma * elements.AVOGADRO)


def photon_limited_flux_mol_m2_s(formula, flux_ratio=1.0, teff_k=TEMPLATE_TEMPERATURES_K[1], active_sun=False):
    """
    The most of `formula` the light can break in a thick air, mol per m2 of ground per second (GEN.982): the unshielded rate per molecule
    (`rate_per_s`) of the column that absorbs the light (`shielded_column_mol_m2`). This is the gross photolysis; what a species' products recombine
    into is not taken off. For the Sun at 1 AU, methane comes to about 78 Tmol a year over the Earth's area.

    Raises:
        NoData: If no mean cross section was read for the species.
    """
    return shielded_column_mol_m2(formula) * rate_per_s(formula, flux_ratio, teff_k, active_sun)


def inverse_square(distance_au):
    """`flux_ratio` of the Sun at `distance_au`."""
    if not distance_au > 0:
        raise ValueError("a distance must be positive")
    return 1.0 / (distance_au * distance_au)
