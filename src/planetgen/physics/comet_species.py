# planetgen/physics/comet_species.py

"""
The species of a comet's ice that only one comet has measured, from the host star's own ratios and one trapping factor each (GEN.1005, k7;
docs/design/composition/comet-species-model.md).

Row k7 of the atmosphere possibility table drew N2, O2, OCS, SO2, Ar, Kr and Xe of a comet as log-normals about the values of 67P (Rubin et al. 2019), with a
width of 0.6 dex that was a design choice. Those seven species have a population of one comet, so no median can be taken from a population. This module writes each
that can be tied to the star's elements as a chain:

    X / H2O  =  (X per CO in the nebula, from the host's own ratio to carbon)  x  (the share of that kept in the ice)  x  (CO / H2O of the comet)

1. **The nebula.** The number of the element per carbon of the host is the proto-solar one (`tuning.STELLAR_SOLAR_LOG_EPS`) times the host's own ratio over it
   (`n_over_c`, `ar_over_c`, `s_over_c`, which are 1 for a host of solar composition: the generator's stars carry [X/Fe] offsets,
   `stellar_abundances.number_ratio`). The carbon in CO is the 1.5 per 10^4 hydrogen of Oberg, Murray-Clay and Bergin 2011 (`comet_ice_budget.NEBULA_PER_1E4_H`) over the
   proto-solar 2.9 per 10^4, a share of 0.52. Nitrogen is all N2 (half the atoms), argon is atomic, and the sulfur of the two sulfur species is shared as the table says.
2. **The ice.** The share of the nebular amount that ends in the ice is one number per species (`TRAPPED_OVER_NEBULAR`), the factor that makes the chain give the value of the one
   comet measured: the trapping of a gas in amorphous water ice depends on the temperature the ice formed at (Bar-Nun: the trapped gas drops six decades between 25 and 100 K,
   and CO at 3.5 percent needs 48 +/- 5 K), a quantity the generator does not carry per comet. [Stand-ins CK1 to CK3, tier T7, tweak L3: tuned to one comet.]
3. **The water.** CO / H2O is the median of the population of comets (row k5: 10^-1.602, `CO_OVER_H2O_MEDIAN`). [Stand-in CK4, tier T3, tweak L2.]

O2, Kr and Xe have no host ratio here (O2 is not tied to the oxygen budget by any source read; krypton and xenon are not in the abundance table) and keep the value of 67P (CK5).

`sigma_dex` replaces the design width of 0.6 dex with 0.4: Bar-Nun's 5 K uncertainty in the formation temperature at a decade of trapped gas per 12 K near the middle of the range. [Stand-in CK6, tier T6, tweak L1.]

Tests against measurement (never inputs): the chain gives the 67P value at a host of solar composition (it is the calibration, and the test says so); a host with twice the nitrogen per carbon gives twice the N2; the
four species measured in populations (CH4, C2H6, H2CO, H2S, rows of k7 that are not stand-ins) have scatter of 0.33 to 0.39 dex, which brackets the width.

Limits. One comet fixes all seven values: a second comet or a laboratory table of trapping by temperature replaces the factors. The sulfur of the two species plus the population's H2S is more
than the nebular sulfur the chain allows for (the sulfur is not conserved across rows). The share of carbon in CO is the nebular share of one disc.
"""

import math

from planetgen import tuning
from planetgen.physics import comet_ice_budget

CO_OVER_H2O_MEDIAN = 10.0 ** -1.602
"""float: The median CO / H2O of the comets of the population (row k5, Robinson et al. 2024 database, per-comet medians). [Stand-in CK4, tier T3, tweak L2.]"""

TRAPPED_OVER_NEBULAR = {"N2": 0.158, "Ar": 0.0145, "OCS": 0.187, "SO2": 0.592}
"""dict: The share of the nebular amount (per CO, at the host's ratio, times CO / H2O) that is in the ice, by species: the factor that gives the value of the one comet measured (67P,
Rubin et al. 2019: N2 8.9e-4, Ar 5.8e-6, OCS 4.1e-4, SO2 1.3e-3 per water). [Stand-ins CK1 to CK3, tier T7, tweak L3.]"""

UNSCALED_OVER_H2O = {"O2": 3.1e-2, "Kr": 4.9e-7, "Xe": 2.4e-7}
"""dict: The species with no host ratio, per water: the values of 67P (Rubin et al. 2019). [Stand-in CK5, tier T7, tweak L3.]"""

SIGMA_DEX = 0.4
"""float: The width, in decades, of every species. [Stand-in CK6, tier T6, tweak L1.]"""

SPECIES = tuple(TRAPPED_OVER_NEBULAR) + tuple(UNSCALED_OVER_H2O)


def carbon_in_co_share():
    """The share of the nebula's carbon that is in CO: 1.5 per 10^4 hydrogen over the proto-solar carbon."""
    return comet_ice_budget.NEBULA_PER_1E4_H["CO"] * 1e-4 / 10.0 ** (tuning.STELLAR_SOLAR_LOG_EPS["C"] - 12.0)


def nebular_per_co(species, n_over_c=1.0, ar_over_c=1.0, s_over_c=1.0):
    """The molecules of `species` per CO in the nebula of a host whose N/C, Ar/C and S/C are `n_over_c`, `ar_over_c`, `s_over_c` times the proto-solar ones."""
    logs = tuning.STELLAR_SOLAR_LOG_EPS
    if species == "N2":
        per_carbon = 0.5 * 10.0 ** (logs["N"] - logs["C"]) * n_over_c
    elif species == "Ar":
        per_carbon = 10.0 ** (logs["Ar"] - logs["C"]) * ar_over_c
    elif species in ("OCS", "SO2"):
        per_carbon = 10.0 ** (logs["S"] - logs["C"]) * s_over_c
    else:
        raise ValueError(f"{species} has no host ratio: {', '.join(TRAPPED_OVER_NEBULAR)} have")
    return per_carbon / carbon_in_co_share()


def median_over_h2o(species, n_over_c=1.0, ar_over_c=1.0, s_over_c=1.0):
    """The median `species` per H2O of a comet round a host with the given ratios over the proto-solar ones."""
    if species in UNSCALED_OVER_H2O:
        return UNSCALED_OVER_H2O[species]
    ratios = (n_over_c, ar_over_c, s_over_c)
    if min(ratios) <= 0.0 or not all(math.isfinite(r) for r in ratios):
        raise ValueError(f"a host ratio over the proto-solar one is above zero, not {ratios}")
    return nebular_per_co(species, *ratios) * TRAPPED_OVER_NEBULAR[species] * CO_OVER_H2O_MEDIAN


def sigma_dex():
    """The width of each species, in decades."""
    return SIGMA_DEX
