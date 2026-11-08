# planetgen/physics/comet_species.py

"""
The species of a comet's ice that only one comet has measured, from the host star's own ratios and one trapping factor each (GEN.1005, k7;
docs/design/composition/comet-species-model.md).

Row k7 of the atmosphere possibility table drew N2, O2, OCS, SO2, Ar, Kr and Xe of a comet as log-normals about the values of 67P [R3835], with a
width of 0.6 dex that was a design choice. Those seven species have a population of one comet, so no median can be taken from a population. This module writes each
that can be tied to the star's elements as a chain:

    X / H2O  =  (X per CO in the nebula, from the host's own ratio to carbon)  x  (the share of that kept in the ice)  x  (CO / H2O of the comet)

1. **The nebula.** The number of the element per carbon of the host is the proto-solar one (`tuning.STELLAR_SOLAR_LOG_EPS`) times the host's own ratio over it
   (`n_over_c`, `ar_over_c`, `s_over_c`, which are 1 for a host of solar composition: the generator's stars carry [X/Fe] offsets,
   `stellar_abundances.number_ratio`). The carbon in CO is the 1.5 per 10^4 hydrogen of [R3322] (`comet_ice_budget.NEBULA_PER_1E4_H`) over the
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
than the nebular sulfur the chain allows for (the sulfur is not conserved across rows). The share of carbon in CO is the nebular share of one disk.
"""

import math

from planetgen import tuning
from planetgen.physics import comet_ice_budget

CO_OVER_H2O_MEDIAN = 10.0 ** -1.602
"""float: The median CO / H2O of the comets of the population (row k5, [R3782] database, per-comet medians). [Stand-in CK4, tier T3, tweak L2.]"""

TRAPPED_OVER_NEBULAR = {"N2": 0.158, "Ar": 0.0145, "OCS": 0.187, "SO2": 0.592}
"""dict: The share of the nebular amount (per CO, at the host's ratio, times CO / H2O) that is in the ice, by species: the factor that gives the value of the one comet measured (67P,
[R3835]: N2 8.9e-4, Ar 5.8e-6, OCS 4.1e-4, SO2 1.3e-3 per water). [Stand-ins CK1 to CK3, tier T7, tweak L3.]"""

UNSCALED_OVER_H2O = {"O2": 3.1e-2, "Kr": 4.9e-7, "Xe": 2.4e-7}
"""dict: The species with no host ratio, per water: the values of 67P [R3835]. [Stand-in CK5, tier T7, tweak L3.]"""

SIGMA_DEX = 0.4
"""float: The width, in decades, of every species. [Stand-in CK6, tier T6, tweak L1.]"""

SPECIES = tuple(TRAPPED_OVER_NEBULAR) + tuple(UNSCALED_OVER_H2O)

# ---------------------------------------------------------------------
# The brackets of two comets (GEN.1021; docs/design/composition/comet-species-model.md section 4)
# ---------------------------------------------------------------------

N2_OVER_CO = (0.0287, 0.09)
"""tuple: N2 per CO of the two comets with the widest range: 67P (best measurement, Rubin and others 2020, read in the GEN.1021 report) and C/2016 R2 (Q(N2)/Q(CO) about 0.09, [R0109] Anderson and others 2022; 0.06 and 0.08 in [R0836]). [Stand-in CB1, tier T3, tweak L0.]"""

AR_OVER_H2O = (5.8e-6, 8.7e-3)
"""tuple: Ar per H2O of the two comets: 67P ([R3835] Rubin and others 2019 table) and Hale-Bopp (Ar/O 0.0058 with Q(O) = 1.5 Q(H2O), [R4215] Stern and others 2000; a coma, not a nucleus). [Stand-in CB1.]"""

OCS_OVER_H2O = (5e-4, 4e-3)
"""tuple: OCS per H2O, the range of the comets measured (0.05 to 0.40 percent; Table 7 of the 46P paper, A&A 648, A49, no comet named; the 67P value 0.07 percent lies inside). [Stand-in CB1.]"""

SO2_OVER_H2O = (3e-4, 2.3e-3)
"""tuple: SO2 per H2O, the range of the comets measured (0.03 to 0.23 percent; the same table; the 67P value 0.127 percent lies inside). [Stand-in CB1.]"""

O2_OVER_H2O = (3.7e-2 - 1.7e-2, 3.8e-2 + 0.85e-2)
"""tuple: O2 per H2O from the lower end of Halley's value (3.7 +/- 1.7 percent, [R3836] Rubin and others 2015) to the upper end of 67P's (3.80 +/- 0.85 percent, [R0386] Bieler and others 2015). [Stand-in CB1.]"""

KR_OVER_AR = 0.089
"""float: Kr per Ar of 67P ([R3834] Rubin and others 2018, 0.089 +/- 0.021), applied to the brackets of argon: krypton follows argon. Only one comet has measured it. [Stand-in CB2, tier T2, tweak L1.]"""

XE_OVER_KR = 0.269
"""float: Xe per Kr of 67P ([R3834], 0.269 +/- 0.064), applied the same way. [Stand-in CB2.]"""

_HOST_RATIO = {"N2": "n_over_c", "Ar": "ar_over_c", "Kr": "ar_over_c", "Xe": "ar_over_c", "OCS": "s_over_c", "SO2": "s_over_c"}
"""dict: Which host ratio over the proto-solar one scales a species' bracket; O2 has none."""


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


def _host_scale(species, n_over_c, ar_over_c, s_over_c):
    ratios = {"n_over_c": n_over_c, "ar_over_c": ar_over_c, "s_over_c": s_over_c}
    if min(ratios.values()) <= 0.0 or not all(math.isfinite(r) for r in ratios.values()):
        raise ValueError(f"a host ratio over the proto-solar one is above zero, not {tuple(ratios.values())}")
    name = _HOST_RATIO.get(species)
    return 1.0 if name is None else ratios[name]


def bracket_over_h2o(species, n_over_c=1.0, ar_over_c=1.0, s_over_c=1.0):
    """
    The `(low, high)` of `species` per H2O between the comets measured, round a host whose N/C, Ar/C and S/C are the given multiples of the proto-solar ones (GEN.1021).

    N2 is the N2 per CO of the two comets times the nebular CO / H2O median (the chain is kept: nitrogen follows CO). Argon is per water (it follows the water ice, not the CO: Hale-Bopp's Ar per CO is 2.7 times
    what the nebula holds), krypton and xenon are argon times the ratios of 67P. OCS and SO2 are the range of the comets, never above the sulfur the nebula holds (a share of 1). O2 has no host ratio.
    """
    if species not in SPECIES:
        raise ValueError(f"{species} is none of {', '.join(SPECIES)}")
    scale = _host_scale(species, n_over_c, ar_over_c, s_over_c)
    if species == "N2":
        low, high = (value * CO_OVER_H2O_MEDIAN for value in N2_OVER_CO)
    elif species in ("Ar", "Kr", "Xe"):
        factor = {"Ar": 1.0, "Kr": KR_OVER_AR, "Xe": KR_OVER_AR * XE_OVER_KR}[species]
        low, high = (value * factor for value in AR_OVER_H2O)
    elif species == "OCS":
        low, high = OCS_OVER_H2O
    elif species == "SO2":
        low, high = SO2_OVER_H2O
    else:
        low, high = O2_OVER_H2O
    low, high = low * scale, high * scale
    if species in ("OCS", "SO2"):
        whole = nebular_per_co(species, s_over_c=s_over_c) * CO_OVER_H2O_MEDIAN       # a share of 1: all the nebular sulfur of the species
        high = min(high, whole)
    return low, high


def low_over_h2o(species, n_over_c=1.0, ar_over_c=1.0, s_over_c=1.0):
    """The low end of `bracket_over_h2o`."""
    return bracket_over_h2o(species, n_over_c, ar_over_c, s_over_c)[0]


def high_over_h2o(species, n_over_c=1.0, ar_over_c=1.0, s_over_c=1.0):
    """The high end of `bracket_over_h2o`."""
    return bracket_over_h2o(species, n_over_c, ar_over_c, s_over_c)[1]


def sigma_dex():
    """The width of each species, in decades."""
    return SIGMA_DEX
