# planetgen/physics/comet_species.py

"""
The species of a comet's ice that few comets have measured, each between the values of the comets that did (GEN.1005, GEN.1021, row k7 of the
atmosphere possibility table; docs/design/composition/comet-species-model.md).

Row k7 drew N2, O2, OCS, SO2, Ar, Kr and Xe of a comet as log-normals about the values of 67P alone, with a trapping share per species tuned to
that one comet (GEN.1005, stand-ins CK1 to CK6). GEN.1021 replaced that calibration with the range between the comets measured, the widest
range where two sources are in question (report `reports/gen1021-comet-species-second-comet-2026-10-30.md`):

1. **The bracket.** `bracket_over_h2o` gives the `(low, high)` of each species per water, and the row draws log-uniform between them. N2 is the
   N2 per CO of 67P and of C/2016 R2 times the median CO / H2O (nitrogen follows CO); argon is per water, 67P to Hale-Bopp (it does not follow
   the CO: Hale-Bopp's argon per CO is 2.7 times what the nebula holds); krypton and xenon follow argon by the ratios of 67P; OCS and SO2 are
   the range of the comets and never above the sulfur the nebula holds; O2 is Halley to 67P. [Stand-ins CB1 and CB2.]
2. **The host.** The ends of N2 scale with the host's nitrogen per carbon over the proto-solar one (`n_over_c`), those of Ar, Kr and Xe with
   `ar_over_c`, those of OCS and SO2 with `s_over_c`; all three are 1 for a host of solar composition. `nebular_per_co` gives the nebular amount
   of a species per CO, from the proto-solar abundances (`tuning.STELLAR_SOLAR_LOG_EPS`) and the carbon in CO, a share of 0.52 (the 1.5 per 10^4
   hydrogen of [R3322], `comet_ice_budget.NEBULA_PER_1E4_H`, over the proto-solar 2.9 per 10^4).
3. **The water.** CO / H2O is the median of the population of comets (row k5: 10^-1.602, `CO_OVER_H2O_MEDIAN`). [Stand-in CK4, tier T3, tweak L2.]

Tests against measurement (never inputs): 67P's own value lies inside every bracket, Hale-Bopp's argon per water and C/2016 R2's N2 per CO lie
inside, a host with twice the nitrogen per carbon doubles the ends of N2, and the sulfur species never pass the sulfur there is.

Limits. Each bracket has two comets (krypton and xenon one); the two long-period comets are coma measurements, not nucleus bulk, and 67P's own argon per water varies
1.4 decades with the season. The sulfur of OCS and SO2 plus the population's H2S is more than the nebular sulfur the chain allows for (the sulfur is not conserved across rows).
"""

import math

from planetgen import tuning
from planetgen.physics import comet_ice_budget

CO_OVER_H2O_MEDIAN = 10.0 ** -1.602
"""float: The median CO / H2O of the comets of the population (row k5, [R3782] database, per-comet medians). [Stand-in CK4, tier T3, tweak L2.]"""

SPECIES = ("N2", "Ar", "OCS", "SO2", "O2", "Kr", "Xe")
"""tuple: The species of row k7 that sit between the comets measured."""

# ---------------------------------------------------------------------
# The brackets between the comets measured (GEN.1021; docs/design/composition/comet-species-model.md section 4)
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
        raise ValueError(f"{species} has no host ratio: N2, Ar, OCS and SO2 have")
    return per_carbon / carbon_in_co_share()


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
