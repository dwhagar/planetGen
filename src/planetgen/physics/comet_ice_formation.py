# planetgen/physics/comet_ice_formation.py

"""
One latent ice formation temperature per comet drives its CO / H2O, N2 and Ar together (GEN.1022, option B of the comet species model;
docs/design/composition/comet-species-model.md section 5).

The comets measured say the species are correlated: the two CO-rich comets (Hale-Bopp, C/2016 R2) are also argon-rich and nitrogen-rich per CO, and
67P is poor in all of them. What they share is the temperature the ice formed at (Bar-Nun: the trapped gas falls six decades between 25 and 100 K, and the CO of
Halley, 3.5 percent of the water, needs ice formed at 48 +/- 5 K). The generator does not store that temperature, so it is the latent of the draw and it is
read off the one quantity the comets measure for every comet, CO / H2O (row k5), whose spread across comets is unchanged:

1. **The temperature** from CO / H2O: `T = 48 K + 12.5 K x log10(0.035 / (CO / H2O))` (`formation_temperature_k`), the laboratory law through the Halley anchor, at 12.5 K per decade (six
   decades in 75 K). A comet at the median CO / H2O of 2.5 percent has T = 49.8 K, and one at Hale-Bopp's 0.2 has 38.5 K.
2. **The position** `xi` (`position`) of that temperature between the two ends of the comets measured: 0 at the median comet's temperature (the warm end, where 67P is) and 1 at
   Hale-Bopp's (the cold end), held at the ends. A comet warmer than the median has the warm-end values.
3. **The species** (`latent_species_over_h2o`). N2 and Ar are interpolated in the logarithm between the two ends of their brackets (`comet_species.bracket_over_h2o`: N2 per CO of 67P to that of
   C/2016 R2, Ar per water of 67P to that of Hale-Bopp) at position `xi`, plus a residual of their own: the part
   of a species that the temperature does not explain (a species varies within one comet: 67P's argon changes by 1.4 decades over its season). N2 is per CO, so it also follows the
   comet's own CO / H2O. Krypton is the comet's argon times the ratio of 67P and xenon its krypton times the ratio of 67P, each with the lab's fractionation (at most 40 percent) as its residual.
4. **What the temperature does not drive** stays an independent bracket of the comets measured (`comet_species.bracket_over_h2o`): OCS and SO2 (less volatile than the CO, trapped whatever the
   temperature) and O2 (the comets agree).

The residual widths are stand-ins (CY2). The mapping of the position onto the temperature ends is one number per end (CY1). The laboratory law is Bar-Nun's two published numbers (CY3);
Ciesla and others 2018 show the trapping depends on the deposition rate and the partial pressure as well, which is why the temperature here is a latent, not a measured quantity.

Tests against measurement (never inputs): the median comet gives the warm end of every species; a comet at Hale-Bopp's CO / H2O gives the cold end; the CO / H2O of every comet drawn is the k5 draw
(the latent does not change its spread); the species are correlated across comets the way the data say (a CO-poor, argon-rich comet is unlikely).
"""

import math

from planetgen.physics import comet_species

HALLEY_CO_OVER_H2O = 0.035
"""float: CO / H2O of Halley (3.5 percent), the laboratory anchor: the ice that traps it formed at `HALLEY_TEMPERATURE_K`. [Bar-Nun, NTRS N90-16610, R0267.]"""

HALLEY_TEMPERATURE_K = 48.0
"""float: The formation temperature the laboratory gives Halley's CO (48 +/- 5 K)."""

K_PER_DECADE = 75.0 / 6.0
"""float: How many kelvin of formation temperature change the trapped gas by a decade: it falls six decades from 25 to 100 K (Bar-Nun)."""

MEDIAN_CO_OVER_H2O = comet_species.CO_OVER_H2O_MEDIAN
"""float: CO / H2O of the median comet (row k5), the warm end of the position."""

COLD_CO_OVER_H2O = 0.2
"""float: CO / H2O of Hale-Bopp (about 0.2, 'notably high CO', Biver and others 1997; recalled, not opened), the cold end of the position. [Stand-in CY1, tier T0, tweak L1.]"""

RESIDUAL_DEX = {"N2": 0.05, "Ar": 0.4}
"""dict: The width, in decades, of what the temperature does not explain: N2 per CO of C/2016 R2 (0.06 to 0.09, a range of 0.18 dex, over the square root of 12) and the argon of 67P over its season
(1.4 decades, over the square root of 12). [Stand-ins CY2, tier T2, tweak L1.]"""

NOBLE_TIE_DEX = 0.15
"""float: The width, in decades, of krypton and xenon about argon: the laboratory's fractionation of at most 40 percent (Ligterink and others 2024). [Stand-in CY2, tier T2, tweak L1.]"""

LATENT_SPECIES = ("N2", "Ar", "Kr", "Xe")
"""tuple: The species the temperature drives."""

INDEPENDENT_SPECIES = ("OCS", "SO2", "O2")
"""tuple: The species it does not."""


def formation_temperature_k(co_over_h2o):
    """The formation temperature, in kelvin, of the ice of a comet with this CO / H2O (the laboratory law through Halley)."""
    if not co_over_h2o > 0.0 or not math.isfinite(co_over_h2o):
        raise ValueError(f"CO / H2O is above zero, not {co_over_h2o}")
    return HALLEY_TEMPERATURE_K + K_PER_DECADE * math.log10(HALLEY_CO_OVER_H2O / co_over_h2o)


def co_over_h2o_at(temperature_k):
    """The CO / H2O that ice formed at `temperature_k` holds: the inverse of `formation_temperature_k`."""
    return HALLEY_CO_OVER_H2O * 10.0 ** ((HALLEY_TEMPERATURE_K - temperature_k) / K_PER_DECADE)


def position(co_over_h2o):
    """0 at the median comet's formation temperature and 1 at Hale-Bopp's, held at the ends: where the comet lies between the warm and the cold end of the comets measured."""
    warm = formation_temperature_k(MEDIAN_CO_OVER_H2O)
    cold = formation_temperature_k(COLD_CO_OVER_H2O)
    return min(1.0, max(0.0, (warm - formation_temperature_k(co_over_h2o)) / (warm - cold)))


def latent_species_over_h2o(co_over_h2o, residuals=None, n_over_c=1.0, ar_over_c=1.0, s_over_c=1.0):
    """
    N2, Ar, Kr and Xe per H2O of a comet whose CO / H2O is `co_over_h2o`, round a host whose N/C, Ar/C and S/C are the given multiples of the proto-solar ones.

    Args:
        co_over_h2o (float): The comet's CO / H2O (the k5 draw); its formation temperature and position follow from it.
        residuals (Mapping): A standard normal for each of `LATENT_SPECIES` (0 for any missing): N2 and Ar about what the temperature gives, Kr about argon times the ratio of 67P and Xe about krypton times the ratio
            of 67P (`NOBLE_TIE_DEX`).

    Returns:
        dict: `{"N2": .., "Ar": .., "Kr": .., "Xe": ..}`.
    """
    z = {name: 0.0 for name in LATENT_SPECIES}
    z.update(residuals or {})
    unknown = set(z) - set(LATENT_SPECIES)
    if unknown:
        raise ValueError(f"{', '.join(sorted(unknown))} has no residual here; the species are {', '.join(LATENT_SPECIES)}")
    xi = position(co_over_h2o)

    def along(species):
        low, high = comet_species.bracket_over_h2o(species, n_over_c, ar_over_c, s_over_c)
        return math.log10(low) + xi * math.log10(high / low)

    n2 = along("N2") + math.log10(co_over_h2o / MEDIAN_CO_OVER_H2O) + RESIDUAL_DEX["N2"] * z["N2"]      # per CO at that position, so per water follows the comet's own CO
    argon = along("Ar") + RESIDUAL_DEX["Ar"] * z["Ar"]
    krypton = argon + math.log10(comet_species.KR_OVER_AR) + NOBLE_TIE_DEX * z["Kr"]
    xenon = krypton + math.log10(comet_species.XE_OVER_KR) + NOBLE_TIE_DEX * z["Xe"]
    return {"N2": 10.0 ** n2, "Ar": 10.0 ** argon, "Kr": 10.0 ** krypton, "Xe": 10.0 ** xenon}
