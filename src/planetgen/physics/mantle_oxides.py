# planetgen/physics/mantle_oxides.py

"""
The seven oxides a rock stores besides SiO2, MgO, CaO, Al2O3 and FeO (GEN.1362, N5 of the 8.1 real-material set)
===============================================================================================================

A real bulk silicate holds Na2O, K2O, TiO2, Cr2O3, MnO, NiO and P2O5 as well. Without them the CIPW norm of a mantle has no albite, orthoclase, ilmenite,
chromite or apatite and the felsic recipe cannot come out of the derived chain. This module gives the seven as functions of what the generator already
has: the host star's atoms per magnesium (GEN.1360 carries Na, K, Ti, Cr, Mn and P with Ni), the volatility of the element (GEN.1361,
`volatile_depletion`) and the share of the element that went to the core.

For one magnesium atom of the body, an element X with star ratio r (atoms of X per Mg):

1. **The bulk.** `r * f`, with f the share of the CI abundance a body keeps (`volatile_depletion.retained_fraction` of the element's T50 at the cut of the body's
   formation zone). Ti, Cr and Ni are at or near 1; Na and K keep about a quarter at the Earth cut, Mn 0.6 and P 0.6.
2. **The core.** Of that, the share `s` is in the metal: a metal and silicate partition `s = D c / (D c + 1 - c)` with c the core's share of the dry
   mass (`core_fraction`) and `D` the coefficient that gives the Earth's own share at the Earth's own core (`EARTH_CORE_SHARE` at c = `EARTH_CORE_FRACTION`;
   [R4608], the core ppm over the bulk ppm times the core mass fraction): Cr 0.55, Mn 0.60, Ni 0.93, P 0.94. A bigger core takes
   a bigger share of each, and Ti, which is lithophile, takes none. Na and K carry the 5 percent the research note caps them at (the Earth's concordance
   table has 20 and 30 percent from one old model). D moves with how oxidised the mantle is: the element goes into the metal as a neutral atom and out of the
   silicate as an oxide of valence n, so D scales as the mantle's FeO to the power minus n/2 (`VALENCE`; Ni, Cr and Mn 2, P 5), set to the Earth's D at the Earth's
   FeO (`EARTH_ROCK_FEO`). A reduced body (Mercury) puts nearly all its P in the metal and a more oxidised one (Mars, FeO 0.147) keeps four times the Earth's P
   in the mantle, which is how the model reaches Mars's phosphorus. The law is the leading term only: it leaves out the pressure, the temperature and the sulfur
   of the real metal, which is why Mars's NiO (0.046 percent against the model's 0.3) is not reached.
3. **The mantle.** What is not in the core is the rock's oxide, `n (1 - s) M(oxide) / (atoms of X in the oxide)` grams per Mg atom.

`oxides(ratios, cut_k, core_fraction)` returns the grams per Mg atom of the seven oxides of the rock and the grams of each element the core holds; `bulk_draw.dry_partition`
adds the first to the rock, the second to the metal, and takes the nickel the mantle keeps out of the metal. The calibration targets are the bulk silicate Earth (tested
within a factor of 1.5 of the 12-oxide table) and Mars (within a factor of 2.5): neither is an input, the Earth's core shares are the anchor of D.
"""

from planetgen.physics import elements, volatile_depletion

OXIDE_OF = {
    "Na": ("Na2O", 2), "K": ("K2O", 2), "Ti": ("TiO2", 1), "Cr": ("Cr2O3", 2), "Mn": ("MnO", 1), "Ni": ("NiO", 1), "P": ("P2O5", 2),
}
"""dict: Element -> (oxide, atoms of the element in one formula unit of it)."""

SYMBOLS = tuple(OXIDE_OF)
"""tuple: The seven elements, in the order the stored columns have them."""

OXIDE_NAMES = tuple(oxide for oxide, _ in OXIDE_OF.values())
"""tuple: The seven oxides, in the order of `SYMBOLS`."""

EARTH_CORE_FRACTION = 0.325
"""float: The Earth's core mass fraction of its whole mass ([R4608]: 32.5 +- 0.3 percent), the c at which `EARTH_CORE_SHARE` was measured."""

EARTH_ROCK_FEO = 0.0805
"""float: The FeO mass fraction of the Earth's bulk silicate (McDonough and Sun 1995: 8.05 percent; the [R4608] concordance has 8.1), the FeO at which
`EARTH_CORE_SHARE` was measured. [Stand-in MO3, with `VALENCE`.]"""

VALENCE = {"Ti": 0, "Cr": 2, "Mn": 2, "Ni": 2, "P": 5, "Na": 0, "K": 0}
"""dict: The valence n of the oxide each element makes, which sets how D falls as the mantle's FeO rises: D goes as FeO^(-n/2) (the exchange of a metal atom for n/2 oxygen).
Cr is taken as 2 (Cr2+ at the oxygen fugacity of core formation; Cr3+ would give 3), the alkalis and Ti as 0 (Na and K have a fixed cap; Ti is lithophile)."""

EARTH_CORE_SHARE = {"Ti": 0.0, "Cr": 0.552, "Mn": 0.601, "Ni": 0.929, "P": 0.942, "Na": 0.05, "K": 0.05}
"""dict: The share of each element that is in the Earth's core (1): the core ppm times 0.325 over the bulk ppm, from the concordance table of Wang, Lineweaver and
Ireland 2018 (Cr 6590 / 3883, Mn 3257 / 1760, Ni 50600 / 17700, P 2774 / 957 ppm). Ti has no core ppm. Na and K are capped at 5 percent, the research note's cap on
the core fractions of the old model (Na 0.20 and K 0.30 there). [Stand-ins MO1 and MO2.]"""

_MOLAR = {symbol: elements.element(symbol).atomic_weight for symbol in SYMBOLS}


def _distribution_coefficient(symbol):
    """The metal over silicate partition `D` that gives `EARTH_CORE_SHARE` at `EARTH_CORE_FRACTION`."""
    share = EARTH_CORE_SHARE[symbol]
    if share <= 0.0:
        return 0.0
    c = EARTH_CORE_FRACTION
    return share / (1.0 - share) * (1.0 - c) / c


DISTRIBUTION = {symbol: _distribution_coefficient(symbol) for symbol in SYMBOLS}
"""dict: Element -> D, the metal over silicate partition ratio by mass (0 for Ti)."""


MIN_ROCK_FEO = 1e-3
"""float: The lowest FeO the redox scaling of D takes (a body with none would put everything siderophile in the metal)."""


def core_share(symbol, core_fraction, rock_feo=EARTH_ROCK_FEO):
    """The share of `symbol` that is in the metal of a body whose core is `core_fraction` of its dry mass and whose rock holds FeO at `rock_feo`: `D c / (D c + 1 - c)`."""
    if not 0.0 <= core_fraction <= 1.0:
        raise ValueError("the core fraction must be in [0, 1]")
    d = DISTRIBUTION[symbol] * (max(rock_feo, MIN_ROCK_FEO) / EARTH_ROCK_FEO) ** (-0.5 * VALENCE[symbol])
    total = d * core_fraction + 1.0 - core_fraction
    return d * core_fraction / total if total > 0.0 else 0.0


def bulk_atoms(ratios, symbol, cut_k=None):
    """Atoms of `symbol` per Mg atom of the whole body: the star's ratio times the share the zone's cut keeps (the Earth's cut when `cut_k` is None)."""
    ratio = getattr(ratios, symbol.lower())
    kept = volatile_depletion.retained_fraction(volatile_depletion.t50_k(symbol), volatile_depletion.CUT_K if cut_k is None else cut_k)
    return ratio * kept


def oxides(ratios, cut_k, core_fraction, rock_feo=EARTH_ROCK_FEO):
    """
    The seven oxides for one Mg atom of a body.

    Args:
        ratios (HostRatios): The star's atoms per Mg for `na`, `k`, `ti`, `cr`, `mn`, `ni` and `p`.
        cut_k (float or None): The volatility cut of the formation zone, K (`volatile_depletion.cut_temperature_k`); None is the Earth's.
        core_fraction (float): The core's share of the dry mass.
        rock_feo (float): The FeO mass fraction of the rock.

    Returns:
        tuple: `(rock, core, nickel_core)`: grams of each oxide of `OXIDE_NAMES` in the rock, grams of each element but nickel in the core, and the grams of nickel in the core.
    """
    rock, core = {}, {}
    for symbol in SYMBOLS:
        atoms = bulk_atoms(ratios, symbol, cut_k)
        share = core_share(symbol, core_fraction, rock_feo)
        oxide, count = OXIDE_OF[symbol]
        rock[oxide] = atoms * (1.0 - share) * (elements.molar_mass(oxide) * 1000.0) / count
        core[symbol] = atoms * share * _MOLAR[symbol]
    nickel = core.pop("Ni")
    return rock, core, nickel
