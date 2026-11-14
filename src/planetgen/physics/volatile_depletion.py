# planetgen/physics/volatile_depletion.py

"""
Volatility depletion of a rocky body against CI chondrites (GEN.1361, N3 of the 8.1 real-material set)
=====================================================================================================

A rocky body keeps less of an element the more volatile it is. Earth's bulk composition set against the CI chondrites (both normalised to
magnesium) shows it for 26 elements: the refractory ones (Ca, Al, Ti, V, Sr, Y, Zr) sit at the CI value, the moderately volatile ones fall off
with the 50 percent condensation temperature (T50, the temperature at which half of the element has condensed from a solar-composition gas at
1e-4 bar) from sulfur (T50 664 K, a tenth kept) up to manganese and phosphorus (about 1150 to 1230 K, 60 to 70 percent kept). The bulk Earth is
[R4608] (READ), the CI chondrites are [R0168] (READ) and the temperatures are [R2726] (READ); `data/volatile_depletion.json` holds the 26 rows.

The law is a hinge in log space:

    log10 f = SLOPE * (T50 - cut)   for T50 below the cut,     f = 1   at or above it,

with f the share of the CI abundance (against magnesium) the body keeps. The least-squares hinge to the 26 rows has SLOPE 0.00165 per K, cut
1316.6 K and an rms of 0.079 dex; the research set quotes the same line as -2.150 + 0.001623 T50 (cut 1324.7 K), which differs from it by under
0.013 dex between 660 and 1300 K. It reproduces the Earth's sodium, potassium and zinc untuned (the table residuals are in the tests).

What this module adds beyond the fit is a lever for a body that did not form where Earth did. The cut is the temperature below which an element is lost,
so a hotter accretion zone loses more and a colder one less: `cut_temperature_k` scales the Earth cut with the blackbody temperature of the formation
distance in the light of the star (`condensate_shares.disk_temperature_k`, the `d^-1/2` law). Mars (1.52 AU, cut 1066 K) then keeps 2.6 times the Earth law's
sodium, where the research note has about twice; beyond about 4 AU (a cut under the T50 of the element) an element is kept whole. The ices and the noble gases
have their own condensation items (N4), not this law, which holds for T50 from `MIN_T50_EARTH_K` (500 K) up whatever the cut.

The Earth fit is the central case and the lever is a stand-in with a counterexample on record: Mercury's K/U (14,311 against the bulk silicate Earth's 9,704,
[R1003]) says the innermost planet kept more potassium than the Earth, where the lever gives it a cut near 2,100 K and a twentieth of the Earth law's potassium.
Mercury's mantle was stripped by a giant impact or formed from material the Sun had not heated as the blackbody does, and the Moon is far poorer than any
distance scaling gives, which is the same impact history: neither is a target of this law. A caller that wants the Earth law for every body passes no cut. The mantle step
of GEN.1362 (`mantle_oxides`, called by `bulk_draw.dry_partition`) reads it with the cut of each body's formation zone.
"""

import json
import math
import os

from planetgen.physics import condensate_shares, elements
from planetgen.physics.elements import DataError, positive_number
from planetgen.util.checks import finite_domain

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "volatile_depletion.json")
"""str: The data file."""


def load(path=DATA_FILE):
    """Reads and checks the file at `path`: `(law, reference zone, rows by symbol)`; a malformed row is refused whole."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise DataError(f"cannot read the volatile depletion table {path}: {error}") from error

    positive = positive_number      # GEN.1667: the shared check

    law = document.get("law") or {}
    law = {"slope_per_k": positive(law.get("slope_per_k"), "the slope"), "cut_k": positive(law.get("cut_k"), "the cut"),
           "min_t50_earth_k": positive(law.get("min_t50_earth_k"), "the lowest T50")}
    reference = document.get("reference") or {}
    zone = (positive(reference.get("distance_au"), "the reference distance"), positive(reference.get("luminosity_w"), "the reference luminosity"))
    rows = {}
    for raw in document.get("elements") or []:
        symbol = raw.get("symbol")
        if not isinstance(symbol, str) or not symbol or symbol in rows:
            raise DataError(f"a volatile depletion row needs a new element symbol, not {symbol!r}")
        ratio = positive(raw.get("earth_over_ci"), f"the Earth over CI of {symbol}")
        logged = raw.get("log10_depletion")
        if not isinstance(logged, (int, float)) or not math.isfinite(logged) or abs(logged - math.log10(ratio)) > 0.02:
            raise DataError(f"the logarithm of {symbol} is not the logarithm of its ratio {ratio}")
        rows[symbol] = (positive(raw.get("t50_k"), f"the T50 of {symbol}"), ratio, float(logged), raw.get("host_phase", ""))
    if len(rows) < 20:
        raise DataError("the volatile depletion table has fewer than 20 elements")
    return law, zone, rows


LAW, REFERENCE_ZONE, ROWS = load()
"""The law's constants, the reference zone `(AU, W)` and the rows `{symbol: (T50 in K, Earth over CI, log10 of it, host phase)}`."""

SLOPE_PER_K = LAW["slope_per_k"]
"""float: log10 of the share kept per K of T50 below the cut, 1/K."""

CUT_K = LAW["cut_k"]
"""float: The Earth fit's cut, K: an element with a T50 at or above it is kept whole."""

MIN_T50_EARTH_K = LAW["min_t50_earth_k"]
"""float: The lowest T50 the Earth fit is used for, K (the coolest fitted element, sulfur, is 664 K); below it lie the ices and the noble gases."""

REFERENCE_BLACKBODY_K = condensate_shares.disk_temperature_k(*REFERENCE_ZONE)
"""float: The blackbody temperature of the zone the Earth fit is made for, K (the Sun's light at 1 AU: 278.3 K)."""


@finite_domain()
def cut_temperature_k(distance_au, luminosity_w):
    """
    The cut for a body formed `distance_au` from a star of `luminosity_w`, K: the Earth cut scaled with the blackbody temperature of the zone,
    `CUT_K * T(d, L) / REFERENCE_BLACKBODY_K`.
    """
    return CUT_K * condensate_shares.disk_temperature_k(distance_au, luminosity_w) / REFERENCE_BLACKBODY_K


@finite_domain()
def log_retention(t50_k, cut_k=CUT_K):
    """
    log10 of the share of an element kept against CI (magnesium-normalised) when the cut is `cut_k` K: `SLOPE_PER_K * (t50_k - cut_k)` below the cut and 0
    from it up. A T50 under `MIN_T50_EARTH_K` is refused: the law is fitted to the moderately volatile elements (the coolest is 664 K) and the ices
    and the noble gases have their own items.
    """
    if not (t50_k > 0.0 and cut_k > 0.0):
        raise ValueError("the T50 and the cut must be positive")
    if t50_k < MIN_T50_EARTH_K:
        raise ValueError(f"a T50 of {t50_k:g} K is under {MIN_T50_EARTH_K:g} K: the law holds for moderately volatile elements, the ices and the "
                         f"noble gases have their own items")
    return min(0.0, SLOPE_PER_K * (t50_k - cut_k))


@finite_domain()
def retained_fraction(t50_k, cut_k=CUT_K):
    """The share of the CI abundance (against magnesium) an element of 50 percent condensation temperature `t50_k` keeps; at most 1."""
    return 10.0 ** log_retention(t50_k, cut_k)


def t50_k(symbol):
    """The 50 percent condensation temperature of one of the table's elements, K."""
    return _row(symbol)[0]


def earth_over_ci(symbol):
    """The Earth over CI ratio (both normalised to magnesium) the table holds for an element: what the law is fitted to."""
    return _row(symbol)[1]


def retained_fractions(cut_k=CUT_K):
    """`{symbol: share kept}` for the table's 26 elements at a cut of `cut_k` K, coolest element first."""
    return {symbol: retained_fraction(row[0], cut_k) for symbol, row in sorted(ROWS.items(), key=lambda item: item[1][0])}


def _row(symbol):
    try:
        return ROWS[symbol]
    except KeyError:
        raise DataError(f"{symbol!r} is none of the table's elements: {', '.join(ROWS)}") from None
