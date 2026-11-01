# planetgen/physics/alpha_tracks.py

"""
Alpha-enhanced stars on solar-scaled tracks (GEN.744, part 1).

The MIST isochrones of `datasets.isochrones` are solar-scaled ([alpha/Fe] = 0); the thick disc and the halo are alpha-enhanced. A
star rich in the alpha elements (O, Mg, Si, Ca, Ti) burns and looks like a solar-scaled star of a higher [Fe/H]: the alpha elements
carry most of the opacity and the CNO catalysts, and the same total metal content is reached with less iron. The rule of [R3878] is

    [M/H] = [Fe/H] + log10(0.638 f_alpha + 0.362),        f_alpha = 10 ** [alpha/Fe],

with 0.638 the alpha elements' share of the metals by number in the solar mixture. A star of the high-alpha branch at [Fe/H] with
[alpha/Fe] = +0.4 therefore reads the grid at [Fe/H] + 0.29. `equivalent_fe_h` gives the [Fe/H] a track lookup takes.

The star's own [alpha/Fe] is the mean of the [O/Fe], [Mg/Fe], [Si/Fe] and [Ca/Fe] of its branch at its [Fe/H]
(`stellar_abundances.mean_xfe`, the hinge fits of [R0349, R3294] that `composition` already uses), so the track
and the composition of one star agree. Nothing in generation calls this yet: the wiring is GEN.744's second step.
"""

import functools
import math

from planetgen import tuning
from planetgen.physics import stellar_abundances

ALPHA_SYMBOLS = ("O", "Mg", "Si", "Ca")
"""tuple: The alpha elements that the branch fits hold; their mean is the star's [alpha/Fe]."""

ALPHA_SHARE = 0.638
"""float: The alpha elements' share of the metals in the solar mixture [R3878]. [Stand-in AL1, tier T2, tweak L0: read
second-hand [R3877], not in a paper opened here.]"""


def alpha_fe(fe_h, branch):
    """The [alpha/Fe] (dex) of a star of this branch (`stellar_abundances.BRANCHES`) at [Fe/H] `fe_h`."""
    return sum(stellar_abundances.mean_xfe(symbol, fe_h, branch) for symbol in ALPHA_SYMBOLS) / len(ALPHA_SYMBOLS)


ZERO_POINT_DEX = alpha_fe(0.0, stellar_abundances.BRANCH_THIN)
"""float: The thin branch's fitted [alpha/Fe] at [Fe/H] = 0 (-0.017 dex). The grid is solar-scaled, so a thin-disc star of solar
[Fe/H] has no alpha excess by definition; `equivalent_fe_h` reads the excess from this zero point, not from 0, so that star sits on the
reference track and only the excess the fits give above it shifts a star (the fits' own offset at the Sun is inside their scatter)."""


def metallicity_shift(alpha_fe_dex):
    """The [M/H] excess over [Fe/H] of a star of this [alpha/Fe]: `log10(0.638 f_alpha + 0.362)` [R3878]."""
    return math.log10(ALPHA_SHARE * 10.0 ** alpha_fe_dex + (1.0 - ALPHA_SHARE))


def equivalent_fe_h(fe_h, branch):
    """The [Fe/H] of the solar-scaled track a star of this branch and [Fe/H] is read from."""
    return fe_h + metallicity_shift(alpha_fe(fe_h, branch) - ZERO_POINT_DEX)


def track_fe_h(fe_h, age_gy):
    """The [Fe/H] of the solar-scaled track a star of this birth [Fe/H] and age reads (GEN.1024): its [Fe/H] raised by its branch's
    alpha excess (`stellar_abundances.branch_of` finds the branch, the same one its composition has)."""
    return equivalent_fe_h(fe_h, stellar_abundances.branch_of(fe_h, age_gy))


@functools.lru_cache(maxsize=1 << 16)
def track_segments(fe_h):
    """
    The pieces of the age axis on which a star of this birth [Fe/H] reads one track, as `(start_gy, end_gy, track [Fe/H])`. A star's branch
    is fixed by its own uniform number and changes with age only at the thick disc's age gate (`tuning.STELLAR_THICK_START_GY`), where the
    high-alpha chance leaves 0, so the star reads the thin branch's track below it and `track_fe_h` at it above, which is the same piece
    when the star is thin. The scatter's windows and bounds (`star_population`) are drawn on these pieces so they mirror `evolve_star`.
    """
    gate = tuning.STELLAR_THICK_START_GY
    young = equivalent_fe_h(fe_h, stellar_abundances.BRANCH_THIN)
    old = track_fe_h(fe_h, gate)
    if old == young:
        return ((0.0, math.inf, young),)
    return ((0.0, gate, young), (gate, math.inf, old))


@functools.lru_cache(maxsize=1 << 16)
def possible_segments(fe_h):
    """
    Every piece of the age axis and track a star of this birth [Fe/H] could read whatever its branch: the thin track below the age gate, and
    above it the thin and the high-alpha track. A bound over a table's [Fe/H] bin holds for any star in it (`star_population`), and a star of
    the bin has its own branch, so the bound is taken over both.
    """
    gate = tuning.STELLAR_THICK_START_GY
    young = equivalent_fe_h(fe_h, stellar_abundances.BRANCH_THIN)
    old = equivalent_fe_h(fe_h, stellar_abundances.BRANCH_HIGH_ALPHA)
    return ((0.0, gate, young), (gate, math.inf, young), (gate, math.inf, old))
