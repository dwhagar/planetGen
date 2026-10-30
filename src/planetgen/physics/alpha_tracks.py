# planetgen/physics/alpha_tracks.py

"""
Alpha-enhanced stars on solar-scaled tracks (GEN.744, part 1).

The MIST isochrones of `datasets.isochrones` are solar-scaled ([alpha/Fe] = 0); the thick disc and the halo are alpha-enhanced. A
star rich in the alpha elements (O, Mg, Si, Ca, Ti) burns and looks like a solar-scaled star of a higher [Fe/H]: the alpha elements
carry most of the opacity and the CNO catalysts, and the same total metal content is reached with less iron. The rule of Salaris,
Chieffi and Straniero 1993 (ApJ 414, 580) is

    [M/H] = [Fe/H] + log10(0.638 f_alpha + 0.362),        f_alpha = 10 ** [alpha/Fe],

with 0.638 the alpha elements' share of the metals by number in the solar mixture. A star of the high-alpha branch at [Fe/H] with
[alpha/Fe] = +0.4 therefore reads the grid at [Fe/H] + 0.29. `equivalent_fe_h` gives the [Fe/H] a track lookup takes.

The star's own [alpha/Fe] is the mean of the [O/Fe], [Mg/Fe], [Si/Fe] and [Ca/Fe] of its branch at its [Fe/H]
(`stellar_abundances.mean_xfe`, the hinge fits of Bensby et al. 2014 and Nissen 2014 that `composition` already uses), so the track
and the composition of one star agree. Nothing in generation calls this yet: the wiring is GEN.744's second step.
"""

import math

from planetgen.physics import stellar_abundances

ALPHA_SYMBOLS = ("O", "Mg", "Si", "Ca")
"""tuple: The alpha elements that the branch fits hold; their mean is the star's [alpha/Fe]."""

ALPHA_SHARE = 0.638
"""float: The alpha elements' share of the metals in the solar mixture (Salaris et al. 1993). [Stand-in AL1, tier T2, tweak L0: read
second-hand (Salaris and Cassisi 2005), not in a paper opened here.]"""


def alpha_fe(fe_h, branch):
    """The [alpha/Fe] (dex) of a star of this branch (`stellar_abundances.BRANCHES`) at [Fe/H] `fe_h`."""
    return sum(stellar_abundances.mean_xfe(symbol, fe_h, branch) for symbol in ALPHA_SYMBOLS) / len(ALPHA_SYMBOLS)


def metallicity_shift(alpha_fe_dex):
    """The [M/H] excess over [Fe/H] of a star of this [alpha/Fe]: `log10(0.638 f_alpha + 0.362)` (Salaris et al. 1993)."""
    return math.log10(ALPHA_SHARE * 10.0 ** alpha_fe_dex + (1.0 - ALPHA_SHARE))


def equivalent_fe_h(fe_h, branch):
    """The [Fe/H] of the solar-scaled track a star of this branch and [Fe/H] is read from."""
    return fe_h + metallicity_shift(alpha_fe(fe_h, branch))
