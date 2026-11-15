# planetgen/physics/heavy_elements.py

"""
The heavy elements of a star and of the rock and metal of its planets (GEN.1371, docs/design/composition/stellar-composition.md)
================================================================================================================================

The s- and r-process split of the elements beyond iron, and where those elements go in a differentiated body. Everything is a pure function of the star's stored [Fe/H] and age (and its
own abundance key) and of the body's stored rock and metal shares, so nothing is stored: the object sheet works it out when it is made (`describe/heavy_rows.py`), the way `radiogenic`
works the thorium of a body from the same star.

**The star** (`birth_heavy`). The elements of `tuning.STELLAR_HEAVY_LOG_EPS` (Sr, Y, Zr, Ba, the rare earths La to Lu, Os, Ir, Pt, Au, Pb) are the A21 solar abundance moved to proto-solar and to the
star's [Fe/H], each shared between the two processes by the r-fraction of the Sun (`tuning.STELLAR_HEAVY_R_FRACTION`, [R5228] as recalled): the r share follows the star's
[Eu/Fe] (the branch's mean and scatter of `stellar_abundances`, the draw the thorium of `radiogenic` uses), the s share follows [s/Fe], the Brewer and Fischer trend of Y
(`tuning.STELLAR_S_PROCESS_XFE`) and a scatter of its own from a stream named by the star (`heavy:{key}`, so no other draw moves):

    A_e = A_sun,e * 10^([Fe/H]) * ( f_r,e * 10^[Eu/Fe] + (1 - f_r,e) * 10^[s/Fe] )

Eu, Th and U are `stellar_abundances.birth_radiogenic`'s and are not repeated here.

**The rock and the metal** (`bulk_ppm`). A refractory element goes into the planet with the rock, so its mass is the star's element-to-magnesium ratio times the magnesium of the rock
(`ROCK_MG_MASS_FRACTION`, stand-in HE3); it is then shared between the silicate and the metal core by its metal-silicate partition coefficient `D` (`PARTITION`, rows HE4 to HE7):

    share in the core = D f_metal / (D f_metal + f_rock)        (equilibrium at the core-forming conditions, all the metal and all the rock in contact)

so the platinum-group metals and gold (D of hundreds) are nearly all in the core and the rare earths and the large-ion lithophiles (D far under one) stay in the mantle. Lead is moderately
siderophile and is not depleted by its volatility here (stand-in HE7: a hot rocky planet that lost its lead is not told apart). The late accretion that restores a mantle's platinum-group
metals (`late_accretion`) is not added: the mantle value is the equilibrium one.

Stand-ins are HE1 to HE7 of section 6.274 of `docs/design/research/gap-filling-plan.md`.
"""

import functools
import math

from planetgen import tuning
from planetgen.physics import elements, stellar_abundances as sa, stellar_composition as sc
from planetgen.util import draw

ELEMENTS = tuple(tuning.STELLAR_HEAVY_LOG_EPS)
"""tuple: The heavy elements of the split, by atomic number."""

GROUPS = {
    "rare earths": ("La", "Ce", "Pr", "Nd", "Sm", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"),
    "platinum group": ("Os", "Ir", "Pt"),
    "gold": ("Au",),
    "large-ion": ("Sr", "Y", "Zr", "Ba"),
    "lead": ("Pb",),
}
"""dict: The groups the sheet reports, by the elements they hold."""

PARTITION = {"rare earths": 0.005, "platinum group": 1000.0, "gold": 1000.0, "large-ion": 0.001, "lead": 3.0}
"""dict: The metal-silicate partition coefficient `D` (mass fraction in the metal over mass fraction in the silicate at equilibrium) of each group at core formation. The platinum group and gold are
the 1,000 that leaves the Earth's mantle 0.3 percent of the chondritic ratio (tweak L2, tuned to that depletion, rows HE4 and HE5); the rare earths and the large-ion lithophiles are lithophile (far under
one, T5, L0, HE6); lead is weakly siderophile (3, T4, L0, HE7)."""

ROCK_MG_MASS_FRACTION = 0.222
"""float: The magnesium mass fraction of the silicate (rock) of a body: the bulk silicate Earth's, 22.2 percent [[R3417] via the research note; second hand]. Stand-in HE3, tier T3, tweak L0."""

_ATOMIC_WEIGHT = {symbol: elements.element(symbol).atomic_weight for symbol in ELEMENTS}


def s_process_xfe(fe_h):
    """The mean [s/Fe] (dex) at [Fe/H] `fe_h`: the Brewer and Fischer line for Y, with [Fe/H] held to the thin disk's fit range (`tuning.STELLAR_FE_H_FIT_RANGE`)."""
    a, b, _sigma = tuning.STELLAR_S_PROCESS_XFE
    low, high = tuning.STELLAR_FE_H_FIT_RANGE[sa.BRANCH_THIN]
    return a + b * min(max(fe_h, low), high)


def birth_heavy(star, branch, eu_normal=0.0, s_normal=0.0):
    """
    The birth abundances of the heavy elements of `star` (a `StellarComposition`): `{"abundances": {symbol: atoms per hydrogen atom}, "mass_fractions": {symbol: mass fraction}, "eu_fe", "s_fe"}`.
    [Eu/Fe] is the branch's mean plus `eu_normal` of its scatter (the same as `birth_radiogenic`'s) and [s/Fe] is `s_process_xfe` plus `s_normal` of `tuning.STELLAR_S_PROCESS_XFE`'s sigma.
    """
    eu_fe = sa.mean_xfe("Eu", star.fe_h, branch) + eu_normal * tuning.STELLAR_XFE_SIGMA[branch]["Eu"]
    s_fe = s_process_xfe(star.fe_h) + s_normal * tuning.STELLAR_S_PROCESS_XFE[2]
    r_scale, s_scale = 10.0 ** eu_fe, 10.0 ** s_fe
    abundances = {}
    for symbol, log_eps in tuning.STELLAR_HEAVY_LOG_EPS.items():
        r_share = tuning.STELLAR_HEAVY_R_FRACTION[symbol]
        abundances[symbol] = 10.0 ** (log_eps + tuning.STELLAR_PROTO_METAL_OFFSET_DEX - 12.0 + star.fe_h) * (r_share * r_scale + (1.0 - r_share) * s_scale)
    hydrogen = elements.element("H").atomic_weight
    mass_fractions = {symbol: star.x * number * _ATOMIC_WEIGHT[symbol] / hydrogen for symbol, number in abundances.items()}
    return {"abundances": abundances, "mass_fractions": mass_fractions, "eu_fe": eu_fe, "s_fe": s_fe}


def _s_normal(key):
    """The scatter draw of [s/Fe] of the star named `key`: its own stream, so no other draw of the star moves."""
    return draw.Stream(f"heavy:{key}").gauss(0.0, 1.0)


def _to_mg(found, star):
    return {symbol: fraction / star.elements["Mg"] for symbol, fraction in found["mass_fractions"].items()}


@functools.cache
def _solar_to_mg():
    """Each heavy element's mass over magnesium's in the solar-scale star: [Fe/H] 0 on the thin branch with every residual and scatter at its mean."""
    star = sa.composition(0.0, sa.abundance_offsets(0.0, sa.BRANCH_THIN, (0.0,) * len(sa.RESIDUAL_ELEMENTS)))
    return _to_mg(birth_heavy(star, sa.BRANCH_THIN), star)


@functools.lru_cache(maxsize=65536)
def star_to_mg(fe_h, age_gy, key):
    """
    Each heavy element's mass over magnesium's in the star of [Fe/H] `fe_h` and age `age_gy` named `key` (`stellar_composition.abundance_key`): `{symbol: ratio}`. The same star always gives the same
    ratios because its draws come from its own streams.
    """
    draws = sc.key_draws(key)
    branch = sa.branch_of(fe_h, age_gy)
    offsets = sa.abundance_offsets(fe_h, branch, draws.normals, draws.u_enhanced, draws.u_carbon)
    star = sa.composition(fe_h, offsets)
    return _to_mg(birth_heavy(star, branch, draws.eu_normal, _s_normal(key)), star)


def offset_dex(fe_h, age_gy, key):
    """Each heavy element's offset from the solar-scale star, dex: log10 of the star's element-to-magnesium ratio over the solar scale's."""
    ratios, solar = star_to_mg(fe_h, age_gy, key), _solar_to_mg()
    return {symbol: math.log10(ratios[symbol] / solar[symbol]) for symbol in ELEMENTS}


def partition_of(symbol):
    """The metal-silicate partition coefficient `D` of the element `symbol` (its group's)."""
    for group, members in GROUPS.items():
        if symbol in members:
            return PARTITION[group]
    raise KeyError(f"{symbol!r} is not a heavy element of the split")


def core_share(d, metal_fraction, rock_fraction):
    """The share of an element in the metal core of a body with the given mass shares of metal and rock, for the partition coefficient `d`: `D f_m / (D f_m + f_r)`. 0 for a body with no metal."""
    if metal_fraction <= 0.0 or rock_fraction <= 0.0:
        return 0.0
    return d * metal_fraction / (d * metal_fraction + rock_fraction)


def bulk_ppm(ratios, rock_fraction, metal_fraction):
    """
    The heavy elements of a body of the given mass shares of rock and metal, parts per million by mass: `{symbol: (silicate ppm, metal ppm)}`. `ratios` is `star_to_mg` (the element over magnesium in
    the star); the body's mass of each is that ratio times the magnesium of its rock, then shared by `core_share`. A body with no metal has all of it in the silicate, one with no rock none anywhere.
    """
    found = {}
    for symbol in ELEMENTS:
        mass = ratios[symbol] * ROCK_MG_MASS_FRACTION * rock_fraction          # per unit mass of the body
        in_core = core_share(partition_of(symbol), metal_fraction, rock_fraction)
        silicate = mass * (1.0 - in_core) / rock_fraction * 1.0e6 if rock_fraction > 0.0 else 0.0
        metal = mass * in_core / metal_fraction * 1.0e6 if metal_fraction > 0.0 else 0.0
        found[symbol] = (silicate, metal)
    return found


def group_ppm(per_element):
    """The sum of each group of `bulk_ppm`'s elements: `{group: (silicate ppm, metal ppm)}`."""
    return {group: (sum(per_element[s][0] for s in members), sum(per_element[s][1] for s in members)) for group, members in GROUPS.items()}
