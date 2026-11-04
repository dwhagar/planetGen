# planetgen/physics/ice_giant_heat.py

"""
The spread of ice-giant internal heat (GEN.1033, registry 6.166)
================================================================

The cooling law of `rogue_surface.giant_internal_flux_w_m2` (GEN.816) runs through Jupiter and Neptune, so it gives a Uranus-mass giant
Neptune's heat, 0.32 W/m2, where Uranus radiates 0.042 (+0.047, -0.042) W/m2 against Neptune's 0.433 +/- 0.046 ([R5096: Table 2], from
[R3466, R5097]). The two planets are alike in mass, radius and composition, so a spread of ten in their heat needs a physical reason.

**The reason.** An ice giant cools adiabatically (its whole interior convects) or it does not. A steep compositional gradient between its
hydrogen-helium envelope and its ice-rich inner part stops the convection there and leaves a thermal boundary layer: the heat of the inner
part leaks through it slowly and the planet looks cold from outside while the inside stays hot ([R5098]: a gradual composition distribution
"acts as a thermal boundary to suppress convection and slow down the internal cooling", enough to give Uranus's luminosity, and stable for
billions of years; [R5101]: even a thin conducting layer of a few kilometres changes the cooling a lot, and Uranus's luminosity is reproduced
only if it has been near equilibrium with the sunlight for long). An adiabatic Uranus cannot be as faint as it is ([R5100]: cooling times of
5.1 Gyr for Uranus and 3.7 for Neptune bracket the age and neither present luminosity is adiabatic). What decides which interior a giant
keeps is its history: a violent impact mixes it ([R5099]: a head-on collision of a massive projectile can mix the interior to an adiabatic
temperature profile, which may explain Neptune's larger flux; an oblique one tilts Uranus and leaves it layered).

**The draw.** An ice giant (a core that found no gas, `ICE_GIANT_MAX_EARTH`) either keeps a stable gradient, with the chance `STRATIFIED_SHARE`,
or is mixed. A mixed giant radiates the cooling law as it is. A layered one radiates a share of it drawn log-uniform from
`SUPPRESSION_RANGE`. The draw is made on a stream of the giant's own (its mass), so every reader of the heat (the one-bar temperature, the
helium rain, the air) gets the same number and the run's stream does not move. Giants above the ice limit, brown dwarfs and rocky bodies
are not touched.

Use:

    share = suppression(mass_earth)     # 1.0 for a gas giant or a mixed ice giant
"""

from planetgen import tuning
from planetgen.util import draw

ICE_GIANT_MAX_EARTH = tuning.GIANT_NO_GAS_CORE_MAX_EARTH
"""float: The heaviest ice giant, Earth masses (a core that found no gas stops growing here, `giant_draw.ICE_GIANT_MAX_EARTH`)."""

STRATIFIED_SHARE = 0.5
"""float: The chance an ice giant keeps a stable compositional gradient. One of the two ice giants measured has one (Uranus; [R5096: Table 2]), so
the share is the rule of succession's (1 + 1) / (2 + 2). The impact history that mixes a giant (`giant_tilt`) is not drawn from this. Registry
row IGH1 (source tier T6, tweak level L0)."""

SUPPRESSION_RANGE = (0.03, 0.30)
"""tuple: The share of the cooling law's flux a layered ice giant radiates, log-uniform. Uranus radiates 0.042 +0.047 -0.042 W/m2 ([R5096: Table 2]) of the
0.32 W/m2 the law gives at its mass, 0.13 of it and from nothing to 0.28; the range is 0.01 to 0.095 W/m2 there, inside that uncertainty. Registry row
IGH2 (source tier T2, tweak level L1)."""

SEED_DIGITS = 9
"""int: Significant digits of the mass (Earth masses) the giant's own stream is seeded from, so the mass the planet holds and the mass its
interior sums to (equal to rounding) give one stream."""


def _share(key):
    with draw.bound(f"ice-giant-heat:{key}"):
        if draw.random() >= STRATIFIED_SHARE:
            return 1.0
        low, high = SUPPRESSION_RANGE
        return low * (high / low) ** draw.random()


def suppression(mass_earth):
    """
    The share of the cooling law's flux a giant of `mass_earth` radiates: 1 for a gas giant, a mixed ice giant or a mass with no giant in it,
    the drawn share for an ice giant that keeps a stable compositional gradient.
    """
    if mass_earth is None or not 0.0 < mass_earth <= ICE_GIANT_MAX_EARTH:
        return 1.0
    return _share(f"{mass_earth:.{SEED_DIGITS}g}")
