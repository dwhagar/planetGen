# planetgen/physics/rogue_pairs.py

"""
Rogue planets that stay a pair (GEN.332, slice D)
==================================================

Free-floating planetary-mass objects are seen in pairs. The JWST survey of the Trapezium [R3467] found 540 candidates of 0.7 to 13 Jupiter masses, 40 of
them in binaries and 2 in triplets, so about 9 percent of the systems of free planets are bound, at projected separations of 28 to 384 AU. The count is disputed ([R5261] and 2025 find
most of the pairs are reddened background sources), so 9 percent is an upper figure for a young cluster and not a settled one. A pair does not last: [[R5262]]
predict about 2 percent of the planet-planet pairs a cluster makes survive 10 Myr, so a pair of an old galaxy is a rare survivor. All three are read second hand (Boss's source rule, 2026-10-24).

**What this module gives.** The share of a young rogue's ejection that leaves it bound to a companion (`bound_share`), the separation of the two (`separation_au`: log-uniform inside the
observed range, held under the radius the galaxy's tide leaves a bound pair, `tidal_limit_au`) and the period they go round each other (`period_yr`, Kepler's third law, circular: no source
gives an eccentricity for wide young pairs, so none is drawn).

**The share.** `OBSERVED_SHARE * SURVIVAL_AT_REFERENCE ** (age / REFERENCE_MYR)`: a pair fades by the same factor in every 10 Myr, so a rogue thrown out a Gyr ago has no companion
(`0.02 ** 100` is 1e-170, which no draw reaches), and the rogues of a young donor are the only ones that do.

**The tide.** In the frame of a donor's circular orbit a body of mass m holds a companion inside `(G m / (4 Omega A)) ** (1/3)` (Oort's constants, `rogue_drift`), about 35,000 AU for two
Jupiters: far outside the observed range, so the cap only matters for a pair of a much lighter mass or the dense inner galaxy's own constants.

Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.291.
"""

import math

from planetgen.physics import constants, rogue_drift
from planetgen.util.checks import finite_domain

OBSERVED_SHARE = 0.09
"""float: The share of the systems of free planets that are bound pairs or triplets in the young Trapezium sample ([R3467], 42 of about 500 systems; also quoted as 8 percent). Source tier T3, tweak level L1."""
SURVIVAL_AT_REFERENCE = 0.02
"""float: The share of a cluster's planet-planet pairs still bound after `REFERENCE_MYR` ([R5262], second hand). Source tier T3, tweak level L0."""
REFERENCE_MYR = 10.0
"""float: The age at which `SURVIVAL_AT_REFERENCE` is given, Myr. Source tier T3, tweak level L0."""
SEPARATION_RANGE_AU = (28.0, 384.0)
"""tuple: The projected separations of the observed pairs, AU [R3467]. A projection is smaller than the true separation by a factor near one, which is not undone. Source tier T3, tweak level L0."""

_KM_S_KPC_TO_PER_S = 1.0e3 / (1.0e3 * constants.PARSEC_M)
"""float: One km/s/kpc in 1/s."""


class RoguePairError(ValueError):
    """An input the model cannot use."""


@finite_domain()
def survival(age_gy):
    """The share of a pair's birth still bound `age_gy` after the ejection: `SURVIVAL_AT_REFERENCE ** (age / REFERENCE_MYR)`."""
    if age_gy < 0.0:
        raise RoguePairError("the age cannot be negative")
    return SURVIVAL_AT_REFERENCE ** (age_gy * 1000.0 / REFERENCE_MYR)


@finite_domain()
def bound_share(age_gy):
    """The chance that a system of free planets `age_gy` after its ejection is a bound pair: the observed young share times `survival`."""
    return OBSERVED_SHARE * survival(age_gy)


@finite_domain()
def tidal_limit_au(total_mass_kg):
    """The radius inside which a pair of `total_mass_kg` stays bound against the galaxy's tide, AU: `(G m / (4 Omega A)) ** (1/3)`."""
    if not total_mass_kg > 0.0:
        raise RoguePairError("the mass must be positive")
    a = rogue_drift.OORT_A_KM_S_KPC * _KM_S_KPC_TO_PER_S
    b = rogue_drift.OORT_B_KM_S_KPC * _KM_S_KPC_TO_PER_S
    return (constants.G * total_mass_kg / (4.0 * (a - b) * a)) ** (1.0 / 3.0) / constants.AU_M


@finite_domain()
def separation_au(draw, total_mass_kg):
    """
    The separation of a pair, AU: log-uniform in `SEPARATION_RANGE_AU` at the uniform `draw`, held to `tidal_limit_au`.

    Returns:
        float: The separation, or `None` when the galaxy's tide would unbind a pair even at the nearest observed separation (no pair is made).
    """
    if not 0.0 <= draw < 1.0:
        raise RoguePairError("the draw must be in [0, 1)")
    low, high = SEPARATION_RANGE_AU
    limit = tidal_limit_au(total_mass_kg)
    if limit <= low:
        return None
    return min(low * (high / low) ** draw, limit)


@finite_domain()
def period_yr(separation_au_, mass_a_kg, mass_b_kg):
    """The period of two bodies on a circular orbit `separation_au_` apart, years: Kepler's third law, `2 pi sqrt(a^3 / (G (m1 + m2)))`."""
    if not separation_au_ > 0.0 or not mass_a_kg > 0.0 or not mass_b_kg > 0.0:
        raise RoguePairError("the separation and the masses must be positive")
    axis = separation_au_ * constants.AU_M
    return 2.0 * math.pi * math.sqrt(axis ** 3 / (constants.G * (mass_a_kg + mass_b_kg))) / constants.SECONDS_PER_YEAR
