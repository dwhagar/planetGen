# planetgen/galaxy/scatter_limits.py

"""
The scatter limits by place (DB.54, docs/design/realistic-galactic-core.md 5.1)
==============================================================================

The galaxy's star scatter draws its line at one luminosity floor and one mass limit (the baseline,
`galaxy_shape.bright_star_min_luminosity_sol` and `bright_star_mass_limit_sol`). Three more settings of the galaxy move
that line by place: a luminosity floor and a mass limit for the core and a luminosity floor for the halo (the halo has
no mass limit of its own). The core's two are `None` until set, and then follow the baseline times
`tuning.CORE_TIER_FACTOR`; the halo's is a number (default `tuning.HALO_MIN_LUMINOSITY_SOL`, never under
`tuning.HALO_MIN_LUMINOSITY_RANGE_SOL[0]`). Which place takes which limit is GEN.333's; this module is the settings.
"""

from planetgen import tuning


def baseline_min_luminosity_sol(level_sol):
    """The baseline luminosity floor: the scatter's level when one has run, else the tuning's default."""
    return tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL if level_sol is None else float(level_sol)


def baseline_mass_limit_sol(limit_sol):
    """The baseline mass limit: the scatter's when it had a mass pass, else the tuning's default."""
    return tuning.PHENOMENON_MIN_MASS_SOLAR if limit_sol is None else float(limit_sol)


def core_min_luminosity_sol(skeleton, level_sol):
    """The core's luminosity floor: the one set, else the baseline (`level_sol`, the scatter's level or `None`) times
    `tuning.CORE_TIER_FACTOR`."""
    if skeleton.core_min_luminosity_sol is not None:
        return float(skeleton.core_min_luminosity_sol)
    return baseline_min_luminosity_sol(level_sol) * tuning.CORE_TIER_FACTOR


def core_mass_limit_sol(skeleton, limit_sol):
    """The core's mass limit: the one set, else the baseline (`limit_sol`, the scatter's or `None`) times
    `tuning.CORE_TIER_FACTOR`."""
    if skeleton.core_mass_limit_sol is not None:
        return float(skeleton.core_mass_limit_sol)
    return baseline_mass_limit_sol(limit_sol) * tuning.CORE_TIER_FACTOR


def halo_min_luminosity_sol(skeleton):
    """The halo's luminosity floor."""
    return float(skeleton.halo_min_luminosity_sol)


def check_halo_min_luminosity(value):
    """Raises `ValueError` naming the reason when `value` is not a number inside the halo floor's range."""
    low, high = tuning.HALO_MIN_LUMINOSITY_RANGE_SOL
    try:
        ok = low <= float(value) <= high
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise ValueError(f"the halo luminosity floor must be from {low:g} to {high:g} solar luminosities, "
                         f"got {value!r}")
