# planetgen/generation/comet_abundance.py

"""
The comet abundance of a galaxy (GEN.1775; Boss, 8:43 PM EST 2026-11-16: "Use physics but provide a scale option like our M Class planet and Rogue planet in generation with a default of Science")
==============================================================================================================================================================================================

A star's own comets are the planetesimals of its outer disk that a giant threw about (`physics/comet_supply`, GEN.1765). On `science`, the default, a star has comets exactly where the generator's own laws
put them: where it holds an outer-disk belt that a giant scattered (`comet_supply.family_weights` is not `None`), and nowhere else. The census of GEN.1765 found nearly every star without one, so a galaxy
made on `science` has few comets.

A galaxy may be told to have more, with one word, the same shape as the Class M abundance (`abundance_tier`) and the Rogue planets level (`topup_stage`):

| Tier | A star with no scattered outer disk holds comets in |
|---|---|
| `science` (the default) | none: only the stars the laws give them |
| `sparse` | 1 in 12 (0.25 of the hand-set chance) |
| `scattered` | 1 in 6 (0.5) |
| `plentiful` | 3 in 10 (1, the chance the generator used before the tier existed, `tuning.SYSTEM_COMET_CHANCE`) |
| `abundant` | 6 in 10 (2) |

A tier is a multiplier on the registered stand-in `tuning.SYSTEM_COMET_CHANCE` (`tuning.COMET_TIER_MULTIPLIERS`), the one hand-set quantity of how many stars hold comets; it never changes what the
physics gives. A star the physics gives comets keeps them on every tier, and how many it holds (`tuning.SYSTEM_COMET_COUNT_RANGE`) and the shape of each are as they were.

**Science is the physics.** On `science` no multiplier acts and a star without a scattered outer disk has no comets, so the comet count of a galaxy is what its disks give.

The tier reaches a system through its config (`SystemConfig.COMET_TIER`, a class attribute set on a config only for a tier other than science), which a galaxy sector sets from the galaxy's stored
choice (`bright_stars.FillContext.apply_population`), as it does the Class M tier.
"""

from planetgen import tuning

TIERS = ("science", "sparse", "scattered", "plentiful", "abundant")
"""tuple[str]: Every tier, from the physics up: the choices of the `comet_abundance` lever (`galaxy.dials`), in its order."""

DEFAULT = TIERS[0]
"""str: The tier of a galaxy that was told none, and of work outside any galaxy."""


def tier_of(subject):
    """
    The tier `subject` is made under: a `SystemConfig`, or anything that carries one as `system_config` (a planet, a star system).
    `DEFAULT` when it carries none.

    Raises:
        ValueError: The config names a tier that does not exist.
    """
    config = getattr(subject, "system_config", subject)
    tier = getattr(config, "COMET_TIER", None) or DEFAULT
    if tier not in TIERS:
        raise ValueError(f"{tier!r} is not one of {', '.join(TIERS)}")
    return tier


def multiplier(subject):
    """The factor on `tuning.SYSTEM_COMET_CHANCE` of the tier `subject` is made under; 0 on `science` (the stand-in does not act)."""
    tier = tier_of(subject)
    return 0.0 if tier == DEFAULT else tuning.COMET_TIER_MULTIPLIERS[tier]


def presence_chance(subject, has_reservoir):
    """
    The chance a star has native comets: 1 for a star with a scattered outer disk (the physics), else the registered stand-in `tuning.SYSTEM_COMET_CHANCE` times the tier's multiplier (0 on `science`).

    Args:
        subject: A `SystemConfig` or what carries one (`tier_of`).
        has_reservoir (bool): Whether the star holds an outer-disk belt a giant scattered (`comet_supply.family_weights(belts) is not None`).
    """
    if has_reservoir:
        return 1.0
    return min(1.0, tuning.SYSTEM_COMET_CHANCE * multiplier(subject))
