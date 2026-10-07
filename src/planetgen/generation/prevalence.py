# planetgen/generation/prevalence.py

"""
Prevalence: how much more or less often a sector or galaxy run's systems
get a feature than they would by chance (GEN.52). The forcing options
(`+comets`, `-habitable_world`) force every system, which a whole sector
can't use (GEN.48); a prevalence is a percentage deviation from the
feature's normal chance instead: +50 gives 1.5 times the usual share of
systems with it, -100 none, and 0 (the default) leaves it alone.

A feature decided by one draw (comets, a binary companion, a wide pair,
a planet's moons, the most orbits a star can hold, and intelligent life:
the population pass's `CIVILIZATION_CHANCE`) has that draw's chance
scaled (`scaled_chance`). A feature that comes out of many draws (a
habitable world, an asteroid belt, a large star, any planets at all) has
no single chance to scale, so its usual share is
measured (`tuning.PREVALENCE_BASE_SHARES`) and each system is forced to
have or lack it often enough to move that share by the percentage
(`resolve`), through the same tri-state flags the forcing options set.
"""

import random

from planetgen import tuning

FEATURES = ("habitable_world", "asteroid_belt", "comets", "large_star", "moons", "max_planets",
            "intelligent_life", "binary_system", "wide_binary", "planets")
"""tuple[str]: Every feature with a prevalence: the forcing options' names."""

MEASURED_FEATURES = tuple(tuning.PREVALENCE_BASE_SHARES)
"""tuple[str]: The features `resolve` forces on or off by their measured
usual share; the rest scale their own draw (`scaled_chance`)."""

MIN_PERCENT = -100.0
"""float: -100% is never; anything lower would be a negative chance."""


def percent(config, feature):
    """`config`'s prevalence for `feature`, in percent (0 when unset)."""
    return float((getattr(config, "PREVALENCE", None) or {}).get(feature) or 0.0)


def scaled_chance(config, feature, chance):
    """
    `chance` (a probability) moved by `config`'s prevalence for `feature`:
    +50 makes it 1.5 times as likely, -100 never, capped at certain.

    Returns:
        float: The probability to draw against.
    """
    return min(1.0, max(0.0, chance * (1.0 + percent(config, feature) / 100.0)))


def resolve(config, rng=random):
    """
    Decides, for one system, the measured features its prevalences move:
    each one still left to chance (`None`) is forced off in `-p`% of
    systems for a negative prevalence, or forced on in just enough of the
    rest for a positive one that the share with it comes to `1 + p/100`
    times its usual share (all of them past 100%). The rules the forcing
    options follow hold for what this decides: a habitable world and a
    belt together need a large star, and either needs planets. Flags
    already forced are left alone. Runs once per config.

    Args:
        config (SystemConfig): The system's config, changed in place.
        rng: The random source.
    """
    if getattr(config, "_prevalence_resolved", False):
        return
    config._prevalence_resolved = True
    decided = {}
    for feature in MEASURED_FEATURES:
        p = percent(config, feature)
        attr = feature.upper()
        if not p or getattr(config, attr) is not None:
            continue
        if p < 0:
            if rng.random() < -p / 100.0:
                decided[attr] = False
        else:
            base = tuning.PREVALENCE_BASE_SHARES[feature]
            if rng.random() < min(1.0, (p / 100.0) * base / (1.0 - base)):
                decided[attr] = True
    if not decided:
        return

    def flag(attr):
        return decided[attr] if attr in decided else getattr(config, attr)

    # A habitable world and a belt together need a large star.
    if flag("HABITABLE_WORLD") is True and flag("ASTEROID_BELT") is True:
        if decided.get("LARGE_STAR") is False:
            del decided["LARGE_STAR"]
        if flag("LARGE_STAR") is False:
            _undo(decided, ("ASTEROID_BELT",) if "ASTEROID_BELT" in decided else ("HABITABLE_WORLD",))
        else:
            decided["LARGE_STAR"] = True
    # Either needs planets.
    if flag("PLANETS") is False and (flag("HABITABLE_WORLD") is True or flag("ASTEROID_BELT") is True):
        if "PLANETS" in decided:
            del decided["PLANETS"]
        else:
            _undo(decided, ("HABITABLE_WORLD", "ASTEROID_BELT"))
    for attr, value in decided.items():
        setattr(config, attr, value)


def _undo(decided, attrs):
    """Drops `resolve`'s decisions on `attrs` (those it made)."""
    for attr in attrs:
        decided.pop(attr, None)
