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


from planetgen import tuning
from planetgen.util import draw

FEATURES = ("habitable_world", "asteroid_belt", "comets", "large_star", "moons", "max_planets",
            "intelligent_life", "binary_system", "wide_binary", "planets")
"""tuple[str]: Every feature with a prevalence: the forcing options' names."""

MEASURED_FEATURES = tuple(tuning.PREVALENCE_BASE_SHARES)
"""tuple[str]: The features `resolve` forces on or off by their measured
usual share; the rest scale their own draw (`scaled_chance`)."""

MIN_PERCENT = -100.0
"""float: -100% is never; anything lower would be a negative chance."""

USUAL_SHARES = {**tuning.PREVALENCE_BASE_SHARES, **tuning.PREVALENCE_DRAW_SHARES,
                "intelligent_life": tuning.CIVILIZATION_CHANCE}
"""dict: Each feature's usual share (0-1) with default options, of systems
or of what `tuning.PREVALENCE_DRAW_SHARES` says (intelligent life: of
worlds whose timeline reaches a technological age)."""


def percent_for_share(feature, share):
    """
    The prevalence (percent from the usual chance) that moves `feature`
    from its usual share (`USUAL_SHARES`) to `share` (0-1): its usual
    share gives 0, none -100, twice it +100.
    """
    return (share / USUAL_SHARES[feature] - 1.0) * 100.0


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


def resolve(config, rng=draw):
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


STAR_MIX = ("single", "close_pair", "wide_pair")
"""tuple[str]: The kinds of system by stars (ADM.45): one star, a close
binary, a wide pair. Their shares of systems always total 100%."""

STAR_MIX_TOTAL_TOLERANCE = 0.01
"""float: How far from 100 (percentage points) a star mix may total, for the
rounding of a typed share."""


def usual_star_mix():
    """The usual share (0-1) of each `STAR_MIX` kind, from the usual binary
    share and the usual share of binaries that are wide pairs."""
    binary = USUAL_SHARES["binary_system"]
    wide = USUAL_SHARES["wide_binary"]
    return {"single": 1.0 - binary, "close_pair": binary * (1.0 - wide), "wide_pair": binary * wide}


def star_mix_prevalences(single, close_pair, wide_pair):
    """
    The prevalences (percent from usual, as `--prevalence` takes them) that
    give systems the shares of stars asked for (ADM.45): `binary_system` is
    the close and wide pairs together, `wide_binary` the wide ones' share of
    those.

    Args:
        single, close_pair, wide_pair (float): Percent of systems, each at
            least 0, totalling 100.

    Returns:
        dict: `binary_system` and `wide_binary` -> percent.

    Raises:
        ValueError: A negative share, or a total that is not 100.
    """
    shares = (single, close_pair, wide_pair)
    if any(share < 0 for share in shares):
        raise ValueError("a share of systems can't be below 0%")
    total = sum(shares)
    if abs(total - 100.0) > STAR_MIX_TOTAL_TOLERANCE:
        verb, amount = ("lower", total - 100.0) if total > 100.0 else ("raise", 100.0 - total)
        raise ValueError(f"the star mix totals {total:.4g}%, not 100%: {verb} one of single, close pair or "
                         f"wide pair by {amount:.4g} points")
    binary = (close_pair + wide_pair) / 100.0
    wide_of_binaries = wide_pair / (close_pair + wide_pair) if binary > 0 else USUAL_SHARES["wide_binary"]
    return {"binary_system": percent_for_share("binary_system", binary),
            "wide_binary": percent_for_share("wide_binary", wide_of_binaries)}
