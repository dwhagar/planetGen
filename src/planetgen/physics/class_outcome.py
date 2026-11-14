# planetgen/physics/class_outcome.py

"""
Whether a class counts as an outcome: the match a body must have for a request or a count (GEN.1315)
====================================================================================================

Every body that some class is eligible for takes a class, fully or fudged (`class_fit`, GEN.1314), and stores the `class_match`: 1 for a full fit, less for a fudged one. The class letter
alone then says little about a poor fit, so what *reads* the class (a request for a class, the count of Class M worlds, the search) reads the match with it: a body fills the outcome
when its match is at least `tuning.CLASS_MATCH_THRESHOLD` (registry row FT5). The outcomes that were letter lists are already conditions: the habitable count and the bar on habitable worlds
read the PHI-6 score and the life of the finished body (`habitability.is_habitable`, `habitable_request.carries_life`, GEN.1387, GEN.1397), so they need no threshold here.
"""

from planetgen import tuning


def fits(match):
    """Whether a stored `class_match` is enough for the class to count as an outcome (`None`, a body with no class, is not)."""
    return match is not None and match >= tuning.CLASS_MATCH_THRESHOLD


def has_class(body, letter):
    """Whether `body` has the class `letter` with a match that counts (a request for that class is met by it, a count of that class includes it)."""
    return getattr(body, "planet_class", None) == letter and fits(getattr(body, "class_match", None))
