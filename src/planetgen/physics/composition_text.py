# planetgen/physics/composition_text.py

"""
What a body is made of, in words, from its stored bulk fractions (GEN.430; docs/design/composition/bulk-composition.md, section on descriptions)
================================================================================================================================================

One writer for the composition phrase of every planet and moon. The fractions are the whole-body shares of
`bulk_composition.BulkFractions`; each material the body holds at least `MIN_SHARE` of is named with the level its share falls in
(`LEVELS`, fixed thresholds), largest first. The same call serves a rocky world, an icy moon and a giant, so the phrase and the stored
fractions cannot disagree (a body's class no longer carries a fixed composition string).

    composition_text.text(body)   # "mostly silicate rock, with a large share of iron-rich metal and a small share of water"
"""

from planetgen.physics import bulk_composition

LEVELS = ((0.5, "mostly"), (0.2, "a large share of"), (0.05, "a moderate share of"), (0.01, "a small share of"))
"""tuple: `(least share, words)`, largest first: the level words of a material by its share of the whole body's mass. A material under the
last threshold is not named."""

MIN_SHARE = LEVELS[-1][0]
"""float: The smallest share of the body's mass that its text names."""

POROUS = 0.1
"""float: A small body with at least this void share of its volume is called loosely packed."""

NAMES = {
    "metal": "iron-rich metal",
    "water": "water",
    "ice": "volatile ices",
    "envelope": "hydrogen-helium gas",
}
"""dict: What each whole-body fraction is called; the rock's name depends on its kind (`rock_name`)."""


def rock_name(carbon_rock):
    """The rock's name: silicate, or silicon carbide and graphite for a carbon-rich body."""
    return "silicon carbide and graphite" if carbon_rock else "silicate rock"


def level_words(share):
    """The words of the level `share` falls in, or `None` below `MIN_SHARE`."""
    for least, words in LEVELS:
        if share >= least:
            return words
    return None


def parts(bulk):
    """
    The named materials of `bulk` as `(name, share, words)`, largest first.

    Args:
        bulk (BulkFractions): A body's bulk composition.
    """
    named = []
    for fraction, share in bulk.shares().items():
        words = level_words(share)
        if words is not None:
            name = rock_name(bulk.carbon_rock) if fraction == "rock" else NAMES[fraction]
            named.append((name, share, words))
    return sorted(named, key=lambda item: -item[1])


def _join(items):
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + f" and {items[-1]}"


def text(body):
    """
    The composition phrase of `body`, from its stored bulk fractions.

    Returns:
        str or None: `None` for a body with no stored composition.
    """
    bulk = bulk_composition.bulk_fractions(body)
    if bulk is None:
        return None
    groups = {}
    for name, _share, words in parts(bulk):
        groups.setdefault(words, []).append(name)
    phrases = [f"{words} {_join(names)}" for words, names in groups.items()]
    phrase = phrases[0] + (f", with {_join(phrases[1:])}" if len(phrases) > 1 else "")
    if bulk.porosity is not None and bulk.porosity >= POROUS:
        phrase += ", loosely packed"
    return phrase
