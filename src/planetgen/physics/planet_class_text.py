# planetgen/physics/planet_class_text.py

"""
The one-line description of each planet class (UX.252, GEN.1212)
================================================================

"a small, round, airless world", "an ocean world, its surface mostly liquid water": the phrase the planet page and the system map
put after "It is". The words are the `phrase` of each row of `physics/data/planet_classes.json` (GEN.1432: one row holds a class's phrase and its page sentence), read when a page is built, so the text is the page's own and
a body keeps no copy of it (GEN.1212 removes the stored `description` column and the text in `tuning.PLANET_CLASSES`).

Use:

    description("O")                # "an ocean world, its surface mostly liquid water"
    description(None)               # the Asteroid type, which has no letter
    description("O", moon=True)     # the same with "planet" written as "moon"; a class text is written generically
"""

import re

from planetgen.physics import planet_classes


class PlanetClassTextError(ValueError):
    """A class asked for does not exist."""


DESCRIPTIONS = {letter: planet_classes.BY_LETTER[letter].phrase for letter in planet_classes.LETTERS}
"""dict: Class letter -> its description phrase (lower case, starting with its article), the `phrase` of its row in `planet_classes.json`. `ASTEROID_PHRASE` is
the phrase of the Asteroid type, which has no letter (GEN.1160)."""

ASTEROID_PHRASE = planet_classes.ASTEROID.phrase


def description(planet_class, moon=False):
    """
    The description phrase of a class; with `moon` the words "planet" and "planets" read "moon" and "moons".

    Raises:
        PlanetClassTextError: There is no such class.
    """
    if planet_class is None:
        return ASTEROID_PHRASE
    try:
        text = DESCRIPTIONS[planet_class]
    except KeyError:
        raise PlanetClassTextError(f"no description for planet class {planet_class!r}") from None
    if moon:
        text = re.sub(r"\bplanets\b", "moons", text)
        text = re.sub(r"\bplanet\b", "moon", text)
    return text
