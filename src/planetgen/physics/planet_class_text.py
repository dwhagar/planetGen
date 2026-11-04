# planetgen/physics/planet_class_text.py

"""
The one-line description of each planet class (UX.252, GEN.1212)
================================================================

"a small icy body", "an ocean world with most of its surface covered in liquid water": the phrase the planet page and the system map
put after "It is". The words are `physics/data/planet_class_text.json`, read when a page is built, so the text is the page's own and
a body keeps no copy of it (GEN.1212 removes the stored `description` column and the text in `tuning.PLANET_CLASSES`).

Use:

    description("O")                # "an ocean world with most of its surface covered in liquid water"
    description(None)               # the Asteroid type, which has no letter
    description("O", moon=True)     # the same with "planet" written as "moon"; a class text is written generically
"""

import json
import os
import re

from planetgen.physics import elements

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "planet_class_text.json")


class PlanetClassTextError(ValueError):
    """The text file is malformed, or a class asked for does not exist."""


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    descriptions = document["descriptions"]
    for letter, text in descriptions.items():
        if len(letter) != 1 or not letter.isupper() or not text.strip():
            raise PlanetClassTextError(f"planet_class_text.json: a class is one capital letter and has words, not {letter!r}: {text!r}")
    if not document["asteroid"].strip():
        raise PlanetClassTextError("planet_class_text.json: the asteroid type has no words")
    return descriptions, document["asteroid"]


DESCRIPTIONS, ASTEROID_PHRASE = _load()
"""dict: Class letter -> its description phrase (lower case, starting with its article). `ASTEROID_PHRASE` is the phrase of the Asteroid
type, which has no letter (GEN.1160)."""


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
