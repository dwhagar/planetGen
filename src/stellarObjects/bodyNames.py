# stellarObjects/bodyNames.py

"""
Derived names for the stars, planets and moons inside one star system.

Only the system itself draws a generated name (and only it goes through
`nameUniqueness`'s database-wide collision handling). Everything inside
the system is named from it:

- A single star shares the system's name: `Voranthis`.
- A binary's two stars each put their own generated word after the system
  name, with no A/B letters: `Voranthis Kelmoor`, `Voranthis Ostra`.
- Planets are numbered in orbit order with roman numerals after the star
  they orbit: `Voranthis I`, `Voranthis II`. A close pair's circumbinary
  planets orbit both stars, so they use the system name; a wide pair's
  planets use their own star's name (`Voranthis Kelmoor I`). Asteroid belts
  take no number.
- Moons add a lowercase letter to their planet's numeral: `Voranthis IIa`,
  `Voranthis IIb`.

Since the system name is unique, every derived name is too, so planets and
moons never need a database-wide uniqueness search. `rename_prefix` keeps
those derived names in step when a system or star is renamed later.
"""

from .names import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from .utils import generate_phoneme_salad_name

_ROMAN_PAIRS = (
    (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"), (100, "C"), (90, "XC"),
    (50, "L"), (40, "XL"), (10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I"),
)


def to_roman(number):
    """`number` (a positive int) as a roman numeral, e.g. `4 -> "IV"`."""
    if number < 1:
        raise ValueError(f"roman numerals start at 1, got {number}")
    parts = []
    for value, numeral in _ROMAN_PAIRS:
        count, number = divmod(number, value)
        parts.append(numeral * count)
    return "".join(parts)


def moon_letters(index):
    """A moon's letter(s) from its 0-based orbit index: `0 -> "a"`,
    `25 -> "z"`, `26 -> "aa"`."""
    letters = ""
    index += 1
    while index:
        index, remainder = divmod(index - 1, 26)
        letters = chr(ord("a") + remainder) + letters
    return letters


def generate_star_word(exclude=()):
    """One generated single-word name for a binary's star, distinct from
    every word in `exclude`."""
    while True:
        word = generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES, allow_split=False)
        if word not in exclude:
            return word


def name_bodies(prefix, bodies):
    """
    Names every planet in `bodies` (one star's orbit-ordered planet/belt
    list) `<prefix> <numeral>`, and each planet's moons
    `<prefix> <numeral><letter>`. Belts are skipped and don't use up a
    number.
    """
    number = 0
    for body in bodies:
        if body.body_type == "a":
            continue
        number += 1
        numeral = to_roman(number)
        body.name = f"{prefix} {numeral}"
        for index, moon in enumerate(body.moons):
            moon.name = f"{prefix} {numeral}{moon_letters(index)}"


def rename_prefix(name, old_prefix, new_prefix):
    """
    `name` with a leading `old_prefix` swapped for `new_prefix`, or `None`
    when `name` doesn't carry it. Matches the whole name or the prefix
    followed by a space, so `Voranthis II` follows a rename of `Voranthis`
    but `Voranthisa` does not. Names someone set by hand usually don't
    start with the old prefix, so they're left alone.
    """
    if name == old_prefix:
        return new_prefix
    if name is not None and name.startswith(old_prefix + " "):
        return new_prefix + name[len(old_prefix):]
    return None
