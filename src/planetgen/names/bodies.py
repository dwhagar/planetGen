# planetgen/names/bodies.py

"""
Derived names for the stars, planets and moons inside one star system.

Only the system itself draws a generated name (and only it goes through
`nameUniqueness`'s database-wide collision handling). Everything inside
the system is named from it:

- A single star shares the system's name: `Voranthis`.
- A close pair's two stars are the system name with a letter (GEN.62):
  `Voranthis A`, `Voranthis B`, shown together as `Voranthis A / B`.
- A wide pair's two stars are two words each (GEN.62): the system name's
  first word, then the star's own word: `Voranthis Kelmoor` and
  `Voranthis Pikkita`. The second star's word is drawn from
  `names.COMPANION_STAR_NAMES` (words for small, little, daughter, son and
  child), so it sounds like a diminutive.
- Planets are numbered in orbit order with roman numerals: `Voranthis I`,
  `Voranthis II`. A single star's and a close pair's planets use the
  system name. A wide pair's primary's planets use the shared first word
  (`Voranthis I`) and its secondary's planets use the secondary's own word
  (`Pikkita I`). Asteroid belts take no number. Planet names can repeat
  across systems, as town names do (Boss, 2026-10-02).
- Moons add a lowercase letter to their planet's numeral: `Voranthis IIa`,
  `Voranthis IIb`.

Since the system name is unique, every derived name is too, so planets and
moons never need a database-wide uniqueness search. `rename_prefix` keeps
those derived names in step when a system or star is renamed later.
"""

from planetgen.names.wordlists import (COMPANION_STAR_NAMES, COMPANION_STAR_PREFIXES, COMPANION_STAR_SUFFIXES, STAR_NAMES,
                    STAR_PREFIXES, STAR_SUFFIXES)
from planetgen.names.wordsalad import generate_phoneme_salad_name

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


COMPANION_STAR_WORD_MAX_LENGTH = 10
"""int: The longest a wide pair's second star word may be, in letters."""


def generate_star_word(exclude=(), companion=False):
    """One generated single-word name for a binary's star, distinct from
    every word in `exclude`. `companion=True` draws a wide pair's second
    star's word from the small/little/child word lists (GEN.62).
    `allow_split=False` only stops a long result being split in two; a base
    name that already has a space (STAR_NAMES' "El Nath") can still come
    through, so a draw with any whitespace is rejected too."""
    if companion:
        # Kept short, as a diminutive is.
        lists = (COMPANION_STAR_NAMES, COMPANION_STAR_PREFIXES, COMPANION_STAR_SUFFIXES)
        max_length = COMPANION_STAR_WORD_MAX_LENGTH
    else:
        lists = (STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)
        max_length = None
    while True:
        word = generate_phoneme_salad_name(*lists, allow_split=False, max_length=max_length,
                                           syllable_fraction=0.5 if companion else 1.0)
        if word not in exclude and not any(c.isspace() for c in word):
            return word


CLOSE_PAIR_LETTERS = ("A", "B")
"""tuple: The letters after the system name for a close pair's primary and
secondary star (GEN.62)."""


def wide_pair_first_word(system_name):
    """The word a wide pair's two stars share, and its primary's planets
    are named for: the system name's first word (GEN.62)."""
    return system_name.split()[0]


def close_pair_label(primary_name, secondary_name):
    """
    The two stars of a close pair named together: `Voranthis A / B` for
    `Voranthis A` and `Voranthis B`, or `<primary> & <secondary>` for names
    that don't share that shape (a pair renamed by hand, or stored before
    GEN.62).
    """
    primary_base, _, primary_letter = primary_name.rpartition(" ")
    secondary_base, _, secondary_letter = secondary_name.rpartition(" ")
    if (primary_base and primary_base == secondary_base
            and (primary_letter, secondary_letter) == CLOSE_PAIR_LETTERS):
        return f"{primary_base} {primary_letter} / {secondary_letter}"
    return f"{primary_name} & {secondary_name}"


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
