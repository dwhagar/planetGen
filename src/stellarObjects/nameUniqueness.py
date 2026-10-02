# stellarObjects/nameUniqueness.py

"""
Name-uniqueness decoration -- pure functions, no database access.

Guarantees, across a whole database, that no two sectors share a name, no
two star systems share a name, and no sector/system pair shares a name.
Enforced as one consistent hierarchy: **sector > system**. Stars, planets
and moons are named from their system (`bodyNames.py` -- `Voranthis II`,
`Voranthis IIa`), so they're unique whenever the system is and never need
a search of their own.

- Two names at the *same* level colliding (sector-vs-sector or
  system-vs-system) is resolved by `resolve_greek_roman_collision`: the
  existing row is decorated with a Greek-letter prefix (`names.GREEK_LETTERS`),
  the new one with the next letter; once all 24 letters for that base name
  are used, decoration switches to a `"Alpha <name> <roman numeral>"`
  scheme (`names.ROMAN_NUMERAL_VALUES`); once that's exhausted too (36
  total occurrences of one base name -- astronomically unlikely), the
  caller must draw an entirely new base name instead.
- A system colliding with a sector (or vice versa) is resolved by
  `resolve_diminutive`: the *system* side (never the sector) gets a
  "small"/"little" word (`names.DIMINUTIVE_PREFIXES`) prefixed onto its
  name, advancing to the next word in the list on a repeat collision for
  the same base name.
- A new star system (or uniquely named phenomenon) name stays within
  `MAX_SYSTEM_NAME_WORDS` (two) words (GEN.46): a one-word base gets the
  Greek tier and one diminutive, a two-word base no decoration, and
  anything that would need a third word draws a fresh name instead.
  Sector names keep every tier.
`strip_decoration` is the inverse of both -- recovering the
underlying base name from an already-decorated one, e.g. for the
one-off backfill script (`src/dedupeNames.py`) to group existing rows.

None of this module talks to the database or knows about `sectors`/
`star_systems` -- see `stellarObjects/_db.py`'s `reserve_sector_name`/
`reserve_system_name` (and their `confirm_*` counterparts) for where
these functions actually get called, against the `sector_name_registry`/
`system_name_registry` tables that track each base name's own progress
through these schemes.
"""

from .names import (
    DIMINUTIVE_PREFIXES, GREEK_LETTERS, ROMAN_NUMERAL_VALUES, ROMAN_NUMERALS_BY_VALUE,
)

MAX_SYSTEM_NAME_WORDS = 2
"""int: The most words a newly decorated star system (or uniquely named
phenomenon) name may have (GEN.46, Boss 2026-10-01: "Name generation
should not produce star names that are more than 2 words long", so that
planet names built on them stay reasonable). A decoration that would go
past it is not used: the caller draws a fresh base name instead. Names
already stored keep whatever they have (Boss, 2026-10-02: no renaming
migration)."""


def word_count(name):
    """How many space-separated words `name` has."""
    return len(name.split())


def fits_word_limit(name, max_words):
    """Whether `name` has at most `max_words` words (`None`: no limit)."""
    return max_words is None or word_count(name) <= max_words


def word_limit_for(base_name, max_words=MAX_SYSTEM_NAME_WORDS):
    """
    The word limit decorations on `base_name` must keep to: `max_words`
    for a generated-shape name, `None` (no limit) for a base already
    longer than that -- a name given by hand ("Lonely Shell Star",
    `--name`), which keeps the old decorations rather than being swapped
    for a random one.
    """
    return max_words if fits_word_limit(base_name, max_words) else None


GREEK_ROMAN_CAPACITY = len(GREEK_LETTERS) + len(ROMAN_NUMERAL_VALUES)
"""int: How many total rows may ever share one base name under
`resolve_greek_roman_collision` (24 Greek letters + 12 roman-numeral
slots) before the caller must fall back to an entirely fresh base name."""

_ROMAN_SUFFIX_TOKENS = set(ROMAN_NUMERALS_BY_VALUE.values())
_GREEK_PREFIX_TOKENS = set(GREEK_LETTERS)
_DIMINUTIVE_PREFIX_TOKENS = set(DIMINUTIVE_PREFIXES)


def resolve_greek_roman_collision(base_name, existing_count, max_words=None):
    """
    Same-level collision resolution for two sectors, or two systems,
    that generated the same `base_name` -- see this module's own
    docstring for the full scheme.

    Args:
        base_name (str): The undecorated name both rows share.
        existing_count (int): How many rows already exist with this base
            name, *before* the row currently being inserted/renamed.
            `0` means no collision at all (this is the first-ever row).
        max_words (int, optional): The most words the new name, and the
            existing row's new name when one is renamed, may have
            (`MAX_SYSTEM_NAME_WORDS` for systems). A decoration that
            would go past it gives `(None, None)`, as an exhausted base
            name does: the caller draws a fresh one. So with a limit of
            2, a one-word base gets the Greek tier only and a two-word
            base no decoration at all. A base already longer than the
            limit (given by hand) is decorated as before.

    Returns:
        tuple: `(new_name, rename)`.
            `new_name` (str or None): The name to give the row currently
                being inserted -- `None` once `GREEK_ROMAN_CAPACITY` (36)
                is reached, meaning every decoration slot for this base
                name is taken and the caller must draw an entirely fresh
                base name and call this function again from scratch
                (`existing_count=0`) against that new name.
            `rename` (tuple or None): `(old_name, new_name)` for an
                *already-existing* row that must be renamed to make room
                for the one being inserted (the sole bare-named row on
                the first collision; that same row again, from its
                Greek-only "Alpha <base>" form to "Alpha <base> I", the
                one time the Greek tier hands off to the roman-numeral
                tier). `None` every other time -- no existing row needs
                touching, the new row's own decoration is enough.
    """
    new_name, rename = _greek_roman_collision(base_name, existing_count)
    if max_words is not None and not fits_word_limit(base_name, max_words):
        max_words = None  # a hand-given name (`word_limit_for`)
    if new_name is not None and not (
        fits_word_limit(new_name, max_words) and (rename is None or fits_word_limit(rename[1], max_words))
    ):
        return None, None
    return new_name, rename


def _greek_roman_collision(base_name, existing_count):
    if existing_count < 0:
        raise ValueError(f"existing_count must be >= 0, got {existing_count!r}")
    n_greek = len(GREEK_LETTERS)
    n_roman = len(ROMAN_NUMERAL_VALUES)
    alpha_name = f"{GREEK_LETTERS[0]} {base_name}"

    if existing_count == 0:
        return base_name, None
    if existing_count == 1:
        # First collision: the sole existing row is still bare -- it
        # becomes Alpha, the row being inserted now becomes Beta.
        return f"{GREEK_LETTERS[1]} {base_name}", (base_name, alpha_name)
    if existing_count < n_greek:
        # existing_count rows already occupy Greek ranks 1..existing_count
        # (GREEK_LETTERS[0..existing_count-1]); the next rank is
        # GREEK_LETTERS[existing_count].
        return f"{GREEK_LETTERS[existing_count]} {base_name}", None
    if existing_count == n_greek:
        # Greek tier exhausted (24 rows, every letter used). Converts the
        # existing "Alpha <base>" row -- the original row, renamed once
        # already -- into the roman-numeral tier's first slot; the row
        # being inserted now takes the second.
        value_i, value_ii = ROMAN_NUMERAL_VALUES[0], ROMAN_NUMERAL_VALUES[1]
        renamed = f"{alpha_name} {ROMAN_NUMERALS_BY_VALUE[value_i]}"
        new_name = f"{alpha_name} {ROMAN_NUMERALS_BY_VALUE[value_ii]}"
        return new_name, (alpha_name, renamed)

    # existing_count == n_greek + 1 must mint roman value index 2 (the
    # n_greek == existing_count branch above already minted indices 0
    # and 1 in that single step), so the offset is +1, not +0.
    roman_index = existing_count - n_greek + 1
    if roman_index < n_roman:
        value = ROMAN_NUMERAL_VALUES[roman_index]
        return f"{alpha_name} {ROMAN_NUMERALS_BY_VALUE[value]}", None

    return None, None


def resolve_diminutive(diminutive_index):
    """
    Cross-level collision resolution between a sector and a system that
    share a base name -- always decorates the *system* side (see this
    module's own docstring for why).

    Args:
        diminutive_index (int or None): The index into
            `names.DIMINUTIVE_PREFIXES` already used for this base name's
            prior cross-level collision(s), or `None` if this is the
            first one.

    Returns:
        tuple: `(prefix, next_index)` -- `prefix` (str or None) is the
            word to prepend to the system's name (`None` once every entry
            in `DIMINUTIVE_PREFIXES` is used -- the caller regenerates the
            newest conflicting name entirely instead, per this feature's
            own "leave the others alone" rule). `next_index` is the value
            to persist for next time (`None` alongside a `None` prefix).
    """
    if diminutive_index is not None and diminutive_index < 0:
        raise ValueError(f"diminutive_index must be None or >= 0, got {diminutive_index!r}")
    next_index = 0 if diminutive_index is None else diminutive_index + 1
    if next_index >= len(DIMINUTIVE_PREFIXES):
        return None, None
    return DIMINUTIVE_PREFIXES[next_index], next_index


def has_diminutive(name):
    """Whether `name` starts with one of `names.DIMINUTIVE_PREFIXES`: a
    system name a sector name can never take, since sectors only ever get
    Greek prefixes."""
    words = name.split(" ")
    return bool(words) and words[0] in _DIMINUTIVE_PREFIX_TOKENS


def strip_decoration(name):
    """
    Recovers the underlying base name from one that may carry any
    combination of this module's own decorations -- the inverse of
    `resolve_greek_roman_collision`/`resolve_diminutive`. Used by `src/dedupeNames.py` to group already-
    stored names by what they'd collide on.

    Strips a trailing roman numeral and leading Greek-letter and
    diminutive prefixes, repeating until none is left, so a stacked name
    ("Little Beta <base>", "<base> IV" with a prefix) comes apart fully.
    Sectors/systems only ever get a Greek prefix and, for systems only,
    also a diminutive prefix or a trailing roman numeral. Safe because the
    three word lists don't overlap and a generated base word is never one
    of them (`utils.is_name_valid` rejects every decoration word).

    Args:
        name (str): A stored `sectors`/`star_systems` name, decorated or
            not.

    Returns:
        str: `name` with every decoration this module could have added
            stripped away.
    """
    # Repeated until nothing changes: `_db.py` stacks decorations -- a
    # diminutive prefix *outside* a Greek one ("Little Beta <base>",
    # `reserve_system_name`), a second diminutive on top of the first
    # ("Petit Little <base>") or a second companion suffix ("<base> Kin
    # Ami"), each time `_rename_existing_*` decorates the same row again.
    tokens = name.split(" ")
    while True:
        before = len(tokens)
        if tokens and tokens[-1] in _ROMAN_SUFFIX_TOKENS:
            tokens = tokens[:-1]
        while tokens and (tokens[0] in _GREEK_PREFIX_TOKENS or tokens[0] in _DIMINUTIVE_PREFIX_TOKENS):
            tokens = tokens[1:]
        if len(tokens) == before:
            return " ".join(tokens)
