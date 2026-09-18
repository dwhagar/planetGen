# stellarObjects/nameUniqueness.py

"""
Name-uniqueness decoration -- pure functions, no database access.

Guarantees, across a whole database, that no two sectors share a name, no
two star systems share a name, no sector/system pair shares a name, and
no planet/moon shares a name with another planet, a moon, a system, or a
sector. Enforced as one consistent hierarchy: **sector > system >
planet/moon**.

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
- A planet or moon colliding with anything -- another planet, a moon, a
  system, or a sector -- is resolved by `resolve_companion`: the
  planet's/moon's own name (always the lower level, regardless of which
  side was generated first) gets a "friend"/"family"/"companion" word
  (`names.COMPANION_SUFFIXES`) appended, again advancing to the next word
  in the list on a repeat.

`strip_decoration` is the inverse of all three -- recovering the
underlying base name from an already-decorated one, e.g. for the
one-off backfill script (`src/dedupeNames.py`) to group existing rows.

None of this module talks to the database or knows about `sectors`/
`star_systems`/`planets`/`moons` -- see `stellarObjects/_db.py`'s
`_reserve_sector_name`/`_reserve_system_name`/`_reserve_body_name` (and
their `_confirm_*` counterparts) for where these functions actually get
called, against the `sector_name_registry`/`system_name_registry`/
`body_name_registry` tables that track each base name's own progress
through these schemes.
"""

from .names import (
    COMPANION_SUFFIXES, DIMINUTIVE_PREFIXES, GREEK_LETTERS,
    ROMAN_NUMERAL_VALUES, ROMAN_NUMERALS_BY_VALUE,
)

GREEK_ROMAN_CAPACITY = len(GREEK_LETTERS) + len(ROMAN_NUMERAL_VALUES)
"""int: How many total rows may ever share one base name under
`resolve_greek_roman_collision` (24 Greek letters + 12 roman-numeral
slots) before the caller must fall back to an entirely fresh base name."""

_ROMAN_SUFFIX_TOKENS = set(ROMAN_NUMERALS_BY_VALUE.values())
_GREEK_PREFIX_TOKENS = set(GREEK_LETTERS)
_DIMINUTIVE_PREFIX_TOKENS = set(DIMINUTIVE_PREFIXES)
_COMPANION_SUFFIX_TOKENS = set(COMPANION_SUFFIXES)


def resolve_greek_roman_collision(base_name, existing_count):
    """
    Same-level collision resolution for two sectors, or two systems,
    that generated the same `base_name` -- see this module's own
    docstring for the full scheme.

    Args:
        base_name (str): The undecorated name both rows share.
        existing_count (int): How many rows already exist with this base
            name, *before* the row currently being inserted/renamed.
            `0` means no collision at all (this is the first-ever row).

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
    next_index = 0 if diminutive_index is None else diminutive_index + 1
    if next_index >= len(DIMINUTIVE_PREFIXES):
        return None, None
    return DIMINUTIVE_PREFIXES[next_index], next_index


def resolve_companion(suffix_index):
    """
    Collision resolution for a planet or moon against anything else that
    must stay unique against it (another planet, a moon, a system, or a
    sector) -- always decorates the planet's/moon's own name, the lowest
    level in the hierarchy. Same shape as `resolve_diminutive`, over
    `names.COMPANION_SUFFIXES` instead, appended as a *suffix* rather
    than a prefix.

    Args:
        suffix_index (int or None): The index into
            `names.COMPANION_SUFFIXES` already used for this base name's
            prior collision(s), or `None` if this is the first one.

    Returns:
        tuple: `(suffix, next_index)` -- see `resolve_diminutive`'s own
            `Returns` for the shape; `suffix` is appended after the name
            rather than prepended before it.
    """
    next_index = 0 if suffix_index is None else suffix_index + 1
    if next_index >= len(COMPANION_SUFFIXES):
        return None, None
    return COMPANION_SUFFIXES[next_index], next_index


def strip_decoration(name):
    """
    Recovers the underlying base name from one that may carry any
    combination of this module's own decorations -- the inverse of
    `resolve_greek_roman_collision`/`resolve_diminutive`/
    `resolve_companion`. Used by `src/dedupeNames.py` to group already-
    stored names by what they'd collide on.

    Checks each decoration independently and unconditionally, in a fixed
    order (trailing roman numeral, then trailing companion suffix, then
    leading Greek letter, then leading diminutive prefix) -- safe for a
    name from any of the four tables, since no real entity ever carries
    more than one of these four vocabularies at once (planets/moons only
    ever get a companion suffix; sectors/systems only ever get a Greek
    prefix and, for systems only, also a diminutive prefix or a trailing
    roman numeral) and the four word lists don't overlap.

    Args:
        name (str): A stored `sectors`/`star_systems`/`planets`/`moons`
            name, decorated or not.

    Returns:
        str: `name` with every decoration this module could have added
            stripped away.
    """
    tokens = name.split(" ")
    if tokens and tokens[-1] in _ROMAN_SUFFIX_TOKENS:
        tokens = tokens[:-1]
    if tokens and tokens[-1] in _COMPANION_SUFFIX_TOKENS:
        tokens = tokens[:-1]
    if tokens and tokens[0] in _GREEK_PREFIX_TOKENS:
        tokens = tokens[1:]
    if tokens and tokens[0] in _DIMINUTIVE_PREFIX_TOKENS:
        tokens = tokens[1:]
    return " ".join(tokens)
