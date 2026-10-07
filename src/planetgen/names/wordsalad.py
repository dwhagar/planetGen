# planetgen/names/wordsalad.py

"""
Word Salad
==========

The syllable and phoneme-salad name generator, its validity checks, and
sector names. GEN.71 replaces it with the phoneme codec.
"""

import random

from planetgen.names.wordlists import (
    BAD_CONSONANTS, COMPANION_SUFFIXES, DICTIONARY_WORDS, DIMINUTIVE_PREFIXES, GREEK_LETTERS, NSFW_WORDS, ROMAN_NUMERALS_BY_VALUE, SECTOR_NAMES, SECTOR_PREFIXES, SECTOR_SUFFIXES, UNIVERSAL_PHONEMES, VOWELS, WORD_SIZE_MEAN,
)


def split_into_syllables(name):
    """
    Splits a word into a list of syllables.

    This is a basic syllable splitting function that helps in the process of
    generating new names by rearranging syllables from existing names.

    Args:
        name (str): The word to be split into syllables.

    Returns:
        list: A list of strings, where each string is a syllable.
    """
    syllables = []
    current_syllable = ""
    for i, char in enumerate(name):
        current_syllable += char
        if char in VOWELS and i < len(name) - 1 and name[i+1] not in VOWELS:
            syllables.append(current_syllable)
            current_syllable = ""
    if current_syllable:
        syllables.append(current_syllable)
    return syllables


_NSFW_PLAIN_WORDS = tuple(w for w in NSFW_WORDS if w and "'" not in w and " " not in w)
_NSFW_SPACED_WORDS = tuple(w for w in NSFW_WORDS if w and ("'" in w or " " in w))


def is_name_valid(name):
    """
    Checks if a generated name is valid based on multiple criteria.

    This function ensures that the generated name is not a common English word,
    does not contain any offensive terms, and follows basic phonetic rules.
    The validation checks are as follows:
    - The name should not exist in the NLTK dictionary of words.
    - The name should not contain any substring from the NSFW (Not Safe For Work) word list,
      also checked with its apostrophes and spaces removed.
    - No word of the name may be a `nameUniqueness` decoration word (`_DECORATION_WORDS`).
    - The name should not contain more than two consecutive vowels.
    - The name should not contain more than two consecutive consonants.
    - The name should not contain any of the defined bad consonant clusters.

    These checks help in generating names that are unique, appropriate, and sound plausible.

    Args:
        name (str): The name to validate. The function expects a lowercase string.

    Returns:
        bool: Returns `True` if the name is valid according to all the rules,
              otherwise returns `False`.
    """
    name_lower = name.lower()
    if name_lower in DICTIONARY_WORDS:
        return False
    if any(token in _DECORATION_WORDS for token in name_lower.split()):
        # A word that is also a nameUniqueness decoration ("Liten",
        # "Ohana", ...) would be stripped off an undecorated name by
        # `nameUniqueness.strip_decoration`, grouping it with another base.
        return False
    # Checked with apostrophes/spaces removed too: an apostrophe spliced in
    # by a phoneme chunk ("pak'i", "k'ike") must not hide an offensive word.
    # A word with no apostrophe or space that appears in any variant also
    # appears in the fully squeezed one, so only the few words that carry
    # one of those need every variant (keeps this hot check one pass).
    squeezed = name_lower.replace("'", "")
    fully_squeezed = squeezed.replace(" ", "")
    if any(word in fully_squeezed for word in _NSFW_PLAIN_WORDS):
        return False
    if _NSFW_SPACED_WORDS:
        variants = (name_lower, squeezed, name_lower.replace(" ", ""))
        if any(word in variant for word in _NSFW_SPACED_WORDS for variant in variants):
            return False
    
    vowel_count = 0
    consonant_count = 0
    for char in name_lower:
        if char in VOWELS:
            vowel_count += 1
            consonant_count = 0
        elif char.isalpha():
            consonant_count += 1
            vowel_count = 0
        else:
            vowel_count = 0
            consonant_count = 0
        if vowel_count > 2 or consonant_count > 2:
            return False

    for cluster in BAD_CONSONANTS:
        if cluster in name_lower:
            return False

    return True


_DECORATION_WORDS = frozenset(
    word.lower() for word in (*GREEK_LETTERS, *DIMINUTIVE_PREFIXES, *COMPANION_SUFFIXES,
                              *ROMAN_NUMERALS_BY_VALUE.values())
)
"""Every word `nameUniqueness` decorates a name with (lowercase) -- see
`is_name_valid`/`generate_phoneme_salad_name`, which keep generated words
from ever being one of them."""


def split_long_word(name):
    """
    Splits a long word into two, capitalizing the second word.

    This function improves the readability of long generated names by splitting
    them into two parts, creating a compound name effect.

    A name can contain an embedded apostrophe (from a base name like
    "Hi'iaka" or a spliced-in `UNIVERSAL_PHONEMES` chunk like "ch'"/"b'a").
    If the naive midpoint split landed right next to one of those, one half
    would end up starting or ending with a bare "'" (e.g. "Amech' Snesis") --
    so the split point is nudged past any apostrophe it would otherwise cut
    beside. If nudging would consume the whole rest of the name, the split
    is skipped and the name is returned unsplit rather than produce an
    empty half.

    Args:
        name (str): The long word to be split.

    Returns:
        str: The split and capitalized name, or the original name if not
             long enough (or if avoiding an apostrophe boundary leaves no
             valid split point).
    """
    if len(name) > WORD_SIZE_MEAN:
        if " " in name:
            # Already more than one word (a base name like "El Nath");
            # splitting again could put a second space beside the first.
            return name
        split_point = len(name) // 2
        while split_point < len(name) and (name[split_point - 1] == "'" or name[split_point] == "'"):
            split_point += 1
        if not (0 < split_point < len(name)):
            return name
        return name[:split_point] + " " + name[split_point:].capitalize()
    return name


UNIVERSAL_PHONEME_CHANCE = 0.4
"""
Odds that `generate_phoneme_salad_name` splices an extra chunk from
`names.UNIVERSAL_PHONEMES` into a generated name's syllable pool. Applies
uniformly to stars, planets, moons, and sectors -- every one of those
funnels through this same function -- rather than needing each type's own
prefix/suffix lists to be extended individually.
"""

MAX_NAME_GENERATION_ATTEMPTS = 10_000
"""
Hard cap on `generate_phoneme_salad_name`'s own retry loop -- see that
function's `Raises` doc. A real generation run finds a valid name within a
handful of attempts; this is only a backstop against `is_name_valid`
rejecting every candidate forever (previously an unconditional `while
True` with no way out at all).
"""


def generate_phoneme_salad_name(name_list, prefix_list, suffix_list, allow_split=True, syllable_fraction=1.0, max_length=None):
    """
    Generates a unique, phonetically pleasing name from a list of base names.

    This function creates new names by taking a base name, shuffling its
    syllables, and adding a prefix and suffix. It includes logic to ensure
    the resulting name is phonetically plausible and passes validation checks.
    With `UNIVERSAL_PHONEME_CHANCE` odds, it also splices in one chunk from
    `names.UNIVERSAL_PHONEMES` -- a cross-linguistic phoneme pool shared by
    every name type -- widening the cultural range names are drawn from
    beyond whatever's in `name_list` itself.

    Args:
        name_list (list): A list of base names to choose from.
        prefix_list (list): A list of possible prefixes.
        suffix_list (list): A list of possible suffixes.
        allow_split (bool): Whether a long result may be split into two
                            space-separated words via `split_long_word`
                            (e.g. `"Xyleth Anore"`). Default `True` for
                            stars/planets/moons, where that reads fine as
                            one name. Callers that combine multiple calls
                            into one already-multi-word name (e.g.
                            `sectorGen.generate_sector_name`, which joins
                            two of these into a two-word sector name) must
                            pass `False` here, or a single call splitting
                            internally would silently make the combined
                            result 3-4 words instead of the intended 2.
        syllable_fraction (float): Fraction (0-1] of the shuffled base
                            name's syllables to actually keep before the
                            prefix/suffix are attached, trimming from the
                            end. Default `1.0` keeps every syllable
                            (unchanged behavior for stars/planets/moons).
                            `sectorGen.generate_sector_name` passes `0.5`
                            here so sector names -- built by joining two
                            of these calls into one two-word name -- come
                            out roughly half as long per word; base names
                            long enough to still shrink always keep at
                            least 1 syllable.
        max_length (int): Optional hard cap, in characters, on the fully
                            assembled name (base syllables + prefix +
                            suffix, before capitalization). `None` (the
                            default) leaves names uncapped. Chopping
                            happens after the suffix is attached, so it's
                            the backstop against `syllable_fraction`
                            alone not being enough -- prefixes, suffixes,
                            and the occasional spliced-in
                            `UNIVERSAL_PHONEMES` chunk are fixed-ish
                            overhead that doesn't shrink with
                            `syllable_fraction`, so a long base name can
                            still produce a longer-than-intended result
                            without this. `sectorGen.generate_sector_name`
                            passes `7` here alongside `syllable_fraction=0.5`
                            to reliably keep each half of a sector name
                            short.

    Returns:
        str: A newly generated, unique name.

    Raises:
        RuntimeError: If `MAX_NAME_GENERATION_ATTEMPTS` consecutive draws
            all fail `is_name_valid` -- confirmed (via a hard-timeout test,
            `test_bughunt_name_exhaustion.py`) to otherwise loop forever
            with no way out. In real generation this is never reached (a
            valid name is essentially always found within a handful of
            attempts); this guard only matters if `is_name_valid`'s own
            filters (or a `name_list`/`prefix_list`/`suffix_list` this
            narrow) were ever misconfigured to reject everything, in which
            case a whole generation run should fail loudly and immediately
            rather than hang with no explanation.
    """
    for _attempt in range(MAX_NAME_GENERATION_ATTEMPTS):
        name = random.choice(name_list)

        syllables = split_into_syllables(name)
        if len(syllables) > 1:
            random.shuffle(syllables)

        if syllable_fraction < 1.0 and len(syllables) > 1:
            keep = max(1, round(len(syllables) * syllable_fraction))
            syllables = syllables[:keep]

        if random.random() < UNIVERSAL_PHONEME_CHANCE:
            syllables.insert(random.randint(0, len(syllables)), random.choice(UNIVERSAL_PHONEMES))

        name = "".join(syllables)

        prefix = random.choice(prefix_list)
        if prefix[-1] in VOWELS and name[0].lower() in VOWELS:
            name = prefix + name[1:]
        elif prefix[-1] not in VOWELS and name[0].lower() not in VOWELS:
            if (prefix[-1] + name[0].lower()) in BAD_CONSONANTS:
                name = prefix + "'" + name
            else:
                name = prefix + name
        else:
            name = prefix + name

        suffix = random.choice(suffix_list)
        if name[-1] in VOWELS and suffix[0].lower() in VOWELS:
            name = name + suffix[1:]
        elif name[-1] not in VOWELS and suffix[0].lower() not in VOWELS:
            if (name[-1] + suffix[0].lower()) in BAD_CONSONANTS:
                name = name + "'" + suffix
            else:
                name = name + suffix
        else:
            name = name + suffix

        if max_length is not None and len(name) > max_length:
            # Chopping at a raw character index can land right after an
            # embedded apostrophe (from a base name like "Hi'iaka" or a
            # spliced-in UNIVERSAL_PHONEMES chunk like "ch'"), leaving the
            # truncated name ending in a bare "'" -- so back off past it.
            cutoff = max_length
            while cutoff > 0 and name[cutoff - 1] == "'":
                cutoff -= 1
            name = name[:cutoff] if cutoff > 0 else name[:max_length]

        name = name.lower()

        if is_name_valid(name):
            if allow_split:
                name = split_long_word(name)
                if any(part.lower() in _DECORATION_WORDS for part in name.split()):
                    # A split half that is itself a decoration word
                    # ("Xxxxx Ohana") -- see `_DECORATION_WORDS`.
                    continue
            name = name[0].upper() + name[1:]
            if "'" in name:
                # Two apostrophes can land adjacent here (e.g. a base name like
                # "Hi'iaka" splits into a syllable starting with "'", and a
                # spliced-in UNIVERSAL_PHONEMES chunk like "ch'" ends with "'";
                # shuffling can place them next to each other), producing an
                # empty string between them once split -- guard against
                # indexing that empty part rather than assuming every part is
                # non-empty.
                parts = name.split("'")
                name = "'".join([part[0].upper() + part[1:] if part else part for part in parts])
            return name

    raise RuntimeError(
        f"generate_phoneme_salad_name: no valid name found in {MAX_NAME_GENERATION_ATTEMPTS} attempts "
        f"(name_list={name_list!r}) -- is_name_valid is rejecting every candidate."
    )


def generate_sector_name():
    """
    Generates a random two-word sector name, each word independently drawn
    from the same phoneme-salad name generator used for star/planet/moon
    names -- using the sector-flavored `SECTOR_NAMES`/`SECTOR_PREFIXES`/
    `SECTOR_SUFFIXES` base lists instead, so generated sectors draw on real
    astronomical regions (galactic arms, superclusters, nebulae) and
    science-fiction sector names rather than reusing star names verbatim.
    No literal "Sector" suffix. `planetgen`'s `sector`/`galaxy` subcommands
    override this entirely via `--name`/`-n`, which hard-sets the whole
    name instead; `planetgen.db.store`'s name-uniqueness machinery
    (`nameUniqueness.py`) also calls this directly, to draw an entirely
    fresh sector name on the rare occasion a collision exhausts every
    decoration this project has for one -- both reasons this lives here,
    in the `planetgen` package, rather than in `planetgen` itself, which
    `store.py` can't import (it would be a backwards/circular dependency --
    `planetgen` already imports `planetgen.db.store`).

    Returns:
        str: A newly generated sector name, e.g. "Voranthis Kelmoor" --
        always exactly two words.
    """
    # allow_split=False: generate_phoneme_salad_name can itself split a
    # long result into two words (e.g. "Xyleth Anore"). Since this
    # function already joins two independent calls into one name, leaving
    # splitting on could silently produce 3-4 words instead of 2.
    # syllable_fraction=0.5 trims each word's base syllables by about
    # half before the prefix/suffix are attached -- many SECTOR_NAMES
    # entries (e.g. "Sagittarius", "Metropolis") are long real place
    # names, and two of them joined together made for unwieldy sector
    # names. max_length=7 backstops that: prefixes, suffixes, and the
    # occasional spliced-in universal phoneme are fixed-ish overhead that
    # doesn't shrink with syllable_fraction, so a long base name could
    # still slip through longer than intended without a hard cap too.
    first_word = generate_phoneme_salad_name(SECTOR_NAMES, SECTOR_PREFIXES, SECTOR_SUFFIXES, allow_split=False, syllable_fraction=0.5, max_length=7)
    second_word = generate_phoneme_salad_name(SECTOR_NAMES, SECTOR_PREFIXES, SECTOR_SUFFIXES, allow_split=False, syllable_fraction=0.5, max_length=7)
    return f"{first_word} {second_word}"
