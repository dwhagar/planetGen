# planetgen/util/plainwords.py

"""
The checkable rules of the plain-words style guide (DOC.83, `docs/design/interface/text-style.md`).

Every sentence the site writes about a planet, a moon, a belt, a star or a system is meant to be short, in plain words, one idea
at a time, with the same class name everywhere. A person reads the guide; this module is the part of it a machine can check,
so a test can read the generated text (`tests/test_plain_words.py`) and the rewrite items (UX.252, UX.253) have something to run
against. Standard library and `physics.life_classes` only; nothing here changes a value, it only reads text.

`check_sentence(text)` returns the rules a sentence breaks and `check_label(label)` the rules a table label breaks, each as
`(rule, detail)`. `sentences(text)` splits a paragraph the way the checks count. `out_of_universe(text)` reads rule 16 (DOC.87): the
words that belong only on the References pages (`/references`), which `exempt(path)` names along with the other pages that keep them.
"""

import re

from planetgen.physics import life_classes

MAX_WORDS = 30
"""int: The most words one sentence may have. The guide aims at 20; a sentence over 30 is two ideas."""

MAX_SEMICOLONS = 1
"""int: The most semicolons one sentence may have: each one after the first is another idea that wants its own sentence."""

BANNED = {
    "approximately": "about",
    "utilize": "use", "utilise": "use",
    "commence": "start",
    "terminate": "end",
    "prior to": "before",
    "subsequent to": "after",
    "in order to": "to",
    "due to the fact that": "because",
    "numerous": "many",
    "sufficient": "enough",
    "a number of": "some",
    "in the vicinity of": "near",
    "the majority of": "most",
    "is able to": "can",
    "has the ability to": "can",
    "in excess of": "more than",
    "e.g.": "for example",
    "i.e.": "that is",
    "etc.": "and so on",
    "stand-in": "(say what the value is, or leave it out)",
    "placeholder": "(say what the value is, or leave it out)",
    "tweak": "(a reader does not need to know how a number was set)",
    "epoch": "(a reader does not need the generator's version)",
}
"""dict: Words and phrases that do not belong in the site's text, each with what to write instead."""

_BANNED = [(re.compile(r"(?<![\w-])" + re.escape(phrase) + r"(?![\w-])", re.IGNORECASE), phrase, plain) for phrase, plain in BANNED.items()]

OUT_OF_UNIVERSE = (
    "speculative", "stand-in", "standin", "tweak", "tweaks", "tweaked", "fitted", "tuned", "tuning", "educated guess",
    "tier", "tiers", "generator", "generated", "seed", "seeds", "simulated", "simulation", "model", "models", "modeled", "modelled",
    "modeling", "modelling", "version", "versions", "placeholder", "epoch",
)
"""tuple: The words that tell a reader the site is a model (rule 16, DOC.87). The site speaks in universe: a working navigation system that
simply knows things. Where a number came from, what is a stand-in, what was fitted to a result and how sure a part is are said only on
the References pages (`/references`); every other page says the thing in universe and, where a reader may want the why, links to the
References entry (`web/helpers.reference_link`). "Level" is not on the list: "sea level" and "ground level" are in universe, and the
provenance sense ("tier and level") always comes with "tier" or "tweak", which are."""

_OUT_OF_UNIVERSE = re.compile(r"(?<![\w-])(" + "|".join(re.escape(word) for word in sorted(OUT_OF_UNIVERSE, key=len, reverse=True)) + r")(?![\w-])",
                              re.IGNORECASE)

EXEMPT_PATHS = ("/references", "/admin")
"""tuple: The page paths where the out-of-universe words are allowed: the References pages, which carry the Out of universe banner, and the
admin pages. (The Debug view is a query option of a page, not a path, and the API reference and `docs/design` are not site pages.)"""

_ABBREVIATIONS = ("e.g.", "i.e.", "etc.", "vs.", "approx.", "Dr.", "St.", "No.")
"""tuple: Where a full stop does not end a sentence."""

_CLASS_NAME = re.compile(r"(?<!uminosity )\bClass (\d+)([a-c]?)\b(?! family)")
"""re.Pattern: A life class named in text ("Class 3a"); a star's luminosity class ("Luminosity class 0") is not one, and "the Class 3 family" names the whole family."""

SMALL_WORDS = {"a", "an", "and", "as", "at", "by", "for", "from", "in", "into", "of", "on", "or", "per", "the", "to", "with"}
"""set: Words a label keeps in lower case wherever they stand."""

PROPER_NAMES = {"Urey", "Earth", "Jupiter", "Sun", "Solar", "Saturn", "Neptune", "Uranus", "Mars", "Venus", "Mercury", "Moon", "Roche", "Hill",
                "Kelvin", "Celsius", "Fahrenheit"}
"""set: Capitalised words a label may carry in the middle of it."""


BRITISH = {
    "colour": "color", "grey": "gray", "centre": "center", "metre": "meter", "litre": "liter", "neighbour": "neighbor",
    "behaviour": "behavior", "favour": "favor", "honour": "honor", "labour": "labor", "vapour": "vapor", "harbour": "harbor",
    "armour": "armor", "flavour": "flavor", "odour": "odor", "fibre": "fiber", "theatre": "theater", "calibre": "caliber",
    "defence": "defense", "offence": "offense", "aluminium": "aluminum", "sulphur": "sulfur", "sulphate": "sulfate",
    "programme": "program", "whilst": "while", "amongst": "among", "judgement": "judgment", "ageing": "aging", "mould": "mold",
    "tyre": "tire", "modelling": "modeling", "travelling": "traveling", "labelled": "labeled", "signalling": "signaling",
    "levelling": "leveling",
}
"""dict: British spellings and the American ones (the site is in American English, Boss 10:33 AM EDT 2026-10-31). A word also matches
with a plural, `-ed`, `-ing`, `-er`, `-ful`, `-less` or `-ly` ending, and `metre` and `litre` with an SI prefix (`kilometre`)."""

ISE_STEMS = ("organ", "real", "recogn", "character", "crystall", "ion", "oxid", "polar", "minim", "maxim", "optim", "util", "summar",
             "categor", "stabil", "steril", "neutral", "special", "emphas", "apolog", "critic", "general", "standard", "normal", "visual",
             "synchron", "hypothes", "ideal", "equal", "final", "local", "capital", "fertil", "mobil", "sensit", "familiar", "mineral",
             "harmon", "central", "formal", "legal", "national", "personal", "social", "solid", "rational", "fossil", "atom", "carbon")
"""tuple: The stems of the verbs that British English spells `-ise` and American English `-ize` (`organise`, `organisation`)."""

BRITISH_NAMES = ("Centre National", "Centre for", "Royal Observatory")
"""tuple: Names that keep their own spelling (an organization's name, a title); they are left out before the spelling is read."""

_SI_PREFIX = r"(?:kilo|milli|centi|micro|nano|deci|mega|giga)?"
_ENDING = r"(?:s|d|ed|ing|er|ers|ful|less|ly)?"
_SPELLINGS = [(re.compile(r"\b" + (_SI_PREFIX if word in ("metre", "litre") else "") + re.escape(word) + _ENDING + r"\b", re.IGNORECASE), word, plain)
              for word, plain in BRITISH.items()]
_SPELLINGS += [(re.compile(r"\b" + stem + r"is(?:e|es|ed|ing|ation|ations)\b", re.IGNORECASE), stem + "is", stem + "iz") for stem in ISE_STEMS]
_SPELLINGS.append((re.compile(r"\banalys(?:e|ed|ing|er|ers)\b", re.IGNORECASE), "analys", "analyz"))
_SPELLINGS.append((re.compile(r"\bcatalys(?:e|ed|ing)\b", re.IGNORECASE), "catalys", "catalyz"))
_QUOTED = re.compile(r'"[^"]*"|\u201c[^\u201d]*\u201d')


def british_spellings(text):
    """
    The British spellings in `text` as `[("spelling", detail)]`, each saying what to write (`'colour': write 'color'`). Text in quotation
    marks and the names of `BRITISH_NAMES` are not read.
    """
    plain = _QUOTED.sub(" ", text)
    for name in BRITISH_NAMES:
        plain = plain.replace(name, " ")
    found = []
    for pattern, british, american in _SPELLINGS:
        for match in pattern.finditer(plain):
            word = match.group(0)
            lower = word.lower()
            fix = lower.replace(british, american, 1)
            fix = fix[:-3] + "ered" if fix.endswith("erd") and british.endswith("re") else fix
            found.append(("spelling", f"'{word}': write '{fix}'"))
    return found


def sentences(text):
    """
    `text` split into its sentences. A full stop, `!` or `?` ends one when a capital letter or a digit follows after a space, and not
    after an abbreviation (`e.g.`) or inside a number (`1.06`).
    """
    protected = text
    for abbreviation in _ABBREVIATIONS:
        protected = protected.replace(abbreviation, abbreviation.replace(".", "\u0000"))
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z0-9\"'(])", protected.strip())
    return [part.replace("\u0000", ".").strip() for part in parts if part.strip()]


def words(sentence):
    """The words of a sentence as the length rule counts them: runs of characters between spaces."""
    return sentence.split()


def check_sentence(sentence):
    """
    The rules one sentence breaks, as `[(rule, detail)]` (empty when it keeps them all). The rules:

    * `length`: more than `MAX_WORDS` words;
    * `ideas`: more than `MAX_SEMICOLONS` semicolons;
    * `banned`: a word or phrase of `BANNED`, with what to write instead;
    * `class-name`: a life class named by its family ("Class 3") where the class has letters ("Class 3a");
    * `spelling`: a British spelling (`colour`, `centre`, `metre`, `analyse`, `organise`); the site is in American English (Boss,
      10:33 AM EDT 2026-10-31). Quoted passages and the names in `BRITISH_NAMES` are left alone.
    """
    found = []
    count = len(words(sentence))
    if count > MAX_WORDS:
        found.append(("length", f"{count} words (the most is {MAX_WORDS})"))
    if sentence.count(";") > MAX_SEMICOLONS:
        found.append(("ideas", f"{sentence.count(';')} semicolons (the most is {MAX_SEMICOLONS})"))
    for pattern, phrase, plain in _BANNED:
        if pattern.search(sentence):
            found.append(("banned", f"{phrase!r}: write {plain}"))
    for digits, letter in _CLASS_NAME.findall(sentence):
        if f"{digits}{letter}" not in life_classes.BY_ID:
            found.append(("class-name", f"'Class {digits}{letter}' is not a class: the classes are {', '.join(life_classes.BY_ID)}"))
    found += british_spellings(sentence)
    return found


def check_label(label):
    """
    The rules a table label breaks, as `[(rule, detail)]`: `spelling` (as for a sentence) and `label-case` when a word after the first starts with a capital letter and is
    not a proper name or an abbreviation (labels are in sentence case: "Distance from the star", not "Distance From The Star").
    """
    found = british_spellings(label)
    for position, word in enumerate(label.split()):
        bare = word.strip("()'’,:")
        bare = bare[:-2] if bare.endswith(("'s", "’s")) else bare
        if position == 0 or not bare or not bare[0].isupper():
            continue
        if bare in PROPER_NAMES or (len(bare) > 1 and bare.upper() == bare) or any(ch.isdigit() for ch in bare):
            continue
        found.append(("label-case", f"{label!r}: write '{bare.lower()}' in lower case"))
    return found


def out_of_universe(text):
    """
    The out-of-universe words in `text` (rule 16, DOC.87) as `[("out-of-universe", detail)]`, one for each different word, in the order
    they first stand. Text in quotation marks is read too: a quoted source is still a hint that the site is a model.
    """
    seen = []
    for match in _OUT_OF_UNIVERSE.finditer(text):
        word = match.group(1).lower()
        if word not in seen:
            seen.append(word)
    return [("out-of-universe", f"'{word}': say it in universe, or link to the References entry") for word in seen]


def exempt(path):
    """True when the page at `path` may use the out-of-universe words: a References page or an admin page (`EXEMPT_PATHS`)."""
    return any(path == prefix or path.startswith(prefix + "/") for prefix in EXEMPT_PATHS)
