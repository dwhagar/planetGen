# planetgen/util/plainwords.py

"""
The checkable rules of the plain-words style guide (DOC.83, `docs/design/interface/text-style.md`).

Every sentence the site writes about a planet, a moon, a belt, a star or a system is meant to be short, in plain words, one idea
at a time, with the same class name everywhere. A person reads the guide; this module is the part of it a machine can check,
so a test can read the generated text (`tests/test_plain_words.py`) and the rewrite items (UX.252, UX.253) have something to run
against. Standard library and `physics.life_classes` only; nothing here changes a value, it only reads text.

`check_sentence(text)` returns the rules a sentence breaks and `check_label(label)` the rules a table label breaks, each as
`(rule, detail)`. `sentences(text)` splits a paragraph the way the checks count.
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

_ABBREVIATIONS = ("e.g.", "i.e.", "etc.", "vs.", "approx.", "Dr.", "St.", "No.")
"""tuple: Where a full stop does not end a sentence."""

_CLASS_NAME = re.compile(r"(?<!uminosity )\bClass (\d+)([a-c]?)\b")
"""re.Pattern: A life class named in text ("Class 3a"); a star's luminosity class ("Luminosity class 0") is not one."""

SMALL_WORDS = {"a", "an", "and", "as", "at", "by", "for", "from", "in", "into", "of", "on", "or", "per", "the", "to", "with"}
"""set: Words a label keeps in lower case wherever they stand."""

PROPER_NAMES = {"Urey", "Earth", "Jupiter", "Sun", "Solar", "Saturn", "Neptune", "Uranus", "Mars", "Venus", "Mercury", "Moon", "Roche", "Hill",
                "Kelvin", "Celsius", "Fahrenheit"}
"""set: Capitalised words a label may carry in the middle of it."""


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
    * `class-name`: a life class named by its family ("Class 3") where the class has letters ("Class 3a").
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
    return found


def check_label(label):
    """
    The rules a table label breaks, as `[(rule, detail)]`: `label-case` when a word after the first starts with a capital letter and is
    not a proper name or an abbreviation (labels are in sentence case: "Distance from the star", not "Distance From The Star").
    """
    found = []
    for position, word in enumerate(label.split()):
        bare = word.strip("()'’,:")
        bare = bare[:-2] if bare.endswith(("'s", "’s")) else bare
        if position == 0 or not bare or not bare[0].isupper():
            continue
        if bare in PROPER_NAMES or (len(bare) > 1 and bare.upper() == bare) or any(ch.isdigit() for ch in bare):
            continue
        found.append(("label-case", f"{label!r}: write '{bare.lower()}' in lower case"))
    return found
