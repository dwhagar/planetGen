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
    "tier", "tiers", "generator", "generated", "random seed", "random seeds", "seeded", "seed value", "seed number", "simulated", "simulation", "model", "models", "modeled", "modelled",
    "modeling", "modelling", "version", "versions", "placeholder", "epoch",
)
"""tuple: The words that tell a reader the site is a model (rule 16, DOC.87). The site speaks in universe: a working navigation system that
simply knows things. Where a number came from, what is a stand-in, what was fitted to a result and how sure a part is are said only on
the References pages (`/references`); every other page says the thing in universe and, where a reader may want the why, links to the
References entry (`web/helpers.reference_link`). "Level" is not on the list: "sea level" and "ground level" are in universe, and the
provenance sense ("tier and level") always comes with "tier" or "tweak", which are. A bare "seed" is not on it either: a seed vault, or the
seeds of a plant, is in universe; only the generator's random seed is out of it."""

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


_CLAUSE_STARTS = frozenset("""which who whom whose that where when while but so because then if though although since as it its they
there this these those he she we you i not also however thus therefore yet whereas unless until after before once""".split())
"""frozenset: Words that begin a clause, not a list item: "red, which is big and cold" is not a list of three."""

_LEADING_PHRASE = frozenset("""in on at by for from with without during under over near within across between through into onto
about above below beyond before after if when while since as once because although though yet then thus however also instead
first second third finally""".split())
"""frozenset: Words that open an introductory phrase ("In the north, rock and ice cover it") whose comma is not a list comma."""

_FINITE_VERBS = frozenset("""is are was were be been being am has have had can could may might will would shall should must do does did
holds hold writes write makes make means mean gives give takes take shows show lies lie sits sit stays stay reaches reach forms form
rises rise falls fall runs run grows grow keeps keep needs need lets let sets set gets get comes come goes go says say tells tell
calls call counts count rates rate""".split())
"""frozenset: Verbs that make a stretch of words a clause, not one item of a list."""

_DETERMINERS = frozenset("a an the its their his her this that some each every no any".split())
"""frozenset: Words an item of a list of noun phrases may open with; a list repeats them ("the day side, the ring and the night side")."""

_AND_OR = re.compile(r"^(?P<first>[^,;:]+?) (?P<join>and|or) (?P<last>[^,;:]+)$")
_MAX_ITEM_WORDS = 4


def _bare(word):
    return word.strip("()\"'“”‘’.:").lower()


def _list_item(words):
    """True when `words` can be one item of a list: a few words, no clause opener, no finite verb."""
    if not words or len(words) > _MAX_ITEM_WORDS:
        return False
    bare = [_bare(w) for w in words]
    return bare[0] not in _CLAUSE_STARTS and bare[0] not in ("and", "or") and not any(w in _FINITE_VERBS for w in bare)


def _tail_item(segment):
    """The words a segment can end a list item with: those after its last verb, so "The surface is rock" ends with "rock"."""
    words = segment.split()
    for index in range(len(words) - 1, -1, -1):
        if _bare(words[index]) in _FINITE_VERBS:
            return words[index + 1:]
    return words[-_MAX_ITEM_WORDS:]


def unpunctuated_lists(text):
    """
    The lists in `text` that have no comma before their last "and" or "or", as `[(before, first, join, last)]`: the item before the comma,
    the item after it, "and" or "or", and the words that follow ("red, green and blue" gives `("red", "green", "and", "blue")`).
    """
    plain = _QUOTED.sub(" ", text)
    found = []
    for sentence in sentences(plain) or [plain]:
        depth, start, segments = 0, 0, []
        for index, ch in enumerate(sentence):
            depth += ch == "("
            depth -= ch == ")"
            if ch == "," and depth <= 0 and not (sentence[index - 1:index].isdigit() and sentence[index + 1:index + 2].isdigit()):
                segments.append(sentence[start:index])
                start = index + 1
        segments.append(sentence[start:])
        for k in range(1, len(segments)):
            match = _AND_OR.match(segments[k].strip().rstrip(".!?"))
            if not match:
                continue
            first_words = match["first"].split()
            last_words = match["last"].split()
            if not _list_item(first_words) or not last_words or _bare(last_words[0]) in _CLAUSE_STARTS:
                continue
            first_word = _bare(first_words[0])
            if first_word in _LEADING_PHRASE or first_word.endswith("ly"):
                continue    # "Clouds of gas and dust, by family and contents": a phrase after the comma, not a list item
            if first_word in ("a", "an", "the") and len(last_words) > 5:
                continue    # "a giant, a deep hydrogen and helium envelope over a core": two items, the second with a pair inside
            if first_word in _DETERMINERS and _bare(last_words[0]) not in _DETERMINERS:
                continue    # "the chronic and acute inspired edges": two adjectives for one noun, not three items
            if len(first_words) > 1 and first_word.endswith("ing"):
                continue    # "living on sulfur and carbon": the pair belongs to the verb, not a list of its own
            if any(_bare(w) in _FINITE_VERBS for w in last_words[:2]):
                continue
            previous = segments[k - 1].split()
            if k == 1 and previous and (_bare(previous[0]) in _LEADING_PHRASE or (len(previous) == 1 and _bare(previous[0]).endswith("ly"))):
                continue
            before = _tail_item(segments[k - 1])
            if not _list_item(before):
                continue
            found.append((" ".join(before), match["first"], match["join"], match["last"]))
    return found


def oxford_commas(text):
    """
    The lists in `text` that have no comma before their last "and" or "or", as `[("oxford-comma", detail)]`: "red, green and blue" is
    flagged and "red, green, and blue" is not (Boss, 2:50 PM EDT 2026-10-31). A list is three or more short items (no verb, no clause
    opener) with a comma between the first two; a clause ("X, which is big and cold"), an introductory phrase ("In the north, rock and
    ice") and a pair ("red and blue") are not lists. Text in quotation marks is not read.
    """
    return [("oxford-comma", f"'{before}, {first} {join} {last}': write a comma before '{join}'") for before, first, join, last in unpunctuated_lists(text)]


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


def check_sentence(sentence, references=False):
    """
    The rules one sentence breaks, as `[(rule, detail)]` (empty when it keeps them all). `references=True` reads a sentence of a References
    page (`/references`, DOC.88), where the process words of `OUT_OF_UNIVERSE` (stand-in, tweak, epoch, placeholder) are allowed and every
    other rule still holds. The rules:

    * `length`: more than `MAX_WORDS` words;
    * `ideas`: more than `MAX_SEMICOLONS` semicolons;
    * `banned`: a word or phrase of `BANNED`, with what to write instead;
    * `class-name`: a life class named by its family ("Class 3") where the class has letters ("Class 3a");
    * `oxford-comma`: a list of three or more with no comma before its last "and" or "or" ("red, green and blue"; Boss, 2:50 PM EDT
      2026-10-31);
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
        if references and phrase in OUT_OF_UNIVERSE:
            continue
        if pattern.search(sentence):
            found.append(("banned", f"{phrase!r}: write {plain}"))
    for digits, letter in _CLASS_NAME.findall(sentence):
        if f"{digits}{letter}" not in life_classes.BY_ID:
            found.append(("class-name", f"'Class {digits}{letter}' is not a class: the classes are {', '.join(life_classes.BY_ID)}"))
    found += british_spellings(sentence)
    found += oxford_commas(sentence)
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
