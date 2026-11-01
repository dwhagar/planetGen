# tests/test_plain_words.py

"""
The plain-words style guide holds for the text the site writes (DOC.83, `docs/design/interface/text-style.md`).

The checkable rules are `planetgen/util/plainwords.py`: a sentence of at most 30 words, at most one semicolon, none of the
banned words, a life class named by its whole id ("Class 3a", not "Class 3"), and table labels in sentence case. This test reads
the text three ways: the prose and the table labels of the Document of ten seeded systems (`describe.build`), the class pages
(`web/lib/classref.py`) and the label of every row of the property sheet (`describe/sheet.py`).

Text that breaks a rule today is listed in `KNOWN` with the item that will rewrite it (UX.252: planets, moons, belts, comets and
their class pages; UX.253: stars, systems and star-class pages); text that is not listed fails. So old text is rewritten once and
new text keeps the guide from the day it is written. A rewrite item deletes its entries and is done when none carry its code.
"""

import re
from pathlib import Path

import pytest

from planetgen.describe import build as describe_build, sheet
from planetgen.describe.quantity import inline_text
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import plainwords
from planetgen.web.lib import classref

from tests.fuzz_support import deterministic_entropy

SEEDS = range(10)
"""range: The seeds of the systems whose Document is read (the `KNOWN` list was made from eighty)."""

OWNERS = ("UX.252", "UX.253")
"""tuple: The items that rewrite the text: planets, moons and belts; stars and systems."""

PLANET_KINDS = ("planet", "moon", "belt", "comet")
PLANET_CLASS_TYPES = ("planet", "rogue-planet", "asteroid-field", "comet")

GUIDE = Path(__file__).resolve().parents[2] / "docs" / "design" / "interface" / "text-style.md"

KNOWN = {
    ('class-name', 'Animals have appeared among'): 'UX.252',  # Animals have appeared among Class 1 life.
    ('class-name', 'Cells with internal structure'): 'UX.252',  # Cells with internal structure have appeared among Class 3 life.
    ('class-name', 'Class # life has'): 'UX.252',  # Class 3 life has grown into colonies, mats or bodies big enough to see.
    ('ideas', 'Size digit d means'): 'UX.252',  # Size digit d means a radius of 10^d to 10^(d+1) AU, from 1 to 4: 1: 10 to 100 AU; 2: 100 to 1,000 AU; 3: 1 × 1
    ('label-case', 'Atmosphere Colour'): 'UX.252',
    ('label-case', 'Distance from the Star'): 'UX.252',
    ('label-case', 'Inner Edge'): 'UX.252',
    ('label-case', 'Internal Pressure'): 'UX.252',
    ('label-case', 'Moon Wobble'): 'UX.252',
    ('label-case', 'Nearest Meteorite Group'): 'UX.252',
    ('label-case', 'Nucleus Diameter'): 'UX.252',
    ('label-case', 'Outer Edge'): 'UX.252',
    ('label-case', 'Surface Pressure'): 'UX.252',
    ('label-case', 'Surface Temperature'): 'UX.252',
    ('spelling', 'Aluminium-26 at formation'): 'UX.252',
    ('spelling', 'Atmosphere Colour'): 'UX.252',
    ('spelling', 'It is the overall'): 'UX.252',
    ('spelling', "Offset from the system's centre, x"): 'UX.253',
    ('spelling', 'Soundings read a gas'): 'UX.252',
    ('length', 'A class is never'): 'UX.252',  # A class is never drawn first: once a planet or moon is complete, its class is the first class in the check ord
    ('length', 'A marginal world, Earth-like'): 'UX.252',  # A marginal world, Earth-like but hard to breathe without a mask and not toxic: air of 30 kPa to 3 bar at 273 t
    ('length', 'Size digit d means'): 'UX.252',  # Size digit d means a radius of 10^d to 10^(d+1) AU, from 1 to 4: 1: 10 to 100 AU; 2: 100 to 1,000 AU; 3: 1 × 1
    ('length', 'Soundings read a gas'): 'UX.252',  # Soundings read a gas giant: a deep envelope of hydrogen and helium over a compact core, with hydrogen squeezed
    ('length', 'Soundings read an ice'): 'UX.252',  # Soundings read an ice giant: a thin envelope of hydrogen and helium over a deep mantle of hot water, ammonia a
    ('length', 'The dilute core (metallic'): 'UX.252',  # The dilute core (metallic fluid) holds 2.8 Earth masses, from 0.98 Mbar down to 1.17 Mbar; it is rock 36%, sup
    ('length', 'The molecular hydrogen and'): 'UX.252',  # The molecular hydrogen and helium envelope (molecular fluid) holds 3.3 Earth masses, from 1 bar down to 0.28 M
}
"""dict: `{(rule, start of the sentence or the whole label): the item that rewrites it}` for the text that breaks a rule today. The
comment is one example. A rewrite item removes its entries; nothing is added without an owner."""


def shape(sentence):
    """The start of a sentence with its numbers blanked: what `KNOWN` names it by."""
    return " ".join(re.sub(r"\d[\d.,]*", "#", sentence).split()[:4])


def _sections(node):
    yield node
    for child in (node.get("children") or []) + (node.get("sections") or []):
        yield from _sections(child)


def _document_text(seed):
    """`[(owner, kind, text)]` for the sentences and labels of one seeded system's Document."""
    with deterministic_entropy(7000 + seed):
        system = StarSystem(system_config=SystemConfig())
    found = []
    for section in _sections(describe_build.describe(system)):
        owner = "UX.252" if section["kind"] in PLANET_KINDS else "UX.253"
        for block in section.get("blocks") or []:
            if block["type"] == "prose":
                found += [(owner, "sentence", inline_text(sentence)) for sentence in block["text"]]
            elif block["type"] == "facts":
                found += [(owner, "label", row["label"]) for row in block["rows"]]
    return found


def _breaks(owner, kind, text):
    """`[(key, owner, text)]` for each rule `text` (a sentence or a label) breaks."""
    if kind == "label":
        return [((rule, text), owner, text) for rule, _detail in plainwords.check_label(text)]
    return [((rule, shape(sentence)), owner, sentence) for sentence in plainwords.sentences(text)
            for rule, _detail in plainwords.check_sentence(sentence)]


def _unlisted(found):
    """The breaks in `found` (`[(owner, kind, text)]`) that `KNOWN` does not list, as readable lines."""
    lines = {}
    for owner, kind, text in found:
        for key, who, example in _breaks(owner, kind, text):
            if key not in KNOWN:
                lines[key] = f"{key[0]} ({who}): {example[:160]!r}"
    return sorted(lines.values())


ADVICE = ("\nWrite it to the guide (docs/design/interface/text-style.md). Only text that already existed is listed in KNOWN, under the "
          "item that will rewrite it; new text is not added to it.")


def test_the_prose_and_labels_of_a_document_keep_the_guide_except_the_listed():
    found = [item for seed in SEEDS for item in _document_text(seed)]
    assert len(found) > 1500, "the corpus is too small to mean anything"
    unlisted = _unlisted(found)
    assert not unlisted, "text that breaks the plain-words guide:\n  " + "\n  ".join(unlisted) + ADVICE


def test_the_class_pages_keep_the_guide_except_the_listed():
    if not classref.preload():
        pytest.skip("a downloaded data set the star classes read is missing")
    found = []
    for slug, kind in classref.catalog().items():
        owner = "UX.252" if slug in PLANET_CLASS_TYPES else "UX.253"
        texts = [kind.get("summary", "")] + list(kind.get("notes", [])) + [entry["summary"] for entry in kind["classes"].values()]
        found += [(owner, "sentence", text) for text in texts if isinstance(text, str)]
        found += [(owner, "label", label) for entry in kind["classes"].values() for label, _value in entry["facts"]]
    assert len(found) > 100
    unlisted = _unlisted(found)
    assert not unlisted, "class page text that breaks the plain-words guide:\n  " + "\n  ".join(unlisted) + ADVICE


def test_the_property_sheet_labels_are_in_sentence_case():
    labels = {spec.label for kind in ("star", "planet", "moon", "belt", "comet") for spec in sheet.specs_for(kind)}
    assert len(labels) > 100
    unlisted = _unlisted([("UX.252", "label", label) for label in sorted(labels)])
    assert not unlisted, "sheet labels that break the plain-words guide:\n  " + "\n  ".join(unlisted) + ADVICE


def test_every_listed_exception_names_its_owner_and_a_rule():
    for (rule, start), owner in KNOWN.items():
        assert owner in OWNERS, (rule, start, owner)
        assert rule in {"length", "ideas", "banned", "class-name", "label-case", "spelling"}, rule


def test_new_text_that_breaks_a_rule_and_is_not_listed_fails():
    unlisted = _unlisted([("UX.252", "sentence", "A brand new sentence uses approximately two words."),
                          ("UX.252", "label", "Brand New Label"),
                          ("UX.252", "sentence", "A brand new sentence that is fine.")])
    assert len(unlisted) == 2 and unlisted[0].startswith("banned") and unlisted[1].startswith("label-case")


@pytest.mark.parametrize("sentence, rule", [
    ("This sentence runs on " + "and on " * 16 + "until it has far too many words in it to be read at a glance.", "length"),
    ("It is hot; it is dry; it is small.", "ideas"),
    ("The star lies approximately 5 pc away.", "banned"),
    ("Dust is blown out in order to clear the disk.", "banned"),
    ("Some are rocky, e.g. Mars.", "banned"),
    ("Class 3 life has grown into mats.", "class-name"),
    ("Class 3a life has grown into mats.", None),
    ("Class 2 life lives in the oceans.", None),
    ("The Class 3 family has two variants.", None),
    ("The colour of the sky is grey.", "spelling"),
    ("The color of the sky is gray.", None),
    ("It is 5 kilometres across.", "spelling"),
    ("The company is organised in the centre of the town.", "spelling"),
    ("We analyse the data and the analyses agree.", "spelling"),
    ("The rise is precise, so exercise the promise.", None),
    ('The book "The Colour Purple" is quoted.', None),
    ("The Centre for Planetary Science wrote it.", None),
    ("Luminosity class 0 stars are hypergiants.", None),
    ("It is a hot world; the sea is dry.", None),
    ("The planet is about 5 pc away.", None),
])
def test_the_checker_flags_each_rule(sentence, rule):
    found = {name for name, _detail in plainwords.check_sentence(sentence)}
    assert found == ({rule} if rule else set())


def test_a_british_spelling_says_the_american_one():
    said = dict((detail.split("'")[1], detail.split("'")[3]) for _rule, detail in plainwords.british_spellings(
        "Colourful kilometres, centred neighbours, organisation, analysed, sulphur, labelled."))
    assert said == {"Colourful": "colorful", "kilometres": "kilometers", "centred": "centered", "neighbours": "neighbors",
                    "organisation": "organization", "analysed": "analyzed", "sulphur": "sulfur", "labelled": "labeled"}


@pytest.mark.parametrize("label, bad", [
    ("Surface temperature", False), ("Distance from the star", False), ("Urey ratio", False), ("Radius (Earth radii)", False),
    ("Surface Temperature", True), ("Distance from the Star", True), ("Primary's Planetary Limit", True),
    ("Atmosphere color", False), ("Atmosphere colour", True),
])
def test_the_checker_reads_a_label_in_sentence_case(label, bad):
    assert bool(plainwords.check_label(label)) == bad


def test_sentences_split_at_a_full_stop_but_not_in_a_number_or_an_abbreviation():
    assert plainwords.sentences("It is 1.06 Mbar deep. That is, e.g. Jupiter, big. Yes!") == [
        "It is 1.06 Mbar deep.", "That is, e.g. Jupiter, big.", "Yes!"]


def test_the_guide_names_every_rule_and_the_two_rewrite_items():
    text = GUIDE.read_text(encoding="utf-8")
    for needle in ("`MAX_WORDS`", "`MAX_SEMICOLONS`", "`BANNED`", "Class 3a", "sentence case", "American", "UX.252", "UX.253", "KNOWN"):
        assert needle in text, needle


def test_the_stored_system_flavour_lines_are_in_plain_words():
    """GEN.1212: `tuning.SYSTEM_FLAVOR` is drawn by index and stored as the line, so its length and order are fixed (the draw
    must not move) and every line, read after "Sensors show", meets the guide."""
    from planetgen import tuning
    lines = tuning.SYSTEM_FLAVOR
    assert len(lines) == 23 and len(set(lines)) == 23
    for line in lines:
        for sentence in plainwords.sentences("Sensors show " + line):
            assert not plainwords.check_sentence(sentence), (sentence, plainwords.check_sentence(sentence))
