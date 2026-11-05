# planetgen/physics/class_code.py

"""
The class code: one letter plus one digit (GEN.1220; design: reports/class-code-system-2026-10-31.md)
=====================================================================================================

A class code is **a letter plus a digit**. The letter (`planet_class`) names a family or a stand-alone kind of body; the digit (`class_variant`, 0 to 9)
names the position inside the family on one scale for every family. Digit 0 is never written, so M0 is written M, and a letter that exists
keeps its meaning with digit 0. Searching or filtering by a letter matches all its digits.

| digit | position | differs from 0 in |
|---|---|---|
| 0 | the ideal (or the only version of a stand-alone class) | |
| 1 | too hot | temperature only |
| 2 | too cold | temperature only |
| 3 | too wet | liquid cover only |
| 4 | too dry | liquid cover only |
| 5 to 8 | two of those at once (hot and wet, hot and dry, cold and wet, cold and dry) | both |
| 9 | another domain (gravity, pressure, light) is the limit | one named domain |

Digits 5 to 9 are reserved: one is switched on only when a registered stand-in (GEN.725) gives its ranges and a shadow report shows enough bodies.
Every place that prints or sorts a class goes through `format_code`, `parse_code` and `sort_key`, so a digit appears everywhere at once.

The letters Boss confirmed (11:59 AM EDT 2026-10-31) for the families the best-fit classifier (GEN.1207) will fill: `LETTERS` below. Q stays retired and W spare.

Use:

    class_code.format_code("F", 2)      # "F2"
    class_code.format_code("M", 0)      # "M"
    class_code.parse_code("F2")         # ("F", 2)
    class_code.parse_code("M")          # ("M", 0)
"""

import re

DIGITS = range(10)
RESERVED_DIGITS = range(5, 10)

POSITION = {
    0: "ideal",
    1: "too hot",
    2: "too cold",
    3: "too wet",
    4: "too dry",
    5: "too hot and too wet",
    6: "too hot and too dry",
    7: "too cold and too wet",
    8: "too cold and too dry",
    9: "limited by another domain",
}
"""dict: What each digit says about the body, on the one scale of every family."""

LETTERS = {
    "M": "temperate world, the oxygen-water ideal",
    "R": "Hothouse world, M's too hot (Boss)",
    "P": "glaciated world, M's too cold",
    "O": "ocean world, M's too wet",
    "H": "desert world, M's too dry",
    "F": "anoxic water family",
    "G": "dense carbon dioxide family",
    "Y": "sulfuric acid family (Boss: a Demon world)",
    "U": "methane-ethane family",
    "X": "ammonia-water family",
    "Z": "hot thin-air world",
}
"""dict: The letters of the families the best-fit classifier names, each with what it stands for (the stand-alone classes keep the letters of `planet_classes`)."""

FAMILIES = {"F": "anoxic water", "G": "dense CO2", "Y": "sulfuric acid", "U": "methane-ethane", "X": "ammonia-water"}
"""dict: The families whose digits 1 to 4 are in use; the oxygen-water family keeps its five positions in five letters (M, R, P, O, H) and so digit 0."""

FIELDS = ("class_variant",)
"""tuple: What a planet or a moon stores of this module, as `planets` and `moons` columns (GEN.1220)."""

_PATTERN = re.compile(r"^([A-Z])([0-9])?$")


class ClassCodeError(ValueError):
    """A class code that is not one letter and an optional digit."""


def format_code(letter, variant=0):
    """The written class code: the letter, with the digit when it is not 0 (`None` for a body without a class)."""
    if letter is None:
        return None
    if variant not in DIGITS and variant is not None:
        raise ClassCodeError(f"class variant {variant!r} is not a digit from 0 to 9")
    return f"{letter}{variant}" if variant else letter


def parse_code(text):
    """`(letter, digit)` of a written class code ('M' is ('M', 0), 'F2' is ('F', 2))."""
    match = _PATTERN.match((text or "").strip().upper())
    if match is None:
        raise ClassCodeError(f"{text!r} is not a class code: one letter and an optional digit")
    return match.group(1), int(match.group(2) or 0)


def sort_key(letter, variant=0):
    """The key that sorts class codes by letter and then by digit."""
    return letter or "", variant or 0
