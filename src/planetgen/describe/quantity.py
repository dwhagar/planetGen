# planetgen/describe/quantity.py

"""
Writes a Document's inline text and quantities as a string (docs/design/architecture/system-description.md,
section 4). A quantity goes through the site's one formatter (`util/format.format_quantity`, so
it reads the same as every page and as `static/unitladder.js`); the page formats differ only in
how a power of ten and emphasis are marked.
"""

import re

from planetgen.physics.quantity import Quantity
from planetgen.util import chem, ladders
from planetgen.util.format import format_number, format_quantity, scientific_text

LADDER_BY_KIND = {ladder.quantity: ladder for ladder in ladders.LADDERS}
"""dict: Each ladder family's `kind`, to its ladder."""

_SUPERSCRIPT_TO_DIGIT = str.maketrans("⁻⁰¹²³⁴⁵⁶⁷⁸⁹", "-0123456789")
_POWER_OF_TEN = re.compile(r"(\d+(?:\.\d+)?) × 10([⁻⁰¹²³⁴-⁹]+)")

FORMATS = ("text", "markdown", "wikitext")
"""tuple: The plain-string formats an inline list is written in."""


def quantity_text(quantity):
    """The plain text of a quantity (`{"kind", "si", ...}`), with Unicode superscripts for a power of ten."""
    kind, value = quantity["kind"], quantity["si"]
    if kind in LADDER_BY_KIND:
        return format_quantity(LADDER_BY_KIND[kind], value, companions=quantity.get("companions", True))
    if value is None:
        return "–"
    if kind == "radius":
        return f"{scientific_text(Quantity('length', value, 'm').km)} km"
    if kind == "luminosity":
        return format_quantity(ladders.POWER, value)
    if kind == "count":
        return f"{int(value):,}"
    if kind == "ratio":
        return format_number(value, ".3f")
    if kind == "angle":
        return f"{format_number(value, '.1f')}°"
    if kind == "dose":
        return f"{format_number(value * 1000.0, '.3g')} mSv a year" if abs(value) < 1.0 else f"{format_number(value, '.3g')} Sv a year"
    if kind == "density":
        return f"{format_number(value, '.3g')} kg/m³"
    raise ValueError(f"no text for quantity kind {kind!r}")


def _powers(text, fmt):
    if fmt == "text":
        return text

    def replace(match):
        exponent = match.group(2).translate(_SUPERSCRIPT_TO_DIGIT)
        if fmt == "markdown":
            return f"{match.group(1)} × 10<sup>{exponent}</sup>"
        return "{{Exp|%s|%s}}" % (match.group(1), exponent)

    return _POWER_OF_TEN.sub(replace, text)


def inline_text(items, fmt="text"):
    """An inline list as one string in `fmt` (`FORMATS`)."""
    if fmt not in FORMATS:
        raise ValueError(f"unknown format {fmt!r} (expected one of {', '.join(FORMATS)})")
    out = []
    for item in items:
        if isinstance(item, str):
            out.append(item)
        elif "qty" in item:
            out.append(_powers(quantity_text(item["qty"]), fmt))
        elif "strong" in item:
            out.append({"text": "{}", "markdown": "**{}**", "wikitext": "'''{}'''"}[fmt].format(item["strong"]))
        elif "chem" in item:
            # Markdown and MediaWiki both take <sub> and <sup>; a terminal gets Unicode.
            out.append(chem.text(item["chem"]) if fmt == "text" else chem.html(item["chem"]))
        else:
            out.append(item["name"])
    return "".join(out)
