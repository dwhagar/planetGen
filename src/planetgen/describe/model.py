# planetgen/describe/model.py

"""
The Document's node constructors and its checker (docs/design/architecture/system-description.md,
section 3).

A Document is plain dicts, lists, strings and numbers, so it is already JSON.
Inline text is a list whose items are a string or one of:

- `{"qty": {"kind", "value", "unit", ...}}`: a quantity, never formatted here. `kind` is
  a ladder family (`planetgen.util.ladders`) or one of `OTHER_KINDS`; `value` is the number in `unit`, a key of the unit
  registry (`planetgen.physics.unit_registry`) or one of `PLAIN_UNITS` (the unit the value was held in, not converted). Extra
  keys: `sol` (the ratio to Sol), `companions` (whether the text shows the units beside the main one),
  `note` (a `ratio`'s unit, scale or noun, written after the number: "times Earth's DNA-damaging dose"; a ratio is never
  written without one).
- `{"ref": "planet:31", "name": "40 Eridani II"}`: a link to another object.
- `{"strong": "Pressure: Red"}`: emphasis, written per format by the renderer.
- `{"chem": "H2O"}`: a chemical formula (UX.149), written per format by `planetgen.util.chem`.
- `{"swatch": {"hex": "#99bfff", "word": "pale blue"}}`: a colour (GEN.420): the word, with the colour as a chip where the format has one.
"""

import re

from planetgen.physics import unit_registry
from planetgen.util import chem as chemistry
from planetgen.util import ladders

LADDER_KINDS = tuple(ladder.quantity for ladder in ladders.LADDERS)
"""tuple: The quantities that have a ladder."""

OTHER_KINDS = ("radius", "luminosity", "dose", "flux", "ratio", "angle", "density", "count")
"""tuple: Quantities with no ladder (yet). `radius` is a body radius: always km in scientific notation."""

KINDS = LADDER_KINDS + OTHER_KINDS

DISPLAY_UNITS = {
    "length": "m", "speed": "m_s", "time": "s", "temperature": "K", "temperature difference": "dK", "fraction": "frac", "pressure": "Pa", "gravity": "m_s2", "mass": "kg",
    "power": "W", "flux": "W_m2", "number density": "m-3", "magnetic field": "T", "radius": "m", "luminosity": "W",
    "angle": "deg", "density": "kg_m3", "dose": "Sv/yr", "ratio": "ratio", "count": "count",
}
"""dict: The unit a kind's text is written from, a key of the unit registry or one of `PLAIN_UNITS`. A quantity held in
another unit of the kind is converted to it (`quantity.display_value`)."""

PLAIN_UNITS = ("ratio", "count", "Sv/yr")
"""tuple: Units with no registry kind: a plain number, a count and a dose rate. They convert to nothing."""

SECTION_KINDS = ("system", "star", "black-hole", "neutron-star", "pair", "planet", "moon", "belt", "comet", "nebula",
                 "quasar", "rogue-planet", "interstellar-comet", "asteroid-field", "supernova-remnant")
BLOCK_TYPES = ("facts", "ratings", "composition", "prose")
PROSE_ROLES = ("overview", "summary", "description", "depth", "interior", "life-class", "life", "timeline", "milestone", "flavour")


def qty(kind, value, unit=None, **extra):
    """A quantity inline item: `value` in `unit` (a registry unit key or one of `PLAIN_UNITS`; the kind's own display
    unit when omitted)."""
    if kind not in KINDS:
        raise ValueError(f"unknown quantity kind {kind!r}")
    if kind == "ratio" and not extra.get("note"):
        raise ValueError("a ratio is never written as a bare number: pass `note`, what it is or its scale")
    unit = unit or DISPLAY_UNITS[kind]
    if unit not in PLAIN_UNITS:
        unit_registry.unit(unit)
    number = None if value is None else (int(value) if kind == "count" else float(value))
    return {"qty": {"kind": kind, "value": number, "unit": unit, **extra}}


def measured(kind, quantity, **extra):
    """
    A quantity inline item (`qty`) for a `Quantity` (`planetgen.physics.quantity`) shown as `kind`: its value and the unit
    it is held in, so no call site multiplies by a unit factor (GEN.348, GEN.574). `None` is an absent value.
    """
    return qty(kind, None) if quantity is None else qty(kind, quantity.value, quantity.unit, **extra)


def ref(ref_text, name):
    """A link inline item."""
    return {"ref": ref_text, "name": name}


def strong(text):
    """An emphasis inline item."""
    return {"strong": text}


def swatch(hex_colour, word):
    """A colour swatch inline item (GEN.420): `hex_colour` is `#rrggbb`, `word` what it is called."""
    return {"swatch": {"hex": hex_colour, "word": word}}


def chem(formula):
    """A chemical formula inline item; `formula` must read as one (`planetgen.util.chem.parse`)."""
    chemistry.parse(formula)
    return {"chem": formula}


def chem_text(words):
    """An inline list from a text that names chemicals (a nebula's contents): each formula `chem.split` finds there
    is a `chem` item. An empty or missing text is an empty list."""
    if not words:
        return []
    return inline(*(chem(part) if isinstance(part, chemistry.Formula) else part for part in chemistry.split(words)))


def inline(*parts):
    """An inline list from strings and inline items (or lists of them), neighbouring strings joined."""
    out = []
    for part in parts:
        for item in (part if isinstance(part, list) else [part]):
            if isinstance(item, str) and out and isinstance(out[-1], str):
                out[-1] += item
            elif item != "" and item is not None:
                out.append(item)
    return out


def fact(key, label, *value):
    """One row of a facts table; the value is inline."""
    return {"key": key, "label": label, "value": inline(*value)}


def facts(caption, rows):
    """A facts block; `rows` are `fact()` rows, `None` rows are dropped."""
    return {"type": "facts", "caption": caption, "rows": [row for row in rows if row is not None]}


def prose(role, *sentences):
    """A prose block: `sentences` are inline lists (or strings); empty ones are dropped."""
    return {"type": "prose", "role": role, "text": [s for s in (inline(sentence) for sentence in sentences) if s]}


def composition(label, parts=None, summary=None):
    """A composition block: `parts` are `{"name", "share"}` (share 0 to 1, or `"level"` a word like
    "high"); `summary` is the same thing as a phrase, for a composition without numbers."""
    return {"type": "composition", "label": label, "parts": list(parts or []), "summary": summary}


def ratings(overall, rows, gear):
    """A ratings block: `overall` is `{"name", "score", "tier"}`, `rows` are `{"domain", "tier", "score",
    "text"}` (`text` inline) and `gear` an inline list."""
    return {"type": "ratings", "overall": overall, "rows": rows, "gear": inline(gear)}


def section(kind, ref_text, title, blocks=(), children=()):
    """A section."""
    return {"kind": kind, "ref": ref_text, "title": title, "blocks": [b for b in blocks if b], "children": list(children)}


def document(kind, ref_text, title, sections):
    """The root."""
    return {"kind": kind, "ref": ref_text, "title": title, "sections": list(sections)}


def walk_inline(node):
    """Every inline item (strings included) anywhere in a Document or part of one."""
    if isinstance(node, dict):
        if "qty" in node or "ref" in node and "name" in node or ("strong" in node or "chem" in node or "swatch" in node) and len(node) == 1:
            yield node
            return
        for value in node.values():
            yield from walk_inline(value)
    elif isinstance(node, list):
        for item in node:
            yield from walk_inline(item)
    elif isinstance(node, str):
        yield node


def sections(document_or_section):
    """Every section under (and including, for a section) the node, depth first."""
    if "sections" in document_or_section:
        for child in document_or_section["sections"]:
            yield from sections(child)
        return
    yield document_or_section
    for child in document_or_section["children"]:
        yield from sections(child)


SPLIT_KEYS = {"star": "stars", "planet": "planets", "moon": "moons", "belt": "belts", "comet": "comets"}
"""dict: The section kinds `split_by_body` gives their own entry, each to the key it is filed under."""


def split_by_body(document):
    """
    A system's Document split into the pieces the web page shows under each row: `overview` (the blocks above the
    list: a binary pair's own facts, the system overview, a close pair's age sentence and any system flavour text, in
    page order) and `stars`, `planets`, `moons`, `belts` and `comets`, each mapping the `id` in a section's `ref`
    (`"planet:12"` -> `"12"`) to that section's own blocks, without its heading or its child sections. A section
    with no `ref` (a body that was never stored) has no entry.
    """
    top = document["sections"][0]
    overview = list(top["blocks"])
    pair = next((child for child in top["children"] if child["kind"] == "pair"), None)
    if pair is not None:
        overview = [pair["blocks"][0], *overview[:1], *pair["blocks"][1:], *overview[1:]]
    pieces = {"overview": overview, **{key: {} for key in SPLIT_KEYS.values()}}
    for piece in sections(top):
        if piece is top or not piece["ref"]:
            continue
        kind, _, db_id = piece["ref"].partition(":")
        if kind in SPLIT_KEYS:
            pieces[SPLIT_KEYS[kind]][db_id] = piece["blocks"]
    return pieces


def score_text(score):
    """A PHI-6 score as the pages write it: "0.82 of 1", never the bare 0.82 (GEN.772)."""
    return f"{score:.2f} of 1"


BARE_NUMBER = re.compile(r"[-+−]?[\d,]*\.?\d+(?:e[-+]?\d+| ?× ?10(?:[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+|<sup>[-−\d]+</sup>))?")
"""re.Pattern: A number alone, as the pages write it ("0.71", "1,200", "1.74 × 10⁶"); a fact whose value matches it has no
unit, scale or noun."""


def is_bare_number(text):
    """Whether `text` is a number alone, with no unit, scale or noun beside it."""
    return BARE_NUMBER.fullmatch(text.strip()) is not None


def _check_no_bare_number(block):
    """Raises ValueError for a fact row whose value reads as a number alone (GEN.772: never a number without units or
    context). A count with a noun for its label is `count`, a quantity has its unit, a ratio its note."""
    from planetgen.describe.quantity import inline_text

    for row in block["rows"]:
        text = inline_text(row["value"])
        if text and is_bare_number(text) and not all(
                isinstance(item, dict) and item.get("qty", {}).get("kind") == "count" for item in row["value"]):
            raise ValueError(f"the fact {row['label']!r} reads {text!r}, a number with no unit, scale or noun")


def check(node):
    """Raises ValueError unless `node` (a Document) has the documented shape: known section kinds, block
    types and quantity kinds, finite-or-None quantities, JSON-safe values."""
    import json
    json.dumps(node, allow_nan=False)
    for sec in sections(node):
        if sec["kind"] not in SECTION_KINDS:
            raise ValueError(f"unknown section kind {sec['kind']!r}")
        for block in sec["blocks"]:
            if block["type"] not in BLOCK_TYPES:
                raise ValueError(f"unknown block type {block['type']!r}")
            if block["type"] == "prose" and block["role"] not in PROSE_ROLES:
                raise ValueError(f"unknown prose role {block['role']!r}")
            if block["type"] == "facts":
                _check_no_bare_number(block)
    for item in walk_inline(node):
        if isinstance(item, dict) and "qty" in item and item["qty"]["kind"] not in KINDS:
            raise ValueError(f"unknown quantity kind {item['qty']['kind']!r}")
        if isinstance(item, dict) and "qty" in item and item["qty"]["unit"] not in PLAIN_UNITS:
            unit_registry.unit(item["qty"]["unit"])
        if isinstance(item, dict) and "swatch" in item:
            swatch_item = item["swatch"]
            if not re.fullmatch(r"#[0-9a-f]{6}", str(swatch_item.get("hex"))) or not swatch_item.get("word"):
                raise ValueError(f"a swatch needs a #rrggbb colour and a word, not {swatch_item!r}")
