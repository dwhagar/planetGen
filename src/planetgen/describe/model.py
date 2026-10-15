# planetgen/describe/model.py

"""
The Document's node constructors and its checker (docs/design/system-description.md,
section 3).

A Document is plain dicts, lists, strings and numbers, so it is already JSON.
Inline text is a list whose items are a string or one of:

- `{"qty": {"kind", "si", ...}}`: a quantity, never formatted here. `kind` is
  a ladder family (`planetgen.util.ladders`) or one of `OTHER_KINDS`; `si` is the
  value in the family's base unit. Extra keys: `sol` (the ratio to Sol), `unit`
  (the unit symbol of a kind with no ladder yet).
- `{"ref": "planet:31", "name": "40 Eridani II"}`: a link to another object.
- `{"strong": "Pressure: Red"}`: emphasis, written per format by the renderer.
"""

from planetgen.util import ladders

LADDER_KINDS = tuple(ladder.quantity for ladder in ladders.LADDERS)
"""tuple: The quantities that have a ladder."""

OTHER_KINDS = ("radius", "luminosity", "dose", "flux", "ratio", "angle", "density", "count")
"""tuple: Quantities with no ladder (yet). `radius` is a body radius: always km in scientific notation."""

KINDS = LADDER_KINDS + OTHER_KINDS

SECTION_KINDS = ("system", "star", "pair", "planet", "moon", "belt", "comet")
BLOCK_TYPES = ("facts", "ratings", "composition", "prose")
PROSE_ROLES = ("overview", "summary", "description", "depth", "life", "timeline", "flavour")


def qty(kind, si, **extra):
    """A quantity inline item."""
    if kind not in KINDS:
        raise ValueError(f"unknown quantity kind {kind!r}")
    value = None if si is None else (int(si) if kind == "count" else float(si))
    return {"qty": {"kind": kind, "si": value, **extra}}


def ref(ref_text, name):
    """A link inline item."""
    return {"ref": ref_text, "name": name}


def strong(text):
    """An emphasis inline item."""
    return {"strong": text}


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
        if "qty" in node or "ref" in node and "name" in node or "strong" in node and len(node) == 1:
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
    for item in walk_inline(node):
        if isinstance(item, dict) and "qty" in item and item["qty"]["kind"] not in KINDS:
            raise ValueError(f"unknown quantity kind {item['qty']['kind']!r}")
