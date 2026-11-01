# planetgen/web/lib/property_help.py

"""
The text of the properties page (`/help/properties`, UX.187): for each physics value of a star, planet or moon page, what it is
called, the unit it is shown in, which bodies have it, and the one-line source the page shows under it. The rows are read from
the property sheet (`describe/sheet.py`) and its source table (`describe/sheet_sources.py`), so this page and the pages it
explains cannot disagree.
"""

from planetgen.db import column_units
from planetgen.describe import sheet
from planetgen.describe.sheet_sources import source_of
from planetgen.physics import unit_registry
from planetgen.web.lib import help_registry

KINDS = (("planet", "Planets"), ("moon", "Moons"), ("star", "Stars"))
"""tuple: The kinds with source lines, in the order the page lists them, and what the page calls them."""

WORDS = "words"
"""str: What the unit column says for a row that is a name and not a number."""


def unit_text(kind, spec):
    """The unit a sheet row is shown in, as the page's unit column writes it: the registry's symbol, a local unit's noun or
    scale, or `WORDS` for a row that is text."""
    if spec.quantity in (sheet.TEXT, sheet.FLAG, sheet.TIER, sheet.RATING, sheet.LIMIT):
        return WORDS
    unit = column_units.COLUMNS[(sheet.TABLES[kind], spec.stored_column())]
    if unit in sheet._SYMBOLS:
        return sheet._SYMBOLS[unit].strip()
    if unit in unit_registry.UNITS:
        return unit_registry.UNITS[unit].symbol
    local = column_units.LOCAL_UNITS.get(unit)
    if local is not None and (local.suffix or local.prefix):
        text = (local.suffix or local.prefix).strip()
        return f"a number {text}" if text.startswith("(") else text
    return unit


def rows():
    """
    The page's rows, one per physics value with a source line, in the sheet's order: `group` (the sheet's heading for it),
    `label`, `unit`, `applies` (the kinds that have it, as the page names them) and `source`. A value two kinds share with
    the same line is one row.
    """
    found = {}
    for kind, name in KINDS:
        for spec in sheet.SPECS[kind]:
            line = source_of(kind, spec.key)
            if line is None:
                continue
            entry = found.setdefault((spec.group, spec.label, line), {
                "group": sheet.GROUPS[spec.group], "label": spec.label, "unit": unit_text(kind, spec), "applies": [], "source": line,
                "order": len(found)})
            entry["applies"].append(name)
    for key, (label, unit) in sheet.CHEMISTRY_ROWS.items():    # UX.254: the star's rows worked out from its stored values
        found[("composition", label, source_of("star", key))] = {
            "group": sheet.GROUPS["composition"], "label": label, "unit": unit, "applies": ["Stars"], "source": source_of("star", key),
            "order": len(found)}
    return [{**entry, "applies": ", ".join(entry["applies"])} for entry in sorted(found.values(), key=lambda e: e["order"])]


def page():
    """The template variables of the properties page: `rows`, grouped by the sheet's headings in `groups`
    (`[{"label", "rows"}]`, in the order of `describe.sheet.GROUPS`)."""
    listed = rows()
    groups = []
    for key, label in sheet.GROUPS.items():
        members = [row for row in listed if row["group"] == label]
        if members:
            groups.append({"label": label, "rows": members})
    return {"groups": groups}


help_registry.register_category("objects", "Objects and classes", "What the objects of the galaxy are, how they are sorted, and where their values come from.", order=10)
help_registry.register_topic("property_help_page", "Properties", "objects",
                             "Where the orbit, heat, spin, field, ring and life values on a star, planet or moon page come from.", order=30)
