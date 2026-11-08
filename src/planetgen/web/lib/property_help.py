# planetgen/web/lib/property_help.py

"""
The rows of the method pages (`/references/method/<area>`, DOC.91; before that the properties page of UX.187): for each physics value of a star, planet, or moon page, what it is
called, the unit it is shown in, which bodies have it, and the one-line source the page shows under it. The rows are read from
the property sheet (`describe/sheet.py`) and its source table (`describe/sheet_sources.py`), so this page and the pages it
explains cannot disagree.
"""

import functools

from planetgen.db import column_units
from planetgen.describe import sheet
from planetgen.describe.sheet_sources import source_of
from planetgen.physics import unit_registry
from planetgen.web.lib import reference_urls

KINDS = (("planet", "Planets"), ("moon", "Moons"), ("star", "Stars"), ("neutron_star", "Neutron stars"))
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


GROUP_AREA = {
    "Identity": "bodies", "Orbit": "orbits", "Physical": "bodies", "Rotation": "orbits", "Composition": "bodies", "Atmosphere": "atmospheres",
    "Climate": "climate", "Water": "climate", "Magnetic field": "rings-magnetism", "Radiation": "radiation", "Habitability": "life", "Life": "life",
    "Zones and limits": "climate", "Stellar activity": "stars", "Galactic orbit": "galaxy", "Disk": "orbits", "Radiogenic heat": "radiation",
    "Gas and light": "appearance", "Compact object": "stars", "Activity": "stars", "Ratings": "life", "Dust": "atmospheres",
}
"""dict: The area of the method pages (`reference_urls.AREAS`) that holds the rows of each heading of the property sheet (DOC.91)."""


def by_area():
    """`{area: [rows]}` for every area of the method pages: the rows of `rows()` whose sheet heading `GROUP_AREA` puts there (a value only stars or neutron stars have is in `stars`), in the sheet's order."""
    found = {area: [] for area in reference_urls.AREAS}
    for row in rows():
        found["stars" if row["applies"] in ("Stars", "Neutron stars") else GROUP_AREA[row["group"]]].append(row)
    return found


@functools.lru_cache(maxsize=1)
def _refs():
    """`{(kind name, label): "area/label"}` for every row of `by_area()`; the sheet does not change while the program runs."""
    found = {}
    for area, rows in by_area().items():
        for row in rows:
            for name in row["applies"].split(", "):
                found.setdefault((name, row["label"]), f"{area}/{row['label']}")
    return found


def method_ref(kind, label):
    """
    The `reference_link("method", ...)` entry (`"area/Label"`) that explains the value `label` on the page of a `kind` ("planet", "moon" or
    "star"), or `None` when no method row has it: the link a body page shows where it used to print the source line (DOC.93).
    """
    return _refs().get((dict(KINDS).get(kind), label))
