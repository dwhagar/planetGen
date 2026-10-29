# planetgen/web/lib/planet_layout.py

"""
The parts of a planet's or moon's page that the redesign (UX.215 to UX.222) takes out of the description and the property
sheet and sets in a new order: the written description's blocks sorted into the section each belongs to, the table of the
air by mole fraction (UX.219), the table of the body's make-up by mass (UX.219) and the civilization paragraph (UX.221).
Imports no Flask.
"""

from html import escape

from planetgen.describe.html import inline_html
from planetgen.physics import gases
from planetgen.util.format import three_figures_or_tiny

PHI6_CAPTIONS = ("Radiation", "Survival equipment")
"""tuple: The captions of the description's facts tables that belong with the PHI-6 rating."""

MASS_KEYS = (("metal_mass_fraction", "Metal (iron and nickel core)"), ("rock_mass_fraction", "Rock (silicate mantle and crust)"),
             ("ice_mass_fraction", "Ice"), ("envelope_mass_fraction", "Gas envelope"))
"""tuple: The stored mass fractions of a body and what each is, in the order the table lists them."""

MASS_SHARED_KEYS = frozenset(key for key, _label in MASS_KEYS)
"""frozenset: The sheet rows the mass table shows, which the Properties list leaves out."""


def split_blocks(blocks):
    """
    The description's blocks of a planet or moon sorted by where the page puts them: `phi6` (the ratings table and the
    radiation and survival equipment tables), `life` (the life and timeline paragraphs) and `rest` (anything else the page
    still shows, such as the sensors' flavour line and a gas giant's depth). The facts table of the body, the class
    sentence and the atmosphere and composition summaries are dropped: the page builds those from the stored columns.
    """
    parts = {"phi6": [], "life": [], "rest": []}
    for block in blocks:
        kind = block["type"]
        if kind == "ratings" or (kind == "facts" and block["caption"] in PHI6_CAPTIONS):
            parts["phi6"].append(block)
        elif kind == "prose" and block["role"] in ("life", "timeline"):
            parts["life"].append(block)
        elif kind == "prose" and block["role"] in ("flavour", "depth"):
            parts["rest"].append(block)
    return parts


def _rows(sheet, group_key):
    return next((group["rows"] for group in sheet["groups"] if group["key"] == group_key), [])


def _number(row):
    if row.get("raw") is not None:
        return row["raw"]
    item = row["value"][0] if row["value"] else None
    return item["qty"]["value"] if isinstance(item, dict) and "qty" in item else None


def gas_table(sheet):
    """
    The air of a body as `[{"name", "fraction", "pressure"}]`, the most abundant gas first: the share of the molecules each
    gas makes (its mole fraction, which is also its share of the volume and of the pressure) as text, and its partial
    pressure as the page's own HTML. Empty for a body with no stored gases.
    """
    held = [(row["key"][2:], row) for row in _rows(sheet, "atmosphere") if row["key"].startswith("p_") and _number(row)]
    total = sum(_number(row) for _gas, row in held)
    table = []
    for formula, row in held:
        share = _number(row) / total * 100.0
        table.append({"name": gases.gas(formula).name.capitalize(),
                      "fraction": "less than 0.001 percent" if share < 0.001 else f"{three_figures_or_tiny(share)} percent",
                      "pressure": inline_html(row["value"])})
    return table


def mass_table(sheet):
    """The make-up of a body by mass as `[{"name", "share"}]`: each stored mass fraction above zero as a percent of its mass
    (the four parts add up to the whole body). Empty when none is stored."""
    rows = {row["key"]: row for row in _rows(sheet, "composition")}
    table = []
    for key, label in MASS_KEYS:
        row = rows.get(key)
        fraction = None if row is None else _number(row)
        if fraction:
            table.append({"name": label, "share": f"{three_figures_or_tiny(fraction * 100.0)} percent of its mass"})
    return table


def civilization_sentences(row, species):
    """The paragraph on the people of a world, as HTML sentences: its dominant species and whether a technological
    civilization or a colony stands on it. Empty for a world with neither, so the page leaves the section out."""
    sentences = []
    if species:
        sentences.append(f'Dominant species: <a href="{escape(species["url"], quote=True)}">{escape(species["name"])}</a>')
    if row is not None and row.get("inhabited"):
        sentences.append("A technological civilization or a colony lives here.")
    return sentences
