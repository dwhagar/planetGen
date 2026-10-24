# planetgen/web/lib/habitability_cells.py

"""
The Habitability cell of a data table (UX.100): a planet's or moon's PHI-6 symbols as the `symbols` cell
`partials/datatable.html` and `static/datatable.js` draw. Imports no Flask.
"""

from planetgen.physics.habitability_world import DOMAINS, TIER_NAMES
from planetgen.web.lib.habitability_symbols import LABELS, OVERALL, overall_tier

COLUMN = "Habitability"
"""str: The column's heading."""

EQUIPMENT_COLUMN = "Equipment"
"""str: The heading of the equipment column (UX.154)."""


def habitability_cell(row):
    """
    The cell for one planet or moon `row` (`phi4` and `tier_<domain>` by name): the overall symbol and the six domain
    symbols, each with its rating, or a muted dash for a body that was not scored. `text` is the same in words.
    """
    overall = overall_tier(row)
    if overall is None:
        return {"text": "—", "muted": True}
    symbols = [{"name": OVERALL, "tier": overall}] + [{"name": name, "tier": row[f"tier_{name}"]} for name in DOMAINS]
    for symbol in symbols:
        symbol["label"] = f"{LABELS[symbol['name']]}: {TIER_NAMES[symbol['tier']]}"
    return {"text": ", ".join(symbol["label"] for symbol in symbols), "symbols": symbols}


def equipment_chips(row):
    """
    The equipment chips of a search `row` (UX.154): `[{text, severity, tier}]`, one per item of its needs list
    (`equipment_items`, required first). `severity` is `required` (drawn solid) or `recommended` (outlined) and `tier` the
    0 to 3 colour of the PHI-6 domain the item comes from. Empty for a body that needs nothing or was not scored.
    """
    chips = []
    for item in row.get("equipment_items") or ():
        tier = row.get(f"tier_{item['domain'].lower()}")
        chips.append({"text": item["item"], "severity": item["severity"], "tier": tier if tier is not None else 0})
    return chips


def equipment_cell(row):
    """The Equipment cell of a data table: the chips, with `text` the same list in words, or a muted dash."""
    chips = equipment_chips(row)
    if not chips:
        return {"text": "—", "muted": True}
    return {"text": row.get("equipment") or "—", "chips": chips}
