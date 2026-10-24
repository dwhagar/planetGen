# planetgen/web/lib/habitability_cells.py

"""
The Habitability cell of a data table (UX.100): a planet's or moon's PHI-6 symbols as the `symbols` cell
`partials/datatable.html` and `static/datatable.js` draw. Imports no Flask.
"""

from planetgen.physics.habitability_world import DOMAINS, TIER_NAMES
from planetgen.web.lib.habitability_symbols import LABELS, OVERALL, overall_tier

COLUMN = "Habitability"
"""str: The column's heading."""


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
