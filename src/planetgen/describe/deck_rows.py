# planetgen/describe/deck_rows.py

"""
The derived rows of the object sheet for a giant's cloud decks and trace-gas ceilings (GEN.1577)
==============================================================================================

`rows(conn, kind, object_id)` reads the stored `planets.cloud_decks` text (`physics/hot_cloud_deck.unpack`) and gives one row for each cloud that forms along the giant's profile (the pressure of its
base and the temperature there, or "deeper than" the search's floor when the profile stays cooler than the cloud's curve all the way down) and one for each trace gas (the most it can hold, in
molecules per hydrogen atom). They are worked out when the sheet is made from the stored text and are not columns of their own. Each row is `(group, key, label, value, raw)`; none for a body with no
stored decks (a rocky planet, a moon, a rogue planet, a giant of a row made before v269).
"""

from planetgen.describe import model
from planetgen.physics import hot_cloud_deck, hot_condensates

CLOUD_KEY_PREFIX = "cloud_base_"
"""str: The start of the key of each cloud's row (`cloud_base_cr`)."""

CEILING_KEY_PREFIX = "trace_ceiling_"
"""str: The start of the key of each trace gas's row (`trace_ceiling_geh4`)."""

PASCAL_PER_BAR = 1.0e5


def is_derived_key(key):
    """Whether `key` is one of the rows this module gives."""
    return key.startswith(CLOUD_KEY_PREFIX) or key.startswith(CEILING_KEY_PREFIX)


def _cloud_row(base):
    where = [model.qty("pressure", base.base_bar * PASCAL_PER_BAR, "Pa"), " at ", model.qty("temperature", base.base_k, "K")]
    value = ["deeper than ", *where] if base.reaches_bottom else where
    return ("atmosphere", CLOUD_KEY_PREFIX + base.species.lower(), f"Base of the {base.species} cloud", value, base.base_bar)


def _ceiling_row(formula, ceiling):
    gas = hot_condensates.HOT_TRACE_GASES[formula]
    value = [model.qty("ratio", ceiling, note="molecules per hydrogen atom at most")]
    return ("atmosphere", CEILING_KEY_PREFIX + formula.lower(), f"Most {gas.name} the gas can hold", value, ceiling)


def rows(conn, kind, object_id):
    """The derived cloud-deck and trace-gas rows of a planet, `(group, key, label, value, raw)` each (see the module text)."""
    if kind != "planet":
        return []
    found = conn.execute("SELECT cloud_decks FROM planets WHERE id = ?", (object_id,)).fetchone()
    record = hot_cloud_deck.unpack(found["cloud_decks"] if found is not None else None)
    if record is None:
        return []
    return ([_cloud_row(base) for base in record.clouds]
            + [_ceiling_row(formula, ceiling) for formula, ceiling in record.ceilings.items() if ceiling is not None and formula in hot_condensates.HOT_TRACE_GASES])
