"""
planetgen/web/formerly.py

The "Formerly" line of an object page (GEN.1140; docs/design/architecture/designation-history.md section 5): what the object used to be
called, newest first, each with why it changed and when. The lines come from the `former_designations` field of the object's JSON
(`galaxy/designations.former_json`); this module only words them.
"""

REASON_TEXT = {
    "collision bump": "renamed to clear a clash",
    "rename": "renamed",
    "regeneration": "rebuilt",
    "parent change": "moved to another parent",
    "admin renumber": "renumbered by an admin",
}
"""dict: The words of each reason a designation is retired (`galaxy/designations.REASONS`)."""


def view(former):
    """
    The rows of the "Formerly" line for the template: `[{designation, what, reason, epoch, date}]`, newest first; `[]` for an object that
    never changed (the line is then not shown).

    Args:
        former (list or None): The `former_designations` of the object's JSON.
    """
    return [{
        "designation": entry["designation"],
        "what": "Name" if entry["kind"] == "name" else "Object ID",
        "reason": REASON_TEXT.get(entry["reason"], entry["reason"]),
        "epoch": entry["retired_epoch"],
        "date": entry.get("retired_on"),
    } for entry in (former or [])]
