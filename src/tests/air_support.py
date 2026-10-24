"""The air of a test body, written the way a test reads best (kPa by species) and stored the way the code reads it."""

from planetgen.physics import atmosphere, gases


def air(**kpa):
    """`{atmosphere.GAS_FIELD: {formula: Pa}}` of `p_<gas>_kpa=value` keywords (a species at 0 is left out)."""
    return {atmosphere.GAS_FIELD: {gases.stored(key[2:-4]).formula: value * 1000.0 for key, value in kpa.items() if value}}


def fold_air(values):
    """`values` (a dict) with its `p_<gas>_kpa` keys moved into the `atmosphere.GAS_FIELD` dict; returns `values`."""
    found = {key: values.pop(key) for key in [key for key in values if key.startswith("p_") and key.endswith("_kpa")]}
    if found:
        values.update(air(**found))
    return values
