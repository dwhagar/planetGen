# stellarObjects/serialization.py

"""
Shared Serialization Helpers
=============================

Small, generic helpers for the allowlist-based `to_dict`/`from_dict` pattern
used across the package's classes (`SystemConfig`, `Star`, `BinaryStarProxy`,
`Planet`, `AsteroidBelt`, `StarSystem`): each class owns an explicit
`SERIALIZABLE_FIELDS` list naming exactly the attributes that make up its
persisted state, and `to_dict`/`from_dict` read/write only those names.

An allowlist (rather than a blanket `vars(obj)`-minus-exclude-list) fails
*safe*: a future new attribute is simply not serialized until someone adds
it to the list, instead of leaking into storage by default. Reconstruction
uses `object.__new__(cls)` plus attribute assignment, bypassing `__init__`
entirely -- `__init__` re-runs full random generation, not just assignment,
so it can't be used to rebuild an object from already-decided data.

Dict keys are the lowercased field name. This matches `SystemConfig`'s
existing `--system-file` JSON convention (its fields are `UPPER_CASE`
attributes) while being a no-op for every other class here, whose
attributes are already lowercase `snake_case`.
"""


def fields_to_dict(obj, fields):
    """
    Builds a plain dict of `obj`'s current values for each name in `fields`.

    Args:
        obj: The object to read attributes from.
        fields (list): Attribute names to include, e.g. a class's
                       `SERIALIZABLE_FIELDS` constant.

    Returns:
        dict: `{field.lower(): getattr(obj, field) for field in fields}`.
    """
    return {field.lower(): getattr(obj, field) for field in fields}


def fields_from_dict(obj, data, fields):
    """
    Assigns values from `data` onto `obj` for each name in `fields`, via
    plain `setattr` (bypassing any property setter -- callers with
    read-only properties backed by a private attribute must list the
    private attribute's own name in `fields`, not the property's).

    A field whose lowercased name is absent from `data` is left untouched,
    matching `SystemConfig.from_dict`'s existing "partial recipe" behavior
    (a hand-authored `--system-file` document may omit keys deliberately
    and rely on the object's own defaults) rather than overwriting it with
    a default `None`.

    Args:
        obj: The object to assign attributes onto.
        data (dict): A dict in the shape `fields_to_dict` produces.
        fields (list): Attribute names to assign, e.g. a class's
                       `SERIALIZABLE_FIELDS` constant.
    """
    for field in fields:
        key = field.lower()
        if key in data:
            setattr(obj, field, data[key])
