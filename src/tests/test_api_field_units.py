# tests/test_api_field_units.py

"""API.33: every number an API request body takes states its unit, in its name (`distance_km`, `radius_ly`) or in its
description ("Unit: count of ..."), over the JSON schema of every body model (`planetgen.api.schemas`)."""

import inspect

from pydantic import BaseModel

from planetgen.api import schemas
from planetgen.db import column_units

NAME_UNITS = tuple(suffix for suffix, _unit in column_units.SUFFIXES) + ("_days",)
KINDS = ("count", "id", "index", "row-defined", "ratio", "flag")
"""tuple: What a number with no physical unit says it is (the column registry's markers, DB.76)."""


def _models():
    for name, model in inspect.getmembers(schemas, inspect.isclass):
        if issubclass(model, BaseModel) and model.__module__ == schemas.__name__:
            yield name, model


def _numeric(node, definitions):
    """Whether a schema node is a number, or a list or mapping of numbers, or one of several types holding one."""
    if "$ref" in node:
        return _numeric(definitions[node["$ref"].split("/")[-1]], definitions)
    if node.get("type") in ("number", "integer"):
        return True
    for key in ("anyOf", "oneOf", "allOf"):
        if any(_numeric(option, definitions) for option in node.get(key, ())):
            return True
    if node.get("type") == "array":
        return _numeric(node.get("items", {}), definitions)
    if node.get("type") == "object" and isinstance(node.get("additionalProperties"), dict):
        return _numeric(node["additionalProperties"], definitions)
    return False


def _stated(name, description):
    if name.endswith(NAME_UNITS):
        return True
    text = description or ""
    if not text.startswith("Unit: "):
        return False
    return len(text) > len("Unit: ") + 1


def test_every_number_in_a_request_body_states_its_unit():
    unstated = []
    for name, model in _models():
        schema = model.model_json_schema()
        definitions = schema.get("$defs", {})
        for field, node in schema.get("properties", {}).items():
            if _numeric(node, definitions) and not _stated(field, node.get("description")):
                unstated.append(f"{name}.{field}")
    assert unstated == [], "a request field that is a number but names no unit: " + ", ".join(unstated)


def test_a_stated_unit_is_a_unit_or_a_kind():
    from planetgen.physics import unit_registry

    for name, model in _models():
        for field, info in model.model_fields.items():
            description = info.description or ""
            if not description.startswith("Unit: "):
                continue
            word = description[len("Unit: "):].replace(",", " ").split()[0]
            known = word in KINDS or word in unit_registry.UNITS or word in {suffix[1:] for suffix, _u in column_units.SUFFIXES} or word in ("days", "ly")
            assert known, f"{name}.{field}: {description!r} starts with {word!r}, neither a unit nor one of {KINDS}"
