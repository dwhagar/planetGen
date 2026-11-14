# tests/test_flavor_kinds.py

"""GEN.1433, GEN.1438: the flavor lines of the kinds of object that are not planets (a star, comet, black hole, neutron star, nebula, supernova remnant, quasar, interstellar comet and asteroid field).
Every kind has at least three lines for a plain object, every line of every kind can be drawn by some object whose attributes meet its gates, and an object draws one line on its own stream, once."""

import json
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.describe import build
from planetgen.generation import life
from planetgen.physics import constants, flavor_tree

KINDS = flavor_tree.MAPPED_SCOPES
PREFIX = {"star": "star", "comet": "comet", "black_hole": "blackhole", "neutron_star": "pulsar", "nebula": "nebula", "remnant": "remnant", "quasar": "quasar",
          "interstellar_comet": "icomet", "asteroid_field": "field"}


def _tree():
    return json.load(open(flavor_tree.DATA_FILE, encoding="utf-8"))["scopes"]


def _satisfying(scope, node_path):
    """An object whose attributes meet every gate on the way down to `node_path` (a list of nodes, root first)."""
    spec = _tree()[scope]["quantities"]
    wanted = {}
    for node in node_path:
        for gate in node.get("gate") or []:
            name = gate["quantity"]
            if "min" in gate:
                wanted[name] = gate["min"] * 1.01 if gate["min"] else 1.0
            elif "lt" in gate:
                wanted[name] = gate["lt"] * 0.5
            elif "equals" in gate:
                wanted[name] = gate["equals"]
            else:
                wanted[name] = gate["not_equals"] + 1
    attrs = {"name": "Test"}
    for name, entry in spec.items():
        value = wanted.get(name, 1.0)
        if entry.get("flag"):
            attrs[entry["attr"]] = bool(value)
        else:
            per = entry.get("per")
            attrs[entry["attr"]] = value * (per if isinstance(per, (int, float)) else getattr(constants, per)) if per else value
    return SimpleNamespace(flavor_key=None, **attrs)


def _paths(scope):
    paths = []

    def visit(node, path):
        path = [*path, node]
        paths.append(path)
        for child in node.get("children") or []:
            visit(child, path)
    visit(_tree()[scope], [])
    return paths


@pytest.mark.parametrize("scope", KINDS)
def test_a_plain_object_of_every_kind_is_offered_at_least_three_lines(scope):
    plain = _satisfying(scope, [_tree()[scope]])
    assert len(flavor_tree.walk(scope, plain)) >= 3


@pytest.mark.parametrize("scope", KINDS)
def test_every_line_of_every_kind_is_offered_to_some_object(scope):
    offered = set()
    for path in _paths(scope):
        offered |= {option.id for option in flavor_tree.walk(scope, _satisfying(scope, path))}
    in_file = {line["id"] for path in _paths(scope) for node in path for line in node.get("lines") or []}
    assert in_file and in_file <= offered


@pytest.mark.parametrize("scope", KINDS)
def test_every_line_of_a_kind_carries_the_prefix_of_its_kind(scope):
    for path in _paths(scope):
        for line in path[-1].get("lines") or []:
            assert line["id"].startswith(PREFIX[scope] + ":"), line["id"]


@pytest.mark.parametrize("scope", KINDS)
def test_an_object_draws_one_line_on_its_own_stream_and_never_twice(scope, monkeypatch):
    monkeypatch.setattr(tuning, "FLAVOR_CHANCE_PLANET", 1.0)
    first, second = _satisfying(scope, [_tree()[scope]]), _satisfying(scope, [_tree()[scope]])
    life.decide_flavor_text_alone(first, scope)
    life.decide_flavor_text_alone(second, scope)
    assert first.flavor_key and first.flavor_key.startswith(PREFIX[scope] + ":") and first.flavor_key == second.flavor_key     # the same object draws the same line
    kept = first.flavor_key
    life.decide_flavor_text_alone(first, scope)
    assert first.flavor_key == kept


@pytest.mark.parametrize("scope", KINDS)
def test_an_object_that_misses_the_chance_has_no_line(scope, monkeypatch):
    monkeypatch.setattr(tuning, "FLAVOR_CHANCE_PLANET", 0.0)
    body = _satisfying(scope, [_tree()[scope]])
    life.decide_flavor_text_alone(body, scope)
    assert body.flavor_key is None


def test_the_description_says_the_line_after_sensors_show_and_says_nothing_without_one():
    line = flavor_tree.text_of("nebula:0")
    blocks = build.flavor_blocks(SimpleNamespace(flavor_key="nebula:0"))
    assert len(blocks) == 1 and "Sensors show " + line in json.dumps(blocks)
    assert build.flavor_blocks(SimpleNamespace(flavor_key=None)) == []
    assert build.flavor_blocks(SimpleNamespace()) == []
