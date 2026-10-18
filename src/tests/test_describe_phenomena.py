# tests/test_describe_phenomena.py

"""
UX.108: the Document of each phenomenon (`planetgen.describe.phenomena`): a nebula, a quasar, a rogue planet, an
interstellar comet, an asteroid field, a supernova remnant, and a black hole or neutron star on its own. Each one
builds a valid Document, says what the old per-class paragraphs said, and renders in every format.
"""

import json
import re

import pytest

from planetgen import describe
from planetgen.describe import html, model
from planetgen.describe.quantity import inline_text
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.quasar import Quasar
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.util import draw

KINDS = {
    Nebula: "nebula", Quasar: "quasar", RoguePlanet: "rogue-planet", InterstellarComet: "interstellar-comet",
    AsteroidField: "asteroid-field", SupernovaRemnant: "supernova-remnant", BlackHole: "black-hole",
    NeutronStar: "neutron-star",
}


def _make(cls, seed):
    draw.set_run_seed(seed)
    return cls(SystemConfig())


def _text(section):
    return " ".join(inline_text(sentence) for block in section["blocks"] if block["type"] == "prose"
                    for sentence in block["text"])


@pytest.mark.parametrize("cls", list(KINDS))
def test_every_phenomenon_builds_a_valid_json_document(cls):
    for seed in range(8):
        thing = _make(cls, seed)
        document = describe.describe_phenomenon(thing)
        model.check(document)
        assert json.loads(json.dumps(document, allow_nan=False)) == document
        assert document["kind"] == "phenomenon" and document["title"] == thing.name
        assert document["sections"][0]["kind"] == KINDS[cls]


@pytest.mark.parametrize("cls", list(KINDS))
def test_every_phenomenon_renders_in_every_format(cls):
    document = describe.describe_phenomenon(_make(cls, 3))
    for fmt in describe.RENDERERS:
        out = describe.render(document, fmt)
        assert document["title"] in out
    assert "[[Category:Phenomena]]" in describe.render(document, "wikitext")
    assert "[[Category:Star Systems]]" not in describe.render(document, "wikitext")
    assert document["title"] in html.render(document)


def test_a_phenomenon_names_the_things_its_old_paragraph_named():
    nebula = _make(Nebula, 1)
    text = _text(describe.describe_phenomenon(nebula)["sections"][0])
    assert f"{nebula.nebula_type} nebula" in text and nebula.formation_cause in text and "galactic centre" in text

    quasar = _make(Quasar, 1)
    text = _text(describe.describe_phenomenon(quasar)["sections"][0])
    assert "active nucleus" in text and ("radio-quiet" in text) == (not quasar.is_radio_loud)

    comet = _make(InterstellarComet, 1)
    text = _text(describe.describe_phenomenon(comet)["sections"][0])
    assert "hyperbolic" in text and ("active coma" in text) == bool(comet.is_active)

    field = _make(AsteroidField, 1)
    section = describe.describe_phenomenon(field)["sections"][0]
    assert any(block["type"] == "composition" and block["parts"] for block in section["blocks"])


def test_a_rogue_planet_keeps_every_surface_sentence():
    for seed in range(20):
        planet = _make(RoguePlanet, seed)
        text = _text(describe.describe_phenomenon(planet)["sections"][0])
        assert "free-floating" in text and "only heat is its own" in text
        assert ("companion body" in text) == bool(planet.has_moons)
        regime = planet.surface_regime
        if regime in ("gas-giant", "brown-dwarf"):
            assert "no solid surface" in text
        elif regime == "bare-rock":
            assert "bare rock surface" in text


def test_a_supernova_remnant_holds_its_compact_remnant_as_a_child():
    kinds = set()
    for seed in range(40):
        remnant = _make(SupernovaRemnant, seed)
        top = describe.describe_phenomenon(remnant)["sections"][0]
        if remnant.compact_remnant is None:
            assert not top["children"]
        else:
            assert [child["kind"] for child in top["children"]] == [
                "black-hole" if isinstance(remnant.compact_remnant, BlackHole) else "neutron-star"]
            assert remnant.compact_remnant.name in _text(top)
        kinds.add(remnant.compact_remnant is None)
    assert kinds == {True, False}


def test_nothing_a_phenomenon_holds_comes_out_unescaped():
    hostile = '<script>alert("x")</script> & "q"'
    for cls in KINDS:
        thing = _make(cls, 2)
        thing.name = hostile
        out = html.render(describe.describe_phenomenon(thing))
        assert "<script" not in out and "&lt;script&gt;" in out
        assert not re.search(r'<(?!/?(?:section|h\d|div|table|caption|thead|tbody|tr|th|td|p|span|strong|a|sup)\b)', out)


def test_something_that_is_not_a_phenomenon_is_refused():
    with pytest.raises(ValueError):
        describe.describe_phenomenon(object())
