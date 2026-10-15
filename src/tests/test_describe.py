# tests/test_describe.py

"""
`planetgen.describe` (UX.104): the Document of a star system.

The Document must hold everything the old per-class text carried (every fact
key, the same overview sentences), be the same whether it is built from a
fresh system or one rebuilt from the database, be plain JSON, and carry no
Markdown, wikitext or HTML.
"""

import json

import pytest

from planetgen.db import store
from planetgen.describe import describe, model
from planetgen.describe import build
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.system import StarSystem
from planetgen.physics import habitability, protection
from planetgen.util import draw, ladders
from planetgen.util.format import format_quantity

_VARIANTS = [
    dict(BINARY_SYSTEM=False, MOONS=True, COMETS=True, ASTEROID_BELT=True),
    dict(BINARY_SYSTEM=True, WIDE_BINARY=False, MOONS=True, COMETS=True),
    dict(BINARY_SYSTEM=True, WIDE_BINARY=True, COMETS=True, MOONS=True),
    dict(BINARY_SYSTEM=False, HABITABLE_WORLD=True, INTELLIGENT_LIFE=True),
    dict(BINARY_SYSTEM=False, STAR_TYPE="M2VII"),
    dict(BINARY_SYSTEM=False, STAR_TYPE="K1III", AGE="old"),
    dict(BINARY_SYSTEM=False, PLANETS=False),
]

_LADDER_BY_KIND = {ladder.quantity: ladder for ladder in ladders.LADDERS}
_MARKUP = ("**", "<sup>", "</", "{{", "}}", "'''", "[[", "]]", "|---", "&amp;")


def _config(**overrides):
    config = SystemConfig()
    for key, value in overrides.items():
        setattr(config, key, value)
    return config


def _systems(count=3):
    for overrides in _VARIANTS:
        draw.set_run_seed(hash(tuple(sorted(overrides.items()))) & 0xFFFF)
        for _ in range(count):
            yield overrides, StarSystem(system_config=_config(**overrides))


def _plain(inline):
    """An inline list as the ladders would write it (only for the kinds that have a ladder)."""
    out = []
    for item in inline:
        if isinstance(item, str):
            out.append(item)
        elif "qty" in item:
            out.append(format_quantity(_LADDER_BY_KIND[item["qty"]["kind"]], item["qty"]["si"]))
        elif "strong" in item:
            out.append(item["strong"])
        else:
            out.append(item["name"])
    return "".join(out)


def _fact_keys(section):
    return {row["key"] for block in section["blocks"] if block["type"] == "facts" for row in block["rows"]}


def _all_sections(document):
    return list(model.sections(document))


@pytest.mark.parametrize("overrides", _VARIANTS)
def test_every_variant_gives_a_checked_json_document(overrides):
    draw.set_run_seed(hash(tuple(sorted(overrides.items()))) & 0xFFFF)
    for _ in range(4):
        document = describe(StarSystem(system_config=_config(**overrides)))
        model.check(document)
        assert json.loads(json.dumps(document, allow_nan=False)) == document


def test_a_system_anchored_by_a_compact_remnant_is_described():
    for remnant_class in (BlackHole, NeutronStar):
        config = SystemConfig()
        document = describe(StarSystem(system_config=config, compact_remnant=remnant_class(config)))
        model.check(document)
        assert any(section["kind"] == "star" for section in model.sections(document))


def test_no_markup_reaches_the_document():
    for _overrides, system in _systems(2):
        for item in model.walk_inline(describe(system)):
            strings = [item] if isinstance(item, str) else [value for value in item.values() if isinstance(value, str)]
            for text in strings:
                assert not any(mark in text for mark in _MARKUP), text


def test_every_old_table_row_has_a_fact():
    """The facts of each section cover the keys the old data table carried."""
    for _overrides, system in _systems():
        document = describe(system)
        by_title = {}
        for section in model.sections(document):
            by_title.setdefault((section["kind"], section["title"]), section)
        for star in system.stars:
            keys = _fact_keys(by_title[("star", star.name)])
            assert set(star.get_table_properties()) <= keys, star.name
        if system.binary_type == "close":
            pair = next(s for s in _all_sections(document) if s["kind"] == "pair")
            assert set(system.star.get_table_properties()) <= _fact_keys(pair)
        if system.binary_type == "wide":
            pair = next(s for s in _all_sections(document) if s["kind"] == "pair")
            assert set(system.wide_binary.get_table_properties()) <= _fact_keys(pair)
        bodies = list(system.planets) + list(system.secondary_planets)
        for body in bodies:
            if body.body_type == "a":
                continue
            assert set(body.get_table_properties()) <= _fact_keys(by_title[("planet", body.name)]), body.name
            for moon in body.moons:
                assert set(moon.get_table_properties()) <= _fact_keys(by_title[("moon", moon.name)]), moon.name


def test_every_body_has_a_section_in_page_order():
    for _overrides, system in _systems():
        document = describe(system)
        titles = [(s["kind"], s["title"]) for s in model.sections(document)]
        expected = [body.name for body in system.planets + system.secondary_planets if body.body_type != "a"]
        planets = [title for kind, title in titles if kind == "planet"]
        assert sorted(planets) == sorted(expected)
        moons = sorted(moon.name for body in system.planets + system.secondary_planets if body.body_type != "a"
                       for moon in body.moons)
        assert sorted(title for kind, title in titles if kind == "moon") == moons
        belts = [body for body in system.planets + system.secondary_planets if body.body_type == "a"]
        assert sum(kind == "belt" for kind, _ in titles) == len(belts)
        assert sum(kind == "comet" for kind, _ in titles) == len(system.comets) + len(system.secondary_comets)


def test_the_overview_says_what_the_old_summary_said():
    for _overrides, system in _systems():
        overview = describe(system)["sections"][0]["blocks"][0]
        assert overview["role"] == "overview"
        assert " ".join(_plain(sentence) for sentence in overview["text"]) == system.summary_paragraph()


def test_the_life_and_flavour_sentences_are_the_old_ones():
    seen = set()
    for _overrides, system in _systems(4):
        document = describe(system)
        for body in system.planets + system.secondary_planets:
            if body.body_type == "a":
                continue
            old = " ".join(body.to_paragraph_list())
            section = next(s for s in model.sections(document) if s["title"] == body.name)
            for block in section["blocks"]:
                if block["type"] == "prose" and block["role"] in ("life", "flavour", "timeline"):
                    for sentence in block["text"]:
                        assert _plain(sentence) in old
                        seen.add(block["role"])
    assert {"life", "timeline"} <= seen


def test_an_airless_body_has_a_zero_pressure_fact_and_no_atmosphere_block():
    found = False
    for _overrides, system in _systems(6):
        for section in model.sections(describe(system)):
            if section["kind"] not in ("planet", "moon"):
                continue
            rows = {row["key"]: row for block in section["blocks"] if block["type"] == "facts" for row in block["rows"]}
            if "pressure" not in rows:
                continue
            pressure = rows["pressure"]["value"][0]["qty"]["si"]
            has_air = any(block["type"] == "composition" and block["label"] == "Atmosphere" for block in section["blocks"])
            assert (pressure > 0) == has_air
            found |= pressure == 0
    assert found


def test_a_gas_giant_carries_its_depth_sentences_with_quantities():
    for _overrides, system in _systems(4):
        for section in model.sections(describe(system)):
            blocks = [b for b in section["blocks"] if b["type"] == "prose" and b["role"] == "depth"]
            if blocks:
                kinds = {item["qty"]["kind"] for item in model.walk_inline(blocks) if isinstance(item, dict) and "qty" in item}
                assert kinds == {"length", "temperature", "pressure"}
                return
    pytest.fail("no gas giant generated")


@pytest.mark.parametrize("overrides", _VARIANTS[:4])
def test_a_stored_system_is_described_like_the_generated_one(mysql_config, overrides):
    """What the page shows comes from rows; it must say what the fresh object says (ids aside)."""
    def normal(node):
        """Without object references (not yet known to a fresh system) and with floats rounded to the
        precision the database column keeps."""
        if isinstance(node, dict):
            return {key: normal(value) for key, value in node.items() if key != "ref"}
        if isinstance(node, list):
            return [normal(item) for item in node]
        if isinstance(node, float):
            return float(f"{node:.9g}")
        return node

    draw.set_run_seed(hash(tuple(sorted(overrides.items()))) & 0xFFFF)
    conn = store.get_connection(mysql_config)
    try:
        for _ in range(3):
            config = _config(**overrides)
            system = StarSystem(system_config=config)
            with conn:
                system_id = store.insert_star_system(conn, system, config)
            loaded = store.load_star_system(conn, system_id)
            assert normal(describe(loaded)) == normal(describe(system))
    finally:
        conn.close()


def test_a_scored_rocky_body_carries_its_six_ratings_as_the_paragraph_says_them():
    seen = 0
    for _overrides, system in _systems(4):
        document = describe(system)
        for body in system.planets + system.secondary_planets:
            if body.body_type == "a":
                continue
            for part in [body, *body.moons]:
                section = next(s for s in model.sections(document) if s["title"] == part.name)
                blocks = [b for b in section["blocks"] if b["type"] == "ratings"]
                scored = part.habitability_world()
                assert bool(blocks) == (scored is not None)
                if not blocks:
                    continue
                block = blocks[0]
                assert [row["domain"] for row in block["rows"]] == list(habitability.DOMAINS)
                bold = lambda text: text
                paragraph = protection.paragraph(*scored[:1], bold, scored[1])
                for row in block["rows"]:
                    assert f"{row['domain']}: {row['tier']}" in paragraph
                    assert _plain(row["text"]) in paragraph
                assert f"{block['overall']['score']:.2f}" in paragraph
                seen += 1
    assert seen
