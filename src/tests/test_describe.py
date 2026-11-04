# tests/test_describe.py

"""
`planetgen.describe` (UX.104): the Document of a star system.

The Document must hold everything the old per-class text carried (every fact
key, the same overview sentences), be the same whether it is built from a
fresh system or one rebuilt from the database, be plain JSON, and carry no
Markdown or HTML.
"""

import json

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.describe import describe, model
from planetgen.describe import build
from planetgen.describe.quantity import display_value, inline_text
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation import stellar_object
from planetgen.generation.system import StarSystem
from planetgen.physics import habitability, protection, regional_temperature
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
            out.append(format_quantity(_LADDER_BY_KIND[item["qty"]["kind"]], display_value(item["qty"])))
        elif "strong" in item:
            out.append(item["strong"])
        elif "swatch" in item:
            out.append(item["swatch"]["word"])
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


@pytest.mark.parametrize("remnant_class, kind", [(BlackHole, "black-hole"), (NeutronStar, "neutron-star")])
def test_a_system_anchored_by_a_compact_remnant_is_described_as_one(remnant_class, kind):
    for seed in range(6):
        draw.set_run_seed(seed)
        config = SystemConfig()
        remnant = remnant_class(config)
        document = describe(StarSystem(system_config=config, compact_remnant=remnant))
        model.check(document)
        section = next(s for s in model.sections(document) if s["kind"] == kind)
        assert _TABLE_KEYS[kind] <= _fact_keys(section)
        text = " ".join(inline_text(sentence) for block in section["blocks"] if block["type"] == "prose"
                        for sentence in block["text"])
        assert "white dwarf" not in text and "main sequence" not in text
        overview = " ".join(inline_text(sentence) for sentence in document["sections"][0]["blocks"][0]["text"])
        assert "heliosphere" not in overview and "0 km" not in overview


def test_a_black_hole_says_how_it_formed_and_whether_it_feeds():
    for seed in range(12):
        draw.set_run_seed(seed)
        config = SystemConfig()
        hole = BlackHole(config)
        section = next(s for s in model.sections(describe(StarSystem(system_config=config, compact_remnant=hole)))
                       if s["kind"] == "black-hole")
        text = " ".join(inline_text(sentence) for block in section["blocks"] if block["type"] == "prose"
                        for sentence in block["text"])
        if hole.mass_class == "supermassive":
            assert "centre of the galaxy" in text
        else:
            assert ("collapsed core" in text) or ("intermediate-mass" in text)
            assert ("accretion disk of infalling matter" in text) == bool(hole.has_accretion_disk)


def test_a_neutron_star_names_its_pulsar_type():
    seen = set()
    for seed in range(60):
        draw.set_run_seed(seed)
        config = SystemConfig()
        config.POPULATION = "young"    # a young population's neutron stars still pulse (GEN.553); an old one's do not
        star = NeutronStar(config)
        section = next(s for s in model.sections(describe(StarSystem(system_config=config, compact_remnant=star)))
                       if s["kind"] == "neutron-star")
        text = " ".join(inline_text(sentence) for block in section["blocks"] if block["type"] == "prose"
                        for sentence in block["text"])
        expected = {"magnetar": "It is a magnetar", "non-pulsing": "no longer sweeps"}.get(
            star.pulsar_type, f"It is a {star.pulsar_type} pulsar")
        assert expected in text
        seen.add(star.pulsar_type)
    assert len(seen) >= 2


def test_no_markup_reaches_the_document():
    for _overrides, system in _systems(2):
        for item in model.walk_inline(describe(system)):
            strings = [item] if isinstance(item, str) else [value for value in item.values() if isinstance(value, str)]
            for text in strings:
                assert not any(mark in text for mark in _MARKUP), text


_TABLE_KEYS = {
    "star": {"type", "radius", "mass", "temp", "lum", "hab", "orbit", "loc"},
    "planet": {"class", "distance", "period", "speed", "radius", "gravity"},
    "close-pair": {"type", "mass", "lum", "hab", "separation", "mutual_orbit", "orbit", "loc"},
    "wide-pair": {"separation", "eccentricity", "mutual_orbit", "wobble", "primary_limit", "secondary_limit"},
    "black-hole": {"type", "mass", "event_horizon", "spin", "disk", "orbit", "loc"},
    "neutron-star": {"type", "mass", "radius", "spin_period", "magnetic_field", "surface_temp", "orbit", "loc"},
}
"""dict: The fact keys each kind of section carries (the rows of the data tables the pages have always shown)."""


def test_every_table_row_has_a_fact():
    """The facts of each section cover the rows the pages' data tables carry."""
    for _overrides, system in _systems():
        document = describe(system)
        by_title = {}
        for section in model.sections(document):
            by_title.setdefault((section["kind"], section["title"]), section)
        for star in system.stars:
            keys = _fact_keys(by_title[("star", star.name)])
            assert _TABLE_KEYS["star"] <= keys, star.name
        if system.binary_type == "close":
            pair = next(s for s in _all_sections(document) if s["kind"] == "pair")
            assert _TABLE_KEYS["close-pair"] <= _fact_keys(pair)
        if system.binary_type == "wide":
            pair = next(s for s in _all_sections(document) if s["kind"] == "pair")
            assert _TABLE_KEYS["wide-pair"] <= _fact_keys(pair)
        bodies = list(system.planets) + list(system.secondary_planets)
        for body in bodies:
            if body.body_type == "a":
                continue
            assert _TABLE_KEYS["planet"] <= _fact_keys(by_title[("planet", body.name)]), body.name
            for moon in body.moons:
                assert _TABLE_KEYS["planet"] <= _fact_keys(by_title[("moon", moon.name)]), moon.name


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


def test_the_overview_keeps_the_old_summary_but_counts_worlds_by_their_rating():
    """GEN.220: the heliosphere and edge sentences close the overview (UX.253 put them in plain words); the habitable line counts overall Green or Blue."""
    for _overrides, system in _systems():
        overview = describe(system)["sections"][0]["blocks"][0]
        assert overview["role"] == "overview"
        text = [_plain(sentence) for sentence in overview["text"]]
        assert any("holds on out to" in sentence for sentence in text[-2:])
        count = build._habitable_count(system)
        line = {0: "There are no habitable worlds in this system.",
                1: "There is 1 habitable world in the system."}.get(
            count, f"There are {count} habitable worlds in the system.")
        assert text.count(line) == 1 and not any("class M" in sentence for sentence in text)


def test_a_world_counts_as_habitable_by_its_overall_tier(monkeypatch):
    for _overrides, system in _systems():
        bodies = [body for planet in system.planets + system.secondary_planets if planet.body_type != "a"
                  for body in (planet, *planet.moons)]
        expected = sum(1 for body in bodies if body.habitability_world() is not None
                       and habitability.tier(protection.rating(*body.habitability_world())["overall"])
                       in ("Blue", "Green"))
        assert build._habitable_count(system) == expected
    monkeypatch.setattr(habitability, "POTENTIALLY_HABITABLE_TIERS", ())
    assert all(build._habitable_count(system) == 0 for _overrides, system in _systems())


def test_the_life_sentences_are_short_and_plain():
    """GEN.221: one sentence for the life class (and the pace), one each for what the life looks like."""
    seen = 0
    for _overrides, system in _systems(4):
        document = describe(system)
        for body in system.planets + system.secondary_planets:
            if body.body_type == "a" or not body.life_class:
                continue
            section = next(s for s in model.sections(document) if s["title"] == body.name)
            life = [_plain(sentence) for block in section["blocks"] if block["type"] == "prose" and block["role"] in ("life-class", "life")
                    for sentence in block["text"]]
            assert life and f"class {body.life_class}{body.life_class_variant or ''} life" in life[0].lower() and "could arise" not in life[0]
            assert all(len(sentence) < 160 and "suitable for the development" not in sentence for sentence in life)
            seen += 1
    assert seen > 0


def test_an_airless_body_has_a_zero_pressure_fact_and_no_atmosphere_block():
    found = False
    for _overrides, system in _systems(6):
        for section in model.sections(describe(system)):
            if section["kind"] not in ("planet", "moon"):
                continue
            rows = {row["key"]: row for block in section["blocks"] if block["type"] == "facts" for row in block["rows"]}
            if "pressure" not in rows:
                continue
            pressure = display_value(rows["pressure"]["value"][0]["qty"])
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
                regions = {field: getattr(part, field, None) for field in regional_temperature.FIELDS}
                result = protection.rating(*scored, regions=regions)
                assert [(row["domain"], row["tier"], row["score"]) for row in block["rows"]] == [
                    (row["domain"], row["tier"], row["score"]) for row in result["rows"]]
                assert [_plain(row["text"]) for row in block["rows"]] == [row["text"] for row in result["rows"]]
                assert block["overall"]["score"] == result["overall"]
                seen += 1
    assert seen


def test_the_star_sentence_is_two_short_sentences_and_makes_no_claim_about_the_universe():
    """GEN.221: the age and lifespan, then what the class means for life (one sentence each)."""
    seen = 0
    for _overrides, system in _systems(6):
        document = describe(system)
        for section in model.sections(document):
            if section["kind"] != "star":
                continue
            star = next(s for s in system.stars if s.name == section["title"])
            summary = [block for block in section["blocks"] if block["type"] == "prose" and block["role"] == "summary"]
            text = " ".join(_plain(sentence) for block in summary for sentence in block["text"])
            assert text.startswith("The star is about ") and "in the universe" not in text and "approximately" not in text
            if star.yerkes_class == "V" and stellar_object.supported_evolutionary_scales(star):
                assert build.LIFESPAN_SENTENCES[stellar_object.supported_evolutionary_scales(star)[-1]] in text
                assert text.count(". ") >= 1
            seen += 1
    assert seen > 0


def test_the_timeline_is_written_from_the_stored_state_with_the_ladders():
    seen = 0
    for _overrides, system in _systems(8):
        document = describe(system)
        for body in system.planets + system.secondary_planets:
            if body.body_type == "a" or not body.has_timeline:
                continue
            section = next(s for s in model.sections(document) if s["title"] == body.name)
            block = next((b for b in section["blocks"] if b["type"] == "prose" and b["role"] == "timeline"), None)
            if block is None:
                continue
            text = " ".join(_plain(sentence) for sentence in block["text"])
            assert "Billion Years" not in text and "Million Years" not in text and "speculative" not in text and "would" not in text
            seen += 1
    assert seen > 0
