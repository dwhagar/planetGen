# tests/test_system_render.py

"""
Schema v29: a system's wikitext/Markdown page is no longer stored but
rendered from its database rows (`stellarObjects/systemRender.py`) -- see
`schema.sql`'s "v29" header note. These tests pin the property that made
dropping the stored copies safe: rendering a system from the database
gives exactly what rendering the freshly generated object gave (which is
what used to be stored), for every kind of system.
"""

import random

import pytest

import checkRenderParity
from stellarObjects import _db
from stellarObjects.compactRemnant import BlackHole, NeutronStar
from stellarObjects.config import SystemConfig
from stellarObjects.doubleStar import BinaryStarProxy
from stellarObjects.evolution import life_stage_from_paragraphs
from stellarObjects.starData import Star
from stellarObjects.systemData import StarSystem
from stellarObjects.systemRender import render_star_system, render_system_sections, render_system_text
from stellarObjects.wideBinary import WideBinaryPair

_VARIANTS = [
    dict(BINARY_SYSTEM=False, MOONS=True, COMETS=True, ASTEROID_BELT=True),
    dict(BINARY_SYSTEM=True, WIDE_BINARY=False, MOONS=True),
    dict(BINARY_SYSTEM=True, WIDE_BINARY=True, COMETS=True),
    dict(BINARY_SYSTEM=False, HABITABLE_WORLD=True, INTELLIGENT_LIFE=True),
    dict(BINARY_SYSTEM=False, STAR_TYPE="M2VII"),
    dict(BINARY_SYSTEM=False, STAR_TYPE="K1III", AGE="old"),
    dict(BINARY_SYSTEM=False, PLANETS=False),
]


def _config(**overrides):
    cfg = SystemConfig()
    for key, value in overrides.items():
        setattr(cfg, key, value)
    return cfg


@pytest.mark.parametrize("overrides", _VARIANTS)
def test_rendering_from_the_database_matches_the_generated_system(mysql_config, overrides):
    random.seed(hash(tuple(sorted(overrides.items()))) & 0xFFFF)
    conn = _db.get_connection(mysql_config)
    try:
        for _ in range(4):
            cfg = _config(**overrides)
            system = StarSystem(system_config=cfg)
            with conn:
                system_id = _db.insert_star_system(conn, system, cfg)
            for fmt in ("wikitext", "markdown"):
                assert render_system_text(conn, system_id, fmt) == render_star_system(system, fmt), fmt
    finally:
        conn.close()


@pytest.mark.parametrize("remnant_class", [BlackHole, NeutronStar])
def test_rendering_a_compact_remnant_anchored_system(mysql_config, remnant_class):
    cfg = SystemConfig()
    system = StarSystem(system_config=cfg, compact_remnant=remnant_class(cfg))
    system_id = _db.save_phenomenon(
        system, cfg, "black-hole" if remnant_class is BlackHole else "neutron-star", config=mysql_config,
    )
    conn = _db.get_connection(mysql_config)
    try:
        for fmt in ("wikitext", "markdown"):
            assert render_system_text(conn, system_id, fmt) == render_star_system(system, fmt)
    finally:
        conn.close()


def _swapped_pair(mysql_config, wide):
    """A binary whose generation-order secondary is the heavier star --
    the case where the pair's own (heavier-first) orientation differs
    from the `stars.role` rows."""
    cfg = _config(BINARY_SYSTEM=True, WIDE_BINARY=wide, PLANETS=False)
    system = StarSystem(system_config=cfg)
    primary, secondary = system.primary_star, system.secondary_star
    heavy, light = max(primary.mass, secondary.mass), min(primary.mass, secondary.mass)
    primary.mass, secondary.mass = light, heavy
    if wide:
        system.wide_binary.primary, system.wide_binary.secondary = secondary, primary
    else:
        system.star._primary, system.star._secondary = secondary, primary
    return system, cfg


@pytest.mark.parametrize("wide", [False, True])
def test_binary_orientation_survives_the_round_trip(mysql_config, wide):
    system, cfg = _swapped_pair(mysql_config, wide)
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        for fmt in ("wikitext", "markdown"):
            assert render_system_text(conn, system_id, fmt) == render_star_system(system, fmt)
    finally:
        conn.close()


def test_from_dict_puts_the_heavier_star_first():
    cfg = SystemConfig()
    light, heavy = Star(cfg), Star(cfg)
    light.mass, heavy.mass = 1.0e30, 2.0e30

    pair = WideBinaryPair.from_dict({}, cfg, light, heavy)
    assert pair.primary is heavy and pair.secondary is light

    proxy_cfg = _config(BINARY_SYSTEM=True, WIDE_BINARY=False, PLANETS=False)
    proxy = StarSystem(system_config=proxy_cfg).star
    data = proxy.to_dict()
    data["primary"], data["secondary"] = data["secondary"], data["primary"]
    reloaded = BinaryStarProxy.from_dict(data, proxy_cfg)
    assert reloaded._primary.mass >= reloaded._secondary.mass


def test_rendering_follows_a_rename(mysql_config):
    cfg = _config(BINARY_SYSTEM=False, PLANETS=False)
    system_id = _db.save_system(StarSystem(system_config=cfg), cfg, config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("UPDATE star_systems SET name = 'Somewhere Else' WHERE id = ?", (system_id,))
        assert render_system_text(conn, system_id, "wikitext").startswith("= Somewhere Else =")
        with pytest.raises(ValueError):
            render_system_text(conn, system_id, "html")
    finally:
        conn.close()


def test_sections_split_the_page_per_body(mysql_config):
    cfg = _config(BINARY_SYSTEM=False, MOONS=True, COMETS=True, ASTEROID_BELT=True)
    system = StarSystem(system_config=cfg)
    system_id = _db.save_system(system, cfg, config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        sections = render_system_sections(conn, system_id)
        markdown = render_system_text(conn, system_id, "markdown")
    finally:
        conn.close()

    # Every piece is a verbatim excerpt of the full page.
    for key in ("stars", "planets", "moons", "belts", "comets"):
        for text in sections[key].values():
            assert text and text in markdown, key
    assert sections["overview"] in markdown
    planets = [p for p in system.planets if p.body_type != "a"]
    assert len(sections["planets"]) == len(planets)
    assert len(sections["moons"]) == sum(len(p.moons) for p in planets)


def test_migrate_v28_to_v29_drops_the_stored_page_text(mysql_config):
    cfg = _config(BINARY_SYSTEM=False, PLANETS=False)
    system_id = _db.save_system(StarSystem(system_config=cfg), cfg, config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        conn.execute("ALTER TABLE star_systems ADD COLUMN wikitext_content LONGTEXT, ADD COLUMN markdown_content LONGTEXT")
        conn.execute("UPDATE star_systems SET wikitext_content = 'old', markdown_content = 'old'")
        conn.execute("DELETE FROM schema_migrations WHERE version IN (29, 30, 31)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (28)")
        conn.commit()
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        columns = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM star_systems").fetchall()}
        assert not columns & set(_db.V29_DROPPED_COLUMNS)
        assert render_system_text(conn, system_id, "markdown").startswith("# ")
    finally:
        conn.close()

    # Already current: a no-op.
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION


def test_life_stage_is_read_back_from_the_timeline_paragraph():
    paragraph = (
        "A speculative evolutionary timeline ... would have been Technological Civilization at "
        "4.50 Billion Years. The planet boasts a fully mature biosphere."
    )
    assert life_stage_from_paragraphs([paragraph]) == "technological_civilization"
    assert life_stage_from_paragraphs(["No significant evolutionary milestones are predicted."]) is None
    assert life_stage_from_paragraphs([]) is None


def test_parity_check_sorts_differences():
    names = ["Onoshur Kin", "Onoshur"]
    stored = "= Onoshur =\n|wobble=1 km from its nominal position\nSame line"
    assert checkRenderParity.classify(stored, stored, names) is None
    assert checkRenderParity.classify(stored, stored.replace("= Onoshur =", "= Onoshur Kin ="), names) == "names"
    assert checkRenderParity.classify(stored, stored.replace("1 km", "2 km"), names) == "live"
    assert checkRenderParity.classify(stored, stored.replace("Same line", "Changed line"), names) == "other"
