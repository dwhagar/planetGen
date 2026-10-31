# tests/test_system_render.py

"""
Schema v29: a system's Markdown page is no longer stored but
rendered from its database rows (`planetgen/db/render.py`) -- see
`schema.sql`'s "v29" header note. These tests pin the property that made
dropping the stored copies safe: rendering a system from the database
gives exactly what rendering the freshly generated object gave (which is
what used to be stored), for every kind of system.
"""


import pytest

from planetgen.util import draw
from planetgen.db import store
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.config import SystemConfig
from planetgen.generation.binary import BinaryStarProxy
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.db.render import describe_system, render_star_system, render_system_text
from planetgen.describe.model import split_by_body
from planetgen.generation.wide_binary import WideBinaryPair

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
    draw.set_run_seed(hash(tuple(sorted(overrides.items()))) & 0xFFFF)
    conn = store.get_connection(mysql_config)
    try:
        for _ in range(4):
            cfg = _config(**overrides)
            system = StarSystem(system_config=cfg)
            with conn:
                system_id = store.insert_star_system(conn, system, cfg)
            for fmt in ("text", "markdown"):
                assert render_system_text(conn, system_id, fmt) == render_star_system(system, fmt), fmt
    finally:
        conn.close()


@pytest.mark.parametrize("remnant_class", [BlackHole, NeutronStar])
def test_rendering_a_compact_remnant_anchored_system(mysql_config, remnant_class):
    cfg = SystemConfig()
    system = StarSystem(system_config=cfg, compact_remnant=remnant_class(cfg))
    system_id = store.save_phenomenon(
        system, cfg, "black-hole" if remnant_class is BlackHole else "neutron-star", config=mysql_config,
    )
    conn = store.get_connection(mysql_config)
    try:
        for fmt in ("text", "markdown"):
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
    # The stability windows are worked out from the masses, which were just swapped by hand, so work them out again (the reload does).
    system.build_ledger()
    return system, cfg


@pytest.mark.parametrize("wide", [False, True])
def test_binary_orientation_survives_the_round_trip(mysql_config, wide):
    system, cfg = _swapped_pair(mysql_config, wide)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        for fmt in ("text", "markdown"):
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
    system_id = store.save_system(StarSystem(system_config=cfg), cfg, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("UPDATE star_systems SET name = 'Somewhere Else' WHERE id = ?", (system_id,))
        assert render_system_text(conn, system_id, "markdown").startswith("# Somewhere Else")
        with pytest.raises(ValueError):
            render_system_text(conn, system_id, "html")
    finally:
        conn.close()


def test_sections_split_the_page_per_body(mysql_config):
    cfg = _config(BINARY_SYSTEM=False, MOONS=True, COMETS=True, ASTEROID_BELT=True)
    system = StarSystem(system_config=cfg)
    system_id = store.save_system(system, cfg, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        sections = split_by_body(describe_system(conn, system_id))
    finally:
        conn.close()

    for key in ("stars", "planets", "moons", "belts", "comets"):
        for blocks in sections[key].values():
            assert blocks, key
    assert sections["overview"]
    planets = [p for p in system.planets if p.body_type != "a"]
    assert len(sections["planets"]) == len(planets)
    assert len(sections["moons"]) == sum(len(p.moons) for p in planets)


def test_rendering_keeps_an_anchored_black_holes_fractional_disk_temperature(mysql_config):
    # The disk temperature is a float; loading it back must not truncate
    # it (9,064.7 K renders as "9,065 K", not "9,064 K"; kept under 5
    # digits so it isn't shown in scientific notation, UX.20).
    cfg = SystemConfig()
    black_hole = BlackHole(cfg)
    black_hole.has_accretion_disk = True
    black_hole.temperature = 9064.7
    system = StarSystem(system_config=cfg, compact_remnant=black_hole)
    system_id = store.save_phenomenon(system, cfg, "black-hole", config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        for fmt in ("text", "markdown"):
            text = render_system_text(conn, system_id, fmt)
            assert "9,065 K" in text
            assert text == render_star_system(system, fmt)
    finally:
        conn.close()


