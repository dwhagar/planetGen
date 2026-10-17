# tests/test_population_streams.py

"""
GEN.306: a population's draws are keyed on the world's unique ID, never on a row id. Row ids come from the id blocks
of however many workers filled the galaxy, so the same galaxy filled another way must get the same species, ages,
technology indices, governments and colors.
"""

import ast
import pathlib

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.population import model as population
from planetgen.population import tech

POPULATION = pathlib.Path(__file__).resolve().parents[1] / "planetgen" / "population"


def _civilized():
    for _ in range(200):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.INTELLIGENT_LIFE = True
        cfg.HABITABLE_WORLD = True
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        if any(population.parse_timeline(getattr(p, "evolutionary_data", None) or []) for p in system.planets):
            return system, cfg
    pytest.fail("could not generate a system with a civilization")


def _plain(config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "M5V"
    cfg.INTELLIGENT_LIFE = False
    cfg.BINARY_SYSTEM = False
    store.save_system(StarSystem(system_config=cfg), cfg, config=config)


def _populate(config, system, cfg, plain_first):
    """Saves `system` after `plain_first` other systems (which take the low row ids), gives its planets fixed
    unique IDs, and runs the population pass until its civilizations have polities. Returns what the pass drew,
    keyed by the homeworld's unique ID."""
    for _ in range(plain_first):
        _plain(config)
    system_id = store.save_system(system, cfg, config=config)
    conn = store.get_connection(config)
    try:
        with conn:
            planets = [row["id"] for row in conn.execute(
                "SELECT id FROM planets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()]
            for number, planet_id in enumerate(planets, start=1):
                conn.execute("UPDATE planets SET uid = ? WHERE id = ?", (number.to_bytes(12, "big"), planet_id))
        population.run_pass(conn)
        with conn:
            conn.execute("UPDATE species SET civilization_age_years = 50000")
        population.run_pass(conn)
        found = {}
        for row in conn.execute(
                "SELECT s.*, pl.uid AS uid, pol.government, pol.color FROM species s "
                "JOIN planets pl ON pl.id = s.homeworld_planet_id LEFT JOIN polities pol ON pol.species_id = s.id"
        ).fetchall():
            found[bytes(row["uid"])] = {
                key: row[key] for key in ("name", "build", "climate", "size", "life_stage", "era", "government", "color",
                                          "tech_level", *(f"{domain}_index" for domain in tech.DOMAINS))}
        return found, [row["id"] for row in conn.execute("SELECT id FROM planets WHERE uid IS NOT NULL").fetchall()]
    finally:
        conn.close()


def test_the_same_world_gets_the_same_species_whatever_its_row_ids(mysql_config, fresh_mysql_config):
    system, cfg = _civilized()
    first, first_ids = _populate(mysql_config, system, cfg, plain_first=0)
    second, second_ids = _populate(fresh_mysql_config, system, cfg, plain_first=3)
    assert first and first == second
    assert first_ids != second_ids, "the second galaxy must have put the world at other row ids"
    assert all(value["government"] and value["color"] and value["tech_level"] is not None for value in first.values())


def test_no_population_stream_is_seeded_from_a_row_id():
    problems = []
    for path in sorted(POPULATION.glob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and getattr(node.func, "attr", None) in ("Stream", "bound"):
                text = ast.unparse(node)
                if "_id" in text or "['id']" in text or '["id"]' in text or " * 7919" in text:
                    problems.append(f"{path.name}:{node.lineno} {text}")
    assert not problems, "\n".join(problems)
