# tests/test_wording_agrees_with_the_index.py

"""
TEST.155: the words and the habitability index say the same thing. For random worlds and the example worlds, and for
the bodies of generated systems as the Document describes them:

- a category's colour, its score and the band its sentence names agree;
- the headline gear and the equipment line are the same list;
- an airless body is never better than Red in Air;
- a body whose text claims life has life in the index;
- no quantity in the describer is written with a raw `:.Nf`.
"""

import dataclasses
import re
from pathlib import Path


import planetgen
from planetgen import describe
from planetgen.describe import model
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import equipment
from planetgen.physics import habitability as hab
from planetgen.physics import protection
from planetgen.util import draw

from tests.test_habitability import EXAMPLES
from tests.test_protection import _worlds

BAND = re.compile(r"inside the (Blue|Green|Yellow|Red) band")


def _all_worlds():
    yield from _worlds(300)
    for world, *_rest in EXAMPLES.values():
        yield world


def test_a_colour_a_score_and_the_band_a_sentence_names_agree():
    for world in _all_worlds():
        result = protection.rating(world)
        assert hab.tier(result["overall"]) in hab.TIERS
        for row in result["rows"]:
            # The stored score is rounded to hundredths, so one on a band edge may sit either side.
            if min(abs(row["score"] - edge) for edge in (0.4, 0.75, 1.0)) > 0.006:
                assert row["tier"] == hab.tier(row["score"]), row
            named = BAND.search(row["text"])
            if named:
                assert named.group(1) == row["tier"], row


def test_the_gear_clause_and_the_equipment_line_name_the_same_items():
    """The headline's gear and the stored equipment line are one list (Boss, 10-18): they come from the same levels."""
    for world in _all_worlds():
        levels = equipment.compute(world)
        clause = protection.rating(world)["gear"]
        line = equipment.line(levels)
        assert bool(clause) == bool(line)
        for name in equipment.required_names(levels) + equipment.recommended_names(levels):
            assert name in clause and name in line, (name, clause, line)


def test_an_airless_body_is_never_better_than_red_in_air():
    airless = [dataclasses.replace(world, pressure_kpa=pressure, gases={})
               for world in _worlds(40) for pressure in (0.0, 1e-9, 1e-6)]
    assert all(not hab.has_atmosphere(world) for world in airless)
    for world in airless:
        row = next(r for r in protection.rating(world)["rows"] if r["domain"] == "Air")
        assert row["tier"] == "Red", row


def _rocky_bodies(system):
    for planet in system.planets + system.secondary_planets:
        if planet.body_type == "a":
            continue
        yield from (planet, *planet.moons)


def _life_text(section):
    return " ".join(describe.quantity.inline_text(sentence) for block in section["blocks"]
                    if block["type"] == "prose" and block["role"] in ("life", "timeline") for sentence in block["text"])


def test_a_generated_body_claims_life_only_when_the_index_gives_it_life():
    claims = 0
    for seed in range(25):
        draw.set_run_seed(seed)
        config = SystemConfig()
        config.MOONS = True
        system = StarSystem(system_config=config)
        document = describe.describe(system)
        by_name = {section["title"]: section for section in model.sections(document)}
        for body in _rocky_bodies(system):
            if body.body_type == "g":
                continue
            text = _life_text(by_name[body.name])
            if not text:
                continue
            claims += 1
            assert body.phi_bio is not None and body.phi_bio > 0, (body.name, body.phi_bio, text)
    assert claims > 0


def test_the_ratings_block_of_a_generated_body_equals_the_stored_scores():
    seen = 0
    for seed in range(12):
        draw.set_run_seed(seed)
        system = StarSystem(system_config=SystemConfig())
        by_name = {section["title"]: section for section in model.sections(describe.describe(system))}
        for body in _rocky_bodies(system):
            scored = body.habitability_world()
            ratings = [b for b in by_name[body.name]["blocks"] if b["type"] == "ratings"]
            if scored is None:
                assert not ratings
                continue
            (block,) = ratings
            stored = hab.phi4(scored[0] if isinstance(scored, tuple) else scored)
            assert block["overall"]["score"] == stored["PHI-4"]
            for row in block["rows"]:
                assert (row["score"], row["tier"]) == stored[row["domain"]]
            seen += 1
    assert seen > 10


def test_no_quantity_in_the_describer_is_written_with_a_raw_format():
    """A score or a share may be written to a fixed number of places; a physical quantity goes through the ladders."""
    root = Path(planetgen.__file__).parent
    allowed = ("score", "overall", "share")
    raw = re.compile(r":[,]?\.\d+[fgeE]|toFixed")
    offenders = []
    for path in [*sorted((root / "describe").glob("*.py")), root / "physics" / "protection.py"]:
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if raw.search(line) and not any(word in line for word in allowed):
                offenders.append(f"{path.name}:{number}: {line.strip()}")
    assert not offenders, "\n".join(offenders)
