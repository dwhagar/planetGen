"""
Regression tests over every system spec in examples/*.json.

Each file is run through the full generation pipeline (star + planets +
life data + string rendering) exactly as `systemGen.py --system-file`
would, and checked for basic invariants. Every file in examples/ is
discovered automatically, so adding a new fixture there adds it to this
suite for free -- see examples/EXAMPLES.md.

Run with: pytest tests/test_examples.py
"""
import glob
import os

import pytest

import systemGen
from stellarObjects.config import SystemConfig
from stellarObjects.systemData import StarSystem
from stellarObjects import program_constants

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "examples")
EXAMPLE_FILES = sorted(glob.glob(os.path.join(EXAMPLES_DIR, "*.json")))

# A handful of independent trials per fixture, since generation is randomized.
TRIALS = 3

# Substrings from STAR_EVOLUTION's evolutionary_constraint_notes that describe
# a *main-sequence* star's habitable-zone lifespan. These should never appear
# in a star's narrative text for a non-main-sequence (Yerkes != 'V') star --
# see get_star_evolutionary_profile / Star.to_paragraph_list.
_MAIN_SEQUENCE_ONLY_NOTE_FRAGMENTS = [
    note["evolutionary_constraint_notes"]
    for note in program_constants.STAR_EVOLUTION.values()
    if "evolutionary_constraint_notes" in note
]


def build_system(path):
    data = systemGen.load_system_file(path)
    cfg = SystemConfig()
    systemGen.apply_system_file(cfg, data)
    return StarSystem(system_config=cfg)


@pytest.mark.parametrize("path", EXAMPLE_FILES, ids=[os.path.basename(p) for p in EXAMPLE_FILES])
def test_example_generates_without_error(path):
    for _ in range(TRIALS):
        build_system(path)


@pytest.mark.parametrize("path", EXAMPLE_FILES, ids=[os.path.basename(p) for p in EXAMPLE_FILES])
def test_example_renders_nonempty_output(path):
    system = build_system(path)
    text = str(system)
    assert text.strip()
    assert system.star.name in text


@pytest.mark.parametrize("path", EXAMPLE_FILES, ids=[os.path.basename(p) for p in EXAMPLE_FILES])
def test_example_star_age_never_exceeds_lifespan(path):
    for _ in range(TRIALS):
        system = build_system(path)
        for star in system.stars:
            if star.lifespan == float('inf'):
                continue
            assert star.age <= star.lifespan * 1.001, (
                f"{path}: {star.name} age {star.age} Gy exceeds lifespan {star.lifespan} Gy"
            )


@pytest.mark.parametrize("path", EXAMPLE_FILES, ids=[os.path.basename(p) for p in EXAMPLE_FILES])
def test_example_narrative_never_misapplies_main_sequence_notes(path):
    """
    Regression guard: a giant/supergiant/subgiant/white-dwarf star's age
    sentence must never quote a main-sequence lifespan note that belongs to
    a different evolutionary story (this produced self-contradicting text
    like "age 9.33 Gy ... main-sequence lifespan 4.3-9.1 Gy" for a white
    dwarf before the fix in Star.to_paragraph_list).
    """
    for _ in range(TRIALS):
        system = build_system(path)
        text = str(system)
        for star in system.stars:
            if star.yerkes_class == "V":
                continue
            for fragment in _MAIN_SEQUENCE_ONLY_NOTE_FRAGMENTS:
                assert fragment not in text, (
                    f"{path}: non-main-sequence star {star.name} ({star.type}) narrative "
                    f"wrongly includes a main-sequence-only note"
                )


@pytest.mark.parametrize("path", EXAMPLE_FILES, ids=[os.path.basename(p) for p in EXAMPLE_FILES])
def test_example_object_counts_are_consistent(path):
    system = build_system(path)
    planet_count, belt_count, moon_count = system.count_objects()
    assert planet_count == system.planet_count
    assert belt_count == system.belt_count
    assert moon_count == system.moon_count
    assert planet_count + belt_count == len(system.planets)
    assert system.hab_count >= system.m_count >= 0
