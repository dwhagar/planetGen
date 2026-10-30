"""
Regression tests over every system spec in examples/systems/*.json.

Each file is run through the full generation pipeline (star + planets +
life data + string rendering) exactly as `planetgen system --system-file`
would, and checked for basic invariants. Every file in examples/systems/
is discovered automatically, so adding a new fixture there adds it to
this suite for free -- see docs/reference/example-systems.md.

Run with: pytest src/tests/test_examples.py
"""
import glob
import os

import pytest

from planetgen.describe import build
from planetgen.generation import run_system
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen import tuning
from tests import render_support

# This file lives at src/tests/, two levels under the repo root (src
# layout), not one -- ".." twice to reach examples/systems/ at the repo
# root.
EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "examples", "systems")
EXAMPLE_FILES = sorted(glob.glob(os.path.join(EXAMPLES_DIR, "*.json")))

# A handful of independent trials per fixture, since generation is randomized.
TRIALS = 3

# The sentences that say what a *main-sequence* star's lifespan allows. These should never appear
# in a star's narrative text for a non-main-sequence (Yerkes != 'V') star --
# see supported_evolutionary_scales and describe.build._star_age_sentence.
_MAIN_SEQUENCE_ONLY_NOTE_FRAGMENTS = list(build.LIFESPAN_SENTENCES.values())


def build_system(path):
    data = run_system.load_system_file(path)
    cfg = SystemConfig()
    run_system.apply_system_file(cfg, data)
    return StarSystem(system_config=cfg)


@pytest.mark.parametrize("path", EXAMPLE_FILES, ids=[os.path.basename(p) for p in EXAMPLE_FILES])
def test_example_renders_nonempty_output(path):
    system = build_system(path)
    text = render_support.text(system)
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
    dwarf before the fix in the star sentence).
    """
    for _ in range(TRIALS):
        system = build_system(path)
        for star in system.stars:
            if star.yerkes_class in ("V", "VI"):  # a subdwarf is a main-sequence star; the note is written for both (describe/build.py)
                continue
            # Only this star's own paragraphs: a main-sequence companion
            # (GEN.54) rightly carries its own note.
            text = render_support.section_prose(system, star.name)
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
    # count_objects() defaults to every star's combined lists -- identical to
    # `len(system.planets)` alone unless this example happened to roll an
    # S-type (wide) binary or a third star (GEN.332).
    assert planet_count + belt_count == len(system.all_planets)
    assert system.habitable_class_count >= system.m_count >= 0
